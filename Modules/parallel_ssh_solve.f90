!================================================
! modules for sea level calculations prepered by rusakov.
! no splitting here. we form full matrice and solve it
!======================================================================
    module parallel_sea_level_no_split
#include <petsc/finclude/petscksp.h>

        use main_basin_pars
        use mpi_parallel_tools
        use basin_grid
        implicit none

        !======================================================================
        !      declarations
        !======================================================================
        PetscInt  istart, iend, its
        PetscReal norm
        Mat       matrix
        Vec       rhs
        Vec       sol
        KSP       ksp

        integer nvar, locnozero
        integer*4, allocatable :: uindx(:,:),vindx(:,:),slhindx(:,:)

        real*8  :: a(7)
        integer :: posi(1), ja(7)

        real*8  :: valsol(3), valrhs(3), val(7)
        integer :: pos(3)

        real*8, pointer :: retarray(:)

        integer, parameter :: nzmax = 10000000,     &   !size of data storage for precondition
                            precond = 0                 !precond = 0 no preconditioner

        PetscInt nnz
        !      integer nnz    ! number of nonzeroes

        integer maxnnz   ! size of a and ja

        !
        !     arrays of indexes. connect two dimensional problem with one dimensional
        !      vector excluding the ground.
        !
        integer  is_formed, matrice_permuted, preconditioner_formed

        integer, parameter::debug = 0,   & !set debug to 1 to see some useful information
                              do_reordering_for_sparsity = 1

        integer, save :: factors(8),    & ! array of pointer for super lu
                           isfactored  = 0

        !======================================================================
        !     end of declarations
        !======================================================================

        contains

        !======================================================================
        !     subroutines
        !======================================================================
        subroutine estimate_matrix_size()
            implicit none
            integer m, n, k, k1, k2, k3
            integer, allocatable :: nozero(:), offset(:)
            integer ierr
            !----------------------------------------------------------------------
            !
            !      estimate matrix size, form numbering of variables.
            !
            !----------------------------------------------------------------------
            allocate(nozero(procs), offset(procs))
            locnozero = 0
            do n = max(ny_start,nnn),min(ny_end,nn)
                do m = max(nx_start,mmm),min(nx_end,mm)
                    if ( lcu(m,n).gt.0.5 ) locnozero = locnozero + 1
                    if ( lcv(m,n).gt.0.5 ) locnozero = locnozero + 1
                    if ( lu (m,n).gt.0.5 ) locnozero = locnozero + 1
                enddo
            enddo
            call mpi_allgather(locnozero,1,mpi_integer,   &
                                    nozero,1,mpi_integer, &
                                      cart_comm, ierr)
            nvar = sum(nozero)
            offset(1) = 0
            do m = 1, procs-1
                offset(m+1) = offset(m) + nozero(m)
            enddo

            !print *, rank, locnozero, offset(rank+1), nvar

            !     construct indexes
            allocate(uindx(bnd_x1:bnd_x2, bnd_y1:bnd_y2))
            allocate(vindx(bnd_x1:bnd_x2, bnd_y1:bnd_y2))
            allocate(slhindx(bnd_x1:bnd_x2, bnd_y1:bnd_y2))
            k = offset(rank + 1)
            do n = max(ny_start,nnn),min(ny_end,nn)
                do m = max(nx_start,mmm),min(nx_end,mm)
                    if ( lcu(m,n).gt.0.5 ) then
                        uindx(m,n) = k
                        k = k + 1
                    else
                        uindx(m,n) = 0
                    endif
                    if ( lcv(m,n).gt.0.5 ) then
                        vindx(m,n) = k
                        k = k + 1
                    else
                        vindx(m,n) = 0
                    endif
                    if ( lu(m,n).gt.0.5 ) then
                        slhindx(m,n) = k
                        k = k + 1
                    else
                        slhindx(m,n) = 0
                    endif
                enddo
            enddo

            !write(*,*) '|pcoord:',p_coord,'|  ', '|rank:',rank,'|  ', &
            !                 '|nx:',nx_start,nx_end, '|  ',                 &
            !                 '|ny:',ny_start,ny_end, '|  '

            !write(*,*) '|pcoord:',p_coord,'|  ', '|rank:',rank,'|  ', &
            !                 '|b_x:',bnd_x1,bnd_x2, '|  ',                  &
            !                 '|b_y:',bnd_y1,bnd_y2, '|  '

            !call petscfinalize(ierr)
            !call mpi_finalize(ierr)
            !stop

            !write(*,*)rank,'nx_str:',uindx(nx_start,ny_start:ny_end)
            !write(*,*)rank,'nx_end:',uindx(nx_end,ny_start:ny_end)
            !write(*,*)rank,'ny_str:',uindx(nx_start:nx_end,ny_start)
            !write(*,*)rank,'ny_end:',uindx(nx_start:nx_end,ny_end)

            call syncborder_int4(uindx, 1)
            call syncborder_int4(vindx, 1)
            call syncborder_int4(slhindx, 1)

            !write(*,*)rank,'bnd_x1:',uindx(bnd_x1,ny_start:ny_end)
            !write(*,*)rank,'bnd_x2:',uindx(bnd_x2,ny_start:ny_end)
            !write(*,*)rank,'bnd_y1:',uindx(nx_start:nx_end,bnd_y1)
            !write(*,*)rank,'bnd_y2:',uindx(nx_start:nx_end,bnd_y2)

            !print *, vindx(max(nx_start,mmm):min(nx_end,mm),          &
            !                    max(ny_start,nnn):min(ny_end,nn))
            !print *, rank, max(nx_start,mmm),min(nx_end,mm),          &
            !                    max(ny_start,nnn),min(ny_end,nn)
            !call petscfinalize(ierr)
            !stop
        end subroutine estimate_matrix_size

        !----------------------------------------------------------------------
        !     form matrice of sea level task
        !
        !    du/dt + rb*u - l*v = mg d sl/dx
        !
        !    du/dt + rb*v + l*u = ng d sl/dx
        !
        !             sl          d            d   n
        ! freesurface*-----  - [  --(hu) +    --(  - * hv ) ] = 0
        !              dt         dx           dy  m
        !
        !   implicit time scheme:
        !
        !           u-uo                             d slh
        !           ----  + rb * u - l * v = m*g*h * -----
        !           tau                              d x
        !
        !
        !           v-vo                             d slh
        !           ----  + rb * v + l * u = n*g*h * -----     (1)
        !           tau                              d y
        !
        !             slh - slho       d          d   n
        ! freesurface*----------  - [  --(u) +    --( - * v ) ] = 0
        !              m*tau           dx         dy  m
        !
        !********************************************************************
        !      remark:   freesurface =0 - rigid lid condition
        !                freesurface =1 - long gravity waves are avaliable
        !----------------------------------------------------------------------
        !
        !     form matrice in compress sparse row format (csr).
        !    use common ordering. the first is u, the second is v , the third is slh.
        !    index m -changed first.
        !
        !  side effect
        !     uncyclize operation on lcu and lcv
        !
        !----------------------------------------------------------------------
        subroutine form_matrice_ssh(tau, rbottom)
            implicit none
            integer :: m, n, k, k1, k2, k3
            real*8  :: tau, slcu, slcv
            real*8  :: rbottom(bnd_x1:bnd_x2, bnd_y1:bnd_y2)
            real*8  :: frictau
            integer :: ileft, iright
            integer :: ierr
            real*8 :: w_time

            call start_timer(w_time)
            call estimate_matrix_size
            call end_timer(w_time)
            if (rank .eq. 0) print *, "Estimate matrix size: ", w_time

            preconditioner_formed = 0
            is_formed = 0
            matrice_permuted = 0

            call start_timer(w_time)
            call MatCreate(PETSc_Comm_World, matrix, ierr)
            call MatSetSizes(matrix, locnozero, locnozero, nvar, nvar, ierr)
            call MatSetFromOptions(matrix, ierr)
            call MatSetup(matrix, ierr)
            call end_timer(w_time)
            if (rank .eq. 0) print *, "Matrix setup: ", w_time

            !call matgetownershiprange(matrix,istart,iend,ierr)
            !print *, rank, istart, iend

            !call petscfinalize(ierr)
            !stop

            call start_timer(w_time)
            do n = ny_start, ny_end
                do m = nx_start, nx_end
                    !----------------------------------------------------------------------
                    !         first equation
                    !----------------------------------------------------------------------
                    !         something like: see reports
                    !
                    !           u-uo                             d slh
                    !           ----  + rb * u - l * v = m*g*h * -----
                    !           tau                              d x
                    !
                    !         coriolis on luh grid
                    !----------------------------------------------------------------------
                    nnz = 0
                    ileft  = m-1
                    iright = m+1
                    if( periodicity_x .ne. 0 ) then
                        if (m.eq.mmm) then
                            ileft  = mm
                            iright = m+1
                        endif
                        if (m.eq.mm) then   !cyclic right boundary
                            ileft  = m-1
                            iright = mmm
                        end if
                    end if !if cyclic condition

                    if ( lcu(m,n).gt.0.5 ) then
                        !diagonal
                        posi(1) = uindx(m, n)
                        nnz = nnz + 1
                        frictau = 1.0d0/dble(tau) + rbottom(m,n)
                        a(nnz)  = frictau
                        ja(nnz) = uindx(m,n)
                        !ia(uindx(m,n)) = nnz

                        slcv=4.0d0
                        !v elements
                        if ( lcv(m,n).gt.0.5 ) then
                            nnz = nnz + 1
                            a(nnz)  = -dble(rlh_s(m,n)/slcv)
                            ja(nnz) = vindx(m,n)
                        endif
                        if ( lcv(m+1,n).gt.0.5.and.lcv(iright,n).gt.0.5 ) then
                            nnz = nnz + 1
                            a(nnz)  = -dble(rlh_s(m,n)/slcv)
                            ja(nnz) = vindx(iright,n)
                        endif
                        if ( lcv(m,n-1).gt.0.5 ) then
                            nnz = nnz + 1
                            a(nnz)  = -dble(rlh_s(m,n-1)/slcv)
                            ja(nnz) = vindx(m,n-1)
                        endif
                        if ( lcv(m+1,n-1).gt.0.5.and.lcv(iright,n-1).gt.0.5 ) then
                            nnz = nnz + 1
                            a(nnz)  = -dble(rlh_s(m,n-1)/slcv)
                            ja(nnz) = vindx(iright,n-1)
                        endif
                        !slh elements
                        if ( lu(m,n).gt.0.5 ) then
                            nnz = nnz + 1
                            !a(nnz)  = dble(grv)/dble(dxt(m,n)*rn)*sweq_imp
                            a(nnz)  = -dble(FreeFallAcc*hhu_rest(m,n))/dble(dxt(m,n))
                            ja(nnz) = slhindx(m,n)
                        endif
                        if ( lu(m+1,n).gt.0.5 .and. lu(iright,n).gt.0.5 ) then
                            nnz = nnz + 1
                            !a(nnz)  =-dble(grv)/dble(dxt(m,n)*rn)*sweq_imp
                            a(nnz)  = dble(FreeFallAcc*hhu_rest(m,n))/dble(dxt(m,n))
                            ja(nnz) = slhindx(iright,n)
                        endif
                        call MatSetValues(matrix, 1, posi, nnz, ja, a, insert_values, ierr)
                    endif
                    !----------------------------------------------------------------------
                    !         second equation
                    !----------------------------------------------------------------------
                    !         something like: see reports
                    !
                    !           v-vo                             d slh
                    !           ----  + rb * v + l * u = n*g*h * -----     (1)
                    !           tau                              d y
                    !
                    !         coriolis on luh grid
                    !----------------------------------------------------------------------
                    nnz = 0
                    ileft  = m-1
                    iright = m+1
                    if(periodicity_x .gt. 0) then
                        if (m.eq.mmm) then !cyclic left boundary
                            ileft  = mm
                            iright = m+1
                        end if
                        if (m.eq.mm) then   !cyclic right boundary
                            ileft  = m-1
                            iright = mmm
                        end if
                    end if !if cyclic condition

                    if ( lcv(m,n).gt.0.5 ) then
                        !diagonal
                        posi(1) = vindx(m, n)
                        nnz = nnz + 1
                        frictau = 1.0d0/dble(tau) + rbottom(m,n)
                        a(nnz)  = frictau
                        ja(nnz) = vindx(m,n)
                        !ia(vindx(m,n)) = nnz

                        slcu=4.0d0
                        !u elements
                        if ( lcu(m,n).gt.0.5 ) then
                            nnz = nnz + 1
                            a(nnz)  = dble(rlh_s(m,n)/slcu)
                            ja(nnz) = uindx(m,n)
                        endif
                        if ( lcu(m-1,n).gt.0.5 ) then
                            if ( lcu(ileft,n).gt.0.5 ) then
                                nnz = nnz + 1
                                a(nnz)  = dble(rlh_s(m-1,n)/slcu)
                                ja(nnz) = uindx(ileft,n)
                            endif
                        endif
                        if ( lcu(m,n+1).gt.0.5 ) then
                            nnz = nnz + 1
                            a(nnz)  = dble(rlh_s(m,n)/slcu)
                            ja(nnz) = uindx(m,n+1)
                        endif
                        if ( lcu(m-1,n+1).gt.0.5 ) then
                            if ( lcu(ileft,n+1).gt.0.5 ) then
                                nnz = nnz + 1
                                a(nnz)  = dble(rlh_s(m-1,n)/slcu)
                                ja(nnz) = uindx(ileft,n+1)
                            end if
                        endif

                        !slh elements
                        if ( lu(m,n).gt.0.5 ) then
                            nnz = nnz + 1
                            !a(nnz)  = dble(grv)/dble(dyt(m,n)*rn)*sweq_imp
                            a(nnz)  = -dble(FreeFallAcc*hhv_rest(m,n))/dble(dyt(m,n))
                            ja(nnz) = slhindx(m,n)
                        endif
                        if ( lu(m,n+1).gt.0.5 ) then
                            nnz = nnz + 1
                            !a(nnz)  =-dble(grv)/dble(dyt(m,n)*rn)*sweq_imp
                            a(nnz)  = dble(FreeFallAcc*hhv_rest(m,n))/dble(dyt(m,n))
                            ja(nnz) = slhindx(m,n+1)
                        endif
                        call MatSetValues(matrix, 1, posi, nnz, ja, a, insert_values, ierr)
                    endif
                    !----------------------------------------------------------------------
                    !         third equation
                    !----------------------------------------------------------------------
                    !         something like: see reports
                    !
                    !
                    !             slh - slho                d           d    h
                    !  freesurface*----------  - m(ij)n ([  --(hu/n) +   --(  - * v ) ] = 0
                    !                  tau                 dx           dy   m
                    !
                    !         coriolis on luh grid
                    !----------------------------------------------------------------------
                    nnz = 0
                    ileft  = m-1
                    iright = m+1
                    if(periodicity_x .gt. 0) then
                        if (m.eq.mmm) then   !cyclic right boundary
                            ileft  = mm
                            iright = m+1
                        end if
                        if (m.eq.mm) then   !cyclic right boundary
                            ileft  = m-1
                            iright = mmm
                        end if
                    end if !if cyclic condition

                    if ( lu(m,n).gt.0.5 ) then
                        !diagonal
                        posi(1) = slhindx(m, n)
                        nnz = nnz + 1
                        a(nnz)  = 1.0d0/dble(tau) !the same have to be in rhs
                        ja(nnz) = slhindx(m,n)
                        !ia(slhindx(m,n)) = nnz

                        !u elements nan
                        if ( lcu(m,n).gt.0.5 ) then
                            nnz = nnz + 1
                            !a(nnz)  = -dble(hhu(m,n))*dble(dyh(m,n))*sweq_imp            &
                            !                        /(dble(dx(m,n))*dble(dy(m,n)*rn))
                            a(nnz)  = dble(dyh(m,n))/(dble(dx(m,n)*dy(m,n)))
                            ja(nnz) = uindx(m,n)
                        endif
                        if ( lcu(m-1,n).gt.0.5.and.lcu(ileft,n).gt.0.5 ) then
                            nnz = nnz + 1
                            !a(nnz)  =  dble(hhu(m-1,n))*dble(dyh(m-1,n))*sweq_imp        &
                            !                       /(dble(dx(m,n))*dble(dy(m,n)*rn))
                            a(nnz)  = -dble(dyh(m-1,n))/(dble(dx(m,n)*dy(m,n)))
                            ja(nnz) = uindx(ileft,n)
                        endif

                        !v elements
                        if ( lcv(m,n).gt.0.5 ) then
                            nnz = nnz + 1
                            !a(nnz)  = -dble(hhv(m,n))*dble(dxh(m,n))*sweq_imp            &
                            !                          /(dble(dx(m,n))*dble(dy(m,n)*rn))
                            a(nnz)  = dble(dxh(m,n))/(dble(dx(m,n)*dy(m,n)))
                            ja(nnz) = vindx(m,n)
                        endif
                        if ( lcv(m,n-1).gt.0.5 ) then
                            nnz = nnz + 1
                            !a(nnz)  = dble(hhv(m,n-1))*dble(dxh(m,n-1))*sweq_imp         &
                            !                          /(dble(dx(m,n))*dble(dy(m,n)*rn))
                            a(nnz)  = -dble(dxh(m,n-1))/(dble(dx(m,n)*dy(m,n)))
                            ja(nnz) = vindx(m,n-1)
                        endif
                        call MatSetValues(matrix, 1, posi, nnz, ja, a, insert_values, ierr)
                    endif
                enddo
            enddo
            call end_timer(w_time)
            if (rank .eq. 0) print *, "MatSetValues: ", w_time

            call start_timer(w_time)
            call MatAssemblyBegin(matrix, mat_final_assembly, ierr)
            call MatAssemblyEnd(matrix, mat_final_assembly, ierr)
            call end_timer(w_time)
            if (rank .eq. 0) print *, "MatAssembly: ", w_time

            !call matview(matrix, petsc_viewer_stdout_world)

            !ia(nvar+1) = nnz+1
            is_formed = 1

            !call petscfinalize(ierr)
            !stop

            call start_timer(w_time)
            call KSPCreate(PETSc_Comm_World, ksp, ierr)
            call KSPSetOperators(ksp, matrix, matrix, ierr)
            call KSPSetFromOptions(ksp, ierr)
            call KSPSetInitialGuessNonZero(ksp, petsc_true, ierr)
            !call kspsettype(ksp, kspgmres)
            call end_timer(w_time)
            if (rank .eq. 0) print *, "KSP create, setting: ", w_time

            !print *, 'iniit ksp'
            !call kspsetup(ksp, ierr)
            !call kspsetinitialguessnonzero(ksp,petsc_true,ierr)

        end subroutine form_matrice_ssh

        !----------------------------------------------------------------------
        !     form rhs
        !     f1 - on lcu grid corresponds to equation for u
        !     f2 - on lcv grid corresponds to equation for v
        !     u0,v0,slh0 are changed on exit ( 1/tau*[u0,v0,slh0] are added)
        !     rhs = [f1,f2,0] + 1/tau*[u0,v0, freesurface * slh0]
        !----------------------------------------------------------------------
        subroutine form_rhs(u0, v0, ssh0, f1, f2, wf, tau)
            implicit none

            real*8 :: f1(bnd_x1:bnd_x2, bnd_y1:bnd_y2), f2(bnd_x1:bnd_x2, bnd_y1:bnd_y2)
            real*8 :: u0(bnd_x1:bnd_x2, bnd_y1:bnd_y2), v0(bnd_x1:bnd_x2, bnd_y1:bnd_y2)
            real*8 :: ssh0(bnd_x1:bnd_x2, bnd_y1:bnd_y2)
            real*8 :: wf(bnd_x1:bnd_x2, bnd_y1:bnd_y2)
            real*8 :: tau
            real*8 :: tau1, shu, shv
            integer :: m, n, ierr

            !tau1 = dble(1.0/tau)
            tau1 = 1.0d0 / tau

            call VecCreateMPI(PETSc_Comm_World, locnozero, nvar, sol, ierr)
            call VecCreateMPI(PETSc_Comm_World, locnozero, nvar, rhs, ierr)

            do n = ny_start, ny_end
                do m = nx_start, nx_end
                    nnz = 0
                    !----------------------------- 1-st equation -----------------------------------
                    if (lcu(m,n).gt.0.5) then
                        nnz = nnz + 1
                        pos(nnz) = uindx(m,n)
                        shu=(  ssh0(m  ,n)*dx(m  ,n)*dy(m  ,n)   &
                            +  ssh0(m+1,n)*dx(m+1,n)*dy(m+1,n) )/dble(2.0*dxt(m,n)*dyh(m,n))
                        valrhs(nnz) = f1(m,n) + u0(m,n)*hhu(m,n)*tau1        &
                            - dfloat(full_free_surface)*FreeFallAcc*shu*(ssh0(m+1,n)-ssh0(m,n))/dxt(m,n)
                        valsol(nnz) = u0(m,n)
                    endif
                    !----------------------------- 2-nd equation -----------------------------------
                    if (lcv(m,n).gt.0.5) then
                        nnz = nnz + 1
                        pos(nnz) = vindx(m,n)
                        shv=(  ssh0(m,n  )*dx(m,n  )*dy(m,n  )   &
                            +   ssh0(m,n+1)*dx(m,n+1)*dy(m,n+1) )/dble(2.0*dxh(m,n)*dyt(m,n))
                        valrhs(nnz) = f2(m,n) + v0(m,n)*hhv(m,n)*tau1        &
                            - dfloat(full_free_surface)*FreeFallAcc*shv*(ssh0(m,n+1)-ssh0(m,n))/dyt(m,n)
                        valsol(nnz) = v0(m,n)
                    endif
                    !----------------------------- 3-rd equation -----------------------------------
                    if (lu(m,n).gt.0.5) then
                        nnz = nnz + 1
                        pos(nnz) = slhindx(m,n)
                        valrhs(nnz) = ssh0(m,n) * tau1 + wf(m,n)/RefDen*dfloat(full_free_surface)
                        valsol(nnz) = ssh0(m,n)
                    endif

                    call VecSetValues(rhs, nnz, pos, valrhs, insert_values, ierr)
                    call VecSetValues(sol, nnz, pos, valsol, insert_values, ierr)
                enddo
            enddo

            call VecAssemblyBegin(rhs, ierr)
            call VecAssemblyEnd(rhs, ierr)

            call VecAssemblyBegin(sol, ierr)
            call VecAssemblyEnd(sol, ierr)

            !call vecview(rhs, petsc_viewer_stdout_world)
            !call vecview(sol, petsc_viewer_stdout_world)

        end subroutine form_rhs
        !--------------------------------------------------------------

        subroutine solve_system(u, v, ssh)
            implicit none
            real*8, intent(in out) :: u(bnd_x1:bnd_x2, bnd_y1:bnd_y2), v(bnd_x1:bnd_x2, bnd_y1:bnd_y2)
            real*8, intent(in out) :: ssh(bnd_x1:bnd_x2, bnd_y1:bnd_y2)

            integer :: ierr
            integer :: n, m, k
            real*8 :: w_time

            call start_timer(w_time)
            call KSPSolve(ksp, rhs, sol, ierr)
            call end_timer(w_time)

            call KSPGetIterationNumber(ksp, its, ierr)
            call KSPGetResidualNorm(ksp, norm, ierr)

            !call kspview(ksp,petsc_viewer_stdout_world,ierr)

            if (rank.eq.0) then
                print *, rank, 'kspsolve', w_time, 'iters', its, 'norm', norm, 'nvar', nvar
            endif
            !call vecview(sol, petsc_viewer_stdout_world)
            !print *, "its:", its, "norm:", norm

            call start_timer(w_time)
            call VecGetArrayf90(sol, retarray, ierr)
            u = 0.0d0; v = 0.0d0; ssh = 0.0d0
            k = 0
            do n = ny_start, ny_end
                do m = nx_start, nx_end
                    if ( lcu(m,n).gt.0.5 ) then
                        u(m,n) = retarray(k + 1)
                        k = k + 1
                    endif
                    if ( lcv(m,n).gt.0.5 ) then
                        v(m,n) = retarray(k + 1)
                        k = k + 1
                    endif
                    if ( lu(m,n).gt.0.5 ) then
                        ssh(m,n) = retarray(k + 1)
                        k = k + 1
                    endif
                enddo
            enddo
            call VecRestoreArrayf90(sol, retarray, ierr)
            call end_timer(w_time)
            if (rank .eq. 0) print *, 'form u, v, ssh:', w_time

            !call start_timer(w_time)
            call syncborder_real8(u, 1)
            call syncborder_real8(v, 1)
            call syncborder_real8(ssh, 1)
            if(periodicity_x/=0) then
                call cyclize8_x(u, nx,ny,1,mmm,mm)
                call cyclize8_x(v, nx,ny,1,mmm,mm)
                call cyclize8_x(ssh, nx,ny,1,mmm,mm)
            end if
            if(periodicity_y/=0) then
                call cyclize8_y(u, nx,ny,1,nnn,nn)
                call cyclize8_y(v, nx,ny,1,nnn,nn)
                call cyclize8_y(ssh, nx,ny,1,nnn,nn)
            end if
            !call end_timer(w_time)
            !if (rank.eq.0) print *, '  sync:', w_time

            call VecDestroy(sol, ierr)
            call VecDestroy(rhs, ierr)

        end subroutine solve_system

        !--------------------------------------------------------------
        !
    end module parallel_sea_level_no_split

    !  additional function
    !-----------------------------------------------------------------------
    !function distdot(n,x,ix,y,iy)
    !    integer n, ix, iy
    !    real*8 distdot, x(*), y(*), ddot
    !    external ddot
    !    distdot = ddot(n,x,ix,y,iy)
    !    return
    !    !-----end-of-distdot
    !    !-----------------------------------------------------------------------
    !end function distdot
    !======================================================================
