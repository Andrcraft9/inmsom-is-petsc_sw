!============================================================
subroutine hh_init(hq, hqp, hu, hup,    &
                   hv, hvp, hh, hhp,    &
                   sh, shp, hq_r, hu_r, hv_r)
    use main_basin_pars
    use mpi_parallel_tools
    use basin_grid
    implicit none

    real(8) hq(bnd_x1:bnd_x2, bnd_y1:bnd_y2), hqp(bnd_x1:bnd_x2, bnd_y1:bnd_y2),    &
            hu(bnd_x1:bnd_x2, bnd_y1:bnd_y2), hup(bnd_x1:bnd_x2, bnd_y1:bnd_y2),    &
            hv(bnd_x1:bnd_x2, bnd_y1:bnd_y2), hvp(bnd_x1:bnd_x2, bnd_y1:bnd_y2),    &
            hh(bnd_x1:bnd_x2, bnd_y1:bnd_y2), hhp(bnd_x1:bnd_x2, bnd_y1:bnd_y2),    &
            hq_r(bnd_x1:bnd_x2, bnd_y1:bnd_y2),hu_r(bnd_x1:bnd_x2, bnd_y1:bnd_y2),    &
            hv_r(bnd_x1:bnd_x2, bnd_y1:bnd_y2)

    real(8) sh(bnd_x1:bnd_x2, bnd_y1:bnd_y2), shp(bnd_x1:bnd_x2, bnd_y1:bnd_y2)
    real(8) slu
    integer m,n

    hq =hq_r + sh *dfloat(full_free_surface)
    hqp=hq_r + shp*dfloat(full_free_surface)

    !$omp parallel do private(m,n,slu)
    do n=ny_start-1,ny_end
        do m=nx_start-1,nx_end

            if(llu(m,n)>0.5) then
                ! interpolating hhq given on T-grid(lu) to hhu given on u-grid(lcu).
                slu = dble(lu(m,n)+lu(m+1,n))
                hu_r(m,n)=(  hq_r(m  ,n)*dx(m  ,n)*dy(m  ,n)*lu(m  ,n)   &
                        +  hq_r(m+1,n)*dx(m+1,n)*dy(m+1,n)*lu(m+1,n) )/slu/dxt(m,n)/dyh(m,n)
                hu(m,n)=(  hq(m  ,n)*dx(m  ,n)*dy(m  ,n)*lu(m  ,n)   &
                      +  hq(m+1,n)*dx(m+1,n)*dy(m+1,n)*lu(m+1,n) )/slu/dxt(m,n)/dyh(m,n)
                hup(m,n)=( hqp(m  ,n)*dx(m  ,n)*dy(m  ,n)*lu(m  ,n)   &
                      + hqp(m+1,n)*dx(m+1,n)*dy(m+1,n)*lu(m+1,n) )/slu/dxt(m,n)/dyh(m,n)
            endif

            if(llv(m,n)>0.5) then
                ! interpolating hhq given on T-grid(lu) to hhv given on v-grid(lcv).
                slu = dble(lu(m,n)+lu(m,n+1))
                hv_r(m,n)=(  hq_r(m,n  )*dx(m,n  )*dy(m,n  )*lu(m,n  )    &
                        +  hq_r(m,n+1)*dx(m,n+1)*dy(m,n+1)*lu(m,n+1) )/slu/dxh(m,n)/dyt(m,n)
                hv(m,n)=(  hq(m,n  )*dx(m,n  )*dy(m,n  )*lu(m,n  )    &
                      +  hq(m,n+1)*dx(m,n+1)*dy(m,n+1)*lu(m,n+1) )/slu/dxh(m,n)/dyt(m,n)
                hvp(m,n)=( hqp(m,n  )*dx(m,n  )*dy(m,n  )*lu(m,n  )    &
                      + hqp(m,n+1)*dx(m,n+1)*dy(m,n+1)*lu(m,n+1) )/slu/dxh(m,n)/dyt(m,n)
            endif

            if(luh(m,n)>0.5) then
                ! interpolating hhq given on T-grid(lu) to hhh given on h-grid(luu).
                slu = dble(lu(m,n)+lu(m+1,n)+lu(m,n+1)+lu(m+1,n+1))
                hh(m,n)=(  hq(m  ,n  )*dx(m  ,n  )*dy(m  ,n  )*lu(m  ,n  )       &
                      +  hq(m+1,n  )*dx(m+1,n  )*dy(m+1,n  )*lu(m+1,n  )       &
                       + hq(m  ,n+1)*dx(m  ,n+1)*dy(m  ,n+1)*lu(m  ,n+1)       &
                      +  hq(m+1,n+1)*dx(m+1,n+1)*dy(m+1,n+1)*lu(m+1,n+1) )/slu/dxb(m,n)/dyb(m,n)
                hhp(m,n)=( hqp(m  ,n  )*dx(m  ,n  )*dy(m  ,n  )*lu(m  ,n  )       &
                      + hqp(m+1,n  )*dx(m+1,n  )*dy(m+1,n  )*lu(m+1,n  )       &
                       +hqp(m  ,n+1)*dx(m  ,n+1)*dy(m  ,n+1)*lu(m  ,n+1)       &
                      + hqp(m+1,n+1)*dx(m+1,n+1)*dy(m+1,n+1)*lu(m+1,n+1) )/slu/dxb(m,n)/dyb(m,n)
            endif

        end do
    end do
    !$omp end parallel do

    call syncborder_real8(hu_r, 1)
    call syncborder_real8(hu, 1)
    call syncborder_real8(hup, 1)
    call syncborder_real8(hv_r, 1)
    call syncborder_real8(hv, 1)
    call syncborder_real8(hvp, 1)
    call syncborder_real8(hh, 1)
    call syncborder_real8(hhp, 1)

    if(periodicity_x/=0) then
        call cyclize_x(hu_r,nx,ny,1,mmm,mm)
        call cyclize_x(hu, nx,ny,1,mmm,mm)
        call cyclize_x(hup,nx,ny,1,mmm,mm)
        call cyclize_x(hv_r,nx,ny,1,mmm,mm)
        call cyclize_x(hv, nx,ny,1,mmm,mm)
        call cyclize_x(hvp,nx,ny,1,mmm,mm)
        call cyclize_x(hh, nx,ny,1,mmm,mm)
        call cyclize_x(hhp,nx,ny,1,mmm,mm)
    end if

    if(periodicity_y/=0) then
        call cyclize_y(hu_r,nx,ny,1,nnn,nn)
        call cyclize_y(hu, nx,ny,1,nnn,nn)
        call cyclize_y(hup,nx,ny,1,nnn,nn)
        call cyclize_y(hv_r,nx,ny,1,nnn,nn)
        call cyclize_y(hv, nx,ny,1,nnn,nn)
        call cyclize_y(hvp,nx,ny,1,nnn,nn)
        call cyclize_y(hh, nx,ny,1,nnn,nn)
        call cyclize_y(hhp,nx,ny,1,nnn,nn)
    end if

endsubroutine hh_init

!============================================================
subroutine hh_update(hq, hqp, hu, hup,    &
                     hv, hvp, hh, hhp,    &
                     sh, hq_r)
    use main_basin_pars
    use mpi_parallel_tools
    use basin_grid
    implicit none

    real(8) hq(bnd_x1:bnd_x2, bnd_y1:bnd_y2), hqp(bnd_x1:bnd_x2, bnd_y1:bnd_y2),    &
            hu(bnd_x1:bnd_x2, bnd_y1:bnd_y2), hup(bnd_x1:bnd_x2, bnd_y1:bnd_y2),    &
            hv(bnd_x1:bnd_x2, bnd_y1:bnd_y2), hvp(bnd_x1:bnd_x2, bnd_y1:bnd_y2),    &
            hh(bnd_x1:bnd_x2, bnd_y1:bnd_y2), hhp(bnd_x1:bnd_x2, bnd_y1:bnd_y2),    &
            hq_r(bnd_x1:bnd_x2, bnd_y1:bnd_y2)

    real(8) sh(bnd_x1:bnd_x2, bnd_y1:bnd_y2)
    real(8) slu
    integer m,n

    hqp=hq
    hup=hu
    hvp=hv
    hhp=hh

    hq =hq_r + sh

    !$omp parallel do private(m,n,slu)
    do n=ny_start-1,ny_end
        do m=nx_start-1,nx_end

        if(llu(m,n)>0.5) then
            ! interpolating hhq given on T-grid(lu) to hhu given on u-grid(lcu).
            slu = dble(lu(m,n)+lu(m+1,n))
            hu(m,n)=( hq(m  ,n)*dx(m  ,n)*dy(m  ,n)*lu(m  ,n)   &
                  + hq(m+1,n)*dx(m+1,n)*dy(m+1,n)*lu(m+1,n) )/slu/dxt(m,n)/dyh(m,n)
        endif

        if(llv(m,n)>0.5) then
            ! interpolating hhq given on T-grid(lu) to hhv given on v-grid(lcv).
            slu = dble(lu(m,n)+lu(m,n+1))
            hv(m,n)=( hq(m,n  )*dx(m,n  )*dy(m,n  )*lu(m,n  )   &
                  + hq(m,n+1)*dx(m,n+1)*dy(m,n+1)*lu(m,n+1) )/slu/dxh(m,n)/dyt(m,n)
        endif

        if(luh(m,n)>0.5) then
            ! interpolating hhq given on T-grid(lu) to hhh given on h-grid(luu).
            slu = dble(lu(m,n)+lu(m+1,n)+lu(m,n+1)+lu(m+1,n+1))
            hh(m,n)=( hq(m  ,n  )*dx(m  ,n  )*dy(m  ,n  )*lu(m  ,n  )       &
                   + hq(m+1,n  )*dx(m+1,n  )*dy(m+1,n  )*lu(m+1,n  )       &
                    +hq(m  ,n+1)*dx(m  ,n+1)*dy(m  ,n+1)*lu(m  ,n+1)       &
                   + hq(m+1,n+1)*dx(m+1,n+1)*dy(m+1,n+1)*lu(m+1,n+1) )/slu/dxb(m,n)/dyb(m,n)
        endif

        end do
    end do
    !$omp end parallel do

    call syncborder_real8(hu, 1)
    call syncborder_real8(hv, 1)
    call syncborder_real8(hh, 1)

    if(periodicity_x/=0) then
        call cyclize_x(hu, nx,ny,1,mmm,mm)
        call cyclize_x(hv, nx,ny,1,mmm,mm)
        call cyclize_x(hh, nx,ny,1,mmm,mm)
    end if

    if(periodicity_y/=0) then
        call cyclize_y(hu, nx,ny,1,nnn,nn)
        call cyclize_y(hv, nx,ny,1,nnn,nn)
        call cyclize_y(hh, nx,ny,1,nnn,nn)
    end if

endsubroutine hh_update
