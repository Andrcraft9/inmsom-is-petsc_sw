!--------------------- SUBROUTINE FOR: -----------------------------------------!
!----------------- Only shallow water solving ----------------------------------!
!-------------------------------------------------------------------------------!
subroutine shallow_water_model_step(tau)
    use main_basin_pars
    use mpi_parallel_tools
    use basin_grid
    use ocean_variables
    use parallel_sea_level_no_split

    use time_integration

    implicit none
    integer :: m, n, k, ierr
    integer*8 :: curstep
    real*8 :: tau, diffslpr
    real*8 :: time_count

    diffslpr = 1.0d+3
    !diffslpr = 0.0d0

!---------------------- Shallow water equ solver -------------------------------

    if (atm_forcing_on == 1) then
        if (rank .eq. 0 ) print *, "ATM MOD NOT READY!"
        call finalize_parallel(ierr)
        call exit(0)
    else
        wf_tot = 0.0d0

        !$omp parallel do private(m,n)
        do n = ny_start, ny_end
              do m = nx_start, nx_end
                   if(lu(m,n)>0.5) then
                       RHSx2d(m,n)= -(diffslpr)*hhu(m,n)/dxt(m,n)/RefDen
                       RHSy2d(m,n)= -(diffslpr)*hhv(m,n)/dyt(m,n)/RefDen
                   endif
              end do
       	end do
       !$omp end parallel do
    endif

    call form_rhs(ubrtr, vbrtr, ssh, RHSx2d, RHSy2d, wf_tot, tau)
    call solve_system(ubrtr, vbrtr, ssh)

    !$omp parallel do private(m,n)
    do n = ny_start-1, ny_end+1
          do m = nx_start-1, nx_end+1
              if(lcu(m,n)>0) then
                  ubrtr(m,n) = ubrtr(m,n) / hhu(m,n)
              endif

              if(lcv(m,n)>0) then
                  vbrtr(m,n) = vbrtr(m,n) / hhv(m,n)
              endif
          enddo
    enddo
    !$omp end parallel do

    !Updating depth functions
    if (full_free_surface>0) then
        call hh_update(hhq, hhu, hhv, hhh, ssh, hhq_rest)
    endif

    ! Check errors
    do n=ny_start,ny_end
      do m=nx_start,nx_end
          if(lu(m,n)>0.5) then
              if(ssh(m,n)<10000.0d0.and.ssh(m,n)>-10000.0d0) then
                  continue
              else
                  write(*,*) rank, 'ERROR!!! In the point m=', m, 'n=', n, 'ssh=', ssh(m,n),   &
                    'step: ', num_step, 'lon: ', geo_lon_t(m, n), 'lat: ', geo_lat_t(m, n)

                  call finalize_parallel(ierr)
                  call exit(0)
              endif
          endif
      enddo
    enddo

endsubroutine shallow_water_model_step
