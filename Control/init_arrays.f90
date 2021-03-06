!array boundary definition for non-mpi arrays
subroutine non_mpi_array_boundary_definition
 use main_basin_pars
 use mpi_parallel_tools
 implicit none

       nx_start=mmm
       nx_end  =mm
       ny_start =nnn
       ny_end  =nn

       bnd_x1=nx_start-2
       bnd_x2=nx_end  +2
       bnd_y1=ny_start-2
       bnd_y2=ny_end  +2

endsubroutine non_mpi_array_boundary_definition
!end of array boundary definition for non-mpi arrays
!-------------------------------------------------------------------------------

!array boundary definition for mpi arrays
subroutine mpi_array_boundary_definition
    use main_basin_pars
    use mpi_parallel_tools
    implicit none

    integer :: ierr, procn, i
    integer :: locn
    integer :: count_threads, num_thread

    period = (/1,1/)
    p_size = (/0,0/)
    ierr = 0

    ! mpi_comm_world
    call mpi_comm_rank(PETSc_Comm_World, rank, ierr)
    call mpi_comm_size(PETSc_Comm_World, procs, ierr)
    call mpi_dims_create(procs, 2, p_size, ierr)
    call mpi_cart_create(PETSc_Comm_World, 2, p_size, period, 0, cart_comm, ierr)
    call mpi_cart_coords(cart_comm, rank, 2, p_coord, ierr)

!-----------------------------------NX------------------------------------------!
    locn = floor(real(nx - 4)/real(p_size(1)))
    nx_start = locn*p_coord(1) + 1 + 2
    if ( p_coord(1) .EQ. p_size(1) - 1 ) then
        locn = (nx - 2) - nx_start + 1
    endif
    nx_end = nx_start + locn - 1
    nx_start = nx_start
!   border area
    bnd_x1 = nx_start - 2
!    if (bnd_x1 < 1) bnd_x1 = 1
    bnd_x2 = nx_end + 2
!    if (bnd_x2 > nx) bnd_x2 = nx

!-----------------------------------NY------------------------------------------!
    locn = floor(real(ny - 4)/real(p_size(2)))
    ny_start = locn*p_coord(2) + 1 + 2
    if ( p_coord(2) .EQ. p_size(2) - 1 ) then
        locn = (ny - 2) - ny_start + 1
    endif
    ny_end = ny_start + locn - 1
    ny_start = ny_start
!   border area
    bnd_y1 = ny_start - 2
!    if (bnd_y1 < 1) bnd_y1 = 1
    bnd_y2 = ny_end + 2
!    if (bnd_y2 > ny) bnd_y2 = ny

    call mpi_comm_size(cart_comm, procn, ierr)
    if (rank .eq. 0) print *, "MPI pocs: ", procn, " Domain decomposition:"
    do i = 0, procn-1
        if (rank .eq. i) then
            print *, "nx ", rank, p_coord, nx_start, nx_end, ny_start, ny_end
            print *, "bnd", rank, p_coord, bnd_x1, bnd_x2, bnd_y1, bnd_y2
        endif
        call mpi_barrier(cart_comm, ierr)
    enddo

    !$omp parallel
    count_threads = omp_get_num_threads()
    num_thread = omp_get_thread_num()
    if (num_thread .eq. 0) print *, "OMP Threads: ", count_threads
    !$omp end parallel

    call mpi_barrier(cart_comm, ierr)

endsubroutine mpi_array_boundary_definition
!-------------------------------------------------------------------------------

!allocation of arrays
subroutine model_grid_allocate
    use main_basin_pars
    use mpi_parallel_tools
    use basin_grid
    implicit none

    allocate(lu(bnd_x1:bnd_x2,bnd_y1:bnd_y2),             &  !mask of t-grid
            lu1(bnd_x1:bnd_x2,bnd_y1:bnd_y2),             &  !mask of t-grid (1 everywhere)
            luu(bnd_x1:bnd_x2,bnd_y1:bnd_y2),             &  !mask of h-grid (0 on boundary)
            luh(bnd_x1:bnd_x2,bnd_y1:bnd_y2),             &  !mask of h-grid (1 on boundary)
            lcu(bnd_x1:bnd_x2,bnd_y1:bnd_y2),             &  !mask of u-grid (0 on boundary)
            lcv(bnd_x1:bnd_x2,bnd_y1:bnd_y2),             &  !mask of v-grid (0 on boundary)
            llu(bnd_x1:bnd_x2,bnd_y1:bnd_y2),             &  !mask of u-grid (0 on boundary)
           llv(bnd_x1:bnd_x2,bnd_y1:bnd_y2) )                !mask of v-grid (0 on boundary)
    allocate   (lbasins(nx,ny))       !integer masks of regional basins
    allocate   (  hhh(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &  !ocean depth on luh (h-points)
                 hhhp(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &  !ocean depth on luh (h-points) at previous step
             hhq_rest(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &  !ocean depth on lu  (t-points) at rest state
             hhu_rest(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &  !ocean depth on lu  (t-points) at rest state
             hhv_rest(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &  !ocean depth on lu  (t-points) at rest state
                  hhq(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &  !ocean depth on lu  (t-points)
                 hhqp(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &  !ocean depth on lu  (t-points) at previous step
                  hhu(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &  !ocean depth on lcu (u-points)
                 hhup(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &  !ocean depth on lcu (u-points) at previous step
                  hhv(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &  !ocean depth on lcv (v-points)
                 hhvp(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &  !ocean depth on lcv (v-points) at previous step
                  rlh_s(bnd_x1:bnd_x2,bnd_y1:bnd_y2),    &  !coriolis-1 parameter on edge (t-centers) points
                  rlh_c(bnd_x1:bnd_x2,bnd_y1:bnd_y2),    &  !coriolis-2 parameter on edge (t-centers) points
                                     z(nz), zw(nz+1),    &  !vertical sigma-levels (t-points and w-points)
                                 hzt(nz+1), dz(nz),      &  !steps between t-levels and w-levels
       dxt(bnd_x1:bnd_x2,bnd_y1:bnd_y2),dyt(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &  !horizontal grid steps between   t-points (in radians or meters)
       dx (bnd_x1:bnd_x2,bnd_y1:bnd_y2),dy (bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &  !horizontal grid steps between u,v-points (in radians or meters)
       dxh(bnd_x1:bnd_x2,bnd_y1:bnd_y2),dyh(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &  !horizontal grid steps between   h-points (in radians or meters)
       dxb(bnd_x1:bnd_x2,bnd_y1:bnd_y2),dyb(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &  !horizontal grid steps between v,u-points (in radians or meters)
             xt(bnd_x1:bnd_x2),yt(bnd_y1:bnd_y2),        &  !horizontal t-grid            x- and y-coordinates (in degrees)
             xu(bnd_x1:bnd_x2),yv(bnd_y1:bnd_y2)    )       !horizontal u-grid and v-grid x- and y-coordinates (in degrees)

    allocate(geo_lon_t(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &    !geographical longitudes of T-points
             geo_lat_t(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &    !geographical latitudes  of T-points
             geo_lon_u(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &    !geographical longitudes of U-points
             geo_lat_u(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &    !geographical latitudes  of U-points
             geo_lon_v(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &    !geographical longitudes of V-points
             geo_lat_v(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &    !geographical latitudes  of V-points
             geo_lon_h(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &    !geographical longitudes of H-points
             geo_lat_h(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &    !geographical latitudes  of H-points
          rotvec_coeff(bnd_x1:bnd_x2,bnd_y1:bnd_y2, 4) )

    lu=0.0; lu1=0.0; luu=0.0; luh=0.0; lcu=0.0; lcv=0.0; llu=0.0; llv=0.0; lbasins=0
    hhh=0.0d0; hhhp=0.0d0
    hhq_rest=0.0d0; hhu_rest=0.0d0; hhv_rest=0.0d0

    hhq=0.0d0; hhqp=0.0d0
    hhu=0.0d0; hhup=0.0d0
    hhv=0.0d0; hhvp=0.0d0
    rlh_s=0.0d0; rlh_c=0.0d0
    z=0.0d0; zw=0.0d0; hzt=0.0d0; dz=0.0d0
    dxt=0.0d0; dyt=0.0d0; dx=0.0d0; dy=0.0d0; dxh=0.0d0; dyh=0.0d0; dxb=0.0d0; dyb=0.0d0
    xt=0.0d0; yt=0.0d0; xu=0.0d0; yv=0.0d0

    geo_lon_t=0.0d0; geo_lat_t=0.0d0; geo_lon_u=0.0d0; geo_lat_u=0.0d0
    geo_lon_v=0.0d0; geo_lat_v=0.0d0; geo_lon_h=0.0d0; geo_lat_h=0.0d0
    rotvec_coeff=0.0d0

endsubroutine model_grid_allocate
!-------------------------------------------------------------------------------

!deallocation of arrays
subroutine model_grid_deallocate
    use basin_grid
    implicit none

    deallocate(rotvec_coeff)
    deallocate(geo_lat_h,geo_lon_h,geo_lat_v,geo_lon_v,geo_lat_u,geo_lon_u,geo_lat_t,geo_lon_t)
    deallocate(yv,xu,yt,xt,dyb,dxb,dyh,dxh,dy,dx,dyt,dxt,dz,hzt,zw,z,rlh_c,rlh_s)
    deallocate(hhvp,hhv,hhup,hhu,hhqp,hhq,hhhp,hhh)
    deallocate(hhq_rest, hhu_rest, hhv_rest)

    deallocate(lbasins)
    deallocate(llv,llu,lcv,lcu,luh,luu,lu1,lu)

endsubroutine model_grid_deallocate
!-------------------------------------------------------------------------------

!allocation of arrays
subroutine ocean_variables_allocate
    use main_basin_pars
    use mpi_parallel_tools
    use ocean_variables
    implicit none

    allocate(ssh(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &  !sea surface height (SSH) at current  time step [m] (internal mode)
             pgrx(bnd_x1:bnd_x2,bnd_y1:bnd_y2),    &  !pressure gradient x-component for RHS
             pgry(bnd_x1:bnd_x2,bnd_y1:bnd_y2),    &  !pressure gradient y-component for RHS
             ubrtr(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &  !barotropic velocity      zonal[m/s] at current  time step [m] (internal mode)
             vbrtr(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &  !barotropic velocity meridional[m/s] at current  time step [m] (internal mode)
             RHSx2d(bnd_x1:bnd_x2,bnd_y1:bnd_y2),  &  !x-component of external force(barotropic)
             RHSy2d(bnd_x1:bnd_x2,bnd_y1:bnd_y2))     !y-component of external force(barotropic)

    allocate (sshp(bnd_x1:bnd_x2,bnd_y1:bnd_y2),        &
              ubrtrp(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &
              vbrtrp(bnd_x1:bnd_x2,bnd_y1:bnd_y2))

    allocate (xxt(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),    &  !auxiliary array 1
              yyt(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz))       !auxiliary array 2

    allocate (tflux_surf(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &       !total surface heat flux [�C*m/s]
              tflux_bot(bnd_x1:bnd_x2,bnd_y1:bnd_y2),       &       !total bottom heat flux [�C*m/s]
              sflux_surf(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &       !total surface salt flux [   m/s]
              sflux_bot(bnd_x1:bnd_x2,bnd_y1:bnd_y2),       &       !total bottom salt flux [   m/s]
          surf_stress_x(bnd_x1:bnd_x2,bnd_y1:bnd_y2),       &       !wind      zonal stress per water density [m^2/s^2]
          surf_stress_y(bnd_x1:bnd_x2,bnd_y1:bnd_y2),       &       !wind meridional stress per water density [m^2/s^2]
           bot_stress_x(bnd_x1:bnd_x2,bnd_y1:bnd_y2),       &       !bottom      zonal stress per water density [m^2/s^2]
           bot_stress_y(bnd_x1:bnd_x2,bnd_y1:bnd_y2),       &       !bottom meridional stress per water density [m^2/s^2]
               divswrad(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),    &       !shortwave radiation divergence coefficients
                   dkft(bnd_x1:bnd_x2,bnd_y1:bnd_y2),       &       !relaxation coefficient for SST, [m/s]
                   dkfs(bnd_x1:bnd_x2,bnd_y1:bnd_y2),       &       !relaxation coefficient for SSS, [m/s]
               sensheat(bnd_x1:bnd_x2,bnd_y1:bnd_y2),       &       !sensible heat flux
                latheat(bnd_x1:bnd_x2,bnd_y1:bnd_y2),       &       !latent heat flux
                 lw_bal(bnd_x1:bnd_x2,bnd_y1:bnd_y2),       &       !longwave radiation balance
                 sw_bal(bnd_x1:bnd_x2,bnd_y1:bnd_y2),       &       !shortwave radiation balance
                 hf_tot(bnd_x1:bnd_x2,bnd_y1:bnd_y2),       &       !total heat flux
                 wf_tot(bnd_x1:bnd_x2,bnd_y1:bnd_y2) )              !total water flux

    allocate (tatm(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &    !Air temperature, [�C]
              qatm(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &    !Air humidity, [kg/kg]
              rain(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &    !rain, [m/s]
              snow(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &    !snow, [m/s]
              wind(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &    !Wind speed module, [m/s]
               lwr(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &    !Downward  longwave radiation, [W/m^2]
               swr(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &    !Downward shortwave radiation, [W/m^2]
              slpr(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &    !Sea level pressure, [Pa]
              uwnd(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &    !Zonal      wind speed, [m/s]
              vwnd(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &    !Meridional wind speed, [m/s]
              taux(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &    !Zonal      wind speed, [m/s]
              tauy(bnd_x1:bnd_x2,bnd_y1:bnd_y2) )       !Meridional wind speed, [m/s]

    allocate(BottomFriction(bnd_x1:bnd_x2,bnd_y1:bnd_y2),        &
                   r_diss(bnd_x1:bnd_x2,bnd_y1:bnd_y2))

    allocate(amuv2d(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &    !depth mean lateral viscosity
             amuv42d(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &    !depth mean lateral viscosity
            r_vort2d(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &    !relative vorticity of depth mean velocity
          stress_t2d(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &    !Horizontal tension tensor component (barotropic)
          stress_s2d(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &    !Horizontal shearing tensor component(barotropic)
    RHSx2d_tran_disp(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &    !dispersion x-component of external force(barotropic)
    RHSy2d_tran_disp(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &    !dispersion x-component of external force(barotropic)
    RHSx2d_diff_disp(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &    !dispersion x-component of external force(barotropic)
    RHSy2d_diff_disp(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &    !dispersion x-component of external force(barotropic)
    RHSx2d_bfc(bnd_x1:bnd_x2,bnd_y1:bnd_y2),           &
    RHSy2d_bfc(bnd_x1:bnd_x2,bnd_y1:bnd_y2))


    ssh = 0.0d0;
    pgrx = 0.0d0; pgry = 0.0d0
    ubrtr = 0.0d0;  vbrtr = 0.0d0
    RHSx2d = 0.0d0; RHSy2d = 0.0d0

    sshp = 0.0d0
    ubrtrp = 0.0d0
    vbrtrp = 0.0d0

    xxt=0.0d0; yyt=0.0d0

    tflux_surf=0.0d0; tflux_bot=0.0d0
    sflux_surf=0.0d0; sflux_bot=0.0d0
    surf_stress_x=0.0d0; surf_stress_y=0.0d0
    bot_stress_x=0.0d0;  bot_stress_y=0.0d0
    divswrad=0.0d0; dkft=0.0d0; dkfs=0.0d0
    sensheat=0.0d0; latheat=0.0d0; lw_bal=0.0d0; sw_bal=0.0d0
    hf_tot=0.0d0; wf_tot=0.0d0

    tatm=0.0d0; qatm=0.0d0; rain=0.0d0; snow=0.0d0; wind=0.0d0
    lwr=0.0d0; swr=0.0d0; slpr=0.0d0; uwnd=0.0d0; vwnd=0.0d0
    taux=0.0d0; tauy=0.0d0

    BottomFriction=0.0d0; r_diss=0.0d0

    BottomFriction=0.0d0; r_diss=0.0d0

    amuv2d=0.0d0; amuv42d=0.0d0; r_vort2d=0.0d0
    stress_t2d=0.0d0; stress_s2d=0.0d0
    RHSx2d_tran_disp=0.0d0; RHSy2d_tran_disp=0.0d0
    RHSx2d_diff_disp=0.0d0; RHSy2d_diff_disp=0.0d0

    RHSx2d_bfc = 0.0d0; RHSy2d_bfc = 0.0d0

endsubroutine ocean_variables_allocate
!-------------------------------------------------------------------------------

!deallocation of arrays
subroutine ocean_variables_deallocate
    use ocean_variables
    implicit none

    deallocate(RHSx2d_bfc, RHSy2d_bfc)
    deallocate(RHSy2d_diff_disp,RHSx2d_diff_disp,RHSy2d_tran_disp,RHSx2d_tran_disp,     &
              stress_s2d,stress_t2d,r_vort2d,amuv42d,amuv2d)
    deallocate(r_diss, BottomFriction)
    deallocate(tauy,taux,vwnd,uwnd,slpr,swr,lwr,wind,snow,rain,qatm,tatm)
    deallocate(wf_tot,hf_tot,sw_bal,lw_bal,latheat,sensheat,dkfs,dkft,           &
               divswrad,bot_stress_y,bot_stress_x,surf_stress_y,surf_stress_x,   &
               sflux_bot,sflux_surf,tflux_bot,tflux_surf)
    deallocate(xxt, yyt)
    deallocate(vbrtrp,ubrtrp,sshp)
    deallocate(RHSy2d,RHSx2d,vbrtr,ubrtr,pgry,pgrx,ssh)

endsubroutine ocean_variables_deallocate
!-------------------------------------------------------------------------------

!================================= ATM DATA ====================================
subroutine atm_arrays_allocate
    use atm_pars
    use atm_forcing
    implicit none

    allocate(xa(nxa),ya(nya))

    allocate( a_hflux(nxa,nya),       &   !heat balance [w/m**2]
           a_swrad(nxa,nya),       &   !sw radiation balance[w/m**2]
           a_wflux(nxa,nya),       &   !precipitation-evaporation[m/s]
           a_stress_x(nxa,nya),    &   !zonal wind stress[pA=n/m**2]
           a_stress_y(nxa,nya),    &   !meridional wind stress[pA=n/m**2]
           a_slpr(nxa,nya),        &   !pressure at sea surface
           a_lwr(nxa,nya),         &   !dw-lw-rad[w/m**2]
           a_swr(nxa,nya),         &   !dw-sw-rad[w/m**2]
           a_rain(nxa,nya),        &   !precipit[m/s]
           a_snow(nxa,nya),        &   !precipit[m/s]
           a_tatm(nxa,nya),        &   !temp of atmosphere[�c]
           a_qatm(nxa,nya),        &   !humidity [g/kg]
           a_uwnd(nxa,nya),        &   !u-wind speed[m/s]
           a_vwnd(nxa,nya)  )          !v-wind speed[m/s]

    xa=0.0d0; ya=0.0d0

    a_hflux=0.0; a_swrad=0.0; a_wflux=0.0; a_stress_x=0.0; a_stress_y=0.0
    a_slpr=0.0;  a_lwr=0.0;   a_swr=0.0; a_rain=0.0; a_snow=0.0
    a_tatm=0.0;  a_qatm=0.0; a_uwnd=0.0; a_vwnd=0.0

    ind_change_heat =0
    ind_change_water=0
    ind_change_stress=0
    ind_change_rad  =0
    ind_change_prec =0
    ind_change_tatm =0
    ind_change_qatm =0
    ind_change_wind =0
    ind_change_slpr =0

    num_rec_heat =0
    num_rec_water=0
    num_rec_stress=0
    num_rec_rad  =0
    num_rec_prec =0
    num_rec_tatm =0
    num_rec_qatm =0
    num_rec_wind =0
    num_rec_slpr =0

endsubroutine atm_arrays_allocate

!-------------------------------------------------------------------------------
subroutine atm_arrays_deallocate
    use atm_pars
    use atm_forcing
    implicit none

    deallocate(a_vwnd,a_uwnd,a_qatm,a_tatm,a_snow,a_rain,a_swr,a_lwr,     &
             a_slpr,a_stress_y,a_stress_x,a_wflux,a_swrad,a_hflux)
    deallocate(ya,xa)

endsubroutine atm_arrays_deallocate

!-------------------------------------------------------------------------------
subroutine atm2oc_allocate
    use mpi_parallel_tools
    use atm2oc_interpol
    implicit none

    allocate(   wght_mtrx_a2o(bnd_x1:bnd_x2,bnd_y1:bnd_y2,4),      &
                i_input_a2o(bnd_x1:bnd_x2,bnd_y1:bnd_y2,4),      &
                j_input_a2o(bnd_x1:bnd_x2,bnd_y1:bnd_y2,4)    )
    wght_mtrx_a2o=0.0d0; i_input_a2o=0; j_input_a2o=0

endsubroutine atm2oc_allocate

!-------------------------------------------------------------------------------
subroutine atm2oc_deallocate
    use atm2oc_interpol
    implicit none

    deallocate(j_input_a2o  , i_input_a2o  , wght_mtrx_a2o  )
endsubroutine atm2oc_deallocate
