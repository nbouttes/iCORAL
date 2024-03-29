! get necessary variables
! call coral model
! Beware: levels must be negative under sea level in the input
! ./model_coral > sortie.txt

      PROGRAM MAIN

      use netcdf
      use ncio

      use coral_mod

      implicit none

      integer compteur_day
      integer compteur_day2
      integer t
      integer n,j,l
      integer js
      integer kendy

!input output corals
      REAL area_coral
      REAL temp
      REAL sal
      REAL phos
      REAL light_surf
      REAL kd
      REAL omega_arag
      REAL topof
      REAL depth 
      REAL :: mass_carb_new
      !REAL :: net_carb
      !REAL :: mass_carb
      REAL d_sl
      REAL tau_bleach_l
      REAL timebleach_l
      REAL temp_too_low

      character*256 outfilename

! Start of routine
      compteur_day=0 ! 1 to 365/360
      compteur_day2=0 ! 1 to 30
      NYR=0 ! current year

! Initialisation
      write(*,*) 'Initialisation'
      call ini_coral

      !file with global values
      open (10,file='Coral_output.txt', status='unknown')
      write (10,'(A6,3A10)') 'NYR  ', ' area_coral ',                    &
             ' prod_coral ', ' mass_coral '


      print*, ' '
      print*, 'start of time loop'
! Temporal loop
      do t=1,tmax,dt ! in days
        print*, ' '
        print*, 't=',t
        compteur_day=compteur_day+1 !day number during the year
        compteur_day2=compteur_day2+1 !day number during the month
        kendy=0 !last day of year

        if (compteur_day2.eq.30) then 
          kmon=1
          compteur_day2=0
        endif

        if (compteur_day.eq.ndays_yr) then 
          kendy=1
          NYR=NYR+1
          print*, 'NYR: ', NYR
          compteur_day=0
        endif

       ! Spatial loop
       do n=1,NOC_CBR
        do j=1,JT
         do l=1,LT
          !print*, 'n,j,l ', n,j,l

          ! Temperature from WOCE
          !temp=25 !test
          temp=temp_woce(l,j,n)

          temp_too_low=temp_too_low_all(l,n)

          !Phophate
          !phos=0.1 !test
          phos=phosphate_woce(l,j,n)

          !Salinity
          !sal=35 !test
          sal=salt_woce(l,j,n)

          ! light at surface
          !light_surf=40*1e6/(24*60*60) !test
          light_surf=par(l,n)*1e6/(24*60*60) !light surf from E/m2/j to muE/m2/s

          !Attenuation coefficient
          !kd=0.041 !test
          kd=Kd_490(l,n)
          !if ((kd.gt.0.9) .or. (kd.lt.0)) print*, 'Kd', kd, l,n
          !if (kd.lt.0) print*, 'Kd', kd, l,n
          if (kd .le.0) then
            !write(*,*) 'kd', i,n, kd
            kd=0.1
            kd_490(l,n)=0.1
          endif
          !if (kd.lt.0) print*, 'Kd', kd, l,n
   

          ! aragonite saturation state      
          !omega_arag=5
          omega_arag=aragonite(l,j,n)
          !print*, 'arag_sat', omega_arag

          !print*, 'Bleaching'
          !for bleaching: computes tau_bleach and timebleach later used in
          !coral to limit production
          !only after the first 30 years when the reference is computed for
          !MMMclim
          !now from the beginning using value from restart for MMMclim-> or
          !not
          ! if (NYR.gt.30) then
          if (NYR.gt.window_MMM) then
            call calc_bleach(l,j,n)
          endif

          !loop on subgrid depth
          do js=1,kmax_hypso

           !In iLOVELCIM : all depth in negative except ZX hence the (-1)*ZX
           !here all depth from the bottom to the surface and all in negative
           !write(*,*) 'test depth ', level_bounds_hypso(js),   &
           ! ZX(j), ZX(j+1)
           if ((level_bounds_hypso(js+1).le.ZX(j+1)) .and.              &
            (level_bounds_hypso(js+1).gt.ZX(j))) then
            !write(*,*) 'test depth ', level_bounds_hypso(js+1),       &
            !          ZX(j), ZX(j+1)
            depth=level_hypso(js)
            !write(*,*) 'depth', depth ! depth is negative

            !area: surface of continental bottom in the grid cell in m2
            !area_coral=SQRO2(i,n)
            !print*, 'area coral (1e3 km2)', i,n, area_coral*1e-9
             area_coral=hypso(l,js,n)
            !area_coral=hypso(l,js,n)*0.1 !10percent of the area for corals accounting for inter-reefal area
            !print*, 'area_coral in 1e3 km2',l,js,n, area_coral*1e-9

            ! corals can exist only if some surface area
            if (area_coral .gt. 1e-12) then
             !print*, 'area_coral in 1e3 km2',l,js,n, area_coral*1e-9

             !topographic factor
             !topof=topoff(l,js,n)
             topof=1 !test
             !if ((topof .le.1) .or. (topof .ge.0)) then
             !!  !print*, 'ok', i,js,n
             !else
             !  print*, 'topof',i,js,n, topof
             !endif

             !call subroutine computing CaCO3 production by corals
             !print*, 'call corals'
             call corals(area_coral,temp,sal,phos,light_surf,omega_arag,depth,        &
              kd,topof,tau_bleach(l,j,n), timebleach(l,j,n), temp_too_low,                 &
              t0_diss_all(l,js,n), prod_before_all(l,js,n))

             !output of corals is net_carb in Pmol/day
             coral_prod(l,j,n)=coral_prod(l,j,n)+net_carb !annual net production in Pmol/year

             !New mass in g/day
             mass_carb_new=net_carb*caco3_molar_mass*1.0e15
             !if (mass_carb_new.ne.0) then
             !   write(*,*) 'mass_carb_new', i,j,n,mass_carb_new
             !endif

             total_prod_coral_an=total_prod_coral_an+net_carb !Production in Pmol/day -> Pmol/year
             !if ( depth .ge. -40) then
             !  total_prod_coral_an_40m=total_prod_coral_an_40m+net_carb !Production in Pmol/day -> Pmol/year
             !endif
             total_mass_coral_an=total_mass_coral_an+mass_carb_new*1e-15!mass_carb en g *1e-15 g->Pg -> Pg/year
             !total_area_coral_an=total_area_coral_an+area_coral
             !!(area_coral*topof*1e-6) !in m2->km2 *1e-6
             coral_cum_mass(l,j,n)=coral_cum_mass(l,j,n)+mass_carb_new !g/year
             !coral_mass_subgrid(i,js,n)=coral_mass_subgrid(i,js,n)             &
             !                           +mass_carb_new
             coral_mass(l,j,n)=coral_mass(l,j,n)+mass_carb_new

             !if ((coral_mass_subgrid(i,js,n).gt.0).and.(KENDY.eq.1)) then ! if last day of year and area with corals
             if ((coral_mass(l,j,n).gt.0).and.(KENDY.eq.1)) then ! if last day of year and area with corals
              total_area_coral_an=total_area_coral_an+(area_coral*topof*1e-6) !in m2 *1e-6 -> in km2
              if (depth .ge. -40) then
               total_area_coral_an_40m=total_area_coral_an_40m              &
                       +(area_coral*topof*1e-6) !in m2 *1e-6 -> in km2
              endif
              coral_area(l,j,n)=coral_area(l,j,n)+area_coral*topof !in m2
              ! write(*,*) 'coral area mbiota in 1e3 km2 ',
              ! coral_area(i,j,n)*1e-9
             endif !coral_mass

            endif !area coral
           endif !levels hypso

          enddo !js

         enddo !l
        enddo !j
       enddo !n


       if (kendy.eq.1) then !end of year

        print*, ' '
        print*, 'CORAL at last day of the year'
        print*, 'total_area_coral_an in 1e3 km2 = ', total_area_coral_an*1e-3
        print*, 'total prod_coral_an in Pmol/year = ' , total_prod_coral_an
        print*, 'total mass_coral_an in Pg/year = ', total_mass_coral_an
        print*, 'coral_CO2=', coral_CO2, 'gC/an'

        temp_too_low_all(:,:)=0.0


        ! write annual global mean in file coral_output.txt
        call out_coral_global

        !Remise a 0 a la fin de l annee
        ! and save 3d output
        total_area_coral_an=0.0
        total_prod_coral_an=0.0
        total_mass_coral_an=0.0
        coral_CO2=0
        do n=1,NOC_CBR
         do j=1,JT
          do l=1,LT
           coral_area_out(n,l,j,NYR)=coral_area(l,j,n)
           coral_prod_out(n,l,j,NYR)=coral_prod(l,j,n)
           tau_bleach_out(n,l,j,NYR)=tau_bleach(l,j,n)
           DHW_out(n,l,j,NYR)=DHW(l,j,n)
          enddo !l
         enddo !j
        enddo !n
        coral_area(:,:,:)=0.0
        coral_prod(:,:,:)=0.0
        DHW_nb(:,:,:)=0.0

       endif !kendy
      enddo !t

!to put in outputs :
!coral_area_out
!coral_prod_out
!tau_bleach_out
!DHW_out

! write output
      !daily outputs
!      if (daily_outputs.eq.1) then
!        outfilename = "results_coral.nc"
!        call nc_create(outfilename,overwrite=.TRUE.,netcdf4=.TRUE.)
!        call nc_write_attr(outfilename,"title","Coral model outputs")
!        call nc_write_attr(outfilename,"institution", &
!                       "Bordeaux")
!        call nc_write_dim(outfilename,"lat",x=-90.d0,dx=0.25d0,nx=imax,units="degrees")
!        call nc_write_dim(outfilename,"lon",x=0.d0,dx=0.25d0,nx=jmax,units="degrees")
!        call nc_write_dim(outfilename,"depth",x=level,units="meters")
!        call nc_write_dim(outfilename,"time",x=1.d0,dx=1.d0,nx=tmax, &
!                      units="days",calendar="365_day", unlimited=.TRUE.)
!        !call nc_write(outfilename,"mass_coral",mass_coral_out(:,:,:,:), dim1="depth",dim2="lon",dim3="lat",dim4="time")
!      endif

      !annual outputs
      outfilename = "results_coral_an.nc"
      print*,'write netcdf files'
      call nc_create(outfilename,overwrite=.TRUE.,netcdf4=.TRUE.)
      call nc_write_attr(outfilename,"title","Coral model outputs")
      call nc_write_attr(outfilename,"institution", &
                       "LSCE")
      !call nc_write_dim(outfilename,"lat",x=-90.d0,dx=0.25d0,nx=imax,units="degrees")
      call nc_write_dim(outfilename,"lon",x=-180.d0,dx=1.0d0,nx=NOC_CBR,units="degrees")
      call nc_write_dim(outfilename,"lat",x=-90.d0,dx=1.0d0,nx=LT,units="degrees")
      call nc_write_dim(outfilename,"depth",x=levels,units="meters")
      call nc_write_dim(outfilename,"time",x=1.d0,dx=1.d0,nx=tmax/ndays_yr, &
                      units="years",calendar="365_day", unlimited=.TRUE.)
      !call nc_write_dim(outfilename,"time",x=1.d0,dx=1.d0,nx=1)
      call nc_write(outfilename,"coral_area",coral_area_out(:,:,:,:), dim1="lon",dim2="lat",dim3="depth",dim4="time") 
      call nc_write(outfilename,"coral_prod",coral_prod_out(:,:,:,:), dim1="lon",dim2="lat",dim3="depth",dim4="time")
      call nc_write(outfilename,"tau_bleach",tau_bleach_out(:,:,:,:), dim1="lon",dim2="lat",dim3="depth",dim4="time")
      call nc_write(outfilename,"DHW",DHW_out(:,:,:,:), dim1="lon",dim2="lat",dim3="depth",dim4="time")


      !close Coral_output.txt
      close(10)

      END

