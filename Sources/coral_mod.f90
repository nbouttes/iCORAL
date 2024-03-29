       module coral_mod
! Based on Kleypas
! Written by N. Bouttes, G. Munhoven, V. Brovkin...
!
!---------------------------------------------------------------------- 
! use section
! declaration section
       implicit none

!------------------------------------------------
! things to modify:
! grid of the model
! beware to modify the resolution in the output netcdf files
      INTEGER LT, JT, NOC_CBR
!tbd      !PARAMETER (imax=60,jmax=120,kmax=10) !loveclim grid (3degres resolution)
!tbd      !PARAMETER (imax=65,jmax=120,kmax=20) !loveclim grid clio
!tbd      !PARAMETER (imax=720,jmax=1440,kmax=4) ! as Jones et al 2015 (0.25 degres resolution) ! kmax=20

      !PARAMETER (LT=720, JT=4, NOC_CBR=1440) ! as Jones et al 2015 (0.25 degres resolution) ! kmax=20
      PARAMETER (LT=180, JT=4, NOC_CBR=360) ! 1 degree resolution
      !PARAMETER (LT=65, JT=20, NOC_CBR=120) ! iLOVECLIM

!depth of our grid
      integer, dimension(JT) :: levels !levels
      integer, dimension(JT+1) :: ZX ! bounds

!tbd---
      !temperature grid (should be the same)
      !INTEGER xaxmax,yaxmax,zaxmax
      !PARAMETER (yaxmax=60,xaxmax=120,zaxmax=10)
      !PARAMETER (yaxmax=720,xaxmax=1440,zaxmax=4) !as Jones et al 2015 (0.25 degres resolution) !zazmax=20 4
      integer, dimension(JT) :: zax ! levels dans WOCE file, should be the same as levels
!tbd end----
! gebco grid for hypsometry different for levels
!      INTEGER kmax_hypso
!      PARAMETER (kmax_hypso=128)
      integer, parameter :: kmax_hypso=150
!time
      INTEGER, parameter :: n_years=2 ! number of years of the simulation
      integer, parameter :: ndays_yr=1 ! number of days in one year 360
      INTEGER tmax
      !PARAMETER (tmax=365*n_years) !max time in days 365 * nb years
      PARAMETER (tmax=ndays_yr*n_years) !max time in days 365 * nb years
      INTEGER, parameter :: dt=1 !time step =1 day
!time sea level
      INTEGER t_sl
      PARAMETER (t_sl=15) !1000 nb of years
!outputs
      INTEGER daily_outputs
      PARAMETER (daily_outputs=0) !1=daily outputs, 0=no

!-------------------------------------------------

!Variables
      integer NYR, kmon
      REAL tau !constante de temps de wash out des coraux
      REAL caco3_molar_mass

      REAL PAR_m
      REAL dens_ocn
      REAL water_z_pp

      REAL net_carb

      !weathering
      REAL C_sed
      REAL C_riv
      REAL C_car_a
      REAL A_riv

      !for bleaching
      REAL tau_bleach_moderat
      REAL tau_bleach_strong
      REAL dhw_thresh_moderat
      REAL dhw_thresh_strong
      !real, allocatable,dimension(:,:) :: bar
      !allocate(bar(nf,nx*ny*nz))

      REAL, dimension(LT,NOC_CBR) :: temp_too_low_all

      !dissolution
      REAL :: lambda_diss 

      real total_area_coral_an
      real total_area_coral_an_40m
      real total_prod_coral_an
      real total_prod_coral_an_40m
      real total_mass_coral_an
      real coral_CO2
      REAL, dimension(LT,JT,NOC_CBR) :: coral_prod
      REAL, dimension(LT,JT,NOC_CBR) :: coral_area
      REAL, dimension(LT,JT,NOC_CBR) :: coral_cum_mass
      !REAL, dimension(LT,kmax_hypso,NOC_CBR) :: coral_mass_subgrid
      REAL, dimension(LT,JT,NOC_CBR) :: coral_mass
      REAL, dimension(LT,JT,NOC_CBR) :: omega_arag3D
      REAL, dimension(LT,JT,NOC_CBR) :: PH
      REAL, dimension(LT,JT,NOC_CBR) :: TM_fix

      !outputs for netcdf files
      REAL, dimension(NOC_CBR,LT,JT,n_years) :: coral_prod_out
      REAL, dimension(NOC_CBR,LT,JT,n_years) :: coral_area_out
      REAL, dimension(NOC_CBR,LT,JT,n_years) :: tau_bleach_out
      REAL, dimension(NOC_CBR,LT,JT,n_years) :: DHW_out


      ! for bleaching
      REAL, dimension(LT,JT,NOC_CBR) :: tau_bleach
      REAL, dimension(LT,JT,NOC_CBR) :: timebleach
      REAL, dimension(LT,JT,NOC_CBR) :: DHW ! degree heating week
      REAL, dimension(LT,JT,NOC_CBR) :: DHW_nb ! nb of weeks when bleaching is activated (moderate)
      REAL, dimension(LT,JT,NOC_CBR) :: MMMclim ! monthly max climatological temperature
!tbd      REAL, dimension(LT,JT,NOC_CBR) :: xsHS 
      REAL, dimension(LT,JT,NOC_CBR) :: temp_mois
      INTEGER, PARAMETER :: window_MMM=30 !30 years: window for running mean of monthly temperature
      INTEGER, PARAMETER :: nb_mois=window_MMM*12
      REAL, dimension(LT,JT,NOC_CBR,nb_mois) :: temp_mois_all
      INTEGER indice_mois
!tbd      INTEGER, PARAMETER :: window_week=17 !number of weeks for DHW
      INTEGER, PARAMETER :: nb_hs=84 !window_week*joursemaine !17 weeks of 5 days = 85 days
      REAL, dimension(LT,JT,NOC_CBR,nb_hs) :: xsHS_all
      INTEGER indice_hs
      ! Temperature with variability for bleaching
      REAL, dimension(LT, JT, NOC_CBR) :: TM_pluswkvar
      INTEGER, SAVE :: k_mbiota_rand
      REAL :: nino3, nino3_var

      !t0 for dissolution when no production
      REAL, dimension (LT, kmax_hypso, NOC_CBR) :: t0_diss_all
      REAL, dimension (LT, kmax_hypso, NOC_CBR) :: prod_before_all
      REAL, dimension(kmax_hypso) :: level_hypso
      REAL, dimension(kmax_hypso+1) :: level_bounds_hypso

      !For input files
      real, dimension (LT,NOC_CBR) :: par ! photosynthetically available radiation
      real, dimension (LT,NOC_CBR) :: Kd_490 ! Attenuation coeffecient
      real, dimension (LT,JT,NOC_CBR) :: aragonite ! saturation state
      real, dimension (LT,JT,NOC_CBR) :: temp_woce 
      real, dimension (LT,JT,NOC_CBR) :: salt_woce
      real, dimension (LT,JT,NOC_CBR) :: phosphate_woce 
      real, dimension (LT,kmax_hypso,NOC_CBR) :: hypso
      real, dimension (LT,NOC_CBR) :: topoff
      !real, dimension (LT) :: latitude_hypso
      !real, dimension (NOC_CBR) :: longitude_hypso

      contains

!----------------------------------------------------------------------
! initialisation
      subroutine ini_coral

      use ncio,      only: nc_read

      ! local
      integer :: n

      !character*256 FILE_NAME 
      !parameter(FILE_NAME = "toto.nc")
      !integer status_nc, ncid

      levels=(/ -35, -25, -15, -5 /) ! depth levels
      ZX=(/ -40,-30,-20,-10,0 /) ! vertical bounds

      !constantes
      tau=4000 !in years
      PAR_m=0.4
      dens_ocn   = 1.03         ! kg/l
      water_z_pp = 25
      !Molar mass of CaCO3 = 100.0869 g/mol
      caco3_molar_mass=100.0869
      !weathering constant in Pmol/day (1e15) (on veut environ 7.5-9 Tmol/year (1e12)).
      !TYER/TDAY=360 day/year
      ! all in form of hco3-
      !C_sed=enfouissement de CaCO3 dans les sediments, inluant les coraux
      !C_sed=7.5*1e-3/(TYER/TDAY) ! C_sed=G_coral a l equilibre
      !C_riv= carbon from carbonate dissolution brought by rivers to the
      !ocean all in HCO3- form
      C_riv=3.9*1e-5  !2*C_sed !Pmol/day
      ! C_sil_a= 0 ! consommation CO2 par alteration silicates
      ! C_vol+C_hyd =0 ! Carbon from volcanism and hydrothermals going
      ! to atmsophere
      ! C_car_a= consommation of CO2 by carbonate alteration
      C_car_a=C_riv /2.
      !Alkalinity
      A_riv=C_riv
      !for bleaching
      !bleaching effect, time scales (yrs)
      tau_bleach_moderat =  5. !20.
      tau_bleach_strong  = 20. !100.
      !DHW index bleaching thresholds, using weeks of 5 days
      !(joursemaine=5)
      dhw_thresh_moderat = 4.!  *7./joursemaine
      dhw_thresh_strong  = 8.!  *7./joursemaine

      !constante for dissolution
      lambda_diss=1/10.

      ! initialisation
      total_area_coral_an=0
      total_area_coral_an_40m=0
      total_prod_coral_an=0
      total_prod_coral_an_40m=0
      total_mass_coral_an=0

      coral_CO2=0

      coral_prod(:,:,:)=0.
      coral_area(:,:,:)=0.
      coral_cum_mass(:,:,:)=0.
      !coral_mass_subgrid(:,:,:)=0.
      coral_mass(:,:,:)=0.
      coral_prod_out(:,:,:,:)=0.
      coral_area_out(:,:,:,:)=0.
      omega_arag3D(:,:,:)=0.
      tau_bleach_out(:,:,:,:)=0.
      DHW_nb(:,:,:)=0. ! nb of weeks with bleaching
      DHW_out(:,:,:,:)=0.
      PH(:,:,:)=0.
      TM_fix(:,:,:)=0.
      !dissolution
      t0_diss_all(:,:,:)=0.0
      prod_before_all(:,:,:)=0.0

      ! for bleaching
     !allocate(tau_bleach(LT,JT,NOC_CBR))
      !if (KLSR.eq.0) then ! only if fresh start, otherwise read in restart
        tau_bleach(:,:,:)=0.0
        timebleach(:,:,:)=0.0
        MMMclim(:,:,:)=0.0 !max of monthly temperature
     !endif
      !indice_mois=0
      indice_mois=nb_mois
      temp_mois(:,:,:)=0.
      DHW(:,:,:)=0.0
      MMMclim(:,:,:)=0.0 !max of monthly temperature
!tbd      xsHS(:,:,:)=0.
      indice_hs=1
!tbd  indice_hs=nb_hs
      xsHS_all(:,:,:,:)=0.
      nino3=0.0
      nino3_var=0.0

      !limite temperature low
      temp_too_low_all(:,:)=0


      !read input files

      !utile uniquement sur ordi Nath pour que ncio soit content, sinon ncio ne marche pas...
      !status_nc =  nf90_create(FILE_NAME, NF90_CLOBBER, ncid)

      !area (surface and hypsometry) from gebco
      print*, 'read hypsometry'
      !!call nc_read_attr("GEBCO/hypsometry.nc", "title", testchar)
      !!write(*,*) "Title: ", trim(testchar)
      call nc_read("../../Input/GEBCO/hypsometry.nc","hypso",hypso) !(:,:,:)
      call nc_read("../../Input/GEBCO/hypsometry.nc","level",level_hypso)
      call nc_read("../../Input/GEBCO/hypsometry.nc","level_bound",level_bounds_hypso)
      !print*, 'latitude', latitude_hypso
      !call nc_read("../../Input/GEBCO/hypsometry.nc","longitude",longitude_hypso)
      !print*, 'longitude', longitude_hypso
      !print*, 'level_hypso', level_hypso
      !print*, 'level_bounds_hypso', level_bounds_hypso
      !d_hypso=level_hypso(2)-level_hypso(1)
      !print*,'resolution vertical hypsometry ', d_hypso

      !topof (topographic factor) computed from gebco
      print*, 'read topof'
      call nc_read("../../Input/GEBCO/topof_regridded.nc","topof",topoff)

      !Temperature from WOCE
      print*, 'read temperature'
      call nc_read("../../Input/TEMP_SAL/TEMP_regridded.nc", "temperature",temp_woce)
      !call nc_read("../../Input/TEMP_SAL/TEMP_regridded.nc", "longitude",xax)
      !call nc_read("../../Input/TEMP_SAL/TEMP_regridded.nc", "latitude",yax)
      call nc_read("../../Input/TEMP_SAL/TEMP_regridded.nc", "levels",zax)
      print*, 'levels WOCE', zax
      !si pas lecture de area (sinon doit etre identique)
      !latitude=yax
      !longitude=xax
!      level=zax ! different de level_hypso levels en m et en valeur negatives (sous la mer), du plus profond vers la surface

      !Salinity from WOCE
      print*, 'read salinity'
      call nc_read("../../Input/TEMP_SAL/SALT_regridded.nc", "salinity",salt_woce)

      !Phosphate from WOCE
      print*, 'read phosphate'
      call nc_read("../../Input/TEMP_SAL/PO4_regridded.nc", "phosphate",phosphate_woce)

      !PAR (insolation)
      print*, 'read insolation (PAR)'
      call nc_read("../../Input/PAR/par_regridded.nc", "par", par)

      !Kd_490 (attenuation of insolation)
      print*, 'read attenuation coefficient (Kd)'
      call nc_read("../../Input/Kd_490/Kd_490_regridded.nc", "Kd_490",Kd_490)

      !Omega aragonite
      print*, 'read aragonite saturation state (Omega)'
      call nc_read("../../Input/Omega/omega_aragonite_regridded.nc", "aragonite",aragonite)

      !TBD add TOPOF

      end subroutine ini_coral

!----------------------------------------------------------------------

!----------------------------------------------------------------------
! coral production
      subroutine corals(area,temp,sal,phos,light_surf,omega_arag,depth, &
            kd,topof,tau_bleach_l, timebleach_l,temp_too_low,           &
            t0_diss, prod_before)

!input output

      REAL area
      REAL temp
!      REAL tomin
!      REAL tomax
      REAL sal
      REAL light_surf
      REAL omega_arag !arag_sat
      REAL depth
      REAL mass_carb !, mass_carb_new
      !REAL d_sl
      REAL kd
      REAL topof
      REAL phos
!      REAL nitr
      REAL tau_bleach_l
      REAL timebleach_l
      REAL temp_too_low
      REAL t0_diss
      REAL prod_before


!local
!     to be moved to init and module
      REAL gmax_coral
      REAL zmax
      REAL Ik !RAD_IK
      REAL Imin
      REAL coral_dens
      REAL pk490
      REAL tmin, tmax
      REAL smin, smax
      REAL nmax, pmax
      REAL coef_ftemp
      REAL arag_k
      REAL dens_ocn
      REAL water_z_pp
      integer i_pp
      REAL tfmin, tfmax
      REAL a_temp, b_temp


!to be deleted      REAL area_prod
      REAL Iz
      !REAL mass_carb_old
      REAL P_carb, D_carb
      REAL g_coral
      REAL RAD_m
      REAL temp_factor


! maximum vertical accumulation rate for corals in m/yr -> per day for us (-> /360) -> m/day
!      gmax_coral=1.04/360 !mm/day
!      gmax_coral=1.0/1000./360 !m/day
!      gmax_coral=(4.0*1E-3)/360 !in m/day (from mm/year)
      gmax_coral=(3.0*1E-3)/360. !in m/day (from mm/year)

! Coefficient for Temperature dependency
      coef_ftemp = 0.2      
! target - caco3 flux of 0.105 Gtc/yr
! saturation light intensity Ik, in W/m2  (conversion between W/m2 and muE m-2 s-1 using 4.6 from Kirk, 1994)
      Ik=350./4.6 !RAD_IK ! from Guypour CLIMBER
      !Ik=50 !200 !50 !?? A voir ?? !pour moi : mE/m2/s
      ! Ik devrait etre en W/m2 ?

! minimum light intensity necessary for reef growth
      Imin=300./4.6 !in W/m2  (conversion between W/m2 and muE m-2 s-1 using 4.6 from Kirk, 1994)

! CaCO3 density, kg/m3 (density of 2.89g/cm3 and porosity 50%)
      coral_dens=1.445*1.0E3 ! in kg/m3
! light extinction coefficient  !!!
! to be deleted      pk490=0.1
! minimum temperature for coral growth (grad. C)
      tmin= 18.1
! maximum temp. for coral growth
      tmax= 31.5
! minimum salinity for coral growth (--)
      smin= 30.0
! maximum salinity for coral growth
      smax= 39.0
!maximum phosphate value for coral growth
      pmax=0.2 !micromol/L
!cc sea level (m)
!cc      isea_lev=250
!c supersaturation parameter
      arag_k=2.86

      dens_ocn   = 1.03         ! kg/l
!temperature min and max for linear production function (as a function
!of temperature)
      tfmin=18
      tfmax=31
      b_temp=1./(tfmax-tfmin)
      a_temp=1-(tfmax*b_temp)

!to be deleted      water_z_pp = 25

!to be deleted      area_prod=0

!!!      sum_prod_coral=0

!zmax
! surace incoming shortwave radiation light_surf in W/m2, multiplied by PAR_m=0.4
! to have PAR stil in W/m2
!      RAD_m=light_surf !*PAR_m ! coupled version: RAD_m=SABST_O(i,n)*PAR_m ! en W/m2?
!      pk490=Kd
      !print*, 'Imin, RAD_m, pk490', Imin, RAD_m, pk490
!      zmax=-log(Imin/RAD_m)/pk490 ! depth is negative
!      print*, 'depth and zmax', depth, zmax

! temperature limitation
      if ((temp.lt.tmin).or.(temp.gt.tmax)                              &
          .or.(sal.lt.smin).or.(sal.gt.smax)                            &
          !.or.(phos.gt.pmax).or.(nitr.gt.nmax)) then
          !.or.(phos.gt.pmax)) then
          .or.(phos.gt.pmax).or.(temp_too_low.gt.0)) then
!          .or.(phos.gt.pmax).or.(temp_too_low.gt.0)                     &
!          .or. (depth .lt. zmax) ) then

            P_carb = 0
            !print*, 'P_carb=0'

            if (temp.gt.tmax) THEN
                                    ! Temperature excess always leads
                                    ! to strong bleaching
              tau_bleach_l = tau_bleach_strong
             timebleach_l = NYR

            endif

      else

! surace incoming shortwave radiation light_surf in W/m2, multiplied by PAR_m=0.4
! to have PAR stil in W/m2
              RAD_m=light_surf*PAR_m ! coupled version: RAD_m=SABST_O(i,n)*PAR_m ! en W/m2?
              !print*, 'light_surf', light_surf
              !print*, 'PAR_m', PAR_m
              pk490=Kd
              !print*, 'pk490', pk490
              !print*, 'depth', depth
              !pk490=0.041 ! test with fixed homogeneous value
              !Iz en W/m2
              Iz=RAD_m*exp(-pk490*(-1*depth)) !depth passe en positif ! luminosite at depth
              g_coral = gmax_coral*tanh(Iz/Ik) !in m/day
              !print*, 'light limitation ',  gmax_coral, tanh(Iz/Ik),    &
              !         g_coral, Ik, Iz
              !g_coral = gmax_coral !test
              if (Iz.le.Imin) then
                 g_coral=0.0
              endif

!also limit in depth if less than 150m no more corals
!              print*, 'depth', depth
              if (depth.le.-150) then
                 g_coral=0.0
              endif
              !print*, 'light limitation', g_coral
              


! limitation by temperature?
! my old version
!        epaisseur=0.004 !0.012 !0.004
!        x_opti=30 !25
!        y_opti=1
!        temp_factor=-epaisseur*(temp-x_opti)**2+y_opti
        !print*, 'temp_factor ', temp_factor
!        !temp_factor=1.0
!        if (temp_factor.lt.0) temp_factor=0.0

! Temperature dependency
!Manon s version
!              if ((temp.lt.24))then
!                g_coral = g_coral*coef_ftemp*(temp-18)
!              else if ((temp.lt.27).and.(temp.gt.24)) then
!                g_coral = g_coral*coef_ftemp*(24-18)
!              else if (temp.gt.27) then
!                g_coral = g_coral*(coef_ftemp*(temp-26)+coef_ftemp*(24-18))
!              endif

!new simpler version: linear
                if (temp .lt. tfmin) then
                   temp_factor=0
                else if (temp .ge. tfmin .and. temp.le.tfmax) then
                  temp_factor=a_temp+b_temp*temp
                else
                  temp_factor=1
                endif
                g_coral = g_coral*temp_factor


! Limitation by recovering from bleaching
              IF (timebleach_l .NE. 0.) THEN
                g_coral = g_coral*(1.-EXP(-(NYR-timebleach_l)           &
                              /tau_bleach_l))
              ENDIF


! limitation by supersaturation
! Langdon & Atkinson, JGR, 2005
              if (omega_arag.gt.1) then
                g_coral = g_coral*(omega_arag-1.)/arag_k
              else
                g_coral = 0.
              endif
              !print*, 'omega arag ', omega_arag, g_coral

! production in Pmol/day/grid cell
! nb units: g_coral in m/day, area in m2 (*1e-6 to be in km2), coral_dens in kg/m3
! caco3_molar_mass in g/mol
! hence P_carb in Pmol/day
              P_carb = g_coral*topof                                    &
                            *area*1.0e-6                                &
                            *1.0e-6*coral_dens/caco3_molar_mass
              !print*, '*topof*area*coral_dens', P_carb
      endif

!Carbonate dissolution if omega < 1, sinon dissolution toute petite or no dissolution
!      arag_sat=5 ! test a commenter
!      if (arag_sat.lt.1) then
!           D_carb=mass_carb*(1.0-(arag_sat)) !ajout non linearité **2? dissolution complete si omega =0 ou pas (1.2 ou 1))?
!      else
!           D_carb=0.0 !0.001*mass_carb !0.0 ! 0 ou une toute petite valeur constante (genre 1percent) ?
!      endif

!      D_carb=0.0 !test

!nb    dissolution of all existing coral if no production
!       if (P_carb .eq. 0) then
!         D_carb=mass_carb/(caco3_molar_mass*1.0e15)
!       else
!         D_carb=0.0
!       endif

!nb    dissolution of half of all existing coral in 10 years if no production
!       if (P_carb .eq. 0) then
!         if (prod_before .eq. 1) then
!            prod_before=0
!            t0_diss=NYR-1
!          endif
!         D_carb=mass_carb/(caco3_molar_mass*1.0e15)                     &
!               * exp(-lambda_diss*(NYR-t0_diss))
!       else
!         D_carb=0.0
!         prod_before=1
!       endif

! For now no dissolution
      D_carb=0.0


! Net production in Pmol/day
      net_carb=P_carb-D_carb

!New mass in g/day
!      mass_carb_new=net_carb*caco3_molar_mass*1.0e15

!net mass
!      mass_carb_old=mass_carb
!      mass_carb=mass_carb_old+mass_carb_new


!      if (mass_carb_new .gt.0) then
!          print*, ' '
!          print*, 'in coral: prod, mass ',                             &
!          mass_carb_new, mass_carb
!          print*, 'light factor: ', tanh(Iz/Ik)
!          print*, 'omega arag factor: ', (omega_arag-1.)/arag_k
!      endif

      end subroutine corals
!----------------------------------------------------------------------

!----------------------------------------------------------------------
      subroutine out_coral_global
! writing of global annual outputs

!nb write in Coral_output.txt
!      write (coral_res_fich,'(i6,5f14.5)')
!      &
      write (10,'(i6,3f14.5,1f15.2)')                                   &
         NYR, total_area_coral_an, total_prod_coral_an,                 &
         total_mass_coral_an
!         total_mass_coral_an, total_area_coral_an_40m,
!         &
!         total_prod_coral_an_40m

      end subroutine
!----------------------------------------------------------------------

!----------------------------------------------------------------------
      subroutine calc_monthly_temp(i,j,n)
! stores montly mean temperature
! over 30 years
! and update max value for grid i,j,n

      integer, intent(in) :: i,j,n
      REAL, dimension(nb_mois) :: temp_mois_temp

!      temp_mois(i,j,n)=temp_mois(i,j,n)+TM(i,j,n)
      temp_mois(i,j,n)=temp_mois(i,j,n)+TM_pluswkvar(i,j,n)
      if (KMON.eq.1) then ! si dernier jour du mois (jour 30)
          temp_mois(i,j,n)=temp_mois(i,j,n)/30. !monthly mean temperature
          ! shift all previous months and fill last one
          temp_mois_temp=temp_mois_all(i,j,n,:)
          temp_mois_all(i,j,n,1:indice_mois-1)=                         &
                    temp_mois_temp(2:indice_mois)
          temp_mois_all(i,j,n,indice_mois)=temp_mois(i,j,n)
          temp_mois(i,j,n)=0.0
          ! max of climatological monthly mean temperature
!          if (NYR.gt.window_MMM) then
           MMMclim(i,j,n)=MAXVAL(temp_mois_all(i,j,n,:))
!          endif
          !write(*,*) 'MMMclim: ', MMMclim(i,j,n)
      endif


      end subroutine calc_monthly_temp
!----------------------------------------------------------------------

!----------------------------------------------------------------------
      subroutine calc_DHW(i,j,n)
      ! computes degree heating weeks DHW (degree/week) following NOAA

!      use loveclim_transfer_mod, only: TM, KWEEK
!      use loveclim_transfer_mod, only: KWEEK

      integer, intent(in) :: i,j,n
      real HS
      real,dimension(nb_hs) :: xsHS_temp


      !first computes Hot Spot every day
!      if (TM(i,j,n).ge.MMMclim(i,j,n)) then
!        HS=TM(i,j,n)-MMMclim(i,j,n)
      if (TM_pluswkvar(i,j,n).ge.MMMclim(i,j,n)) then
        HS=TM_pluswkvar(i,j,n)-MMMclim(i,j,n)
      else
        HS=0.0
      endif
      !then computes excess Hot Spot if >1 degree and keep value in
      !matrix
      if (HS.ge.1) then
      !!!old version: tbd
        !xsHS(i,j,n)=xsHS(i,j,n)+HS
      !else
      !  xsHS(i,j,n)=0.0
      !!!end old version
        !first fill all values, then shift and replace last value
        if (indice_hs.lt.nb_hs) then !nb_hs=84 days
          xsHS_all(i,j,n,indice_HS)=HS
        else
          xsHS_temp(:)=xsHS_all(i,j,n,:)
          xsHS_all(i,j,n,1:indice_hs-1)=xsHS_temp(2:indice_hs)
          xsHS_all(i,j,n,indice_hs)=HS
        endif
      endif

!!! old version : tbd
!      if (KWEEK.eq.1) then ! if end of week
!        !average over the week to have the weekly excess hot spot
!        xsHS(i,j,n)=xsHS(i,j,n)/5 ! week of 5 days
!
!        !shifts and replaces last index value
!        xsHS_temp(:)=xsHS_all(i,j,n,:)
!        xsHS_all(i,j,n,1:indice_hs-1)=xsHS_temp(2:indice_hs)
!        xsHS_all(i,j,n,indice_hs)=xsHS(i,j,n)
!        !then computes Degree Heating Week=sum of excess
!        !temperature over 17 weeks of 5 days=85 days, equivalent to 12 weeks of 7 days=84days
!        DHW(i,j,n)=SUM(xsHS_all(i,j,n,:))
!        if (DHW(i,j,n).ge.dhw_thresh_moderat) then ! si plus que 4 degree for moderate bleaching
!          DHW_nb(i,j,n)=DHW(i,j,n)+1 ! compte le nombre de weeks qui declenche bleaching
!        endif
!        xsHS(i,j,n)=0.0 ! set back to 0
!        !write(*,*) 'DHW: ', DHW(i,j,n)
!      endif
!!! end old version


        !then computes Degree Heating Week=sum of excess
        !temperature over 84 days / 7 (in degree per week)
        DHW(i,j,n)=SUM(xsHS_all(i,j,n,:)/7.)
        if (DHW(i,j,n).ge.dhw_thresh_moderat) then ! si plus que 4 degree/week for moderate bleaching
          DHW_nb(i,j,n)=DHW_nb(i,j,n)+1 ! compte le nombre de weeks qui declenche bleaching
        endif
!tbd        xsHS(i,j,n)=0.0 ! set back to 0
      end subroutine calc_DHW
!----------------------------------------------------------------------

!----------------------------------------------------------------------
      subroutine calc_bleach(i,j,n)

      !For bleaching, from Guy Munhoven based on NOAA
      !computes timebleach and tau_bleach used in corals

      integer, intent(in) :: i,j,n

                                    ! Degree-Heating Weeks control
                                    ! ----------------------------
          IF (DHW(i,j,n) .GE. dhw_thresh_strong) THEN
                                    ! Strong bleaching event detected!
            tau_bleach(i,j,n) = tau_bleach_strong
            timebleach(i,j,n) = NYR

          ELSEIF (DHW(i,j,n) .GE. dhw_thresh_moderat) THEN
                                    ! Moderate bleaching event detected
            IF (tau_bleach(i,j,n) .LT. tau_bleach_moderat) THEN
                                    ! Corals in healthy state (i.e., not
                                    ! recovering from a previous
                                    ! bleaching event)
                                    ! simply set the bleaching
                                    ! parameters
              tau_bleach(i,j,n) = tau_bleach_moderat
              timebleach(i,j,n) = NYR

            ELSEIF (tau_bleach(i,j,n) .LT. tau_bleach_strong) THEN
                                    ! Corals are recovering from a
                                    ! previous
                                    ! moderate event.
              IF ((NYR - timebleach(i,j,n)) .LE. 2.) THEN
                                    ! If this event was less than 2
                                    ! years ago,

                tau_bleach(i,j,n) = tau_bleach_strong
                timebleach(i,j,n) = NYR

              ELSE
                                    ! This event was more than 2 years
                                    ! ago,
                                    ! leave in moderate state, but reset
                                    ! bleaching date.
                timebleach(i,j,n) = NYR

              ENDIF

            ELSE
                                    ! Corals are recovering from a
                                    ! previous
                                    ! strong event. Leave in recovery
                                    ! from
                                    ! strong but reset the bleaching
                                    ! date.
              timebleach(i,j,n) = NYR

            ENDIF

          ELSE
                                    ! No bleaching event triggered by
                                    ! DHW
                                    ! thresholds. Check if bleaching
                                    ! control
                                    ! can possibly be reset (i.e.,
                                    ! previous
                                    ! bleaching events were sufficiently
                                    ! long ago -- typically more than
                                    ! 4*tau).
            IF (timebleach(i,j,n) .NE. 0.) THEN
                                    ! Coral is currently in a recovery
                                    ! phase:
              IF ((NYR-timebleach(i,j,n)) .GT. 4.*tau_bleach(i,j,n))    &
                 THEN                  ! It was triggered more than 4*tau ago:
                timebleach(i,j,n) = 0.      ! Reset bleaching time
                tau_bleach(i,j,n) = 0.      ! Reset recovery time-scale
              ENDIF

            ENDIF

          ENDIF



                                    ! Report moderate and strong hot
                                    ! events,
                                    ! but limit print-out of hot events
                                    ! to
                                    ! temperatures above 25 degC to
                                    ! ignore
                                    ! numerous cold "hot events"
!          IF (TMMM(i,n).gt.25.) THEN
!            IF (DHW_curr(i,n).ge.dhw_thresh_strong) THEN
!              write(1000,*) 'Strong Hot event ', i, n, NYR, NJUL,
!     @                    DHW_curr(i,n)
!            ELSEIF (DHW_curr(i,n).ge.dhw_thresh_moderat) THEN
!              write(1000,*) 'Moderate Hot event ', i, n, NYR, NJUL,
!     @                    DHW_curr(i,n)
!            ELSE
!              CONTINUE
!            ENDIF
!          ENDIF



      end subroutine calc_bleach
!----------------------------------------------------------------------



       end module coral_mod
