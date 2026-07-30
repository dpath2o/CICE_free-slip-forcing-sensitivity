#ifdef ncdf
#define USE_NETCDF
#endif
!=======================================================================
!
! Reads and interpolates forcing data for atmosphere and ocean quantities.
!
! authors: Elizabeth C. Hunke and William H. Lipscomb, LANL
!
! 2004 WHL: Block structure added
! 2005 WHL: ECMWF option added
! 2006 ECH: LY option added
! 2006 WHL: Module name changed from ice_flux_in
! 2006 ECH: Fixed bugs, rearranged routines, edited comments, etc.
!           Added NCAR ocean forcing file
!           Converted to free source form (F90)
! 2007: netcdf version of read_data added by Alison McLaren, Met Office
! 2024 DPA: added hourly-ERA5 and daily/weekly/monthly ocean forcing
! 2025 DPA: revised daily/weekly/monthly ocean forcing to be inline with CICE princples
!

module ice_forcing
  use ice_kinds_mod
  use ice_boundary     , only: ice_HaloUpdate
  use ice_blocks       , only: nx_block, ny_block
  use ice_domain       , only: halo_info
  use ice_domain_size  , only: ncat, max_blocks, nx_global, ny_global
  use ice_communicate  , only: my_task, master_task
  use ice_calendar     , only: istep, istep1, msec, mday, mmonth, myear, yday, daycal, daymo, days_per_year, compute_days_between
  use ice_fileunits    , only: nu_diag, nu_forcing
  use ice_exit         , only: abort_ice
  use ice_read_write   , only: ice_open, ice_read, ice_get_ncvarsize, ice_read_vec_nc, ice_open_nc, ice_read_nc, ice_close_nc
  use ice_timers       , only: ice_timer_start, ice_timer_stop, timer_readwrite, timer_bound, timer_forcing
  use ice_arrays_column, only: oceanmixed_ice, restore_bgc
  use ice_constants    , only: c0, c1, c2, c3, c4, c5, c8, c10, c12, c15, c20, c180, c360, c365, c1000, c3600
  use ice_constants    , only: p001, p01, p1, p2, p25, p5, p6
  use ice_constants    , only: cm_to_m
  use ice_constants    , only: field_loc_center, field_type_scalar, field_type_vector, field_loc_NEcorner
  use icepack_intfc    , only: icepack_warnings_flush, icepack_warnings_aborted
  use icepack_intfc    , only: icepack_sea_freezing_temperature
  use icepack_intfc    , only: icepack_init_wave, icepack_init_parameters
  use icepack_intfc    , only: icepack_query_tracer_indices, icepack_query_parameters
  ! define subroutines that other modules will use
  implicit none
  private
  public :: alloc_forcing, init_forcing_atmo, init_forcing_ocn, init_snowtable, &
       init_tides_metadata, alloc_tides_fields, init_tides_fields, &
       compute_tides_at_time, &
       get_forcing_atmo, get_forcing_ocn, get_wave_spec, &
       read_clim_data, read_clim_data_nc, read_data_nc_point, &
       interpolate_data, interp_coeff_monthly, interp_coeff
  ! create variable namespace and allocate types
  ! year information
  integer (kind=int_kind), public :: &
       ycycle          , & ! number of years in forcing cycle, set by namelist
       fyear_init      , & ! first year of data in forcing cycle, set by namelist
       fyear           , & ! current year in forcing cycle, varying during the run
       fyear_final         ! last year in cycle, computed at init
  ! TIDES
  integer (kind=int_kind), public :: &
       tide_nconst = 0, &
       tide_nj = 0, &
       tide_ni = 0
  logical (kind=log_kind), public :: &
       tide_metadata_loaded = .false.
  real (kind=dbl_kind), dimension(:), allocatable, public :: &
       tide_omega, &
       tide_phase, &
       tide_alpha, &
       tide_amplitude
  character(len=256), public :: &
       tide_constituent_order = ''
  logical (kind=log_kind), public :: &
       tide_fields_loaded = .false.
  real (kind=dbl_kind), dimension(:,:,:,:), allocatable, public :: &
       tide_hRe, &
       tide_hIm, &
       tide_URe, &
       tide_UIm, &
       tide_VRe, &
       tide_VIm
  real (kind=dbl_kind), dimension(:,:,:), allocatable, public :: &
       tide_wct, &
       tide_mask, &
       tide_depth_fac, & ! 0..1 reliability taper from CATS wct and CICE bathymetry
       tide_h_eff       ! effective depth for transport-to-current conversion [m]
  ! Tidal-current diagnostics.
  !
  ! These fields are diagnostic only.  They separate the harmonic tide
  ! pathway from the background ocean forcing and allow offline analysis of
  ! sub-daily current variability without writing every dynamics substep.
  !
  ! raw  : current reconstructed from harmonics before final cap/limiter use
  ! eff  : final perturbation actually added to uocn/vocn
  ! day* : within-model-day accumulators reset when myear/mmonth/mday changes
  real (kind=dbl_kind), dimension(:,:,:), allocatable, public :: &
       tide_u_raw,          &
       tide_v_raw,          &
       tide_speed_raw,      &
       tide_u_eff,          &
       tide_v_eff,          &
       tide_speed_eff,      &
       tide_speed_daymax,   &
       tide_speed_dayrms,   &
       tide_speed_daysum,   &
       tide_speed2_daysum,  &
       tide_n_subday,       &
       tide_n_curr_over_fi, &
       tide_n_capped
  real (kind=dbl_kind), public :: &
       tide_curr_diag_thresh = 5.0e-4_dbl_kind   ! m/s; diagnostic only
  integer (kind=int_kind), public :: &
       tide_diag_myear  = -9999, &
       tide_diag_mmonth = -9999, &
       tide_diag_mday   = -9999
  ! input data file names
  character (char_len_long) :: &
       F_uwind, &
       F_vwind, &
       F_wind, &
       F_strax, &
       F_stray, &
       F_tair, &
       F_humid, &
       F_rhoa, &
       F_fsw, &
       F_flw, &
       F_rain, &
       F_sst, &
       F_sss, &
       F_sublim, &
       F_snow, &
       F_ERA5, &
       F_ocn, &
       F_AFIM
  ! input data file names
  character (char_len_long), dimension(:), allocatable, public :: topmelt_file, botmelt_file
  ! interpolation coefficients
  real (kind=dbl_kind), public :: c1intp, c2intp
  ! old record number (save between steps)
  integer (kind=int_kind) :: oldrecnum = 0, oldrecnum4X = 0
  ! cloud fraction
  real (kind=dbl_kind), dimension(:,:,:), allocatable, public :: cldf
  ! 3-d field values at 2 temporal data points
  real (kind=dbl_kind), dimension(:,:,:,:), allocatable, public :: &
       fsw_data,    &
       cldf_data,   &
       fsnow_data,  &
       Tair_data,   &
       uatm_data,   &
       vatm_data,   &
       wind_data,   &
       strax_data,  &
       stray_data,  &
       Qa_data,     &
       rhoa_data,   &
       flw_data,    &
       ocn_data,    &
       sst_data,    &
       sss_data,    &
       uocn_data,   &
       vocn_data,   &
       sublim_data, &
       frain_data,   &
       pair_data,    &
       snowfall_data,&
       blh_data,     &
       windgust_data,&
       uatm100_data, &
       vatm100_data
  ! 4-d (categorical) field values at 2 temporal data points
  real (kind=dbl_kind), dimension(:,:,:,:,:), allocatable, public :: topmelt_data, botmelt_data
  ! data formats and types
  character(char_len), public :: &
       atm_data_format,  & ! 'bin'=binary or 'nc'=netcdf
       ocn_data_format,  & ! 'bin'=binary or 'nc'=netcdf
       atm_data_type,    & ! 'default', 'monthly', 'ncar', 'box2001', 'hadgem', 'oned', 'calm', 'uniform', 'ERA5', 'JRA55_***'
       atm_data_version, & ! date of atm_forcing file creation
       bgc_data_type,    & ! 'default', 'clim'
       ocn_data_type,    & ! 'default', 'clim', 'ncar', 'oned', 'calm', 'box2001', 'hadgem_sst' or 'hadgem_sst_uvocn', 'uniform', 'AFIM'
       era5_mod_var,     & ! <--------- ONLY VALID IN ERA5_data --------->
                           !    'none'
                           !    'ttlpcp_ant_coast' or 'precip_ant_coast'
                           !    Future options may include:
                           !    'u','v','wspd','tair','qair','sw','lw'
       tide_data_type,   & ! 'none' or 'harmonic'
       tide_data_format, & ! 'none' or 'CICE_TMD3'
       ice_data_type,    & ! 'latsst', 'box2001', 'boxslotcyl', etc
       ice_data_conc,    & ! 'p5','p8','p9','c1','parabolic', 'box2001', etc
       ice_data_dist,    & ! 'box2001','gauss', 'uniform', etc
       precip_units        ! 'mm_per_month', 'mm_per_sec', 'mks','m_per_sec'
  ! modified atmospheric forcing field scale factor
  real (kind=dbl_kind), public :: era5_mod_fac ! <--------- ONLY VALID IN ERA5_data --------->
                                               ! a multiplicative for scaling atmospheric data
                                               ! DEFAULT: 1; multiply era5_mod_var by whatever this value is set to
  ! rotate wind/stress to computational grid from true north directed
  logical (kind=log_kind), public :: rotate_wind
  ! directory information
  character(char_len_long), public :: &
       atm_data_dir,    & ! top directory for atmospheric data
       ocn_data_dir,    & ! top directory for ocean data
       wave_spec_dir,   & ! dir name for wave spectrum
       wave_spec_file,  & ! file name for wave spectrum
       oceanmixed_file, & ! file name for ocean forcing data
       tide_data_file     ! harmonic tide file on CICE T-grid
  ! number of fields to search for in forcing file
  integer (kind=int_kind), parameter :: nfld = 8, ndfld = 4
 ! number of months and days in year for ocean forcing
  integer (kind=int_kind), parameter :: m_per_yr = 12, wk_per_yr = 52
  ! as in the dummy atm (latm)
  real (kind=dbl_kind), parameter, public :: frcvdr = 0.28_dbl_kind, & ! frac of incoming sw in vis direct band
                                             frcvdf = 0.24_dbl_kind, & ! frac of incoming sw in vis diffuse band
                                             frcidr = 0.31_dbl_kind, & ! frac of incoming sw in near IR direct band
                                             frcidf = 0.17_dbl_kind    ! frac of incoming sw in near IR diffuse band
  ! ocn data for 12 months
  real (kind=dbl_kind), dimension (:,:,:,:,:), allocatable, public :: ocn_frc_m, ocn_frc_d, ocn_frc_w
  ! restore sst if true
  logical (kind=log_kind), public :: restore_ocn
  ! restoring time scale (days)
  integer (kind=int_kind), public :: trestore
  ! AFIM_data() subroutine specific
  real(kind=dbl_kind), public :: hmix_0  ! initial (or fixed if standalone) mixed layer depth
  ! restoring time scale (sec)
  real (kind=dbl_kind), public :: trest
  ! tidal switches
  logical(kind=log_kind), public :: &
       tide_use_currents, & ! apply tidal u/v from harmonics
       tide_use_ssh,      & ! apply tidal eta -> ss_tltx/ss_tlty
       tide_use_bathymetry_limit = .true. ! use CICE bathymetry with CATS wct to taper currents
  real (kind=dbl_kind), public :: &
       tide_curr_fac       = 1.0_dbl_kind,  & ! multiply tidal currents by this factor
       tide_speed_cap      = 2.0_dbl_kind,  & ! emergency cap on tidal current vector speed [m/s]
       tide_wct_min        = 20.0_dbl_kind, & ! zero tide currents at/below this reliable depth [m]
       tide_wct_full       = 80.0_dbl_kind, & ! full tide currents at/above this reliable depth [m]
       tide_h_eff_min      = 50.0_dbl_kind, & ! minimum effective depth for U/h conversion [m]
       tide_depth_mismatch = 4.0_dbl_kind,  & ! diagnostic threshold for CATS/CICE depth mismatch ratio
       tide_ramp_days      = 20.0_dbl_kind    ! half-cosine ramp-up time from branch start [days]
  ! prints forcing debugging output if true
  logical (kind=log_kind), public :: debug_forcing
  ! jday time vector from atm forcing files
  real (dbl_kind), dimension(:), allocatable, public :: jday_atm
  ! Number of atm forcing timesteps
  integer (kind=int_kind), public :: Njday_atm
  ! filename for snow lookup table
  character (len=char_len_long), public :: snw_filename
  ! snow table 
  character (char_len), public :: snw_rhos_fname , &  ! snow table 1d rhos field name 
                                  snw_Tgrd_fname , &  ! snow table 1d Tgrd field name
                                  snw_T_fname    , &  ! snow table 1d T field name
                                  snw_tau_fname  , &  ! snow table 3d tau field nam
                                  snw_kappa_fname, &  ! snow table 3d kappa field name
                                  snw_drdt0_fname     ! snow table 3d drdt0 field name
  ! PRIVATE:
  ! default mixed layer depth in m
  real (dbl_kind), parameter :: mixed_layer_depth_default = 60
  ! local debug flag
  logical (kind=log_kind), parameter :: local_debug = .false.

  !=======================================================================

contains

  !=======================================================================
  ! Allocate space for all variables
  subroutine alloc_forcing
    integer (int_kind) :: ierr
    character(len=*), parameter :: subname = '(alloc_forcing)'
    if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'
    allocate ( cldf        (nx_block,ny_block,max_blocks), & ! cloud fraction
         fsw_data    (nx_block,ny_block,2,max_blocks), & ! field values at 2 temporal data points
         cldf_data   (nx_block,ny_block,2,max_blocks), &
         fsnow_data  (nx_block,ny_block,2,max_blocks), &
         Tair_data   (nx_block,ny_block,2,max_blocks), &
         uatm_data   (nx_block,ny_block,2,max_blocks), &
         vatm_data   (nx_block,ny_block,2,max_blocks), &
         wind_data   (nx_block,ny_block,2,max_blocks), &
         strax_data  (nx_block,ny_block,2,max_blocks), &
         stray_data  (nx_block,ny_block,2,max_blocks), &
         Qa_data     (nx_block,ny_block,2,max_blocks), &
         rhoa_data   (nx_block,ny_block,2,max_blocks), &
         flw_data    (nx_block,ny_block,2,max_blocks), &
         ocn_data    (nx_block,ny_block,2,max_blocks), &
         sst_data    (nx_block,ny_block,2,max_blocks), &
         sss_data    (nx_block,ny_block,2,max_blocks), &
         uocn_data   (nx_block,ny_block,2,max_blocks), &
         vocn_data   (nx_block,ny_block,2,max_blocks), &
         sublim_data (nx_block,ny_block,2,max_blocks), &
         frain_data  (nx_block,ny_block,2,max_blocks), &
         pair_data   (nx_block,ny_block,2,max_blocks), &
         snowfall_data(nx_block,ny_block,2,max_blocks), &
         blh_data    (nx_block,ny_block,2,max_blocks), &
         windgust_data(nx_block,ny_block,2,max_blocks), &
         uatm100_data(nx_block,ny_block,2,max_blocks), &
         vatm100_data(nx_block,ny_block,2,max_blocks), &
         topmelt_data(nx_block,ny_block,2,max_blocks,ncat), &
         botmelt_data(nx_block,ny_block,2,max_blocks,ncat), &
         ocn_frc_m   (nx_block,ny_block,max_blocks,nfld,m_per_yr), &
         ocn_frc_d   (nx_block,ny_block,max_blocks,nfld,days_per_year), &
         ocn_frc_w   (nx_block,ny_block,max_blocks,nfld,wk_per_yr), &
         topmelt_file(ncat), &
         botmelt_file(ncat), &
         stat=ierr)
    if (ierr/=0) call abort_ice('(alloc_forcing): Out of Memory')
    ! initialize this, not set in box2001 (and some other forcings?)
    cldf          = c0
    pair_data     = c0
    snowfall_data = c0
    frain_data    = c0
    blh_data      = c0
    windgust_data = c0
    uatm100_data  = c0
    vatm100_data  = c0
  end subroutine alloc_forcing

  !=======================================================================
  subroutine init_forcing_atmo
    ! Determine the current and final year of the forcing cycle based on
    ! namelist input; initialize the atmospheric forcing data filenames.
    use ice_calendar, only: use_leap_years
    integer (kind=int_kind)     :: modadj   ! adjustment for mod function
    character(len=*), parameter :: subname = '(init_forcing_atmo)'
    if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'
    if (ycycle == 0) then
       fyear       = myear
       fyear_final = myear
    else
       modadj      = abs((min(0,myear-fyear_init)/ycycle+1)*ycycle)
       fyear       = fyear_init + mod(myear - fyear_init + modadj, ycycle)
       fyear_final = fyear_init + ycycle - 1
    endif
    if (my_task == master_task) then
       write(nu_diag,*) subname,' fyear, fyear_init, fyear_final : ',fyear,fyear_init,fyear_final
       write(nu_diag,*) subname,' atm_data_type                  : ',trim(atm_data_type)
    endif
    if (my_task == master_task) then
       write (nu_diag,*) ' Initial forcing data year = ',fyear
       write (nu_diag,*) ' Final   forcing data year = ',fyear_final
    endif
    if (trim(atm_data_type) == 'hadgem' .and. trim(precip_units) /= 'mks') then
       if (my_task == master_task) then
          write (nu_diag,*) 'WARNING: HadGEM atmospheric data chosen with wrong precip_units'
          write (nu_diag,*) 'WARNING:   Changing precip_units to mks (i.e. kg/m2 s).'
       endif
       call abort_ice(error_message=subname//' HadGEM precip_units error', file=__FILE__, line=__LINE__)
    endif
    if (use_leap_years .and. (trim(atm_data_type) /= 'JRA55_gx1' .and. &
         trim(atm_data_type) /= 'JRA55_gx3' .and. &
         trim(atm_data_type) /= 'JRA55_tx1' .and. &
         trim(atm_data_type) /= 'hycom' .and. &
         trim(atm_data_type) /= 'box2001' .and. &
         trim(atm_data_type) /= 'ERA5')) then
       write(nu_diag,*) 'use_leap_years option is currently only supported for: JRA55, hycom, box2001, and ERA5'
       call abort_ice(error_message=subname, file=__FILE__, line=__LINE__)
    endif
    !-------------------------------------------------------------------
    ! Get filenames for input forcing data
    !-------------------------------------------------------------------
    ! default forcing values from init_flux_atm
    if (trim(atm_data_type) == 'ncar') then
       call NCAR_files(fyear)
    elseif (trim(atm_data_type) == 'JRA55_gx1') then
       call JRA55_gx1_files(fyear)
    elseif (trim(atm_data_type) == 'JRA55_gx3') then
       call JRA55_gx3_files(fyear)
    elseif (trim(atm_data_type) == 'JRA55_tx1') then
       call JRA55_tx1_files(fyear)
    elseif (trim(atm_data_type) == 'hadgem') then
       call hadgem_files(fyear)
    elseif (trim(atm_data_type) == 'monthly') then
       call monthly_files(fyear)
    elseif (trim(atm_data_type) == 'oned') then
       call oned_files
    elseif (trim(atm_data_type) == 'ISPOL') then
       call ISPOL_files
    elseif (trim(atm_data_type) == 'box2001') then
       call box2001_data_atm
    elseif (trim(atm_data_type) == 'uniform_northeast') then
       call uniform_data_atm('NE')
    elseif (trim(atm_data_type) == 'uniform_north') then
       call uniform_data_atm('N')
    elseif (trim(atm_data_type) == 'uniform_east') then
       call uniform_data_atm('E')
    elseif (trim(atm_data_type) == 'uniform_south') then
       call uniform_data_atm('S')
    elseif (trim(atm_data_type) == 'uniform_west') then
       call uniform_data_atm('W')
    elseif (trim(atm_data_type) == 'calm') then
       call uniform_data_atm('N',c0) ! direction does not matter when c0
    elseif (trim(atm_data_type) == 'hycom') then
       call hycom_atm_files
    elseif (trim(atm_data_type) == 'ERA5') then
       call ERA5_monthly_files(fyear, mmonth)
    elseif (trim(atm_data_type) == 'default') then
       ! don't need to do anything more
    else
       call abort_ice (error_message=subname//' ERROR atm_data_type unknown = '// trim(atm_data_type), file=__FILE__, line=__LINE__)
    endif
  end subroutine init_forcing_atmo

  !=======================================================================
  subroutine init_forcing_ocn(dt)
    ! Set sea surface salinity and freezing temperature to annual mean value using a 12-month climatology.
    ! Read sst data for current month, and adjust sst based on freezing temperature.
    ! No interpolation in time.
    ! SST is subsequently prognosed if CICE is run with a mixed layer ocean (oceanmixed_ice = T), and can be restored to data (restore_ocn = T).
    use ice_blocks, only     : nx_block, ny_block
    use ice_domain, only     : nblocks
    use ice_domain_size, only: max_blocks 
    use ice_flux, only       : sss, sst, Tf, uocn, vocn, hmix, frzmlt, qdp, ss_tltx, ss_tlty
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan
    real (kind=dbl_kind), intent(in) :: dt            ! this is namelist input parameter
    integer (kind=int_kind)          :: i, j, iblk, & ! horizontal indices
                                        k         , & ! month index
                                        fid       , & ! file id for netCDF file
                                        nbits
    logical (kind=log_kind)          :: diag
    real (kind=dbl_kind)             :: secday
    character (char_len)             :: fieldname
    real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: work1
    integer (kind=int_kind) :: modadj   ! adjustment for mod function
    character(len=*), parameter :: subname = '(init_forcing_ocn)'
    ! double precision data
    nbits = 64
    if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'
    call icepack_query_parameters(secday_out=secday)
    call icepack_warnings_flush(nu_diag)
    if (icepack_warnings_aborted()) call abort_ice(error_message=subname, file=__FILE__, line=__LINE__)
    call alloc_forcing()
    sst_data (:,:,:,:) = c0
    sss_data (:,:,:,:) = c0
    uocn_data(:,:,:,:) = c0
    vocn_data(:,:,:,:) = c0
    if (restore_ocn .or. restore_bgc) then
       ! If restore ocean values is set to 0 in namelist then use ocean forcing data at each thermodynamic step
       if (trestore == 0) then 
          trest = dt
       else
          trest = real(trestore,kind=dbl_kind) * secday ! seconds
       endif
    endif
    !-------------------------------------------------------------------
    ! Sea surface salinity (SSS)
    ! initialize to annual climatology created from monthly data
    !-------------------------------------------------------------------
    if (trim(ocn_data_type) == 'clim') then
       F_sss = trim(ocn_data_dir)//'/sss.mm.100x116.da' ! gx3 only
       if (my_task == master_task) then
          write (nu_diag,*) ' '
          write (nu_diag,*) 'SSS climatology computed from:'
          write (nu_diag,*) trim(F_sss)
       endif
       if (my_task == master_task) call ice_open(nu_forcing, F_sss, nbits)
       sss(:,:,:) = c0
       do k = 1,12            ! loop over 12 months
          call ice_read (nu_forcing, k, work1, 'rda8', debug_forcing, field_loc_center, field_type_scalar)
          !$OMP PARALLEL DO PRIVATE(iblk,i,j)
          do iblk = 1, nblocks
             do j = 1, ny_block
                do i = 1, nx_block
                   sss(i,j,iblk) = sss(i,j,iblk) + work1(i,j,iblk)
                enddo
             enddo
          enddo !$OMP END PARALLEL DO
       enddo                  ! k
       !$OMP PARALLEL DO PRIVATE(iblk,i,j)
       do iblk = 1, nblocks
          do j = 1, ny_block
             do i = 1, nx_block
                sss(i,j,iblk) = sss(i,j,iblk) / c12   ! annual average
                sss(i,j,iblk) = max(sss(i,j,iblk),c0)
             enddo
          enddo
       enddo !$OMP END PARALLEL DO
       call ocn_freezing_temperature
       if (my_task == master_task) close(nu_forcing)
       !-------------------------------------------------------------------
       ! Sea surface temperature (SST)
       ! initialize to data for current month
       !-------------------------------------------------------------------
       if (nx_global == 320) then ! gx1
          F_sst = trim(ocn_data_dir)//'/sst_clim_hurrell.dat'
       else                   ! gx3
          F_sst = trim(ocn_data_dir)//'/sst.mm.100x116.da'
       endif
       if (my_task == master_task) then
          write (nu_diag,*) ' '
          write (nu_diag,*) 'Initial SST file:', trim(F_sst)
       endif
       if (my_task == master_task) call ice_open(nu_forcing, F_sst, nbits)
       call ice_read(nu_forcing, mmonth, sst, 'rda8', debug_forcing, field_loc_center, field_type_scalar)
       if (my_task == master_task) close(nu_forcing)
       ! Make sure sst is not less than freezing temperature Tf
       !$OMP PARALLEL DO PRIVATE(iblk,i,j)
       do iblk = 1, nblocks
          do j = 1, ny_block
             do i = 1, nx_block
                sst(i,j,iblk) = max(sst(i,j,iblk),Tf(i,j,iblk))
             enddo
          enddo
       enddo !$OMP END PARALLEL DO
    elseif (trim(ocn_data_type) == 'hadgem_sst' .or. trim(ocn_data_type) == 'hadgem_sst_uvocn') then
       diag = .true.   ! write diagnostic information
       F_sst = trim (ocn_data_dir)//'/MONTHLY/sst.1997.nc'
       if (my_task == master_task) then
          write (nu_diag,*) ' '
          write (nu_diag,*) 'Initial SST file:', trim(F_sst)
          call ice_open_nc(F_sst,fid)
       endif
       fieldname='sst'
       call ice_read_nc(fid,mmonth,fieldname,sst,diag)
       if (my_task == master_task) call ice_close_nc(fid)
       ! Make sure sst is not less than freezing temperature Tf
       !$OMP PARALLEL DO PRIVATE(iblk,i,j)
       do iblk = 1, nblocks
          do j = 1, ny_block
             do i = 1, nx_block
                sst(i,j,iblk) = max(sst(i,j,iblk),Tf(i,j,iblk))
             enddo
          enddo
       enddo !$OMP END PARALLEL DO
    elseif (trim(ocn_data_type) == 'ncar') then
       call ocn_data_ncar_init
       ! call ocn_data_ncar_init_3D
    elseif (trim(ocn_data_type) == 'hycom') then
       call ocn_data_hycom_init
    elseif (trim(ocn_data_type) == 'box2001') then
       call box2001_data_ocn
    ! uniform forcing options
    elseif (trim(ocn_data_type) == 'uniform_northeast') then
       call uniform_data_ocn('NE',p1)
    elseif (trim(ocn_data_type) == 'uniform_east') then
       call uniform_data_ocn('E',p1)
    elseif (trim(ocn_data_type) == 'uniform_north') then
       call uniform_data_ocn('N',p1)
    elseif (trim(ocn_data_type) == 'calm') then
       call uniform_data_ocn('N',c0) ! directon does not matter for c0
    elseif (trim(ocn_data_type) == 'default') then
       ! don't need to do anything more
    elseif (trim(ocn_data_type) == 'AFIM') then
       ! Initialize forcing year
       if (ycycle == 0) then
          fyear       = myear
          fyear_final = myear
       else
          modadj      = abs((min(0,myear-fyear_init)/ycycle+1)*ycycle)
          fyear       = fyear_init + mod(myear - fyear_init + modadj, ycycle)
          fyear_final = fyear_init + ycycle - 1
       endif
       if (debug_forcing .or. local_debug) then
         if (my_task == master_task) then
            write(nu_diag,*) subname, ' AFIM ocean forcing initialization'
            write(nu_diag,*) subname, ' fyear, fyear_init, fyear_final: ', fyear, fyear_init, fyear_final
            if (restore_ocn) write(nu_diag,*) subname, ' SST restoring timescale (days): ', trestore
         endif
       endif
       ! Set up forcing file for current year
       call AFIM_files(fyear)
       if (debug_forcing .or. local_debug) then
         if (my_task == master_task) then
            write(nu_diag,*) subname, ' Reading initial ocean state from: ', trim(F_AFIM)
         endif
       endif
       ! Read SSS first (we need it for Tf calculation)
       call ice_open_nc(F_AFIM, fid)
       call ice_read_nc(fid, 1, 'sss', sss, debug_forcing, field_loc=field_loc_center, field_type=field_type_scalar)
       call ice_read_nc(fid, 1, 'u', uocn, debug_forcing, field_loc=field_loc_center, field_type=field_type_scalar)
       call ice_read_nc(fid, 1, 'v', vocn, debug_forcing, field_loc=field_loc_center, field_type=field_type_scalar)
       call ice_close_nc(fid)
       ! Set default values for SSS, uocn, vocn where needed
       do iblk = 1, nblocks
          do j = 1, ny_block
             do i = 1, nx_block
                ! SSS bounds
                if (sss(i,j,iblk) < 0.0_dbl_kind .or. sss(i,j,iblk) > 50.0_dbl_kind .or. sss(i,j,iblk) /= sss(i,j,iblk)) then
                   sss(i,j,iblk) = 34.0_dbl_kind
                endif
                ! Velocity bounds
                if (abs(uocn(i,j,iblk)) > 10.0_dbl_kind .or. uocn(i,j,iblk) /= uocn(i,j,iblk)) then
                   uocn(i,j,iblk) = c0
                endif
                if (abs(vocn(i,j,iblk)) > 10.0_dbl_kind .or. vocn(i,j,iblk) /= vocn(i,j,iblk)) then
                   vocn(i,j,iblk) = c0
                endif
             enddo
          enddo
       enddo
       ! TIDES
       if (trim(tide_data_type) == 'harmonic') then
          call init_tides_metadata
          call alloc_tides_fields
          call init_tides_fields
       endif
       ! Compute freezing temperature from SSS
       call ocn_freezing_temperature
       ! Initialize SST to Tf (freezing point) for ALL points
       ! This ensures SST always has a valid value before restoring is applied
       do iblk = 1, nblocks
          do j = 1, ny_block
             do i = 1, nx_block
                sst(i,j,iblk) = Tf(i,j,iblk)
             enddo
          enddo
       enddo
       !----------------------------------------------------------------
       ! initialize mixed layer variables for oceanmixed_ice!
       ! Without this, the ocean mixed layer calculation produces NaN.
       !----------------------------------------------------------------
       hmix   (:,:,:) = hmix_0   ! ocean mixed layer depth (60m default)
       frzmlt (:,:,:) = c0       ! freezing/melting potential (W/m^2)
       qdp    (:,:,:) = c0       ! deep ocean heat flux (W/m^2)
       ss_tltx(:,:,:) = c0       ! sea surface tilt (m/m)
       ss_tlty(:,:,:) = c0
       if (debug_forcing .or. local_debug) then
         if (my_task == master_task) then
            write(nu_diag,*) subname, ' AFIM initialization complete'
            write(nu_diag,*) subname, ' SST initialized to Tf, hmix = ', mixed_layer_depth_default
         endif
       endif
    else
       call abort_ice (error_message=subname//' ERROR ocn_data_type unknown = '// trim(ocn_data_type), file=__FILE__, line=__LINE__)
    endif
  end subroutine init_forcing_ocn

  !=======================================================================
  !=======================================================================
  ! Initialise metadata for harmonic tidal forcing.
  !
  ! This routine reads only lightweight/global metadata from the tide
  ! coefficient file:
  !
  !   - number of tidal constituents
  !   - tide-file horizontal dimensions
  !   - constituent angular frequencies and phase metadata
  !   - optional global attribute describing constituent order
  !
  ! The large spatial coefficient fields are NOT read here. They are read
  ! later by init_tides_fields().
  !
  ! Design notes:
  !   * Only master_task opens the NetCDF file directly.
  !   * Metadata are then broadcast to all MPI tasks.
  !   * The tide grid is required to match the active CICE grid exactly
  !     because this prototype file has already been regridded to the CICE
  !     T-grid. No runtime horizontal interpolation is performed here.
  !   * The harmonic phase clock is tied to the CATS/TMD reference epoch
  !     used when constructing the coefficient file.
  !
  ! Interpretation:
  !   Successful completion of this routine means that CICE can identify
  !   the tidal constituents and grid geometry. It does NOT yet mean that
  !   tidal currents are dynamically active; that requires init_tides_fields()
  !   and compute_tides_at_time().
  !=======================================================================
  subroutine init_tides_metadata
#ifdef USE_NETCDF
    use netcdf, only: nf90_noerr, nf90_global, nf90_inq_dimid, nf90_inquire_dimension, &
                      nf90_inq_varid, nf90_get_var, nf90_get_att, nf90_strerror
#endif
    use ice_broadcast, only: broadcast_scalar
    integer (kind=int_kind) :: ncid
    integer (kind=int_kind) :: dimid
    integer (kind=int_kind) :: varid
    integer (kind=int_kind) :: ierr
    integer (kind=int_kind) :: n
    character(len=*), parameter :: subname = '(init_tides_metadata)'
    if (trim(tide_data_type) /= 'harmonic') return
    if (tide_metadata_loaded) return
#ifndef USE_NETCDF
    call abort_ice(error_message=subname//' ERROR tide metadata requires USE_NETCDF', &
         file=__FILE__, line=__LINE__)
#else
    if (trim(tide_data_format) /= 'CICE_TMD3') then
       call abort_ice(error_message=subname//' ERROR unsupported tide_data_format = '//trim(tide_data_format), &
            file=__FILE__, line=__LINE__)
    endif
    if (trim(tide_data_file) == 'unknown_tide_file') then
       call abort_ice(error_message=subname//' ERROR tide_data_file not set', &
            file=__FILE__, line=__LINE__)
    endif
    ! Set safe defaults on all ranks before master_task reads the file.
    ! These values are overwritten after successful metadata read and broadcast.
    ! If the read fails before broadcast, these defaults avoid undefined state.
    tide_nconst = 0
    tide_nj     = 0
    tide_ni     = 0
    tide_constituent_order = 'unknown'
    if (my_task == master_task) then
       ! File-path diagnostics.
       ! These are deliberately verbose because failed NetCDF opens on Gadi
       ! are often caused by trailing blanks, unexpected default strings, or
       ! namelist path edits not propagating as expected.
       write(nu_diag,*) subname//' pre-open len(tide_data_file)      = ', len(tide_data_file)
       write(nu_diag,*) subname//' pre-open len_trim(tide_data_file) = ', len_trim(tide_data_file)
       write(nu_diag,*) subname//' pre-open raw tide_data_file       = >', tide_data_file, '<'
       write(nu_diag,*) subname//' pre-open trim(tide_data_file)     = >', trim(tide_data_file), '<'
       ! Read dimensions from the pre-regridded CICE_TMD3 file.
       ! Expected file layout:
       !   constituents = number of harmonic constituents
       !   nj           = CICE ny_global
       !   ni           = CICE nx_global
       call ice_open_nc(trim(tide_data_file), ncid)
       ierr = nf90_inq_dimid(ncid, 'constituents', dimid)
       if (ierr /= nf90_noerr) then
          call abort_ice(error_message=subname//' ERROR missing dim constituents: '//trim(nf90_strerror(ierr)), &
               file=__FILE__, line=__LINE__)
       endif
       ierr = nf90_inquire_dimension(ncid, dimid, len=tide_nconst)
       if (ierr /= nf90_noerr) then
          call abort_ice(error_message=subname//' ERROR reading dim constituents: '//trim(nf90_strerror(ierr)), &
               file=__FILE__, line=__LINE__)
       endif
       ierr = nf90_inq_dimid(ncid, 'nj', dimid)
       if (ierr /= nf90_noerr) then
          call abort_ice(error_message=subname//' ERROR missing dim nj: '//trim(nf90_strerror(ierr)), &
               file=__FILE__, line=__LINE__)
       endif
       ierr = nf90_inquire_dimension(ncid, dimid, len=tide_nj)
       if (ierr /= nf90_noerr) then
          call abort_ice(error_message=subname//' ERROR reading dim nj: '//trim(nf90_strerror(ierr)), &
               file=__FILE__, line=__LINE__)
       endif
       ierr = nf90_inq_dimid(ncid, 'ni', dimid)
       if (ierr /= nf90_noerr) then
          call abort_ice(error_message=subname//' ERROR missing dim ni: '//trim(nf90_strerror(ierr)), &
               file=__FILE__, line=__LINE__)
       endif
       ierr = nf90_inquire_dimension(ncid, dimid, len=tide_ni)
       if (ierr /= nf90_noerr) then
          call abort_ice(error_message=subname//' ERROR reading dim ni: '//trim(nf90_strerror(ierr)), &
               file=__FILE__, line=__LINE__)
       endif
    endif
    call broadcast_scalar(tide_nconst, master_task)
    call broadcast_scalar(tide_nj,     master_task)
    call broadcast_scalar(tide_ni,     master_task)
    ! The tide file must already be on the active CICE T-grid.
    ! A mismatch here means the wrong grid file was supplied, or the model
    ! is being run with a different domain/grid than the tide preprocessing
    ! assumed. Abort rather than silently applying spatially inconsistent
    ! forcing.
    if (tide_nj /= ny_global .or. tide_ni /= nx_global) then
       if (my_task == master_task) then
          write(nu_diag,*) subname//' ERROR tide file dims do not match CICE grid'
          write(nu_diag,*) subname//' tide_nj, tide_ni = ', tide_nj, tide_ni
          write(nu_diag,*) subname//' ny_global, nx_global = ', ny_global, nx_global
       endif
       call abort_ice(error_message=subname//' tide grid mismatch', &
            file=__FILE__, line=__LINE__)
    endif
    if (allocated(tide_omega))     deallocate(tide_omega)
    if (allocated(tide_phase))     deallocate(tide_phase)
    if (allocated(tide_alpha))     deallocate(tide_alpha)
    if (allocated(tide_amplitude)) deallocate(tide_amplitude)
    ! Allocate one-dimensional constituent metadata arrays.
    ! These are small and replicated on every rank after broadcast.
    !
    ! tide_omega     : angular frequency [rad s-1]
    ! tide_phase     : phase offset used in arg = omega*t + phase [rad]
    ! tide_alpha     : nodal/equilibrium correction metadata, if supplied
    ! tide_amplitude : constituent amplitude metadata, mainly diagnostic here
    allocate(tide_omega(tide_nconst))
    allocate(tide_phase(tide_nconst))
    allocate(tide_alpha(tide_nconst))
    allocate(tide_amplitude(tide_nconst))
    tide_omega     = c0
    tide_phase     = c0
    tide_alpha     = c0
    tide_amplitude = c0
    if (my_task == master_task) then
       ierr = nf90_inq_varid(ncid, 'omega', varid)
       if (ierr /= nf90_noerr) then
          call abort_ice(error_message=subname//' ERROR missing var omega: '//trim(nf90_strerror(ierr)), &
               file=__FILE__, line=__LINE__)
       endif
       ierr = nf90_get_var(ncid, varid, tide_omega)
       if (ierr /= nf90_noerr) then
          call abort_ice(error_message=subname//' ERROR reading omega: '//trim(nf90_strerror(ierr)), &
               file=__FILE__, line=__LINE__)
       endif
       ierr = nf90_inq_varid(ncid, 'phase', varid)
       if (ierr /= nf90_noerr) then
          call abort_ice(error_message=subname//' ERROR missing var phase: '//trim(nf90_strerror(ierr)), &
               file=__FILE__, line=__LINE__)
       endif
       ierr = nf90_get_var(ncid, varid, tide_phase)
       if (ierr /= nf90_noerr) then
          call abort_ice(error_message=subname//' ERROR reading phase: '//trim(nf90_strerror(ierr)), &
               file=__FILE__, line=__LINE__)
       endif
       ierr = nf90_inq_varid(ncid, 'alpha', varid)
       if (ierr /= nf90_noerr) then
          call abort_ice(error_message=subname//' ERROR missing var alpha: '//trim(nf90_strerror(ierr)), &
               file=__FILE__, line=__LINE__)
       endif
       ierr = nf90_get_var(ncid, varid, tide_alpha)
       if (ierr /= nf90_noerr) then
          call abort_ice(error_message=subname//' ERROR reading alpha: '//trim(nf90_strerror(ierr)), &
               file=__FILE__, line=__LINE__)
       endif
       ierr = nf90_inq_varid(ncid, 'amplitude', varid)
       if (ierr /= nf90_noerr) then
          call abort_ice(error_message=subname//' ERROR missing var amplitude: '//trim(nf90_strerror(ierr)), &
               file=__FILE__, line=__LINE__)
       endif
       ierr = nf90_get_var(ncid, varid, tide_amplitude)
       if (ierr /= nf90_noerr) then
          call abort_ice(error_message=subname//' ERROR reading amplitude: '//trim(nf90_strerror(ierr)), &
               file=__FILE__, line=__LINE__)
       endif
       ierr = nf90_get_att(ncid, nf90_global, 'constituent_order', tide_constituent_order)
       if (ierr /= nf90_noerr) tide_constituent_order = 'unknown'
       call ice_close_nc(ncid)
    endif
    do n = 1, tide_nconst
       call broadcast_scalar(tide_omega(n),     master_task)
       call broadcast_scalar(tide_phase(n),     master_task)
       call broadcast_scalar(tide_alpha(n),     master_task)
       call broadcast_scalar(tide_amplitude(n), master_task)
    enddo
    call broadcast_scalar(tide_constituent_order, master_task)
    tide_metadata_loaded = .true.
    ! Echo the tide configuration to ice_diag.d.
    ! This makes each experiment self-documenting, which is critical when
    ! comparing tide_curr_fac, current caps, depth tapers, and ramp settings.
    !
    ! These values come from ice_init.F90 / forcing_nml, but are reported here
    ! because this is where the tide forcing pathway is actually initialised.
    if (my_task == master_task) then
       write(nu_diag,*) subname//' loaded tide metadata from: ', trim(tide_data_file)
       write(nu_diag,*) subname//' tide namelist settings:'
       write(nu_diag,*) subname//'   tide_data_type             = ', trim(tide_data_type)
       write(nu_diag,*) subname//'   tide_data_format           = ', trim(tide_data_format)
       write(nu_diag,*) subname//'   tide_data_file             = ', trim(tide_data_file)
       write(nu_diag,*) subname//'   tide_use_currents          = ', tide_use_currents
       write(nu_diag,*) subname//'   tide_use_ssh               = ', tide_use_ssh
       write(nu_diag,*) subname//'   tide_use_bathymetry_limit  = ', tide_use_bathymetry_limit
       write(nu_diag,*) subname//'   tide_curr_fac              = ', tide_curr_fac
       write(nu_diag,*) subname//'   tide_speed_cap [m/s]       = ', tide_speed_cap
       write(nu_diag,*) subname//'   tide_wct_min [m]           = ', tide_wct_min
       write(nu_diag,*) subname//'   tide_wct_full [m]          = ', tide_wct_full
       write(nu_diag,*) subname//'   tide_h_eff_min [m]         = ', tide_h_eff_min
       write(nu_diag,*) subname//'   tide_depth_mismatch        = ', tide_depth_mismatch
       write(nu_diag,*) subname//'   tide_ramp_days [days]      = ', tide_ramp_days
       if (tide_use_currents) then
          write(nu_diag,*) subname//'   tidal currents are applied to ocean forcing'
       else
          write(nu_diag,*) subname//'   tidal currents are NOT applied'
       endif
       if (tide_use_ssh) then
          write(nu_diag,*) subname//'   tidal SSH forcing is enabled'
       else
          write(nu_diag,*) subname//'   tidal SSH forcing is disabled'
       endif
       if (tide_use_bathymetry_limit) then
          write(nu_diag,*) subname//'   bathymetry-aware current limiter is enabled'
       else
          write(nu_diag,*) subname//'   bathymetry-aware current limiter is disabled'
       endif
       write(nu_diag,*) subname//' tide metadata:'
       write(nu_diag,*) subname//'   tide_nconst               = ', tide_nconst
       write(nu_diag,*) subname//'   tide grid                 = ', tide_nj, tide_ni
       write(nu_diag,*) subname//'   constituent_order         = ', trim(tide_constituent_order)
       if (tide_nconst > 0) then
          write(nu_diag,*) subname//'   omega(1) [rad/s]         = ', tide_omega(1)
          write(nu_diag,*) subname//'   phase(1) [rad]           = ', tide_phase(1)
          write(nu_diag,*) subname//'   alpha(1)                 = ', tide_alpha(1)
          write(nu_diag,*) subname//'   amplitude(1)             = ', tide_amplitude(1)
       endif
    endif
#endif
  end subroutine init_tides_metadata
  !=======================================================================
  !=======================================================================
  ! Allocate distributed tidal coefficient and limiter fields.
  !
  ! Spatial harmonic fields are stored on CICE distributed blocks:
  !
  !   tide_hRe/tide_hIm : complex SSH harmonic coefficients
  !   tide_URe/tide_UIm : complex zonal/barotropic transport coefficients
  !   tide_VRe/tide_VIm : complex meridional/barotropic transport coefficients
  !   tide_wct          : water-column thickness from the tide product
  !   tide_mask         : valid-ocean mask from the tide product
  !   tide_depth_fac    : 0..1 reliability taper used before applying tides
  !   tide_h_eff        : effective denominator used for transport/current
  !                       conversion
  !
  ! Notes:
  !   * These arrays are allocated once per run after metadata are known.
  !   * They are zero-initialised so failed/partial reads do not leave stale
  !     values in memory.
  !   * The harmonic coefficient arrays retain a constituent dimension;
  !     compute_tides_at_time() reconstructs instantaneous tides by summing
  !     over that dimension.
  !=======================================================================
  subroutine alloc_tides_fields
    integer (kind=int_kind) :: ierr
    character(len=*), parameter :: subname = '(alloc_tides_fields)'
    if (trim(tide_data_type) /= 'harmonic') return
    if (.not. tide_metadata_loaded) then
       call abort_ice(error_message=subname//' ERROR tide metadata not loaded', &
            file=__FILE__, line=__LINE__)
    endif
    if (tide_fields_loaded) return
    if (.not. allocated(tide_hRe)) then
       allocate(tide_hRe      (nx_block,ny_block,max_blocks,tide_nconst), &
                tide_hIm      (nx_block,ny_block,max_blocks,tide_nconst), &
                tide_URe      (nx_block,ny_block,max_blocks,tide_nconst), &
                tide_UIm      (nx_block,ny_block,max_blocks,tide_nconst), &
                tide_VRe      (nx_block,ny_block,max_blocks,tide_nconst), &
                tide_VIm      (nx_block,ny_block,max_blocks,tide_nconst), &
                tide_wct      (nx_block,ny_block,max_blocks) , &
                tide_mask     (nx_block,ny_block,max_blocks) , &
                tide_depth_fac(nx_block,ny_block,max_blocks) , &
                tide_h_eff    (nx_block,ny_block,max_blocks) , &
                tide_u_raw         (nx_block,ny_block,max_blocks) , &
                tide_v_raw         (nx_block,ny_block,max_blocks) , &
                tide_speed_raw     (nx_block,ny_block,max_blocks) , &
                tide_u_eff         (nx_block,ny_block,max_blocks) , &
                tide_v_eff         (nx_block,ny_block,max_blocks) , &
                tide_speed_eff     (nx_block,ny_block,max_blocks) , &
                tide_speed_daymax  (nx_block,ny_block,max_blocks) , &
                tide_speed_dayrms  (nx_block,ny_block,max_blocks) , &
                tide_speed_daysum  (nx_block,ny_block,max_blocks) , &
                tide_speed2_daysum (nx_block,ny_block,max_blocks) , &
                tide_n_subday      (nx_block,ny_block,max_blocks) , &
                tide_n_curr_over_fi(nx_block,ny_block,max_blocks) , &
                tide_n_capped      (nx_block,ny_block,max_blocks) , &
                 stat = ierr )
       if (ierr /= 0) then
          call abort_ice(error_message=subname//' ERROR allocating tide fields', &
               file=__FILE__, line=__LINE__)
       endif
    endif
    tide_hRe            = c0
    tide_hIm            = c0
    tide_URe            = c0
    tide_UIm            = c0
    tide_VRe            = c0
    tide_VIm            = c0
    tide_wct            = c0
    tide_mask           = c0
    tide_depth_fac      = c0
    tide_h_eff          = c0
    tide_u_raw          = c0
    tide_v_raw          = c0
    tide_speed_raw      = c0
    tide_u_eff          = c0
    tide_v_eff          = c0
    tide_speed_eff      = c0
    tide_speed_daymax   = c0
    tide_speed_dayrms   = c0
    tide_speed_daysum   = c0
    tide_speed2_daysum  = c0
    tide_n_subday       = c0
    tide_n_curr_over_fi = c0
    tide_n_capped       = c0
    tide_diag_myear     = -9999
    tide_diag_mmonth    = -9999
    tide_diag_mday      = -9999
  end subroutine alloc_tides_fields
  !=======================================================================
  !=======================================================================
  ! Read and scatter harmonic tide coefficient fields.
  !
  ! This routine reads the full CICE_TMD3 tide coefficient file and scatters
  ! global arrays onto the CICE decomposition. It also precomputes the
  ! bathymetry-aware reliability limiter used later when converting tidal
  ! transports to currents.
  !
  ! Input file expectations:
  !   static 2D fields:
  !      wct       : CATS/TMD water-column thickness [m]
  !      cats_mask : valid tide-ocean mask
  !
  !   harmonic 3D fields:
  !      hRe/hIm : real/imaginary SSH coefficients
  !      URe/UIm : real/imaginary transport coefficients
  !      VRe/VIm : real/imaginary transport coefficients
  !
  ! Important distinction:
  !   init_tides_fields() does not compute tides at a model time. It only
  !   loads the coefficient fields required for later harmonic prediction.
  !
  ! Physical interpretation:
  !   The transport coefficients are retained as the primary tidal product.
  !   The bathymetry/depth limiter is used later to avoid converting reliable
  !   transports into unphysical currents in unresolved shallow or mismatched
  !   coastal cells.
  !=======================================================================
  subroutine init_tides_fields
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use ice_gather_scatter, only: scatter_global
  use ice_domain, only: distrb_info
  use ice_grid, only: bathymetry, tmask
  use ice_global_reductions, only: global_minval, global_maxval
#ifdef USE_NETCDF
  use netcdf, only: nf90_noerr, nf90_inq_varid, nf90_get_var, nf90_strerror
#endif
  integer (kind=int_kind) :: fid
  real (kind=dbl_kind) :: bathy_min_global, bathy_max_global
  character(len=*), parameter :: subname = '(init_tides_fields)'
  if (trim(tide_data_type) /= 'harmonic') return
  if (.not. tide_metadata_loaded) then
     call abort_ice(error_message=subname//' ERROR tide metadata not loaded', &
          file=__FILE__, line=__LINE__)
  endif
  if (.not. allocated(tide_hRe)) then
     call abort_ice(error_message=subname//' ERROR tide fields not allocated', &
          file=__FILE__, line=__LINE__)
  endif
  if (tide_fields_loaded) return
#ifndef USE_NETCDF
  call abort_ice(error_message=subname//' ERROR direct tide reads require USE_NETCDF', &
       file=__FILE__, line=__LINE__)
#else
  if (my_task == master_task) then
     write(nu_diag,*) subname//' len_trim(tide_data_file) = ', len_trim(tide_data_file)
     write(nu_diag,*) subname//' tide_data_file = ', trim(tide_data_file)
  endif
  fid = -1
  if (my_task == master_task) then
     call ice_open_nc(trim(tide_data_file), fid)
  endif
  ! Static 2D fields from prototype file.
  !
  ! wct is the water-column thickness associated with the tide product.
  ! It is the internally consistent denominator for converting CATS/TMD
  ! barotropic transports to approximate depth-mean currents.
  !
  ! cats_mask marks cells where the tide product regards the coefficient
  ! fields as valid ocean. It is combined later with the CICE tmask and,
  ! optionally, the CICE bathymetry field.
  call read_static_2d(fid, 'wct',       tide_wct)
  call read_static_2d(fid, 'cats_mask', tide_mask)
  ! Precompute depth/reliability factors once.
  ! This avoids repeating bathymetry/wct/mask checks every timestep and
  ! makes the current-limiting logic deterministic for the whole run.
  call init_tide_depth_limiter
  ! Harmonic 3D fields from prototype file.
  !
  ! The file stores complex harmonic coefficients as separate real and
  ! imaginary arrays. At runtime, compute_tides_at_time() reconstructs the
  ! instantaneous tide using:
  !
  !   X(t) = XRe*cos(omega*t + phase) - XIm*sin(omega*t + phase)
  !
  ! for each constituent and then sums over constituents.
  call read_harmonic_3d(fid, 'hRe', tide_hRe)
  call read_harmonic_3d(fid, 'hIm', tide_hIm)
  call read_harmonic_3d(fid, 'URe', tide_URe)
  call read_harmonic_3d(fid, 'UIm', tide_UIm)
  call read_harmonic_3d(fid, 'VRe', tide_VRe)
  call read_harmonic_3d(fid, 'VIm', tide_VIm)
  if (my_task == master_task) then
     call ice_close_nc(fid)
  endif
  tide_fields_loaded = .true.
  ! These min/max values are local to master_task's owned blocks after
  ! scattering, not global diagnostics. They are useful as smoke tests
  ! that the fields are nonzero and finite, but should not be used as
  ! domain-wide extrema.
  if (my_task == master_task) then
     write(nu_diag,*) subname//' loaded tide coefficient fields from: ', trim(tide_data_file)
     write(nu_diag,*) subname//' master raw/scattered diagnostics complete'
     write(nu_diag,*) subname//' tide_wct       local-block min/max = ', minval(tide_wct),       maxval(tide_wct)
     write(nu_diag,*) subname//' tide_mask      local-block min/max = ', minval(tide_mask),      maxval(tide_mask)
     write(nu_diag,*) subname//' tide_depth_fac local-block min/max = ', minval(tide_depth_fac), maxval(tide_depth_fac)
     write(nu_diag,*) subname//' tide_h_eff     local-block min/max = ', minval(tide_h_eff),     maxval(tide_h_eff)
     write(nu_diag,*) subname//' tide_hRe(1) local-block min/max = ', minval(tide_hRe(:,:,:,1)), maxval(tide_hRe(:,:,:,1))
     write(nu_diag,*) subname//' tide_URe(1) local-block min/max = ', minval(tide_URe(:,:,:,1)), maxval(tide_URe(:,:,:,1))
     write(nu_diag,*) subname//' tide_VRe(1) local-block min/max = ', minval(tide_VRe(:,:,:,1)), maxval(tide_VRe(:,:,:,1))
  endif
#endif
contains
  !---------------------------------------------------------------------
  ! Precompute a bathymetry-aware reliability limiter for tidal currents.
  !
  ! The tide model provides barotropic transports. To apply these as CICE
  ! ocean-current perturbations, the transports must be converted to
  ! approximate currents using a water-column thickness.
  !
  ! Direct U/h conversion can produce unrealistically large currents where:
  !   - the tide water-column thickness is very small,
  !   - the regridded CATS/TMD coast does not exactly match the CICE coast,
  !   - CICE bathymetry and CATS/TMD wct disagree strongly,
  !   - narrow coastal/shelf features are unresolved at 1/4 degree.
  !
  ! This limiter addresses those issues by computing:
  !
  !   tide_depth_fac : 0..1 smooth taper applied to tide amplitude
  !   tide_h_eff     : effective minimum depth used for U/h conversion
  !
  ! The emergency speed cap in compute_tides_at_time() remains as a final
  ! guardrail, but should only catch rare pathological remnants after this
  ! smoother bathymetry/wct treatment.
  subroutine init_tide_depth_limiter
    integer (kind=int_kind) :: i, j, iblk
    integer (kind=int_kind) :: n_valid, n_tapered, n_zeroed, n_mismatch
    real (kind=dbl_kind) :: h_cats, h_cice, h_reliable
    real (kind=dbl_kind) :: xdep, depth_fac, mismatch_ratio
    real (kind=dbl_kind), parameter :: eps_depth = 1.0e-12_dbl_kind
    ! Check whether CICE bathymetry is available on ocean cells.
    ! If tide_use_bathymetry_limit is true, bathymetry is part of the
    ! physical limiter and a missing/zero field should be treated as a
    ! configuration error rather than silently ignored.
    bathy_min_global = global_minval(bathymetry, distrb_info, tmask)
    bathy_max_global = global_maxval(bathymetry, distrb_info, tmask)
    if (tide_use_bathymetry_limit .and. bathy_max_global <= c0) then
       if (my_task == master_task) then
          write(nu_diag,*) subname//' ERROR tide_use_bathymetry_limit is true but bathymetry appears unavailable'
          write(nu_diag,*) subname//' bathymetry global min/max = ', bathy_min_global, bathy_max_global
       endif
       call abort_ice(error_message=subname//' bathymetry unavailable for tide limiter', &
            file=__FILE__, line=__LINE__)
    endif
    if (tide_wct_full <= tide_wct_min) then
       call abort_ice(error_message=subname//' ERROR tide_wct_full must be greater than tide_wct_min', &
            file=__FILE__, line=__LINE__)
    endif
    if (tide_h_eff_min <= c0) then
       call abort_ice(error_message=subname//' ERROR tide_h_eff_min must be positive', &
            file=__FILE__, line=__LINE__)
    endif
    tide_depth_fac = c0
    tide_h_eff     = c0
    ! integer counts
    n_valid        = 0
    n_tapered      = 0
    n_zeroed       = 0
    n_mismatch     = 0
    do iblk = 1, max_blocks
       do j = 1, ny_block
          do i = 1, nx_block
             h_cats = max(c0, tide_wct(i,j,iblk))
             if (tide_use_bathymetry_limit) then
                h_cice     = max(c0, bathymetry(i,j,iblk))
                ! Conservative reliability depth.
                ! h_cats comes from the tide product and is internally
                ! consistent with the CATS/TMD transport coefficients.
                ! h_cice comes from the active CICE grid/bathymetry and
                ! reflects the model cell in which the current will be applied.
                ! Taking min(h_cats,h_cice) means either product can flag a
                ! cell as shallow/unreliable. This is intentionally conservative
                ! near the Antarctic coast, where grid and mask mismatches are
                ! most likely.
                h_reliable = min(h_cats, h_cice)
             else
                h_cice = h_cats
                h_reliable = h_cats
             endif
             ! A cell receives tidal-current forcing only if:
             !   1. it is ocean in the CICE grid       : tmask
             !   2. it is valid in the tide product    : tide_mask > 0.5
             !   3. the conservative reliable depth is deeper than the
             !      minimum allowed depth              : h_reliable > tide_wct_min
             ! Otherwise all tidal forcing is zeroed in that cell.
             if (tmask(i,j,iblk) .and. tide_mask(i,j,iblk) > p5 .and. h_reliable > tide_wct_min) then
                ! Non-dimensional depth coordinate for the taper:
                !   xdep = 0 at tide_wct_min
                !   xdep = 1 at tide_wct_full
                ! Cells shallower than tide_wct_min receive no tide-current
                ! forcing. Cells deeper than tide_wct_full receive full
                ! forcing, subject to tide_curr_fac and the emergency cap.
                xdep = (h_reliable - tide_wct_min) / max(eps_depth, tide_wct_full - tide_wct_min)
                xdep = max(c0, min(c1, xdep))
                ! Smooth 0..1 taper using "smootherstep":
                !   depth_fac = x^3 * (10 - 15x + 6x^2)
                ! This gives zero slope at both endpoints, avoiding abrupt
                ! spatial jumps in tidal current amplitude near the shallow
                ! cutoff. That is preferable to a hard on/off mask at the
                ! continental shelf/coastal boundary.
                depth_fac                = xdep**3 * (10.0_dbl_kind - 15.0_dbl_kind*xdep + 6.0_dbl_kind*xdep**2)
                tide_depth_fac(i,j,iblk) = depth_fac
                ! Effective depth used in transport-to-current conversion.
                ! We keep h_cats as the primary denominator because the
                ! transport coefficients and wct come from the same tide
                ! product. However, we impose tide_h_eff_min so that a small
                ! or noisy wct does not create extreme currents through U/h.
                ! Note: tide_depth_fac handles reliability/tapering;
                ! tide_h_eff handles the numerical magnitude of U/h.
                tide_h_eff(i,j,iblk) = max(h_cats, tide_h_eff_min)
                n_valid = n_valid + 1
                if (depth_fac < c1) n_tapered = n_tapered + 1
                ! Diagnostic only: count cells where CATS/TMD wct and CICE
                ! bathymetry differ by more than tide_depth_mismatch.
                ! These cells are not automatically zeroed by this test;
                ! the mismatch count is intended to help interpret later
                ! diagnostics and identify problematic coastal regions.
                if (tide_use_bathymetry_limit .and. h_cats > c0 .and. h_cice > c0) then
                   mismatch_ratio = max(h_cats, h_cice) / max(eps_depth, min(h_cats, h_cice))
                   if (mismatch_ratio > tide_depth_mismatch) n_mismatch = n_mismatch + 1
                endif
             else
                tide_depth_fac(i,j,iblk) = c0
                tide_h_eff(i,j,iblk)     = c0
                n_zeroed = n_zeroed + 1
             endif
          enddo
       enddo
    enddo
    ! Report limiter settings and local counts:
    ! Counts are local to master_task's block ownership unless explicitly
    ! reduced elsewhere. Treat them as useful run fingerprints, not global
    ! cell totals. For publication-quality diagnostics, add global reductions.
    if (my_task == master_task) then
       write(nu_diag,*) subname//' tide depth limiter settings:'
       write(nu_diag,*) subname//'   tide_use_bathymetry_limit = ', tide_use_bathymetry_limit
       write(nu_diag,*) subname//'   tide_wct_min/full         = ', tide_wct_min, tide_wct_full
       write(nu_diag,*) subname//'   tide_h_eff_min            = ', tide_h_eff_min
       write(nu_diag,*) subname//'   tide_depth_mismatch       = ', tide_depth_mismatch
       write(nu_diag,*) subname//'   bathymetry global min/max  = ', bathy_min_global, bathy_max_global
       write(nu_diag,*) subname//' tide depth limiter local counts valid/tapered/zeroed/mismatch = ', &
                         n_valid, n_tapered, n_zeroed, n_mismatch
    endif

  end subroutine init_tide_depth_limiter
  !---------------------------------------------------------------------
  ! Replace NaN/Inf values read from NetCDF with zero.
  !
  ! This is a defensive read-time sanitation step. It prevents isolated
  ! missing or invalid coefficient values from poisoning the harmonic sum
  ! and producing NaN currents/stresses later in the dynamics.
  !
  ! A zero here means "do not apply tidal contribution from this invalid
  ! coefficient value"; it does not fix upstream data quality issues.
  subroutine sanitize2d(a)
    real (kind=dbl_kind), intent(inout) :: a(:,:)
    where (.not. ieee_is_finite(a)) a = c0
  end subroutine sanitize2d
  !---------------------------------------------------------------------
  subroutine sanitize3d(a)
    real (kind=dbl_kind), intent(inout) :: a(:,:,:)
    where (.not. ieee_is_finite(a)) a = c0
  end subroutine sanitize3d
  !---------------------------------------------------------------------
  ! Read a static 2D tide field on master_task and scatter to CICE blocks.
  !
  ! File/global convention:
  !   The preprocessed CICE_TMD3 file is expected to be on the CICE T-grid
  !   with dimensions compatible with nx_global, ny_global.
  !
  ! Parallel convention:
  !   master_task reads the global field, sanitises it, then scatter_global()
  !   distributes it to the local block decomposition used by CICE.
  !
  ! These fields have no constituent dimension.
  subroutine read_static_2d(fid, varname, dest)
    integer (kind=int_kind), intent(in) :: fid
    character(len=*), intent(in) :: varname
    real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks), intent(out) :: dest
    integer (kind=int_kind) :: ierr, varid
    real (kind=dbl_kind), allocatable :: g2(:,:)   ! Fortran/global order: (ni,nj)
    ! Only master_task needs a full global work array.
    ! Non-master ranks allocate a tiny dummy array so the scatter interface
    ! can be called uniformly on all tasks.
    if (my_task == master_task) then
       allocate(g2(nx_global,ny_global))
    else
       allocate(g2(1,1))
    endif
    g2 = c0
    if (my_task == master_task) then
       ierr = nf90_inq_varid(fid, trim(varname), varid)
       if (ierr /= nf90_noerr) then
          call abort_ice(error_message=subname//' ERROR missing var '//trim(varname)//': '//trim(nf90_strerror(ierr)), &
               file=__FILE__, line=__LINE__)
       endif
       ierr = nf90_get_var(fid, varid, g2)
       if (ierr /= nf90_noerr) then
          call abort_ice(error_message=subname//' ERROR reading '//trim(varname)//': '//trim(nf90_strerror(ierr)), &
               file=__FILE__, line=__LINE__)
       endif
       call sanitize2d(g2)
       write(nu_diag,*) subname//' raw ', trim(varname), ' min/max = ', minval(g2), maxval(g2)
    endif
    ! Scatter from global Fortran-order array to local CICE blocks.
    ! These tide fields are T-grid scalar fields, so use field_loc_center
    ! and field_type_scalar.
    call scatter_global(dest, g2, master_task, distrb_info, &
         field_loc_center, field_type_scalar)
    deallocate(g2)
  end subroutine read_static_2d
  !---------------------------------------------------------------------
  ! Read a harmonic 3D tide field constituent-by-constituent and scatter.
  !
  ! The NetCDF variable is conceptually:
  !
  !   var(constituents, nj, ni)
  !
  ! but the Fortran NetCDF interface presents the storage order such that
  ! the slab read below uses:
  !
  !   start = (/1, 1, n/)
  !   count = (/nx_global, ny_global, 1/)
  !
  ! Each constituent is read into a 2D global work array, sanitised, and
  ! scattered into dest(:,:,:,n).
  !
  ! This avoids holding an additional full global 3D field in memory.
  subroutine read_harmonic_3d(fid, varname, dest)
    integer (kind=int_kind), intent(in) :: fid
    character(len=*), intent(in) :: varname
    real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks,tide_nconst), intent(out) :: dest
    integer (kind=int_kind) :: ierr, varid, n
    real (kind=dbl_kind), allocatable :: g2(:,:)   ! Fortran/global order: (ni,nj)
    if (my_task == master_task) then
       allocate(g2(nx_global,ny_global))
    else
       allocate(g2(1,1))
    endif
    g2 = c0
    if (my_task == master_task) then
       ierr = nf90_inq_varid(fid, trim(varname), varid)
       if (ierr /= nf90_noerr) then
          call abort_ice(error_message=subname//' ERROR missing var '//trim(varname)//': '//trim(nf90_strerror(ierr)), &
               file=__FILE__, line=__LINE__)
       endif
    endif
    ! Print only constituent 1 as a compact smoke test.
    ! Full constituent-by-constituent diagnostics would be verbose
    ! and expensive in ice_diag.d.
    do n = 1, tide_nconst
       g2 = c0
       if (my_task == master_task) then
          ! File variable is hRe(constituents,nj,ni)
          ! Fortran netCDF view is effectively (ni,nj,constituents)
          ierr = nf90_get_var(fid, varid, g2, &
                              start=(/1,1,n/), &
                              count=(/nx_global,ny_global,1/))
          if (ierr /= nf90_noerr) then
             write(nu_diag,*) subname//' ERROR reading ', trim(varname), ' constituent n = ', n
             write(nu_diag,*) subname//' attempted start = (/ 1, 1, ', n, ' /)'
             write(nu_diag,*) subname//' attempted count = (/ ', nx_global, ', ', ny_global, ', 1 /)'
             call abort_ice(error_message=subname//' ERROR reading '//trim(varname)//': '//trim(nf90_strerror(ierr)), &
                  file=__FILE__, line=__LINE__)
          endif
          call sanitize2d(g2)
          if (n == 1) then
             write(nu_diag,*) subname//' raw ', trim(varname), '(1) min/max = ', minval(g2), maxval(g2)
          endif
       endif
       call scatter_global(dest(:,:,:,n), g2, master_task, distrb_info, &
            field_loc_center, field_type_scalar)
    enddo
    deallocate(g2)
  end subroutine read_harmonic_3d
  end subroutine init_tides_fields
  !=======================================================================
  subroutine reset_tide_current_diagnostics_if_needed

    character(len=*), parameter :: subname = '(reset_tide_current_diagnostics_if_needed)'

    if (.not. allocated(tide_speed_eff)) return

    if (tide_diag_myear  /= myear  .or. &
        tide_diag_mmonth /= mmonth .or. &
        tide_diag_mday   /= mday) then

       tide_speed_daymax   = c0
       tide_speed_dayrms   = c0
       tide_speed_daysum   = c0
       tide_speed2_daysum  = c0
       tide_n_subday       = c0
       tide_n_curr_over_fi = c0
       tide_n_capped       = c0

       tide_diag_myear  = myear
       tide_diag_mmonth = mmonth
       tide_diag_mday   = mday

       if (debug_forcing .and. my_task == master_task) then
          write(nu_diag,*) subname//' reset diagnostics for date = ', &
               tide_diag_myear, tide_diag_mmonth, tide_diag_mday
       endif
    endif

  end subroutine reset_tide_current_diagnostics_if_needed

  !=======================================================================
  subroutine update_tide_current_diagnostics(i, j, iblk, u_raw, v_raw, u_eff, v_eff)

    integer (kind=int_kind), intent(in) :: i, j, iblk
    real    (kind=dbl_kind), intent(in) :: u_raw, v_raw
    real    (kind=dbl_kind), intent(in) :: u_eff, v_eff

    real (kind=dbl_kind) :: speed_raw
    real (kind=dbl_kind) :: speed_eff

    if (.not. allocated(tide_speed_eff)) return

    speed_raw = sqrt(u_raw*u_raw + v_raw*v_raw)
    speed_eff = sqrt(u_eff*u_eff + v_eff*v_eff)

    tide_u_raw    (i,j,iblk) = u_raw
    tide_v_raw    (i,j,iblk) = v_raw
    tide_speed_raw(i,j,iblk) = speed_raw

    tide_u_eff    (i,j,iblk) = u_eff
    tide_v_eff    (i,j,iblk) = v_eff
    tide_speed_eff(i,j,iblk) = speed_eff

    tide_speed_daymax (i,j,iblk) = max(tide_speed_daymax(i,j,iblk), speed_eff)
    tide_speed_daysum (i,j,iblk) = tide_speed_daysum(i,j,iblk)  + speed_eff
    tide_speed2_daysum(i,j,iblk) = tide_speed2_daysum(i,j,iblk) + speed_eff*speed_eff
    tide_n_subday     (i,j,iblk) = tide_n_subday(i,j,iblk)      + c1

    tide_speed_dayrms(i,j,iblk) = sqrt( &
         tide_speed2_daysum(i,j,iblk) / max(tide_n_subday(i,j,iblk), c1) )

    if (speed_eff > tide_curr_diag_thresh) then
       tide_n_curr_over_fi(i,j,iblk) = tide_n_curr_over_fi(i,j,iblk) + c1
    endif

  end subroutine update_tide_current_diagnostics
  !=======================================================================
  ! Compute instantaneous tidal SSH, transports, and current perturbations.
  !
  ! This routine is called during forcing updates. It reconstructs the
  ! instantaneous harmonic tide at the current model time and converts
  ! barotropic transports into approximate depth-mean tidal currents.
  !
  ! Output fields:
  !   eta_tide  : tidal sea-surface height perturbation after limiter/ramp
  !   Utide_tr  : tidal U transport after limiter/ramp/cap
  !   Vtide_tr  : tidal V transport after limiter/ramp/cap
  !   utide_cur : tidal u-current perturbation applied to ocean forcing
  !   vtide_cur : tidal v-current perturbation applied to ocean forcing
  !
  ! Physical pathway:
  !   harmonic coefficients
  !        -> instantaneous tidal SSH/transports
  !        -> bathymetry/wct reliability limiter
  !        -> transport-to-current conversion
  !        -> branch-start amplitude ramp
  !        -> emergency current cap
  !        -> ocean-current perturbation seen by CICE dynamics
  !
  ! Important distinction:
  !   tsec_ref is the absolute harmonic phase clock relative to the CATS/TMD
  !   reference epoch. It must NOT be reset at restart.
  !
  !   tide_ramp_t0_sec is the branch/run start time used only for gradually
  !   increasing the forcing amplitude. It deliberately starts from the first
  !   call in this branch run.
  !=======================================================================
  subroutine compute_tides_at_time(eta_tide, Utide_tr, Vtide_tr, utide_cur, vtide_cur)
    integer (kind=int_kind) :: n, i, j, iblk
    integer (kind=int_kind) :: days_since_ref
    real    (kind=dbl_kind) :: tsec_ref
    real    (kind=dbl_kind) :: arg, carg, sarg
    real    (kind=dbl_kind) :: ramp
    real    (kind=dbl_kind) :: spd, sf
    integer (kind=int_kind) :: n_capped
    real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks), intent(out) :: &
         eta_tide, Utide_tr, Vtide_tr, utide_cur, vtide_cur
    logical(kind=log_kind), save :: tide_ramp_initialized = .false.
    real(kind=dbl_kind),    save :: tide_ramp_t0_sec = c0
    real(kind=dbl_kind)          :: tsec_run, x
    real(kind=dbl_kind), parameter :: pi_local = 3.1415926535897932384626433832795_dbl_kind
    character(len=*), parameter :: subname = '(compute_tides_at_time)'
    if (trim(tide_data_type) /= 'harmonic') then
       eta_tide  = c0
       Utide_tr  = c0
       Vtide_tr  = c0
       utide_cur = c0
       vtide_cur = c0
       return
    endif
    if (.not. tide_metadata_loaded) then
       call abort_ice(error_message=subname//' ERROR tide metadata not loaded', &
            file=__FILE__, line=__LINE__)
    endif
    if (.not. tide_fields_loaded) then
       call abort_ice(error_message=subname//' ERROR tide fields not loaded', &
            file=__FILE__, line=__LINE__)
    endif
    if (allocated(tide_speed_eff)) then
       if (tide_diag_myear  /= myear  .or. &
            tide_diag_mmonth /= mmonth .or. &
            tide_diag_mday   /= mday) then

          tide_speed_daymax   = c0
          tide_speed_dayrms   = c0
          tide_speed_daysum   = c0
          tide_speed2_daysum  = c0
          tide_n_subday       = c0
          tide_n_curr_over_fi = c0
          tide_n_capped       = c0

          tide_diag_myear  = myear
          tide_diag_mmonth = mmonth
          tide_diag_mday   = mday
       endif
    endif
    ! Absolute harmonic phase clock.
    ! CATS/TMD phase metadata are referenced to 1992-01-01 00:00:00.
    ! This clock is used only to evaluate cos(omega*t + phase) and must
    ! remain tied to the tidal reference epoch, not to the CICE restart date.
    days_since_ref = compute_days_between(1992, 1, 1, myear, mmonth, mday)
    tsec_ref       = real(days_since_ref, kind=dbl_kind) * 86400.0_dbl_kind + real(msec, kind=dbl_kind)
    eta_tide       = c0
    Utide_tr       = c0
    Vtide_tr       = c0
    utide_cur      = c0
    vtide_cur      = c0
    call reset_tide_current_diagnostics_if_needed
    ! Reconstruct instantaneous harmonic tide by summing all constituents.
    ! For each complex coefficient X = XRe + i XIm, the real-valued tide is:
    !   X(t) = XRe*cos(arg) - XIm*sin(arg)
    ! where arg = omega*t + phase.
    ! eta_tide is SSH-like. Utide_tr and Vtide_tr are barotropic transports,
    ! not currents yet.
    do n = 1, tide_nconst
       arg  = tide_omega(n) * tsec_ref + tide_phase(n)
       carg = cos(arg)
       sarg = sin(arg)
       eta_tide(:,:,:) = eta_tide(:,:,:) + tide_hRe(:,:,:,n) * carg - tide_hIm(:,:,:,n) * sarg
       Utide_tr(:,:,:) = Utide_tr(:,:,:) + tide_URe(:,:,:,n) * carg - tide_UIm(:,:,:,n) * sarg
       Vtide_tr(:,:,:) = Vtide_tr(:,:,:) + tide_VRe(:,:,:,n) * carg - tide_VIm(:,:,:,n) * sarg
    enddo
    ! Convert transports to approximate tidal currents.
    !
    ! The tide product supplies barotropic transports. For CICE ocean forcing
    ! we need current perturbations, so we divide by tide_h_eff.
    !
    !   utide_cur = tide_curr_fac * tide_depth_fac * Utide_tr / tide_h_eff
    !   vtide_cur = tide_curr_fac * tide_depth_fac * Vtide_tr / tide_h_eff
    !
    ! tide_depth_fac:
    !   0..1 reliability taper precomputed from tide_wct, CICE bathymetry,
    !   CICE tmask, and cats_mask.
    !
    ! tide_h_eff:
    !   effective water-column thickness. It is based on CATS/TMD wct but
    !   bounded below by tide_h_eff_min to avoid U/h blow-ups.
    !
    ! Applying tide_depth_fac to eta_tide, Utide_tr, and Vtide_tr keeps the
    ! diagnostic output consistent with the actual applied current forcing.
    where (tide_depth_fac > c0 .and. tide_h_eff > c0)
       utide_cur = tide_curr_fac * tide_depth_fac * (Utide_tr / tide_h_eff)
       vtide_cur = tide_curr_fac * tide_depth_fac * (Vtide_tr / tide_h_eff)
       eta_tide  = tide_depth_fac * eta_tide
       Utide_tr  = tide_depth_fac * Utide_tr
       Vtide_tr  = tide_depth_fac * Vtide_tr
    elsewhere
       eta_tide  = c0
       Utide_tr  = c0
       Vtide_tr  = c0
       utide_cur = c0
       vtide_cur = c0
    end where
    ! Save converted current before branch-start ramp and emergency cap.
    tide_u_raw     = utide_cur
    tide_v_raw     = vtide_cur
    tide_speed_raw = sqrt(utide_cur*utide_cur + vtide_cur*vtide_cur)
    ! Smooth amplitude ramp from the start of this branch run.
    ! Do not use tsec_ref directly for the ramp, because tsec_ref is measured
    ! from the 1992 tide reference epoch. If used directly, any 1993/1994
    ! restart would begin with ramp=1 and the ramp would do nothing.
    ! The first call defines tide_ramp_t0_sec. The ramp clock is therefore:
    !   tsec_run = tsec_ref - tide_ramp_t0_sec
    ! The half-cosine ramp:
    !   ramp = 0.5 * (1 - cos(pi*x))
    ! has zero slope at the start and end, reducing the chance of shocking
    ! the momentum balance when tides are activated.
    if (.not. tide_ramp_initialized) then
       tide_ramp_t0_sec      = tsec_ref
       tide_ramp_initialized = .true.
    endif
    tsec_run = max(c0, tsec_ref - tide_ramp_t0_sec)
    if (tide_ramp_days > c0) then
       x = min(c1, tsec_run / (tide_ramp_days * 86400.0_dbl_kind))
       ! Half-cosine ramp: zero slope at start and end.
       ramp = p5 * (c1 - cos(pi_local * x))
    else
       ramp = c1
    endif
    eta_tide  = ramp * eta_tide
    Utide_tr  = ramp * Utide_tr
    Vtide_tr  = ramp * Vtide_tr
    utide_cur = ramp * utide_cur
    vtide_cur = ramp * vtide_cur
    ! Save pre-cap tidal current diagnostics.
    if (allocated(tide_speed_eff)) then
       tide_u_raw     = utide_cur
       tide_v_raw     = vtide_cur
       tide_speed_raw = sqrt(tide_u_raw*tide_u_raw + tide_v_raw*tide_v_raw)
    endif
    ! Emergency cap on tidal current vector speed.
    ! This is not intended to be the primary physical limiter. The main
    ! physical/numerical control is the bathymetry/wct taper above.
    ! The cap catches any remaining pathological cells after the smoother
    ! limiter and effective-depth conversion. If n_capped is large, the
    ! experiment is still being controlled by the cap and the limiter settings
    ! should be revisited.
    ! n_capped is local to this MPI task. The debug print on master_task is
    ! therefore a smoke test, not a global capped-cell count.
    n_capped = 0
    if (tide_speed_cap > c0) then
       do iblk = 1, max_blocks
          do j = 1, ny_block
             do i = 1, nx_block
                spd = sqrt(utide_cur(i,j,iblk)**2 + vtide_cur(i,j,iblk)**2)
                if (spd > tide_speed_cap .and. spd > c0) then
                   sf                  = tide_speed_cap / spd
                   utide_cur(i,j,iblk) = utide_cur(i,j,iblk) * sf
                   vtide_cur(i,j,iblk) = vtide_cur(i,j,iblk) * sf
                   Utide_tr(i,j,iblk)  = Utide_tr(i,j,iblk)  * sf
                   Vtide_tr(i,j,iblk)  = Vtide_tr(i,j,iblk)  * sf
                   n_capped            = n_capped + 1
                   if (allocated(tide_n_capped)) then
                      tide_n_capped(i,j,iblk) = tide_n_capped(i,j,iblk) + c1
                   endif
                endif
             enddo
          enddo
       enddo
    endif
    ! Save final applied tidal current diagnostics.
    if (allocated(tide_speed_eff)) then
       tide_u_eff     = utide_cur
       tide_v_eff     = vtide_cur
       tide_speed_eff = sqrt(tide_u_eff*tide_u_eff + tide_v_eff*tide_v_eff)
       tide_speed_daymax  = max(tide_speed_daymax, tide_speed_eff)
       tide_speed_daysum  = tide_speed_daysum  + tide_speed_eff
       tide_speed2_daysum = tide_speed2_daysum + tide_speed_eff*tide_speed_eff
       tide_n_subday      = tide_n_subday      + c1
       tide_speed_dayrms = sqrt(tide_speed2_daysum / max(tide_n_subday, c1))
       where (tide_speed_eff > tide_curr_diag_thresh)
          tide_n_curr_over_fi = tide_n_curr_over_fi + c1
       end where
    endif
    if (debug_forcing .and. my_task == master_task) then
       write(nu_diag,*) subname//' date/time         = ', myear, mmonth, mday, msec
       write(nu_diag,*) subname//' tsec_ref          = ', tsec_ref
       write(nu_diag,*) subname//' ramp              = ', ramp
       write(nu_diag,*) subname//' tide_curr_fac,    = ', tide_curr_fac
       write(nu_diag,*) subname//' tide_speed_cap    = ', tide_speed_cap
       write(nu_diag,*) subname//' n_capped local    = ', n_capped
       write(nu_diag,*) subname//' eta_tide  min/max = ', minval(eta_tide),  maxval(eta_tide)
       write(nu_diag,*) subname//' Utide_tr  min/max = ', minval(Utide_tr),  maxval(Utide_tr)
       write(nu_diag,*) subname//' Vtide_tr  min/max = ', minval(Vtide_tr),  maxval(Vtide_tr)
       write(nu_diag,*) subname//' utide_cur min/max = ', minval(utide_cur), maxval(utide_cur)
       write(nu_diag,*) subname//' vtide_cur min/max = ', minval(vtide_cur), maxval(vtide_cur)
    endif
  end subroutine compute_tides_at_time
  !=======================================================================

  !=======================================================================
  subroutine ocn_freezing_temperature
    ! Compute ocean freezing temperature Tf based on tfrz_option
    ! 'minus1p8'         Tf = -1.8 C (default)
    ! 'linear_salt'      Tf = -depressT * sss
    ! 'mushy'            Tf conforms with mushy layer thermo (ktherm=2)
    use ice_blocks, only: nx_block, ny_block
    use ice_domain, only: nblocks
    use ice_flux, only: sss, Tf
    ! local variables
    integer (kind=int_kind) :: i, j, iblk  ! horizontal indices
    character(len=*), parameter :: subname = '(ocn_freezing_temperature)'
    if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'
    !$OMP PARALLEL DO PRIVATE(iblk,i,j)
    do iblk = 1, nblocks
       do j = 1, ny_block
          do i = 1, nx_block
             Tf(i,j,iblk) = icepack_sea_freezing_temperature(sss(i,j,iblk))
          enddo
       enddo
    enddo !$OMP END PARALLEL DO
    call icepack_warnings_flush(nu_diag)
    if (icepack_warnings_aborted()) call abort_ice(error_message=subname, file=__FILE__, line=__LINE__)
  end subroutine ocn_freezing_temperature

  !=======================================================================
  subroutine get_forcing_atmo
    ! Get atmospheric forcing data and interpolate as necessary
    use ice_blocks, only: block, get_block
    use ice_domain, only: nblocks, blocks_ice
    use ice_flux, only: Tair, fsw, flw, frain, fsnow, Qa, rhoa, uatm, vatm, strax, stray, zlvl, wind, swvdr, swvdf, swidr, swidf, potT, sst
    use ice_state, only: aice, trcr
    use ice_grid, only: ANGLET, hm
    ! local variables
    integer (kind=int_kind) :: iblk, &            ! block index
                               ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
                               modadj, &          ! adjustment to make mod a postive number
                               fyear_old, &       ! fyear setting on last timestep
                               nt_Tsfc
    ! block information for current block
    type (block) :: this_block
    character(len=*), parameter :: subname = '(get_forcing_atmo)'
    if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'
    call ice_timer_start(timer_forcing)
    fyear_old = fyear
    modadj    = abs((min(0,myear-fyear_init)/ycycle+1)*ycycle)
    fyear     = fyear_init + mod(myear-fyear_init+modadj,ycycle)
    if (trim(atm_data_type) /= 'default' .and. (istep <= 1 .or. fyear /= fyear_old)) then
       if (my_task == master_task) then
          write (nu_diag,*) ' Set current forcing data year = ',fyear
       endif
    endif
    call icepack_query_tracer_indices(nt_Tsfc_out=nt_Tsfc)
    call icepack_warnings_flush(nu_diag)
    if (icepack_warnings_aborted()) call abort_ice(error_message=subname, file=__FILE__, line=__LINE__)
    !-------------------------------------------------------------------
    ! Read and interpolate atmospheric data
    !-------------------------------------------------------------------
    if (local_debug .and. my_task == master_task) then
       write(nu_diag,*) subname,' fyear         : ',fyear
       write(nu_diag,*) subname,' atm_data_type : ',trim(atm_data_type)
       write(nu_diag,*) subname,' msec          : ',int(msec)
    endif
    if (trim(atm_data_type) == 'ncar') then
       call ncar_data
    elseif (trim(atm_data_type) == 'JRA55_gx1') then
       call JRA55_data
    elseif (trim(atm_data_type) == 'JRA55_gx3') then
       call JRA55_data
    elseif (trim(atm_data_type) == 'JRA55_tx1') then
       call JRA55_data
    elseif (trim(atm_data_type) == 'hadgem') then
       call hadgem_data
    elseif (trim(atm_data_type) == 'monthly') then
       call monthly_data
    elseif (trim(atm_data_type) == 'oned') then
       call oned_data
    elseif (trim(atm_data_type) == 'box2001') then
       call box2001_data_atm
    elseif (trim(atm_data_type) == 'uniform_northeast') then
       call uniform_data_atm('NE')
    elseif (trim(atm_data_type) == 'uniform_north') then
       call uniform_data_atm('N')
    elseif (trim(atm_data_type) == 'uniform_east') then
       call uniform_data_atm('E')
    elseif (trim(atm_data_type) == 'uniform_south') then
       call uniform_data_atm('S')
    elseif (trim(atm_data_type) == 'uniform_west') then
       call uniform_data_atm('W')
    elseif (trim(atm_data_type) == 'calm') then
       call uniform_data_atm('N',c0) ! direction does not matter when c0
    elseif (trim(atm_data_type) == 'hycom') then
       call hycom_atm_data
    elseif (trim(atm_data_type) == 'ERA5') then
       call ERA5_data
    else    ! default values set in init_flux
       return
    endif
    !-------------------------------------------------------------------
    ! Convert forcing data to fields needed by ice model
    !-------------------------------------------------------------------
    !$OMP PARALLEL DO PRIVATE(iblk,ilo,ihi,jlo,jhi,this_block)
    do iblk = 1, nblocks
       this_block = get_block(blocks_ice(iblk),iblk)
       ilo = this_block%ilo
       ihi = this_block%ihi
       jlo = this_block%jlo
       jhi = this_block%jhi
       call prepare_atm_forcing (nx_block, ny_block, ilo, ihi, jlo, jhi, &
                                 hm    (:,:,iblk),   &
                                 Tair  (:,:,iblk),   &
                                 fsw   (:,:,iblk),   &
                                 cldf  (:,:,iblk),   &
                                 flw   (:,:,iblk),   &
                                 frain (:,:,iblk),   &
                                 fsnow (:,:,iblk),   &
                                 Qa    (:,:,iblk),   &
                                 rhoa  (:,:,iblk),   &
                                 uatm  (:,:,iblk),   &
                                 vatm  (:,:,iblk),   &
                                 strax (:,:,iblk),   &
                                 stray (:,:,iblk),   &
                                 zlvl  (:,:,iblk),   &
                                 wind  (:,:,iblk),   &
                                 swvdr (:,:,iblk),   &
                                 swvdf (:,:,iblk),   &
                                 swidr (:,:,iblk),   &
                                 swidf (:,:,iblk),   &
                                 potT  (:,:,iblk),   &
                                 ANGLET(:,:,iblk),   &
                                 trcr  (:,:,nt_Tsfc,iblk), &
                                 sst   (:,:,iblk),   &
                                 aice  (:,:,iblk) )
    enddo ! iblk
    !$OMP END PARALLEL DO
    call ice_timer_start(timer_bound)
    call ice_HaloUpdate (swvdr, halo_info, field_loc_center, field_type_scalar)
    call ice_HaloUpdate (swvdf, halo_info, field_loc_center, field_type_scalar)
    call ice_HaloUpdate (swidr, halo_info, field_loc_center, field_type_scalar)
    call ice_HaloUpdate (swidf, halo_info, field_loc_center, field_type_scalar)
    call ice_timer_stop(timer_bound)
    call ice_timer_stop(timer_forcing)
  end subroutine get_forcing_atmo

  !=======================================================================
  subroutine get_forcing_ocn(dt)
    ! Read and interpolate annual climatologies of SSS and SST.
    ! Restore model SST to data if desired.
    ! Interpolate ocean fields to U grid if necessary.
    ! time step
    use ice_blocks, only             :  block, get_block
    use ice_domain, only             :  nblocks, blocks_ice, distrb_info
    use ice_grid, only               :  tmask, umask
    use ice_global_reductions, only  :  global_minval, global_maxval
    use ice_flux, only               :  sss, sst, Tf, uocn, vocn, ss_tltx, ss_tlty
    integer (kind=int_kind)          :: i, j, iblk, &      ! block index
                                        ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
                                        modadj, &          ! adjustment to make mod a postive number
                                        fyear_old          ! fyear setting on last timestep
    type (block)                     :: this_block
    real (kind=dbl_kind), intent(in) :: dt
    real (kind=dbl_kind)             :: vmin, vmax, dayfrac, weekfrac, monthfrac, secday
    real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks) :: &
         eta_tide, Utide_tr, Vtide_tr, utide_cur, vtide_cur
    character(len=*), parameter      :: subname = '(get_forcing_ocn)'
    call icepack_query_parameters(secday_out=secday)
    if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,' fdbg start'
    call ice_timer_start(timer_forcing)
    fyear_old = fyear
    modadj    = abs((min(0,myear-fyear_init)/ycycle+1)*ycycle)
    fyear     = fyear_init + mod(myear-fyear_init+modadj,ycycle)
    if (trim(atm_data_type) /= 'default' .and. (istep <= 1 .or. fyear /= fyear_old)) then
       if (my_task == master_task) then
          write (nu_diag,*) ' Set current forcing data year = ',fyear
       endif
    endif
    if (local_debug .and. my_task == master_task) then
       write(nu_diag,*) subname,' fyear         : ',fyear
       write(nu_diag,*) subname,' ocn_data_type : ',trim(ocn_data_type)
    endif
    if (trim(ocn_data_type) == 'clim') then
       call ocn_data_clim(dt)
    elseif (trim(ocn_data_type) == 'ncar' .or. trim(ocn_data_type) == 'ISPOL') then
       call ocn_data_ncar(dt)
    elseif (trim(ocn_data_type) == 'hadgem_sst' .or. trim(ocn_data_type) == 'hadgem_sst_uvocn') then
       call ocn_data_hadgem(dt)
    elseif (trim(ocn_data_type) == 'oned') then
       call ocn_data_oned
    elseif (trim(ocn_data_type) == 'hycom') then
       !MHRI: NOT IMPLEMENTED YET
       ! call ocn_data_hycom(dt)
    elseif (trim(ocn_data_type) == 'box2001') then
       call box2001_data_ocn
    elseif (trim(ocn_data_type) == 'uniform_northeast') then
       ! tcraig, not time varying
       call uniform_data_ocn('NE',p1)
    elseif (trim(ocn_data_type) == 'uniform_east') then
       call uniform_data_ocn('E',p1)
    elseif (trim(ocn_data_type) == 'uniform_north') then
       call uniform_data_ocn('N',p1)
    elseif (trim(ocn_data_type) == 'calm') then
       call uniform_data_ocn('N',c0) ! directon does not matter for c0
    elseif (trim(ocn_data_type) == 'AFIM') then
       call AFIM_data(dt)
       if (trim(tide_data_type) == 'harmonic') then
          call compute_tides_at_time(eta_tide, Utide_tr, Vtide_tr, utide_cur, vtide_cur)
          if (tide_use_currents) then
             uocn(:,:,:) = uocn(:,:,:) + utide_cur(:,:,:)
             vocn(:,:,:) = vocn(:,:,:) + vtide_cur(:,:,:)
          endif
          ! v1: keep SSH coupling disabled until currents are stable
          if (tide_use_ssh) then
             ss_tltx(:,:,:) = c0
             ss_tlty(:,:,:) = c0
          endif
       endif
    endif
    call ice_timer_stop(timer_forcing)
  end subroutine get_forcing_ocn

  !=======================================================================
  subroutine read_data (flag, recd, yr, ixm, ixx, ixp, maxrec, data_file, field_data, field_loc, field_type)
    ! If data is at the beginning of a one-year record, get data from
    !  the previous year.
    ! If data is at the end of a one-year record, get data from the
    !  following year.
    ! If no earlier data exists (beginning of fyear_init), then
    !  (1) For monthly data, get data from the end of fyear_final.
    !  (2) For more frequent data, let the ixm value equal the
    !      first value of the year.
    ! If no later data exists (end of fyear_final), then
    !  (1) For monthly data, get data from the beginning of fyear_init.
    !  (2) For more frequent data, let the ixp value
    !      equal the last value of the year.
    ! In other words, we assume persistence when daily or 6-hourly
    !   data is missing, and we assume periodicity when monthly data
    !   is missing.
    use ice_diagnostics, only: debug_model_step
    logical (kind=log_kind), intent(in) :: flag
    integer (kind=int_kind), intent(in) :: recd         , & ! baseline record number
                                           yr           , & ! year of forcing data
                                           ixm, ixx, ixp, & ! record numbers of 3 data values relative to recd
                                           maxrec           ! maximum record value
    ! 2 values needed for interpolation
    real (kind=dbl_kind), dimension(nx_block,ny_block,2,max_blocks), intent(inout) :: field_data
    integer (kind=int_kind), intent(in) :: field_loc, &     ! location of field on staggered grid
                                           field_type       ! type of field (scalar, vector, angle)
    ! data file to be read
    character (char_len_long) :: data_file
    integer (kind=int_kind) :: nbits , & ! = 32 for single precision, 64 for double
                               nrec  , & ! record number to read
                               n2, n4, & ! like ixm and ixp, but adjusted at beginning and end of data
                               arg       ! value of time argument in field_data
    character(len=*), parameter :: subname = '(read_data)'
    if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'
    call ice_timer_start(timer_readwrite)  ! reading/writing
    nbits = 64              ! double precision data
    if (istep1 > debug_model_step) debug_forcing = .true.  !! debugging
    if (my_task==master_task .and. (debug_forcing)) then
       write(nu_diag,*) '  ', trim(data_file)
    endif
    if (flag) then
       !-----------------------------------------------------------------
       ! Initialize record counters
       ! (n2, n4 will change only at the very beginning or end of
       !  a forcing cycle.)
       !-----------------------------------------------------------------
       n2 = ixm
       n4 = ixp
       arg = 0
       !-----------------------------------------------------------------
       ! read data
       !-----------------------------------------------------------------
       if (ixm /= -99) then
          ! currently in first half of data interval
          if (ixx <= 1) then
             if (yr > fyear_init) then ! get data from previous year
                call file_year_atm (data_file, yr-1)
             else             ! yr = fyear_init, no prior data exists
                if (maxrec > 12) then ! extrapolate from first record
                   if (ixx == 1) n2 = ixx
                else          ! go to end of fyear_final
                   call file_year_atm (data_file, fyear_final)
                endif
             endif            ! yr > fyear_init
          endif               ! ixx <= 1
          call ice_open (nu_forcing, data_file, nbits)
          arg = 1
          nrec = recd + n2
          call ice_read (nu_forcing, nrec, field_data(:,:,arg,:), 'rda8', debug_forcing, field_loc, field_type)
          if (ixx==1 .and. my_task == master_task) close(nu_forcing)
       endif                  ! ixm ne -99
       ! always read ixx data from data file for current year
       call file_year_atm (data_file, yr)
       call ice_open (nu_forcing, data_file, nbits)
       arg = arg + 1
       nrec = recd + ixx
       call ice_read (nu_forcing, nrec, field_data(:,:,arg,:), 'rda8', debug_forcing, field_loc, field_type)
       if (ixp /= -99) then
          ! currently in latter half of data interval
          if (ixx==maxrec) then
             if (yr < fyear_final) then ! get data from following year
                if (my_task == master_task) close(nu_forcing)
                call file_year_atm (data_file, yr+1)
                call ice_open (nu_forcing, data_file, nbits)
             else             ! yr = fyear_final, no more data exists
                if (maxrec > 12) then ! extrapolate from ixx
                   n4 = ixx
                else          ! go to beginning of fyear_init
                   if (my_task == master_task) close(nu_forcing)
                   call file_year_atm (data_file, fyear_init)
                   call ice_open (nu_forcing, data_file, nbits)
                endif
             endif            ! yr < fyear_final
          endif               ! ixx = maxrec
          arg = arg + 1
          nrec = recd + n4
          call ice_read (nu_forcing, nrec, field_data(:,:,arg,:), 'rda8', debug_forcing, field_loc, field_type)
       endif                  ! ixp /= -99
       if (my_task == master_task) close(nu_forcing)
    endif                     ! flag
    call ice_timer_stop(timer_readwrite)  ! reading/writing
  end subroutine read_data

  !=======================================================================
  subroutine read_data_nc (flag, recd, yr, ixm, ixx, ixp, maxrec, data_file, fieldname, field_data, field_loc, field_type)
    ! If data is at the beginning of a one-year record, get data from
    !  the previous year.
    ! If data is at the end of a one-year record, get data from the
    !  following year.
    ! If no earlier data exists (beginning of fyear_init), then
    !  (1) For monthly data, get data from the end of fyear_final.
    !  (2) For more frequent data, let the ixm value equal the
    !      first value of the year.
    ! If no later data exists (end of fyear_final), then
    !  (1) For monthly data, get data from the beginning of fyear_init.
    !  (2) For more frequent data, let the ixp value
    !      equal the last value of the year.
    ! In other words, we assume persistence when daily or 6-hourly
    !   data is missing, and we assume periodicity when monthly data
    !   is missing.
    !
    ! Adapted by Alison McLaren, Met Office from read_data
    use ice_diagnostics, only: debug_model_step
    logical (kind=log_kind), intent(in) :: flag
    integer (kind=int_kind), intent(in) :: recd         , & ! baseline record number
                                           yr           , & ! year of forcing data
                                           ixm, ixx, ixp, & ! record numbers of 3 data values relative to recd
                                           maxrec           ! maximum record value
    ! data file to be read
    character (char_len_long) :: data_file
    ! field name in netCDF file
    character (char_len), intent(in) :: fieldname
    integer (kind=int_kind), intent(in) :: field_loc, & ! location of field on staggered grid
                                           field_type   ! type of field (scalar, vector, angle)
    ! 2 values needed for interpolation
    real (kind=dbl_kind), dimension(nx_block,ny_block,2,max_blocks), intent(out) :: field_data
    integer (kind=int_kind) :: nrec  , & ! record number to read
                               n2, n4, & ! like ixm and ixp, but adjusted at beginning and end of data
                               arg   , & ! value of time argument in field_data
                               fid       ! file id for netCDF routines
    character(len=*), parameter :: subname = '(read_data_nc)'
    if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'
    call ice_timer_start(timer_readwrite)  ! reading/writing
    if (istep1 > debug_model_step) debug_forcing = .true.  !! debugging
    if (my_task==master_task .and. (debug_forcing)) then
       write(nu_diag,*) '  ', trim(data_file)
    endif
    if (flag) then
       !-----------------------------------------------------------------
       ! Initialize record counters
       ! (n2, n4 will change only at the very beginning or end of a forcing cycle.)
       !-----------------------------------------------------------------
       n2 = ixm
       n4 = ixp
       arg = 0
       !-----------------------------------------------------------------
       ! read data
       !-----------------------------------------------------------------
       if (ixm /= -99) then
          ! currently in first half of data interval
          if (ixx <= 1) then
             if (yr > fyear_init) then ! get data from previous year
                call file_year_atm (data_file, yr-1)
             else             ! yr = fyear_init, no prior data exists
                if (maxrec > 12) then ! extrapolate from first record
                   if (ixx == 1) n2 = ixx
                else          ! go to end of fyear_final
                   call file_year_atm (data_file, fyear_final)
                endif
             endif            ! yr > fyear_init
          endif               ! ixx <= 1
          call ice_open_nc (data_file, fid)
          arg = 1
          nrec = recd + n2
          call ice_read_nc(fid, nrec, fieldname, field_data(:,:,arg,:), debug_forcing, field_loc, field_type)
          if (ixx==1) call ice_close_nc(fid)
       endif                  ! ixm ne -99
       ! always read ixx data from data file for current year
       call file_year_atm (data_file, yr)
       call ice_open_nc (data_file, fid)
       arg = arg + 1
       nrec = recd + ixx
       call ice_read_nc(fid, nrec, fieldname, field_data(:,:,arg,:), debug_forcing, field_loc, field_type)
       if (ixp /= -99) then
          ! currently in latter half of data interval
          if (ixx==maxrec) then
             if (yr < fyear_final) then ! get data from following year
                call ice_close_nc(fid)
                call file_year_atm (data_file, yr+1)
                call ice_open_nc (data_file, fid)
             else             ! yr = fyear_final, no more data exists
                if (maxrec > 12) then ! extrapolate from ixx
                   n4 = ixx
                else          ! go to beginning of fyear_init
                   call ice_close_nc(fid)
                   call file_year_atm (data_file, fyear_init)
                   call ice_open_nc (data_file, fid)
                endif
             endif            ! yr < fyear_final
          endif               ! ixx = maxrec
          arg = arg + 1
          nrec = recd + n4
          call ice_read_nc(fid, nrec, fieldname, field_data(:,:,arg,:), debug_forcing, field_loc, field_type)
       endif                  ! ixp /= -99
       call ice_close_nc(fid)
    endif                     ! flag
    call ice_timer_stop(timer_readwrite)  ! reading/writing
  end subroutine read_data_nc

  !=======================================================================
  subroutine read_data_nc_hycom (flag, recd, data_file, fieldname, field_data, field_loc, field_type)
    !  Data is assumed to cover the entire time period of simulation.
    !  It is not bounded by start of year nor end of year
    !  Data must be accesible both before and after (or on) the point in time
    !    Assume increasing timeaxis within the forcing files, but they do not
    !      have to be equal spaced. Read time vector from "MT" in "init_hycom"
    !
    ! Adapted by Mads Hvid Ribergaard, DMI from read_data_nc
    use ice_diagnostics, only: debug_model_step
    use ice_timers, only: ice_timer_start, ice_timer_stop, timer_readwrite
    logical (kind=log_kind), intent(in) :: flag
    ! baseline record number
    integer (kind=int_kind), intent(in) :: recd
    ! data file to be read
    character (char_len_long) :: data_file
    ! field name in netCDF file
    character (char_len), intent(in) :: fieldname
    integer (kind=int_kind), intent(in) :: field_loc, & ! location of field on staggered grid
                                           field_type   ! type of field (scalar, vector, angle)
    ! 2 values needed for interpolation
    real (kind=dbl_kind), dimension(nx_block,ny_block,2,max_blocks), intent(out) :: field_data
    ! file id for netCDF routines
    integer (kind=int_kind) :: fid
    character(len=*), parameter :: subname = '(read_data_nc_hycom)'
    if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'
    call ice_timer_start(timer_readwrite)  ! reading/writing
    if (istep1 > debug_model_step) debug_forcing = .true.  !! debugging
    if (my_task==master_task .and. (debug_forcing)) then
       write(nu_diag,*) '  ', trim(data_file)
    endif
    if (flag) then
       call ice_open_nc (data_file, fid)
       !-----------------------------------------------------------------
       ! read data
       !-----------------------------------------------------------------
       call ice_read_nc(fid, recd , fieldname, field_data(:,:,1,:), debug_forcing, field_loc, field_type)
       call ice_read_nc(fid, recd+1, fieldname, field_data(:,:,2,:), debug_forcing, field_loc, field_type)
       call ice_close_nc(fid)
    endif                     ! flag
    call ice_timer_stop(timer_readwrite)  ! reading/writing
  end subroutine read_data_nc_hycom

  !=======================================================================
  subroutine read_clim_data (readflag, recd, ixm, ixx, ixp, data_file, field_data, field_loc, field_type)
    ! Read data needed for interpolation, as in read_data.
    ! Assume a one-year cycle of climatological data, so that there is
    !  no need to get data from other years or to extrapolate data beyond
    !  the forcing time period.
    use ice_diagnostics, only: debug_model_step
    logical (kind=log_kind),intent(in) :: readflag
    integer (kind=int_kind), intent(in) :: recd            , & ! baseline record number
                                           ixm,ixx,ixp         ! record numbers of 3 data values relative to recd
    character (char_len_long), intent(in) ::  data_file
    integer (kind=int_kind), intent(in) :: field_loc, &      ! location of field on staggered grid
                                           field_type        ! type of field (scalar, vector, angle)
    ! 2 values needed for interpolation
    real (kind=dbl_kind), dimension(nx_block,ny_block,2,max_blocks), intent(inout) :: field_data
    integer (kind=int_kind) :: nbits, & ! = 32 for single precision, 64 for double
                               nrec , & ! record number to read
                               arg      ! value of time argument in field_data
    character(len=*), parameter :: subname = '(read_clim_data)'
    if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'
    call ice_timer_start(timer_readwrite)  ! reading/writing
    nbits = 64                ! double precision data
    if (istep1 > debug_model_step) debug_forcing = .true.  !! debugging
    if (my_task==master_task .and. (debug_forcing)) write(nu_diag,*) '  ', trim(data_file)
    if (readflag) then
       !-----------------------------------------------------------------
       ! read data
       !-----------------------------------------------------------------
       call ice_open (nu_forcing, data_file, nbits)
       arg = 0
       if (ixm /= -99) then
          arg = 1
          nrec = recd + ixm
          call ice_read (nu_forcing, nrec, field_data(:,:,arg,:), 'rda8', debug_forcing, field_loc, field_type)
       endif
       arg = arg + 1
       nrec = recd + ixx
       call ice_read (nu_forcing, nrec, field_data(:,:,arg,:), 'rda8', debug_forcing, field_loc, field_type)
       if (ixp /= -99) then
          arg = arg + 1
          nrec = recd + ixp
          call ice_read (nu_forcing, nrec, field_data(:,:,arg,:), 'rda8', debug_forcing, field_loc, field_type)
       endif
       if (my_task == master_task) close (nu_forcing)
    endif                     ! readflag
    call ice_timer_stop(timer_readwrite)  ! reading/writing
  end subroutine read_clim_data

  !=======================================================================
  subroutine read_clim_data_nc (readflag, recd, ixm, ixx, ixp, data_file, fieldname, field_data, field_loc, field_type)
    ! Read data needed for interpolation, as in read_data.
    ! Assume a one-year cycle of climatological data, so that there is
    !  no need to get data from other years or to extrapolate data beyond
    !  the forcing time period.
    use ice_diagnostics, only: debug_model_step
    logical (kind=log_kind),intent(in) :: readflag
    integer (kind=int_kind), intent(in) :: recd            , & ! baseline record number
                                           ixm,ixx,ixp         ! record numbers of 3 data values relative to recd
    character (char_len_long), intent(in) ::  data_file
    character (char_len), intent(in) :: fieldname              ! field name in netCDF file
    integer (kind=int_kind), intent(in) :: field_loc, &        ! location of field on staggered grid
                                           field_type          ! type of field (scalar, vector, angle)
    ! 2 values needed for interpolation
    real (kind=dbl_kind), dimension(nx_block,ny_block,2,max_blocks), intent(out) :: field_data
    integer (kind=int_kind) :: nrec, & ! record number to read
                               arg , & ! value of time argument in field_data
                               fid     ! file id for netCDF routines
    character(len=*), parameter :: subname = '(read_clim_data_nc)'
    if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'
    call ice_timer_start(timer_readwrite)  ! reading/writing
    if (istep1 > debug_model_step) debug_forcing = .true.  !! debugging
    if (my_task==master_task .and. (debug_forcing)) write(nu_diag,*) '  ', trim(data_file)
    if (readflag) then
       !-----------------------------------------------------------------
       ! read data
       !-----------------------------------------------------------------
       call ice_open_nc (data_file, fid)
       arg = 0
       if (ixm /= -99) then
          arg = 1
          nrec = recd + ixm
          call ice_read_nc(fid, nrec, fieldname, field_data(:,:,arg,:), debug_forcing, field_loc, field_type)
       endif
       arg = arg + 1
       nrec = recd + ixx
       call ice_read_nc(fid, nrec, fieldname, field_data(:,:,arg,:), debug_forcing, field_loc, field_type)
       if (ixp /= -99) then
          arg = arg + 1
          nrec = recd + ixp
          call ice_read_nc(fid, nrec, fieldname, field_data(:,:,arg,:), debug_forcing, field_loc, field_type)
       endif
       if (my_task == master_task) call ice_close_nc (fid)
    endif                     ! readflag
    call ice_timer_stop(timer_readwrite)  ! reading/writing
  end subroutine read_clim_data_nc

  !=======================================================================
  subroutine interp_coeff_monthly (recslot)
    ! Compute coefficients for interpolating monthly data to current time step.
    ! slot (1 or 2) for current record
    integer (kind=int_kind), intent(in) :: recslot 
    real (kind=dbl_kind) :: secday, & ! seconds in day
                            tt    , & ! days elapsed in current year
                            t1, t2    ! days elapsed at month midpoint
    ! month mid-points
    real (kind=dbl_kind) :: daymid(0:13)
    character(len=*), parameter :: subname = '(interp_coeff_monthly)'
    if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'
    call icepack_query_parameters(secday_out=secday)
    call icepack_warnings_flush(nu_diag)
    if (icepack_warnings_aborted()) call abort_ice(error_message=subname, file=__FILE__, line=__LINE__)
    daymid(1:13) = 14._dbl_kind   ! time frame ends 0 sec into day 15
    daymid(0)    = 14._dbl_kind - daymo(12)  ! Dec 15, 0 sec
    ! compute days since Jan 1, 00h, yday is the day counter for the year
    tt = real(yday-1,kind=dbl_kind) + real(msec,kind=dbl_kind)/secday
    ! Find neighboring times
    if (recslot==2) then      ! first half of month
       t2 = daycal(mmonth) + daymid(mmonth)   ! midpoint, current month
       if (mmonth == 1) then
          t1 = daymid(0)                 ! Dec 15 (0 sec)
       else
          t1 = daycal(mmonth-1) + daymid(mmonth-1) ! midpoint, previous month
       endif
    else ! second half of month
       t1 = daycal(mmonth) + daymid(mmonth)    ! midpoint, current month
       t2 = daycal(mmonth+1) + daymid(mmonth+1)! day 15 of next month (0 sec)
    endif
    if (tt < t1 .or. tt > t2) then
       write(nu_diag,*) subname,' ERROR in tt',tt,t1,t2
       call abort_ice (error_message=subname//' ERROR in tt', file=__FILE__, line=__LINE__)
    endif
    ! Compute coefficients
    c1intp = (t2 - tt) / (t2 - t1)
    c2intp =  c1 - c1intp
  end subroutine interp_coeff_monthly

  !=======================================================================
  subroutine interp_coeff (recnum, recslot, secint, dataloc)
    ! Compute coefficients for interpolating data to current time step.
    ! Works for any data interval that divides evenly into a
    !  year (daily, 6-hourly, etc.)
    ! Use interp_coef_monthly for monthly data.
    integer (kind=int_kind), intent(in) :: recnum , & ! record number for current data value
                                           recslot, & ! spline slot for current record
                                           dataloc    ! = 1 for data located in middle of time interval
                                                      ! = 2 for date located at end of time interval
    ! seconds in data interval
    real (kind=dbl_kind), intent(in) :: secint
    ! seconds in a day
    real (kind=dbl_kind) :: secday
    real (kind=dbl_kind) :: tt    , & ! seconds elapsed in current year
                            t1, t2, & ! seconds elapsed at data points
                            rcnum     ! recnum => dbl_kind
    character(len=*), parameter :: subname = '(interp_coeff)'
    if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'
    call icepack_query_parameters(secday_out=secday)
    call icepack_warnings_flush(nu_diag)
    if (icepack_warnings_aborted()) call abort_ice(error_message=subname, file=__FILE__, line=__LINE__)
    ! compute seconds since Jan 1, 00h, yday is the day counter for the year
    tt = real(yday-1,kind=dbl_kind)*secday + real(msec,kind=dbl_kind)
    ! Find neighboring times
    rcnum = real(recnum,kind=dbl_kind)
    if (recslot==2) then           ! current record goes in slot 2
       if (dataloc==1) then        ! data located at middle of interval
          t2 = (rcnum-p5)*secint
       else                        !  data located at end of interval
          t2 = rcnum*secint
       endif
       t1 = t2 - secint            !  - 1 interval
    else                           ! recslot = 1
       if (dataloc==1) then        ! data located at middle of interval
          t1 = (rcnum-p5)*secint
       else
          t1 = rcnum*secint        ! data located at end of interval
       endif
       t2 = t1 + secint            !  + 1 interval
    endif
    ! Compute coefficients
    c1intp =  abs((t2 - tt) / (t2 - t1))
    c2intp =  c1 - c1intp
    if (local_debug .and. my_task == master_task) then
       write(nu_diag,*) subname,'fdbg yday,sec = ',yday,msec
       write(nu_diag,*) subname,'fdbg tt = ',tt
       write(nu_diag,*) subname,'fdbg c12intp = ',c1intp,c2intp
    endif
  end subroutine interp_coeff

  !=======================================================================
  subroutine interp_coeff2 (tt, t1, t2)
    ! Compute coefficients for interpolating data to current time step.
    ! Works for any data interval using decimal daynumbers
    real (kind=dbl_kind), intent(in) :: tt    , &  ! current decimal daynumber
                                        t1, t2     ! first+last decimal daynumber
    character(len=*), parameter :: subname = '(interp_coeff2)'
    if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'
    ! Compute coefficients
    c1intp =  abs((t2 - tt) / (t2 - t1))
    c2intp =  c1 - c1intp
  end subroutine interp_coeff2

  !=======================================================================
  subroutine interpolate_data (field_data, field)
    ! Linear interpolation
    ! author: Elizabeth C. Hunke, LANL
    use ice_domain, only: nblocks
    ! 2 values used for interpolation
    real (kind=dbl_kind), dimension(nx_block,ny_block,2,max_blocks), intent(in) :: field_data
    ! interpolated field
    real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks), intent(out) :: field
    integer (kind=int_kind) :: i,j, iblk
    character(len=*), parameter :: subname = '(interpolate data)'
    if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'
    !$OMP PARALLEL DO PRIVATE(iblk,i,j)
    do iblk = 1, nblocks
       do j = 1, ny_block
          do i = 1, nx_block
             field(i,j,iblk) = c1intp * field_data(i,j,1,iblk) + c2intp * field_data(i,j,2,iblk)
          enddo
       enddo
    enddo !$OMP END PARALLEL DO
  end subroutine interpolate_data

  !=======================================================================
  subroutine file_year_atm(data_file, yr)
    character(char_len_long), intent(inout) :: data_file
    integer(kind=int_kind),   intent(in)    :: yr
    character(char_len_long)                :: tmpname
    integer(kind=int_kind)                  :: i
    character(len=*), parameter             :: subname = '(file_year_atm)'
    i = index(data_file, '_', back=.true.) + 1
    tmpname = data_file
    if (i > 1) then
       write(data_file, '(a,i4.4,a)') trim(tmpname(1:i-1)), yr, '.nc'
    else
       call abort_ice(error_message=subname//' ERROR: could not find underscore in base name', file=__FILE__, line=__LINE__)
    endif
  end subroutine file_year_atm 

  !=======================================================================
  subroutine file_year_ocn(data_file, fyr)
    character(char_len_long), intent(inout) :: data_file
    integer(kind=int_kind),   intent(in)    :: fyr
    character(char_len_long)                :: tmpname
    integer(kind=int_kind)                  :: i
    character(len=*),          parameter    :: subname = '(file_year_ocn)'
    i = index(data_file, '_', back=.true.) + 1
    tmpname = data_file
    write(data_file, '(a,i4.4,a)') trim(tmpname(1:i-1)), fyr, '.nc'
  end subroutine file_year_ocn

  !=======================================================================
  subroutine file_year_append(data_file, year)
    character(char_len_long), intent(inout) :: data_file
    integer(kind=int_kind),   intent(in)    :: year
    character(char_len_long)                :: tmpname
    integer(kind=int_kind)                  :: i
    character(len=*), parameter             :: subname = '(file_year_append)'
    i = index(data_file, '_', back=.true.) + 1
    tmpname = data_file
    if (i > 1) then
       write(data_file, '(a,i4.4,a)') trim(tmpname(1:i-1)), year, '.nc'
    else
       call abort_ice(error_message=subname//' ERROR: could not find underscore in base name', file=__FILE__, line=__LINE__)
    endif
  end subroutine file_year_append

  !=======================================================================
  subroutine prepare_atm_forcing (nx_block, ny_block, ilo, ihi, jlo, jhi, hm, Tair, fsw, cldf, flw, frain, fsnow, Qa, rhoa, uatm, vatm, &
                              strax, stray, zlvl, wind, swvdr, swvdf, swidr, swidf, potT, ANGLET, Tsfc, sst, aice)
    integer (kind=int_kind), intent(in) :: nx_block, ny_block, & ! block dimensions
                                           ilo,ihi,jlo,jhi       ! beginning and end of physical domain
    real (kind=dbl_kind), dimension(nx_block,ny_block), intent(in) :: ANGLET  , & ! ANGLE converted to T-cells
                                                                      Tsfc    , & ! ice skin temperature
                                                                      sst     , & ! sea surface temperature
                                                                      aice    , & ! ice area fraction
                                                                      hm          ! land mask
    real (kind=dbl_kind), dimension(nx_block,ny_block), intent(inout) :: fsw     , & ! incoming shortwave radiation (W/m^2)
                                                                         cldf    , & ! cloud fraction
                                                                         frain   , & ! rainfall rate (kg/m^2 s)
                                                                         fsnow   , & ! snowfall rate (kg/m^2 s)
                                                                         Tair    , & ! air temperature  (K)
                                                                         Qa      , & ! specific humidity (kg/kg)
                                                                         rhoa    , & ! air density (kg/m^3)
                                                                         uatm    , & ! wind velocity components (m/s)
                                                                         vatm    , &
                                                                         strax   , & ! wind stress components (N/m^2)
                                                                         stray   , &
                                                                         zlvl    , & ! atm level height (m)
                                                                         wind    , & ! wind speed (m/s)
                                                                         flw     , & ! incoming longwave radiation (W/m^2)
                                                                         swvdr   , & ! sw down, visible, direct  (W/m^2)
                                                                         swvdf   , & ! sw down, visible, diffuse (W/m^2)
                                                                         swidr   , & ! sw down, near IR, direct  (W/m^2)
                                                                         swidf   , & ! sw down, near IR, diffuse (W/m^2)
                                                                         potT        ! air potential temperature  (K)
    integer (kind=int_kind) :: i, j
    real (kind=dbl_kind) :: workx, worky, precip_factor, zlvl0, secday, Tffresh, puny
    logical (kind=log_kind) :: calc_strair
    character(len=*), parameter :: subname = '(prepare_atm_forcing)'
    if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'
    call icepack_query_parameters(Tffresh_out=Tffresh, puny_out=puny)
    call icepack_query_parameters(secday_out=secday)
    call icepack_query_parameters(calc_strair_out=calc_strair)
    call icepack_warnings_flush(nu_diag)
    if (icepack_warnings_aborted()) call abort_ice(error_message=subname, file=__FILE__, line=__LINE__)
    do j = jlo, jhi
       do i = ilo, ihi
          zlvl0 = c10 ! default
          !-----------------------------------------------------------------
          ! make sure interpolated values are physically realistic
          !-----------------------------------------------------------------
          cldf (i,j) = max(min(cldf(i,j),c1),c0)
          fsw  (i,j) = max(fsw(i,j),c0)
          fsnow(i,j) = max(fsnow(i,j),c0)
          rhoa (i,j) = max(rhoa(i,j),c0)
          Qa   (i,j) = max(Qa(i,j),c0)
          !        if (rhoa(i,j) .lt. puny) rhoa(i,j) = 1.3_dbl_kind
          !        if (Tair(i,j) .lt. puny) Tair(i,j) = Tffresh
          !        if (Qa(i,j) .lt. puny) Qa(i,j) = 0.0035_dbl_kind
       enddo ! i
    enddo  ! j
    !-----------------------------------------------------------------
    ! calculations specific to datasets
    !-----------------------------------------------------------------
    if (trim(atm_data_type) == 'ncar') then
       ! precip is in mm/month
       zlvl0 = c10
       do j = jlo, jhi
          do i = ilo, ihi
             ! correct known biases in NCAR data (as in CESM latm)
             Qa (i,j) = Qa (i,j) * 0.94_dbl_kind
             fsw(i,j) = fsw(i,j) * 0.92_dbl_kind
             ! downward longwave as in Parkinson and Washington (1979)
             call longwave_parkinson_washington(Tair(i,j), cldf(i,j), flw(i,j))
          enddo
       enddo
    elseif (trim(atm_data_type) == 'oned') then  ! rectangular grid
       ! precip is in kg/m^2/s
       zlvl0 = c10
       do j = jlo, jhi
          do i = ilo, ihi
             !-----------------------------------------------------------------
             ! compute downward longwave as in Parkinson and Washington (1979)
             !-----------------------------------------------------------------
             ! downward longwave as in Parkinson and Washington (1979)
             call longwave_parkinson_washington(Tair(i,j), cldf(i,j), flw(i,j))
             ! longwave based on Rosati and Miyakoda, JPO 18, p. 1607 (1988)
             ! call longwave_rosati_miyakoda(cldf(i,j), Tsfc(i,j), aice(i,j), sst(i,j), Qa(i,j), Tair(i,j), hm(i,j), flw(i,j))
          enddo
       enddo
    endif ! atm_data_type
    !-----------------------------------------------------------------
    ! Compute other fields needed by model
    !-----------------------------------------------------------------
    ! convert precipitation units to kg/m^2 s
    if (trim(precip_units) == 'mm_per_month') then
       precip_factor = c12/(secday*real(days_per_year,kind=dbl_kind))
    elseif (trim(precip_units) == 'mm_per_day') then
       precip_factor = c1/secday
    elseif (trim(precip_units) == 'mm_per_sec' .or. trim(precip_units) == 'mks') then
       precip_factor = c1    ! mm/sec = kg/m^2 s
    elseif (trim(precip_units) == 'm_per_sec') then
       precip_factor = c1000
    endif
    !------------------------------------------------------
    ! non-specific dataset calculations
    !------------------------------------------------------
    do j = jlo, jhi
       do i = ilo, ihi
          zlvl(i,j) = zlvl0
          potT(i,j) = Tair(i,j)
          ! divide shortwave into spectral bands
          swvdr(i,j) = fsw(i,j)*frcvdr        ! visible direct
          swvdf(i,j) = fsw(i,j)*frcvdf        ! visible diffuse
          swidr(i,j) = fsw(i,j)*frcidr        ! near IR direct
          swidf(i,j) = fsw(i,j)*frcidf        ! near IR diffuse
          ! convert precipitation units to kg/m^2 s
          fsnow(i,j) = fsnow(i,j) * precip_factor
       enddo  ! i
    enddo ! j
    ! determine whether precip is rain or snow
    ! HadGEM forcing provides separate snowfall and rainfall rather
    ! than total precipitation
    if (trim(atm_data_type) /= 'hadgem') then
       do j = jlo, jhi
          do i = ilo, ihi
             frain(i,j) = c0
             if (Tair(i,j) >= Tffresh) then
                frain(i,j) = fsnow(i,j)
                fsnow(i,j) = c0
             endif
          enddo ! i
       enddo ! j
    endif
    !--------------------------------------------------------------------
    ! calculate wind rotation
    !--------------------------------------------------------------------
    if (calc_strair) then
       if (rotate_wind) then
          do j = jlo, jhi
             do i = ilo, ihi
                wind(i,j) = sqrt(uatm(i,j)**2 + vatm(i,j)**2)
                !-----------------------------------------------------------------
                ! Rotate zonal/meridional vectors to local coordinates.
                ! Velocity comes in on T grid, but is oriented geographically ---
                ! need to rotate to pop-grid FIRST using ANGLET
                ! then interpolate to the U-cell centers  (otherwise we
                ! interpolate across the pole).
                ! Use ANGLET which is on the T grid !
                ! Atmo variables are needed in T cell centers in subroutine
                ! atmo_boundary_layer, and are interpolated to the U grid later as
                ! necessary.
                !-----------------------------------------------------------------
                workx      = uatm(i,j) ! wind velocity, m/s
                worky      = vatm(i,j)
                ! convert to POP grid, note: uatm, vatm, wind are on the T-grid here
                uatm (i,j) = workx*cos(ANGLET(i,j)) + worky*sin(ANGLET(i,j))
                vatm (i,j) = worky*cos(ANGLET(i,j)) - workx*sin(ANGLET(i,j))
             enddo ! i
          enddo ! j
       else ! not rotated
          do j = jlo, jhi
             do i = ilo, ihi
                wind(i,j) = sqrt(uatm(i,j)**2 + vatm(i,j)**2)
             enddo ! i
          enddo ! j
       endif ! rotated
    else ! strax, stray, wind are read from files
       if (rotate_wind) then
          do j = jlo, jhi
             do i = ilo, ihi
                workx      = strax(i,j) ! wind stress
                worky      = stray(i,j)
                ! convert to POP grid, note: strax, stray, wind are on the T-grid here
                strax(i,j) = workx*cos(ANGLET(i,j)) + worky*sin(ANGLET(i,j))
                stray(i,j) = worky*cos(ANGLET(i,j)) - workx*sin(ANGLET(i,j))
             enddo ! i
          enddo ! j
       else ! not rotated
          ! wind (speed) is already read from file, so all is in place
       endif ! rotated
    endif ! calc_strair
  end subroutine prepare_atm_forcing

  !=======================================================================
  subroutine longwave_parkinson_washington(Tair, cldf, flw)
    ! compute downward longwave as in Parkinson and Washington (1979)
    ! (for now)
    ! Parkinson, C. L. and W. M. Washington (1979),
    ! Large-scale numerical-model of sea ice,
    ! JGR, 84, 311-337, doi:10.1029/JC084iC01p00311
    real(kind=dbl_kind), intent(in) :: Tair , & ! air temperature  (K)
                                       cldf     ! cloud fraction
    ! incoming longwave radiation (W/m^2)
    real(kind=dbl_kind), intent(out) :: flw
    real(kind=dbl_kind) :: Tffresh, stefan_boltzmann
    character(len=*), parameter :: subname = '(longwave_parkinson_washington)'
    if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'
    call icepack_query_parameters(Tffresh_out=Tffresh, stefan_boltzmann_out=stefan_boltzmann)
    call icepack_warnings_flush(nu_diag)
    if (icepack_warnings_aborted()) call abort_ice(error_message=subname, file=__FILE__, line=__LINE__)
    flw = stefan_boltzmann*Tair**4 * (c1 - 0.261_dbl_kind * exp(-7.77e-4_dbl_kind*(Tffresh - Tair)**2)) * (c1 + 0.275_dbl_kind*cldf)
  end subroutine longwave_parkinson_washington

  !=======================================================================
  subroutine longwave_rosati_miyakoda(cldf, Tsfc, aice, sst, Qa, Tair, hm, flw)
    ! based on
    ! Rosati, A. and K. Miyakoda (1988),
    ! A general-circulation model for upper ocean simulation,
    ! J. Physical Oceanography, 18, 1601-1626,
    ! doi:10.1175/1520-0485(1988)018<1601:AGCMFU>2.0.CO;2
    real(kind=dbl_kind), intent(in) :: cldf , & ! cloud fraction
                                       Tsfc , & ! ice skin temperature
                                       aice , & ! ice area fraction
                                       sst  , & ! sea surface temperature
                                       Qa   , & ! specific humidity (kg/kg)
                                       Tair , & ! air temperature  (K)
                                       hm       ! land mask
    ! incoming longwave radiation (W/m^2)
    real(kind=dbl_kind), intent(out) :: flw
    real(kind=dbl_kind) :: fcc , & ! cloudiness modification
                           sstk, & ! ice/ocean surface temperature (K)
                           rtea, & ! square root of the vapour pressure
                           ptem, & ! potential air temperature (K)
                           qlwm
    real(kind=dbl_kind) :: Tffresh, stefan_boltzmann, emissivity
    character(len=*), parameter :: subname = '(longwave_rosati_miyakoda)'
    if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'
    call icepack_query_parameters(Tffresh_out=Tffresh, stefan_boltzmann_out=stefan_boltzmann, emissivity_out=emissivity)
    call icepack_warnings_flush(nu_diag)
    if (icepack_warnings_aborted()) call abort_ice(error_message=subname, file=__FILE__, line=__LINE__)
    fcc  = c1 - 0.8_dbl_kind * cldf
    sstk = (Tsfc * aice + sst * (c1 - aice)) + Tffresh
    rtea = sqrt(c1000*Qa /  (0.622_dbl_kind+0.378_dbl_kind*Qa))
    ptem = Tair ! get this from stability?
    qlwm = ptem * ptem * ptem  * ( ptem*(0.39_dbl_kind-0.05_dbl_kind*rtea)*fcc + c4*(sstk-ptem) )
    flw = emissivity*stefan_boltzmann * ( sstk**4 - qlwm )
    flw = flw * hm ! land mask
  end subroutine longwave_rosati_miyakoda

  !=======================================================================
  subroutine ncar_files (yr)
    ! Construct filenames based on the LANL naming conventions for NCAR data.
    ! Edit for other directory structures or filenames.
    ! Note: The year number in these filenames does not matter, because
    !       subroutine file_year_atm will insert the correct year.
    ! current forcing year
    integer (kind=int_kind), intent(in) :: yr
    character(len=*), parameter :: subname = '(ncar_files)'
    if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'
    F_fsw   = trim(atm_data_dir)//'/MONTHLY/swdn.1996.dat'
    call file_year_atm(F_fsw,yr)
    F_flw   = trim(atm_data_dir)//'/MONTHLY/cldf.1996.dat'
    call file_year_atm(F_flw,yr)
    F_rain  = trim(atm_data_dir)//'/MONTHLY/prec.1996.dat'
    call file_year_atm(F_rain,yr)
    F_uwind = trim(atm_data_dir)//'/4XDAILY/u_10.1996.dat'
    call file_year_atm(F_uwind,yr)
    F_vwind = trim(atm_data_dir)//'/4XDAILY/v_10.1996.dat'
    call file_year_atm(F_vwind,yr)
    F_tair  = trim(atm_data_dir)//'/4XDAILY/t_10.1996.dat'
    call file_year_atm(F_tair,yr)
    F_humid = trim(atm_data_dir)//'/4XDAILY/q_10.1996.dat'
    call file_year_atm(F_humid,yr)
    F_rhoa  = trim(atm_data_dir)//'/4XDAILY/dn10.1996.dat'
    call file_year_atm(F_rhoa,yr)
    if (my_task == master_task) then
       write (nu_diag,*) ' '
       write (nu_diag,*) 'Forcing data year =', fyear
       write (nu_diag,*) 'Atmospheric data files:'
       write (nu_diag,*) trim(F_fsw)
       write (nu_diag,*) trim(F_flw)
       write (nu_diag,*) trim(F_rain)
       write (nu_diag,*) trim(F_uwind)
       write (nu_diag,*) trim(F_vwind)
       write (nu_diag,*) trim(F_tair)
       write (nu_diag,*) trim(F_humid)
       write (nu_diag,*) trim(F_rhoa)
    endif                     ! master_task
  end subroutine ncar_files

  !=======================================================================
  subroutine ncar_data
    use ice_flux, only: fsw, fsnow, Tair, uatm, vatm, rhoa, Qa
    integer (kind=int_kind) :: ixm,ixx,ixp, & ! record numbers for neighboring months
                               recnum     , & ! record number
                               maxrec     , & ! maximum record number
                               recslot    , & ! spline slot for current record
                               dataloc    , & ! = 1 for data located in middle of time interval
                                              ! = 2 for date located at end of time interval
                               midmonth       ! middle day of month
    real (kind=dbl_kind) :: secday, &         ! number of seconds in day
                            sec6hr            ! number of seconds in 6 hours
    logical (kind=log_kind) :: readm, read6
    character(len=*), parameter :: subname = '(ncar_data)'
    if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'
    call icepack_query_parameters(secday_out=secday)
    call icepack_warnings_flush(nu_diag)
    if (icepack_warnings_aborted()) call abort_ice(error_message=subname, file=__FILE__, line=__LINE__)
    !-------------------------------------------------------------------
    ! monthly data
    ! Assume that monthly data values are located in the middle of the
    ! month.
    !-------------------------------------------------------------------
    midmonth = 15  ! data is given on 15th of every month
    ! midmonth = fix(p5 * real(daymo(mmonth)))  ! exact middle
    ! Compute record numbers for surrounding months
    maxrec = 12
    ixm  = mod(mmonth+maxrec-2,maxrec) + 1
    ixp  = mod(mmonth,         maxrec) + 1
    if (mday >= midmonth) ixm = -99  ! other two points will be used
    if (mday <  midmonth) ixp = -99
    ! Determine whether interpolation will use values 1:2 or 2:3
    ! recslot = 2 means we use values 1:2, with the current value (2)
    !  in the second slot    ! recslot = 1 means we use values 2:3, with the current value (2)
    !  in the first slot
    recslot = 1                             ! latter half of month
    if (mday < midmonth) recslot = 2        ! first half of month
    ! Find interpolation coefficients
    call interp_coeff_monthly (recslot)
    ! Read 2 monthly values
    readm = .false.
    if (istep==1 .or. (mday==midmonth .and. msec==0)) readm = .true.
    if (trim(atm_data_format) == 'bin') then
       call read_data (readm, 0, fyear, ixm, mmonth, ixp, maxrec, F_fsw, fsw_data, field_loc_center, field_type_scalar)
       call read_data (readm, 0, fyear, ixm, mmonth, ixp, maxrec, F_flw, cldf_data, field_loc_center, field_type_scalar)
       call read_data (readm, 0, fyear, ixm, mmonth, ixp, maxrec, F_rain, fsnow_data, field_loc_center, field_type_scalar)
    else
       call abort_ice (error_message=subname//'nonbinary atm_data_format unavailable', file=__FILE__, line=__LINE__)
    endif
    ! Interpolate to current time step
    call interpolate_data (fsw_data,   fsw)
    call interpolate_data (cldf_data,  cldf)
    call interpolate_data (fsnow_data, fsnow)
    !-------------------------------------------------------------------
    ! 6-hourly data
    ! Assume that the 6-hourly value is located at the end of the
    !  6-hour period.  This is the convention for NCEP reanalysis data.
    !  E.g. record 1 gives conditions at 6 am GMT on 1 January.
    !-------------------------------------------------------------------
    dataloc = 2               ! data located at end of interval
    sec6hr = secday/c4        ! seconds in 6 hours
    maxrec = 1460             ! 365*4
    ! current record number
    recnum = 4*int(yday) - 3 + int(real(msec,kind=dbl_kind)/sec6hr)
    ! Compute record numbers for surrounding data
    ixm = mod(recnum+maxrec-2,maxrec) + 1
    ixx = mod(recnum-1,       maxrec) + 1
    ! Compute interpolation coefficients
    ! If data is located at the end of the time interval, then the
    !  data value for the current record always goes in slot 2.
    recslot = 2
    ixp = -99
    call interp_coeff (recnum, recslot, sec6hr, dataloc)
    read6 = .false.
    if (istep==1 .or. oldrecnum /= recnum) read6 = .true.
    if (trim(atm_data_format) == 'bin') then
       call read_data (read6, 0, fyear, ixm, ixx, ixp, maxrec, F_tair, Tair_data, field_loc_center, field_type_scalar)
       call read_data (read6, 0, fyear, ixm, ixx, ixp, maxrec, F_uwind, uatm_data, field_loc_center, field_type_vector)
       call read_data (read6, 0, fyear, ixm, ixx, ixp, maxrec, F_vwind, vatm_data, field_loc_center, field_type_vector)
       call read_data (read6, 0, fyear, ixm, ixx, ixp, maxrec, F_rhoa, rhoa_data, field_loc_center, field_type_scalar)
       call read_data (read6, 0, fyear, ixm, ixx, ixp, maxrec, F_humid, Qa_data, field_loc_center, field_type_scalar)
    else
       call abort_ice (error_message=subname//'nonbinary atm_data_format unavailable', file=__FILE__, line=__LINE__)
    endif
    call interpolate_data (Tair_data, Tair)
    call interpolate_data (uatm_data, uatm)
    call interpolate_data (vatm_data, vatm)
    call interpolate_data (rhoa_data, rhoa)
    call interpolate_data (Qa_data,   Qa)
    ! Save record number for next time step
    oldrecnum = recnum
  end subroutine ncar_data

  !=======================================================================
  subroutine JRA55_gx1_files(yr)
    ! current forcing year
    integer (kind=int_kind), intent(in) :: yr
    character(len=*), parameter :: subname = '(JRA55_gx1_files)'
    if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'
    F_uwind = trim(atm_data_dir)//'/8XDAILY/JRA55_03hr_forcing_2005.nc'
    call file_year_atm(F_uwind,yr)
    if (my_task == master_task) then
       write (nu_diag,*) ' '
       write (nu_diag,*) 'Atmospheric data files:'
       write (nu_diag,*) trim(F_uwind)
    endif
  end subroutine JRA55_gx1_files

  !=======================================================================
  subroutine JRA55_tx1_files(yr)
    ! current forcing year
    integer (kind=int_kind), intent(in) :: yr
    character(len=*), parameter :: subname = '(JRA55_tx1_files)'
    if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'
    F_uwind = trim(atm_data_dir)//'/8XDAILY/JRA55_03hr_forcing_2005.nc'
    call file_year_atm(F_uwind,yr)
    if (my_task == master_task) then
       write (nu_diag,*) ' '
       write (nu_diag,*) 'Atmospheric data files:'
       write (nu_diag,*) trim(F_uwind)
    endif
  end subroutine JRA55_tx1_files

  !=======================================================================
  subroutine JRA55_gx3_files(yr)
    ! current forcing year
    integer (kind=int_kind), intent(in) :: yr
    character(len=*), parameter :: subname = '(JRA55_gx3_files)'
    if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'
    F_uwind = trim(atm_data_dir)//'/8XDAILY/JRA55_gx3_03hr_forcing_2005.nc'
    call file_year_atm(F_uwind,yr)
    if (my_task == master_task) then
       write (nu_diag,*) ' '
       write (nu_diag,*) 'Atmospheric data files:'
       write (nu_diag,*) trim(F_uwind)
    endif
  end subroutine JRA55_gx3_files

  !=======================================================================
  subroutine init_forcing_file(base_dir, prefix, year, file_out, label)
    character(len=*), intent(in)  :: base_dir, prefix, label
    integer(kind=int_kind), intent(in) :: year
    character(char_len_long), intent(out) :: file_out
    file_out = trim(base_dir)//trim(prefix)//'_'
    if (my_task == master_task) write(nu_diag, *) trim(label)//'_files was passed year as: ', year
    call file_year_append(file_out, year)
    if (my_task == master_task) write(nu_diag, *) trim(label)//'_files has file defined as : ', trim(file_out)
  end subroutine init_forcing_file

  !=======================================================================
  subroutine ERA5_files(fyr)
    integer (kind=int_kind), intent(in) :: fyr  !
    call init_forcing_file(atm_data_dir, 'era5_for_cice6', fyr, F_ERA5, 'ERA5')
  end subroutine ERA5_files

  !=======================================================================
  subroutine AFIM_files(fyr)
    integer (kind=int_kind), intent(in) :: fyr  !
    call init_forcing_file(ocn_data_dir, 'ORAS', fyr, F_AFIM, 'AFIM')
  end subroutine AFIM_files

  !=======================================================================
  subroutine JRA55_data
    use ice_blocks, only: block, get_block
    use ice_global_reductions, only: global_minval, global_maxval
    use ice_domain, only: nblocks, distrb_info
    use ice_flux, only: fsnow, Tair, uatm, vatm, Qa, fsw, flw
    use ice_grid, only: hm, tmask, umask
    use ice_state, only: aice
    use ice_calendar, only: days_per_year
    integer (kind=int_kind) :: ncid        , & ! netcdf file id
         i, j, n1    , &
         lfyear      , & ! local year value
         recnum      , & ! record number
         maxrec      , & ! maximum record number
         iblk            ! block index
    integer (kind=int_kind), save :: frec_info(2,2) = -99    ! remember prior values to reduce reading first dim is yr, recnum second dim is data1 data2
    real (kind=dbl_kind) :: sec3hr          , & ! number of seconds in 3 hours
         secday          , & ! number of seconds in day
         eps, tt         , & ! for interpolation coefficients
         Tffresh, vmin, vmax
    character(len=64) :: fieldname !netcdf field name
    character (char_len_long) :: F_uwind_old
    character(len=*), parameter :: subname = '(JRA55_data)'
    if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'
    call icepack_query_parameters(Tffresh_out=Tffresh)
    call icepack_query_parameters(secday_out=secday)
    call icepack_warnings_flush(nu_diag)
    if (icepack_warnings_aborted()) call abort_ice(error_message=subname, file=__FILE__, line=__LINE__)
    sec3hr = secday/c8        ! seconds in 3 hours
    maxrec = days_per_year * 8
    if (local_debug .and. my_task == master_task) then
       write(nu_diag,*) subname,'fdbg dpy, maxrec = ',days_per_year,maxrec
    endif
    !-------------------------------------------------------------------
    ! 3-hourly data
    ! states are instantaneous, 1st record is 00z Jan 1
    ! fluxes are 3 hour averages, 1st record is 00z-03z Jan 1
    ! interpolate states, do not interpolate fluxes
    !-------------------------------------------------------------------
    ! File is NETCDF with winds in NORTH and EAST direction
    ! file variable names are:
    ! glbrad   (shortwave W/m^2), 3 hr average
    ! dlwsfc   (longwave W/m^2), 3 hr average
    ! wndewd   (eastward wind m/s), instantaneous
    ! wndnwd   (northward wind m/s), instantaneous
    ! airtmp   (air temperature K), instantaneous
    ! spchmd   (specific humidity kg/kg), instantaneous
    ! ttlpcp   (precipitation kg/m s-1), 3 hr average
    !-------------------------------------------------------------------
    F_uwind_old = F_uwind
    if (F_uwind /= F_uwind_old .and. my_task == master_task) then
       write(nu_diag,*) subname,' reading forcing file = ',trim(F_uwind)
    endif
    call ice_open_nc(F_uwind,ncid)
    do n1 = 1, 2
       lfyear = fyear
       call file_year_atm(F_uwind,lfyear)
       if (n1 == 1) then
          recnum = 8*int(yday) - 7 + int(real(msec,kind=dbl_kind)/sec3hr)
          if (my_task == master_task .and. (recnum <= 2 .or. recnum >= maxrec-1)) then
             write(nu_diag,*) subname,' reading forcing file 1st ts = ',trim(F_uwind)
          endif
       elseif (n1 == 2) then
          recnum = 8*int(yday) - 7 + int(real(msec,kind=dbl_kind)/sec3hr) + 1
          if (recnum > maxrec) then
             lfyear = fyear + 1  ! next year
             if (lfyear > fyear_final) lfyear = fyear_init
             recnum = 1
             call file_year_atm(F_uwind,lfyear)
             if (my_task == master_task) then
                write(nu_diag,*) subname,' reading forcing file 2nd ts = ',trim(F_uwind)
             endif
             call ice_close_nc(ncid)
             call ice_open_nc(F_uwind,ncid)
          endif
       endif
       if (local_debug .and. my_task == master_task) then
          write(nu_diag,*) subname,'fdbg read recnum = ',recnum,n1
       endif
       ! to reduce reading, check whether it's the same data as last read
       if (lfyear /= frec_info(1,n1) .or. recnum /= frec_info(2,n1)) then
          ! check whether we can copy values from 2 to 1, should be faster than reading
          ! can only do this from 2 to 1 or 1 to 2 without setting up a temporary
          ! it's more likely that the values from data2 when time advances are needed in data1
          ! compare n1=1 year/record with data from last timestep at n1=2
          if (n1 == 1 .and. lfyear == frec_info(1,2) .and. recnum == frec_info(2,2)) then
             Tair_data(:,:,1,:) =  Tair_data(:,:,2,:)
             uatm_data(:,:,1,:) =  uatm_data(:,:,2,:)
             vatm_data(:,:,1,:) =  vatm_data(:,:,2,:)
             Qa_data(:,:,1,:) =    Qa_data(:,:,2,:)
             fsw_data(:,:,1,:) =   fsw_data(:,:,2,:)
             flw_data(:,:,1,:) =   flw_data(:,:,2,:)
             fsnow_data(:,:,1,:) = fsnow_data(:,:,2,:)
          else
             fieldname = 'airtmp'
             call ice_read_nc(ncid, recnum, fieldname, Tair_data (:,:,n1,:), local_debug, field_loc=field_loc_center, field_type=field_type_scalar)
             fieldname = 'wndewd'
             call ice_read_nc(ncid, recnum, fieldname, uatm_data (:,:,n1,:), local_debug, field_loc=field_loc_center, field_type=field_type_scalar)
             fieldname = 'wndnwd'
             call ice_read_nc(ncid, recnum, fieldname, vatm_data (:,:,n1,:), local_debug, field_loc=field_loc_center, field_type=field_type_scalar)
             fieldname = 'spchmd'
             call ice_read_nc(ncid, recnum, fieldname, Qa_data   (:,:,n1,:), local_debug, field_loc=field_loc_center, field_type=field_type_scalar)
             fieldname = 'glbrad'
             call ice_read_nc(ncid, recnum, fieldname, fsw_data  (:,:,n1,:), local_debug, field_loc=field_loc_center, field_type=field_type_scalar)
             fieldname = 'dlwsfc'
             call ice_read_nc(ncid, recnum, fieldname, flw_data  (:,:,n1,:), local_debug, field_loc=field_loc_center, field_type=field_type_scalar)
             fieldname = 'ttlpcp'
             call ice_read_nc(ncid, recnum, fieldname, fsnow_data(:,:,n1,:), local_debug, field_loc=field_loc_center, field_type=field_type_scalar)
          endif  ! copy data from n1=2 from last timestep to n1=1
       endif  ! input data is same as last timestep
       frec_info(1,n1) = lfyear
       frec_info(2,n1) = recnum
    enddo  ! n1
    call ice_close_nc(ncid)
    ! reset F_uwind to original year
    call file_year_atm(F_uwind,fyear)
    ! Compute interpolation coefficients
    eps = 1.0e-6
    tt = real(mod(msec,nint(sec3hr)),kind=dbl_kind)
    c2intp = tt / sec3hr
    if (c2intp < c0 .and. c2intp > c0-eps) c2intp = c0
    if (c2intp > c1 .and. c2intp < c1+eps) c2intp = c1
    c1intp = 1.0_dbl_kind - c2intp
    if (c2intp < c0 .or. c2intp > c1) then
       write(nu_diag,*) subname,' ERROR: c2intp = ',c2intp
       call abort_ice (error_message=subname//' ERROR: c2intp out of range', file=__FILE__, line=__LINE__)
    endif
    if (local_debug .and. my_task == master_task) then
       write(nu_diag,*) subname,'fdbg c12intp = ',c1intp,c2intp
    endif
    ! Interpolate
    call interpolate_data (Tair_data, Tair)
    call interpolate_data (uatm_data, uatm)
    call interpolate_data (vatm_data, vatm)
    call interpolate_data (Qa_data, Qa)
    fsw(:,:,:) = fsw_data(:,:,1,:)
    flw(:,:,:) = flw_data(:,:,1,:)
    fsnow(:,:,:) = fsnow_data(:,:,1,:)
    !$OMP PARALLEL DO PRIVATE(iblk,i,j)
    do iblk = 1, nblocks
       ! limit summer Tair values where ice is present
       do j = 1, ny_block
          do i = 1, nx_block
             if (aice(i,j,iblk) > p1) Tair(i,j,iblk) = min(Tair(i,j,iblk), Tffresh+p1)
          enddo
       enddo
       do j = 1, ny_block
          do i = 1, nx_block
             Qa  (i,j,iblk) = Qa  (i,j,iblk) * hm(i,j,iblk)
             Tair(i,j,iblk) = Tair(i,j,iblk) * hm(i,j,iblk)
             uatm(i,j,iblk) = uatm(i,j,iblk) * hm(i,j,iblk)
             vatm(i,j,iblk) = vatm(i,j,iblk) * hm(i,j,iblk)
             fsw (i,j,iblk) = fsw (i,j,iblk) * hm(i,j,iblk)
             flw (i,j,iblk) = flw (i,j,iblk) * hm(i,j,iblk)
             fsnow(i,j,iblk) = fsnow (i,j,iblk) * hm(i,j,iblk)
          enddo
       enddo
    enddo  ! iblk !$OMP END PARALLEL DO
    if (debug_forcing .or. local_debug) then
       if (my_task.eq.master_task) write (nu_diag,*) subname,'fdbg JRA55_bulk_data'
       vmin = global_minval(fsw,distrb_info,tmask)
       vmax = global_maxval(fsw,distrb_info,tmask)
       if (my_task.eq.master_task) write (nu_diag,*) subname,'fdbg fsw',vmin,vmax
       vmin = global_minval(flw,distrb_info,tmask)
       vmax = global_maxval(flw,distrb_info,tmask)
       if (my_task.eq.master_task) write (nu_diag,*) subname,'fdbg flw',vmin,vmax
       vmin =global_minval(fsnow,distrb_info,tmask)
       vmax =global_maxval(fsnow,distrb_info,tmask)
       if (my_task.eq.master_task) write (nu_diag,*) subname,'fdbg fsnow',vmin,vmax
       vmin = global_minval(Tair,distrb_info,tmask)
       vmax = global_maxval(Tair,distrb_info,tmask)
       if (my_task.eq.master_task) write (nu_diag,*) subname,'fdbg Tair',vmin,vmax
       vmin = global_minval(uatm,distrb_info,umask)
       vmax = global_maxval(uatm,distrb_info,umask)
       if (my_task.eq.master_task) write (nu_diag,*) subname,'fdbg uatm',vmin,vmax
       vmin = global_minval(vatm,distrb_info,umask)
       vmax = global_maxval(vatm,distrb_info,umask)
       if (my_task.eq.master_task) write (nu_diag,*) subname,'fdbg vatm',vmin,vmax
       vmin = global_minval(Qa,distrb_info,tmask)
       vmax = global_maxval(Qa,distrb_info,tmask)
       if (my_task.eq.master_task) write (nu_diag,*) subname,'fdbg Qa',vmin,vmax
    endif  ! debug_forcing
  end subroutine JRA55_data

  !=======================================================================
  !=======================================================================
  ! ERA5: subroutines to read and manipulate ERA5 forcing data.
  !
  ! ERA5 forcing is preprocessed offline from native ERA5 single-level
  ! reanalysis fields on Gadi to the CICE T grid. The preprocessing uses
  ! xESMF patch regridding with nearest-source-to-destination extrapolation
  ! and precomputed ESMF weights. Monthly regridded files are retained
  ! as the primary forcing product:
  !
  !    era5_for_cice6_YYYY_MM.nc 
  !
  ! The offline ERA5-to-CICE variable mapping is:
  !
  !    ERA5 2t       -> airtmp   (2-m air temperature, K)
  !    ERA5 msdwlwrf -> dlwsfc   (downward longwave flux, W/m^2)
  !    ERA5 msdwswrf -> glbrad   (downward shortwave flux, W/m^2)
  !    ERA5 mtpr     -> ttlpcp   (total precipitation rate, kg/m^2/s)
  !    ERA5 10u      -> wndewd   (eastward 10-m wind, m/s)
  !    ERA5 10v      -> wndnwd   (northward 10-m wind, m/s)
  !    ERA5 2d + sp  -> spchmd   (specific humidity, kg/kg)
  !    ERA5 sp       -> pair     (surface pressure, Pa)
  !    ERA5 msr      -> snowfall (snowfall rate, kg/m^2/s)
  !    ERA5 mtpr-msr -> rainfall (rainfall rate, kg/m^2/s)
  !    ERA5 blh      -> blh      (boundary-layer height, m)
  !    ERA5 10fg     -> windgust (10-m wind gust, m/s)
  !    ERA5 100u     -> wnd100ewd (eastward 100-m wind, m/s)
  !    ERA5 100v     -> wnd100nwd (northward 100-m wind, m/s)
  !
  ! Specific humidity is computed offline from 2-m dewpoint temperature
  ! and surface pressure using
  !
  !    E = a1 * exp(a3 * (d2m - T0) / (d2m - a4))
  !    q = (Rdry/Rvap) * E / (sp - (1 - Rdry/Rvap) * E)
  !
  ! where
  !
  !    Rdry = 287.0597 J kg^-1 K^-1
  !    Rvap = 461.5250 J kg^-1 K^-1
  !    a1   = 611.21 Pa
  !    a3   = 17.502
  !    a4   = 32.19 K
  !    T0   = 273.16 K
  !
  ! Within CICE, states are treated as instantaneous fields:
  !
  !    airtmp, wndewd, wndnwd, spchmd
  !
  ! Fluxes are treated as hourly averages:
  !
  !    glbrad, dlwsfc, ttlpcp
  !
  ! The first state record is YYYY-01-01-00:00Z. The first flux record
  ! represents the hourly average over YYYY-12-31-23:00Z to
  ! YYYY-01-01-00:00Z. State fields are interpolated in time; flux fields
  ! are not interpolated and are applied as record values.
  !
  ! File is NetCDF with winds in geographic east/north directions. Wind
  ! rotation to the model grid is handled later in the forcing workflow.
  !
  ! Forcing terms currently consist of:
  !
  !    glbrad   (shortwave W/m^2), 1 hr average
  !    dlwsfc   (longwave W/m^2), 1 hr average
  !    wndewd   (eastward wind m/s), instantaneous
  !    wndnwd   (northward wind m/s), instantaneous
  !    airtmp   (air temperature K), instantaneous
  !    spchmd   (specific humidity kg/kg), instantaneous
  !    ttlpcp   (precipitation kg/m^2/s), 1 hr average
  !
  !=======================================================================
  !=======================================================================
  subroutine ERA5_monthly_files(yr, mon)
    ! Build monthly ERA5 forcing filename.
    !
    ! Expected file pattern:
    !   ${atm_data_dir}/era5_for_cice6_YYYY_MM.nc
    !
    ! Example:
    !   /g/data/gv90/da1339/afim_input/ERA5/0p25/bilinear/monthly_cice6/era5_for_cice6_1994_10.nc
    integer(kind=int_kind), intent(in) :: yr
    integer(kind=int_kind), intent(in) :: mon
    character(len=4) :: cyear
    character(len=2) :: cmonth
    character(len=*), parameter :: subname = '(ERA5_monthly_files)'
    if (mon < 1 .or. mon > 12) then
       write(nu_diag,*) subname, ' ERROR invalid month = ', mon
       call abort_ice(error_message=subname//' invalid month', file=__FILE__, line=__LINE__)
    endif
    write(cyear, '(i4.4)') yr
    write(cmonth,'(i2.2)') mon
    F_ERA5 = trim(atm_data_dir)//'/era5_for_cice6_'//cyear//'_'//cmonth//'.nc'
    if (debug_forcing .or. local_debug) then
       if (my_task == master_task) then
          write(nu_diag,*) subname, ' ERA5 monthly file = ', trim(F_ERA5)
       endif
    endif
  end subroutine ERA5_monthly_files

  !=======================================================================
  logical(kind=log_kind) function ERA5_has_var(ncid, fieldname)
    ! Return true if fieldname exists in an already-open NetCDF file.
#ifdef USE_NETCDF
    use netcdf, only : nf90_noerr, nf90_inq_varid
#endif
    integer(kind=int_kind), intent(in) :: ncid
    character(len=*), intent(in)       :: fieldname
#ifdef USE_NETCDF
    integer :: ierr
    integer :: varid
    ierr = nf90_inq_varid(ncid, trim(fieldname), varid)
    ERA5_has_var = (ierr == nf90_noerr)
#else
    ERA5_has_var = .false.
#endif
  end function ERA5_has_var

  !=======================================================================
  subroutine ERA5_read_optional(ncid, recnum, fieldname, field_data, nslot)
    ! Read an optional ERA5 variable if present; otherwise zero its slot.
    !
    ! This lets older ERA5 files without boundary-layer diagnostics remain
    ! readable while new files can carry pair/snowfall/rainfall/blh/gust/100m wind.
    integer(kind=int_kind), intent(in) :: ncid
    integer(kind=int_kind), intent(in) :: recnum
    integer(kind=int_kind), intent(in) :: nslot
    character(len=*), intent(in)       :: fieldname
    real(kind=dbl_kind), dimension(nx_block,ny_block,2,max_blocks), intent(inout) :: field_data
    character(len=*), parameter :: subname = '(ERA5_read_optional)'
    if (ERA5_has_var(ncid, fieldname)) then
       call ice_read_nc(ncid, recnum, fieldname, field_data(:,:,nslot,:), &
            local_debug, field_loc=field_loc_center, field_type=field_type_scalar)
    else
       field_data(:,:,nslot,:) = c0
       if (debug_forcing .or. local_debug) then
          if (my_task == master_task) then
             write(nu_diag,*) subname, ' optional variable missing: ', trim(fieldname)
          endif
       endif
    endif
  end subroutine ERA5_read_optional

  !=======================================================================
  subroutine scale_era5_ant_coastal_precip(precip)
    ! Scale an ERA5 precipitation field over Antarctic coastal ocean cells.
    !
    ! This is an intentionally simple forcing-sensitivity perturbation:
    !   - applies only when era5_mod_var selects an Antarctic coastal
    !     precipitation perturbation
    !   - uses ocean cells adjacent to land
    !   - restricts to Southern Ocean / Antarctic latitudes
    !   - applies multiplicative era5_mod_fac
    !
    ! The field supplied through precip determines the experiment:
    !   - ttlpcp  for total-precipitation perturbations
    !   - fsnow   for snowfall-only perturbations
    !   - frain   for rainfall-only perturbations, if ever needed
    !
    ! For the Stage 3 experiment:
    !   fsnow <- era5_mod_fac * fsnow
    ! over Antarctic coastal ocean cells only.
    use ice_domain, only : nblocks
    use ice_grid,   only : hm, TLAT
    real(kind=dbl_kind), dimension(nx_block,ny_block,max_blocks), intent(inout) :: precip
    integer(kind=int_kind) :: i, j, iblk
    real(kind=dbl_kind)    :: pi, lat_cutoff
    character(len=*), parameter :: subname = '(scale_era5_ant_coastal_precip)'
    pi = c4 * atan(c1)
    ! TLAT is stored in radians. Use -50 deg as a deliberately blunt
    ! Southern Ocean / Antarctic coastal mask for this sensitivity test.
    lat_cutoff = -50.0_dbl_kind * pi / c180
    if (debug_forcing .or. local_debug) then
       if (my_task == master_task) then
          write(nu_diag,*) subname, ' applying era5_mod_fac = ', era5_mod_fac
          write(nu_diag,*) subname, ' lat_cutoff_deg      = ', -50.0_dbl_kind
       endif
    endif
    !$OMP PARALLEL DO PRIVATE(iblk,j,i)
    do iblk = 1, nblocks
       do j = 1, ny_block
          do i = 1, nx_block
             if (is_ant_fastice_forcing_cell(i,j,iblk,lat_cutoff)) then
                precip(i,j,iblk) = era5_mod_fac * precip(i,j,iblk)
             endif
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
  end subroutine scale_era5_ant_coastal_precip

  !=======================================================================
  logical(kind=log_kind) function is_ant_fastice_forcing_cell(i, j, iblk, lat_cutoff)
    ! Antarctic fast-ice-relevant ocean selector.
    !
    ! A cell is selected if:
    !   1. it is ocean on the T grid,
    !   2. it is south of lat_cutoff,
    !   3. it is either:
    !        a. directly adjacent to land/boundary on the T grid, or
    !        b. within the nonzero lateral-drag/form-factor footprint.
    !
    ! This keeps the crude coastal fringe while also targeting the
    ! high-resolution coastline / grounded-iceberg geometry represented in F2.
    use ice_grid, only : hm, TLAT, F2N, F2E
    integer(kind=int_kind), intent(in) :: i, j, iblk
    real(kind=dbl_kind),    intent(in) :: lat_cutoff
    logical(kind=log_kind) :: ocean_cell
    logical(kind=log_kind) :: adjacent_land
    logical(kind=log_kind) :: ant_cell
    logical(kind=log_kind) :: f2_cell
    real(kind=dbl_kind) :: F2_mag
    real(kind=dbl_kind), parameter :: F2_min = 1.0e-12_dbl_kind
    ocean_cell    = hm(i,j,iblk) > p5
    ant_cell      = TLAT(i,j,iblk) <= lat_cutoff
    adjacent_land = .false.
    if (i > 1       ) adjacent_land = adjacent_land .or. hm(i-1,j,iblk) <= p5
    if (i < nx_block) adjacent_land = adjacent_land .or. hm(i+1,j,iblk) <= p5
    if (j > 1       ) adjacent_land = adjacent_land .or. hm(i,j-1,iblk) <= p5
    if (j < ny_block) adjacent_land = adjacent_land .or. hm(i,j+1,iblk) <= p5
    f2_cell = abs(F2E(i,j,iblk)) > F2_min .or. abs(F2N(i,j,iblk)) > F2_min .or. &
         (i > 1 .and. abs(F2E(i-1,j,iblk)) > F2_min) .or. &
         (j > 1 .and. abs(F2N(i,j-1,iblk)) > F2_min)
    ! F2_mag  = sqrt(F2N(i,j,iblk)**2 + F2E(i,j,iblk)**2)
    ! f2_cell = F2_mag > F2_min
    is_ant_fastice_forcing_cell = ocean_cell .and. ant_cell .and. (adjacent_land .or. f2_cell)
  end function is_ant_fastice_forcing_cell

  !=======================================================================
  subroutine ERA5_debug_minmax(label, field)
    use ice_global_reductions, only : global_minval, global_maxval
    use ice_domain, only            : distrb_info
    use ice_grid, only              : tmask
    character(len=*), intent(in) :: label
    real(kind=dbl_kind), dimension(nx_block,ny_block,max_blocks), intent(in) :: field
    real(kind=dbl_kind) :: vmin, vmax
    if (debug_forcing .or. local_debug) then
       vmin = global_minval(field, distrb_info, tmask)
       vmax = global_maxval(field, distrb_info, tmask)
       if (my_task == master_task) write(nu_diag,*) '(ERA5_debug_minmax) ', trim(label), vmin, vmax
    endif
  end subroutine ERA5_debug_minmax

  !=======================================================================
  subroutine ERA5_data
    ! Reads monthly hourly ERA5 forcing files on the CICE T grid.
    !
    ! Stage 1 behaviour:
    !   - monthly file read-in
    !   - old forcing semantics preserved
    !   - ttlpcp is still mapped to fsnow
    !   - snowfall/rainfall/pair/blh/gust/100m winds are read if present,
    !     but not yet used to alter physics
    !
    ! Expected file pattern:
    !   ${atm_data_dir}/era5_for_cice6_YYYY_MM.nc

    use ice_blocks, only            : block, get_block
    use ice_global_reductions, only : global_minval, global_maxval
    use ice_domain, only            : nblocks, distrb_info
    use ice_flux, only              : fsnow, frain, Tair, uatm, vatm, Qa, fsw, flw
    use ice_grid, only              : hm, tmask, umask, TLAT
    use ice_state, only             : aice
    use ice_calendar, only          : mday, mmonth, myear

    integer(kind=int_kind) :: ncid
    integer(kind=int_kind) :: i, j, n1
    integer(kind=int_kind) :: lfyear
    integer(kind=int_kind) :: lmonth
    integer(kind=int_kind) :: recnum
    integer(kind=int_kind) :: rec_to_read
    integer(kind=int_kind) :: maxrec
    integer(kind=int_kind) :: iblk
    integer(kind=int_kind) :: modadj

    ! remember prior values to reduce reading:
    ! dim 1: year, month, recnum
    ! dim 2: data slot 1/2
    integer(kind=int_kind), save :: frec_info(3,2) = -99

    real(kind=dbl_kind) :: sec1hr
    real(kind=dbl_kind) :: secday
    real(kind=dbl_kind) :: eps, tt
    real(kind=dbl_kind) :: Tffresh
    real(kind=dbl_kind) :: vmin, vmax

    character(len=64)           :: fieldname
    character(len=*), parameter :: subname = '(ERA5_data)'

    call icepack_query_parameters(Tffresh_out=Tffresh)
    call icepack_query_parameters(secday_out=secday)
    call icepack_warnings_flush(nu_diag)
    if (icepack_warnings_aborted()) call abort_ice(error_message=subname, file=__FILE__, line=__LINE__)

    sec1hr = secday / 24.0_dbl_kind

    !-------------------------------------------------------------------
    ! 1. Work out forcing-cycle year and current month.
    !-------------------------------------------------------------------
    if (ycycle == 0) then
       lfyear = myear
    else
       modadj = abs((min(0, myear - fyear_init) / ycycle + 1) * ycycle)
       lfyear = fyear_init + mod(myear - fyear_init + modadj, ycycle)
    endif

    lmonth = mmonth

    ! Monthly ERA5 files contain one hourly record per hour in this month.
    maxrec = 24 * daymo(lmonth)

    recnum = 24 * (mday - 1) + int(real(msec,kind=dbl_kind) / sec1hr) + 1

    if (recnum < 1) recnum = 1
    if (recnum > maxrec) recnum = maxrec

    call ERA5_monthly_files(lfyear, lmonth)

    if (debug_forcing .or. local_debug) then
       if (my_task == master_task) then
          write(nu_diag,*) subname, ' myear, lfyear, mmonth : ', myear, lfyear, lmonth
          write(nu_diag,*) subname, ' monthly maxrec        : ', maxrec
          write(nu_diag,*) subname, ' current recnum        : ', recnum
          write(nu_diag,*) subname, ' ERA5 file             : ', trim(F_ERA5)
       endif
    endif

    !-------------------------------------------------------------------
    ! 2. Open monthly file and read current/next records.
    !
    ! State variables are interpolated using slots 1 and 2.
    ! Flux variables remain slot-1 values below, preserving the current
    ! convention that hourly fluxes are applied as record values.
    !
    ! At the final record of the month, persist the final record for slot 2.
    ! A later stage can add next-month look-ahead if needed.
    !-------------------------------------------------------------------
    call ice_open_nc(F_ERA5, ncid)

    do n1 = 1, 2

       if (n1 == 1) then
          rec_to_read = recnum
       else
          rec_to_read = min(recnum + 1, maxrec)
       endif

       if (lfyear /= frec_info(1,n1) .or. &
           lmonth /= frec_info(2,n1) .or. &
           rec_to_read /= frec_info(3,n1)) then

          ! If slot 1 now equals last step's slot 2, copy instead of reread.
          if (n1 == 1 .and. &
               lfyear == frec_info(1,2) .and. &
               lmonth == frec_info(2,2) .and. &
               rec_to_read == frec_info(3,2)) then

             Tair_data    (:,:,1,:) = Tair_data    (:,:,2,:)
             uatm_data    (:,:,1,:) = uatm_data    (:,:,2,:)
             vatm_data    (:,:,1,:) = vatm_data    (:,:,2,:)
             Qa_data      (:,:,1,:) = Qa_data      (:,:,2,:)
             fsw_data     (:,:,1,:) = fsw_data     (:,:,2,:)
             flw_data     (:,:,1,:) = flw_data     (:,:,2,:)
             fsnow_data   (:,:,1,:) = fsnow_data   (:,:,2,:)
             pair_data    (:,:,1,:) = pair_data    (:,:,2,:)
             snowfall_data(:,:,1,:) = snowfall_data(:,:,2,:)
             frain_data   (:,:,1,:) = frain_data   (:,:,2,:)
             blh_data     (:,:,1,:) = blh_data     (:,:,2,:)
             windgust_data(:,:,1,:) = windgust_data(:,:,2,:)
             uatm100_data (:,:,1,:) = uatm100_data (:,:,2,:)
             vatm100_data (:,:,1,:) = vatm100_data (:,:,2,:)

          else

             ! Required variables: abort if missing.
             fieldname = 'airtmp'
             call ice_read_nc(ncid, rec_to_read, fieldname, Tair_data(:,:,n1,:), local_debug, field_loc=field_loc_center, field_type=field_type_scalar)
             ! call ERA5_debug_minmax('raw airtmp', Tair_data(:,:,n1,:))

             fieldname = 'wndewd'
             call ice_read_nc(ncid, rec_to_read, fieldname, uatm_data(:,:,n1,:), local_debug, field_loc=field_loc_center, field_type=field_type_scalar)
             ! call ERA5_debug_minmax('raw eastward wind', uatm_data(:,:,n1,:))

             fieldname = 'wndnwd'
             call ice_read_nc(ncid, rec_to_read, fieldname, vatm_data(:,:,n1,:), local_debug, field_loc=field_loc_center, field_type=field_type_scalar)
             ! call ERA5_debug_minmax('raw northward wind', vatm_data(:,:,n1,:))

             fieldname = 'spchmd'
             call ice_read_nc(ncid, rec_to_read, fieldname, Qa_data(:,:,n1,:)  , local_debug, field_loc=field_loc_center, field_type=field_type_scalar)
             ! call ERA5_debug_minmax('raw specific humidity', Qa_data(:,:,n1,:))

             fieldname = 'glbrad'
             call ice_read_nc(ncid, rec_to_read, fieldname, fsw_data(:,:,n1,:) , local_debug, field_loc=field_loc_center, field_type=field_type_scalar)
             ! call ERA5_debug_minmax('raw downward shortwave radiation', fsw_data(:,:,n1,:))

             fieldname = 'dlwsfc'
             call ice_read_nc(ncid, rec_to_read, fieldname, flw_data(:,:,n1,:) , local_debug, field_loc=field_loc_center, field_type=field_type_scalar)
             ! call ERA5_debug_minmax('raw downward longwave radiation', flw_data(:,:,n1,:))

             fieldname = 'ttlpcp'
             call ice_read_nc(ncid, rec_to_read, fieldname, fsnow_data(:,:,n1,:), local_debug, field_loc=field_loc_center, field_type=field_type_scalar)
             ! call ERA5_debug_minmax('raw total preceipitation', fsnow_data(:,:,n1,:))

             fieldname = 'pair'
             call ice_read_nc(ncid, rec_to_read, fieldname, pair_data(:,:,n1,:), local_debug, field_loc=field_loc_center, field_type=field_type_scalar)
             ! call ERA5_debug_minmax('raw surface pressure', pair_data(:,:,n1,:))

             fieldname = 'snowfall'
             call ice_read_nc(ncid, rec_to_read, fieldname, snowfall_data(:,:,n1,:), local_debug, field_loc=field_loc_center, field_type=field_type_scalar)
             ! call ERA5_debug_minmax('raw snowfall', snowfall_data(:,:,n1,:))

             fieldname = 'rainfall'
             call ice_read_nc(ncid, rec_to_read, fieldname, frain_data(:,:,n1,:), local_debug, field_loc=field_loc_center, field_type=field_type_scalar)
             ! call ERA5_debug_minmax('raw rainfall', frain_data(:,:,n1,:))

             fieldname = 'blh'
             call ice_read_nc(ncid, rec_to_read, fieldname, blh_data(:,:,n1,:), local_debug, field_loc=field_loc_center, field_type=field_type_scalar)
             ! call ERA5_debug_minmax('raw boundary layer height', blh_data(:,:,n1,:))

             fieldname = 'windgust'
             call ice_read_nc(ncid, rec_to_read, fieldname, windgust_data(:,:,n1,:), local_debug, field_loc=field_loc_center, field_type=field_type_scalar)
             ! call ERA5_debug_minmax('raw surface wind gust magnitude', windgust_data(:,:,n1,:))

             fieldname = 'wnd100ewd'
             call ice_read_nc(ncid, rec_to_read, fieldname, uatm100_data(:,:,n1,:), local_debug, field_loc=field_loc_center, field_type=field_type_scalar)
             ! call ERA5_debug_minmax('raw 100 metre eastward wind', uatm100_data(:,:,n1,:))

             fieldname = 'wnd100nwd'
             call ice_read_nc(ncid, rec_to_read, fieldname, vatm100_data(:,:,n1,:), local_debug, field_loc=field_loc_center, field_type=field_type_scalar)
             ! call ERA5_debug_minmax('raw 100 metre northward wind', vatm100_data(:,:,n1,:))

          endif

          frec_info(1,n1) = lfyear
          frec_info(2,n1) = lmonth
          frec_info(3,n1) = rec_to_read

       endif

    enddo

    call ice_close_nc(ncid)

    !-------------------------------------------------------------------
    ! 3. Interpolation coefficients within the current hour.
    !-------------------------------------------------------------------
    eps    = 1.0e-6_dbl_kind
    tt     = real(mod(msec, nint(sec1hr)), kind=dbl_kind)
    c2intp = tt / sec1hr

    if (c2intp < c0 .and. c2intp > c0 - eps) c2intp = c0
    if (c2intp > c1 .and. c2intp < c1 + eps) c2intp = c1

    c1intp = c1 - c2intp

    if (c2intp < c0 .or. c2intp > c1) then
       write(nu_diag,*) subname, ' ERROR: c2intp = ', c2intp
       call abort_ice(error_message=subname//' ERROR: c2intp out of range', file=__FILE__, line=__LINE__)
    endif

    if (debug_forcing .or. local_debug) then
       if (my_task == master_task) write(nu_diag,*) subname, ' c1,c2 = ', c1intp, c2intp
    endif

    !-------------------------------------------------------------------
    ! 4. Apply current Stage 2 semantics.
    !
    ! States are interpolated.
    ! Fluxes are applied as current hourly record values.
    !
    ! Stage 2 uses ERA5 precipitation phase directly:
    !   fsnow <- snowfall
    !   frain <- rainfall
    !
    ! ttlpcp is still read into fsnow_data as a diagnostic/check field,
    ! but it is no longer passed directly into CICE as snowfall.
    !-------------------------------------------------------------------
    call interpolate_data(Tair_data, Tair)
    call interpolate_data(uatm_data, uatm)
    call interpolate_data(vatm_data, vatm)
    call interpolate_data(Qa_data,   Qa)

    fsw  (:,:,:) = fsw_data      (:,:,1,:)
    flw  (:,:,:) = flw_data      (:,:,1,:)
    fsnow(:,:,:) = snowfall_data (:,:,1,:)
    frain(:,:,:) = frain_data    (:,:,1,:)

    ! Defensive clipping. These should already be non-negative after shuga
    ! preprocessing, but this protects CICE from tiny numerical artifacts.
    fsnow(:,:,:) = max(fsnow(:,:,:), c0)
    frain(:,:,:) = max(frain(:,:,:), c0)

    if (trim(era5_mod_var) == 'snow_ant_coast' .or. &
         trim(era5_mod_var) == 'snowfall_ant_coast') then
       call scale_era5_ant_coastal_precip(fsnow)
    endif

    !-------------------------------------------------------------------
    ! 5. Existing masking and summer Tair limiter.
    !-------------------------------------------------------------------
    !$OMP PARALLEL DO PRIVATE(iblk,i,j)
    do iblk = 1, nblocks

       do j = 1, ny_block
          do i = 1, nx_block
             if (aice(i,j,iblk) > p1) Tair(i,j,iblk) = min(Tair(i,j,iblk), Tffresh + p1)
          enddo
       enddo

       do j = 1, ny_block
          do i = 1, nx_block
             Qa   (i,j,iblk) = Qa   (i,j,iblk) * hm(i,j,iblk)
             Tair (i,j,iblk) = Tair (i,j,iblk) * hm(i,j,iblk)
             uatm (i,j,iblk) = uatm (i,j,iblk) * hm(i,j,iblk)
             vatm (i,j,iblk) = vatm (i,j,iblk) * hm(i,j,iblk)
             fsw  (i,j,iblk) = fsw  (i,j,iblk) * hm(i,j,iblk)
             flw  (i,j,iblk) = flw  (i,j,iblk) * hm(i,j,iblk)
             fsnow(i,j,iblk) = fsnow(i,j,iblk) * hm(i,j,iblk)
             frain(i,j,iblk) = frain(i,j,iblk) * hm(i,j,iblk)
          enddo
       enddo

    enddo
    !$OMP END PARALLEL DO

    !-------------------------------------------------------------------
    ! 6. Diagnostics.
    !-------------------------------------------------------------------
    if (debug_forcing .or. local_debug) then

       if (my_task.eq.master_task) write(nu_diag,*) subname, ' ERA5 monthly forcing record = ', recnum

       vmin = global_minval(fsw,distrb_info,tmask)
       vmax = global_maxval(fsw,distrb_info,tmask)
       if (my_task.eq.master_task) write(nu_diag,*) subname, '  fsw', vmin, vmax

       vmin = global_minval(flw,distrb_info,tmask)
       vmax = global_maxval(flw,distrb_info,tmask)
       if (my_task.eq.master_task) write(nu_diag,*) subname, '  flw', vmin, vmax

       vmin = global_minval(fsnow,distrb_info,tmask)
       vmax = global_maxval(fsnow,distrb_info,tmask)
       if (my_task.eq.master_task) write(nu_diag,*) subname, '  fsnow', vmin, vmax

       vmin = global_minval(frain,distrb_info,tmask)
       vmax = global_maxval(frain,distrb_info,tmask)
       if (my_task.eq.master_task) write(nu_diag,*) subname, '  frain', vmin, vmax

       vmin = global_minval(Tair,distrb_info,tmask)
       vmax = global_maxval(Tair,distrb_info,tmask)
       if (my_task.eq.master_task) write(nu_diag,*) subname, '  Tair', vmin, vmax

       vmin = global_minval(uatm,distrb_info,umask)
       vmax = global_maxval(uatm,distrb_info,umask)
       if (my_task.eq.master_task) write(nu_diag,*) subname, '  uatm', vmin, vmax

       vmin = global_minval(vatm,distrb_info,umask)
       vmax = global_maxval(vatm,distrb_info,umask)
       if (my_task.eq.master_task) write(nu_diag,*) subname, '  vatm', vmin, vmax

       vmin = global_minval(Qa,distrb_info,tmask)
       vmax = global_maxval(Qa,distrb_info,tmask)
       if (my_task.eq.master_task) write(nu_diag,*) subname, '  Qa', vmin, vmax

       vmin = global_minval(pair_data(:,:,1,:),distrb_info,tmask)
       vmax = global_maxval(pair_data(:,:,1,:),distrb_info,tmask)
       if (my_task.eq.master_task) write(nu_diag,*) subname, '  pair_data(slot1)', vmin, vmax

       vmin = global_minval(snowfall_data(:,:,1,:),distrb_info,tmask)
       vmax = global_maxval(snowfall_data(:,:,1,:),distrb_info,tmask)
       if (my_task.eq.master_task) write(nu_diag,*) subname, '  snowfall_data(slot1)', vmin, vmax

       vmin = global_minval(frain_data(:,:,1,:),distrb_info,tmask)
       vmax = global_maxval(frain_data(:,:,1,:),distrb_info,tmask)
       if (my_task.eq.master_task) write(nu_diag,*) subname, '  rainfall/frain_data(slot1)', vmin, vmax

       vmin = global_minval(blh_data(:,:,1,:),distrb_info,tmask)
       vmax = global_maxval(blh_data(:,:,1,:),distrb_info,tmask)
       if (my_task.eq.master_task) write(nu_diag,*) subname, '  blh_data(slot1)', vmin, vmax

       vmin = global_minval(windgust_data(:,:,1,:),distrb_info,tmask)
       vmax = global_maxval(windgust_data(:,:,1,:),distrb_info,tmask)
       if (my_task.eq.master_task) write(nu_diag,*) subname, '  windgust_data(slot1)', vmin, vmax

    endif

  end subroutine ERA5_data

   !=======================================================================
   subroutine AFIM_data(dt)
      ! Reads daily ocean reanalysis ("sst", "sss", "u", and "v") from yearly
      ! NetCDF files and applies to CICE model fields with optional SST restoring.
      !
      ! Features:
      ! - Daily forcing with linear interpolation between days
      ! - SST restoring controlled by trestore (days) from ice_in namelist
      ! - NaN handling for reanalysis data (replaces NaN with sensible flags)
      ! - Year cycling via fyear_init and ycycle namelist parameters
      !
      ! Required namelist settings (forcing_nml in ice_in):
      !   ocn_data_type = 'AFIM'
      !   ocn_data_dir  = '/path/to/daily/sfc/'
      !   restore_ocn   = .true.    ! enable SST restoring
      !   trestore      = 1         ! restoring timescale in days
      !   fyear_init    = 1993      ! first year of forcing data
      !   ycycle        = 4         ! number of years in forcing cycle
      !
      ! Expected file naming: {ocn_data_dir}/ORAS_{YYYY}.nc
      ! Expected variables: sst(time,nj,ni), sss(time,nj,ni), u(time,nj,ni), v(time,nj,ni)
      !
      ! Author: DP@H2O, University of Tasmania
      !
      use ice_blocks, only            : nx_block, ny_block
      use ice_read_write, only        : ice_read_nc, ice_open_nc, ice_close_nc
      use ice_global_reductions, only : global_minval, global_maxval
      use ice_domain, only            : distrb_info, nblocks
      use ice_domain_size, only       : max_blocks
      use ice_state, only             : aice
      use ice_flux, only              : sss, sst, uocn, vocn, Tf
      use ice_grid, only              : hm, tmask, umask
      use ice_calendar, only          : days_per_year, mmonth, myear, yday, msec
      use ice_constants, only         : c0, c1, p5
      use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
      implicit none
      real(kind=dbl_kind), intent(in) :: dt
      ! Local variables
      integer(kind=int_kind)          :: fid, i, j, iblk
      integer(kind=int_kind)          :: recnum, maxrec, fyear, modadj
      real(kind=dbl_kind)             :: secday, vmin, vmax
      real(kind=dbl_kind)             :: val
      ! Temporary arrays for current and previous time records
      real(kind=dbl_kind), dimension(nx_block, ny_block, max_blocks) :: &
         sst_curr, sst_prev, sss_curr, sss_prev, &
         uocn_curr, uocn_prev, vocn_curr, vocn_prev
      logical, save       :: first_call = .true.
      real(kind=dbl_kind) :: sst_interp, sss_interp, u_interp, v_interp
      ! Thresholds for sanity checks
      real(kind=dbl_kind), parameter  :: sst_thresh      = 54.0_dbl_kind
      real(kind=dbl_kind), parameter  :: sss_high_thresh = 47.0_dbl_kind
      real(kind=dbl_kind), parameter  :: sss_low_thresh  =  5.0_dbl_kind
      real(kind=dbl_kind), parameter  :: ispd_thresh     =  2.5_dbl_kind
      ! NaN fallback values
      real(kind=dbl_kind), parameter  :: sst_fallback = -1.8_dbl_kind
      real(kind=dbl_kind), parameter  :: sss_fallback = 34.0_dbl_kind
      real(kind=dbl_kind), parameter  :: uv_fallback  =  0.0_dbl_kind
      ! Very large number to detect Inf
      real(kind=dbl_kind), parameter  :: huge_val = 1.0e30_dbl_kind
      integer(kind=int_kind), save    :: last_fyear = -999999
      character(len=*), parameter     :: subname = '(AFIM_data)'
      !-------------------------------------------------------------------
      ! 1) what time is it?
      !-------------------------------------------------------------------
      call icepack_query_parameters(secday_out=secday)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, file=__FILE__, line=__LINE__)
      if (ycycle == 0) then
         fyear = myear
      else
         modadj = abs((min(0, myear - fyear_init) / ycycle + 1) * ycycle)
         fyear  = fyear_init + mod(myear - fyear_init + modadj, ycycle)
      endif
      !-------------------------------------------------------------------
      ! 2) update forcing file if new year
      !-------------------------------------------------------------------
      if (fyear /= last_fyear) then
         call AFIM_files(fyear)
         last_fyear = fyear
         if (debug_forcing .or. local_debug) then
            if (my_task == master_task) then
               write(nu_diag,*) subname, ' myear, fyear, mmonth:', myear, fyear, mmonth
            endif 
         endif
      endif
      !-------------------------------------------------------------------
      ! 3) record number ... assumes forcing is daily
      !-------------------------------------------------------------------
      maxrec = days_per_year
      if (mod(fyear, 4) == 0) maxrec = 366
      if (mod(fyear, 100) == 0 .and. mod(fyear, 400) /= 0) maxrec = 365
      recnum = int(yday)
      if (recnum < 1) recnum = 1
      if (recnum > maxrec) recnum = maxrec
      call interp_coeff(recnum, 2, secday, 2)
      if (debug_forcing .or. local_debug) then
         if (my_task == master_task) then
            write(nu_diag,*) subname, ' yday, recnum, maxrec:', yday, recnum, maxrec
            write(nu_diag,*) subname, ' c1intp, c2intp      :', c1intp, c2intp
         endif
      endif 
      !-------------------------------------------------------------------
      ! 4) read-in forcing fields ... previous or current time step
      !-------------------------------------------------------------------
      call ice_open_nc(F_AFIM, fid)
      call ice_read_nc(fid, recnum, 'sst', sst_curr, debug_forcing, field_loc=field_loc_center, field_type=field_type_scalar)
      call ice_read_nc(fid, recnum, 'sss', sss_curr, debug_forcing, field_loc=field_loc_center, field_type=field_type_scalar)
      call ice_read_nc(fid, recnum, 'u', uocn_curr, debug_forcing,  field_loc=field_loc_center, field_type=field_type_scalar)
      call ice_read_nc(fid, recnum, 'v', vocn_curr, debug_forcing, field_loc=field_loc_center, field_type=field_type_scalar)
      if (recnum > 1) then
         call ice_read_nc(fid, recnum-1, 'sst', sst_prev, debug_forcing, field_loc=field_loc_center, field_type=field_type_scalar)
         call ice_read_nc(fid, recnum-1, 'sss', sss_prev, debug_forcing, field_loc=field_loc_center, field_type=field_type_scalar)
         call ice_read_nc(fid, recnum-1, 'u', uocn_prev, debug_forcing,  field_loc=field_loc_center, field_type=field_type_scalar)
         call ice_read_nc(fid, recnum-1, 'v', vocn_prev, debug_forcing,  field_loc=field_loc_center, field_type=field_type_scalar)
      else
         sst_prev  = sst_curr
         sss_prev  = sss_curr
         uocn_prev = uocn_curr
         vocn_prev = vocn_curr
      endif
      call ice_close_nc(fid)
      !-------------------------------------------------------------------
      ! 5) deal with NaN/Inf in temporary arrays using ieee_is_nan
      !-------------------------------------------------------------------
      do iblk = 1, nblocks
         do j = 1, ny_block
            do i = 1, nx_block
               ! Current record
               if (ieee_is_nan(sst_curr(i,j,iblk)) .or. abs(sst_curr(i,j,iblk)) > huge_val) then
                  sst_curr(i,j,iblk) = sst_fallback
               endif
               if (ieee_is_nan(sss_curr(i,j,iblk)) .or. abs(sss_curr(i,j,iblk)) > huge_val) then
                  sss_curr(i,j,iblk) = sss_fallback
               endif
               if (ieee_is_nan(uocn_curr(i,j,iblk)) .or. abs(uocn_curr(i,j,iblk)) > huge_val) then
                  uocn_curr(i,j,iblk) = uv_fallback
               endif
               if (ieee_is_nan(vocn_curr(i,j,iblk)) .or. abs(vocn_curr(i,j,iblk)) > huge_val) then
                  vocn_curr(i,j,iblk) = uv_fallback
               endif
               ! Previous record
               if (ieee_is_nan(sst_prev(i,j,iblk)) .or. abs(sst_prev(i,j,iblk)) > huge_val) then
                  sst_prev(i,j,iblk) = sst_fallback
               endif
               if (ieee_is_nan(sss_prev(i,j,iblk)) .or. abs(sss_prev(i,j,iblk)) > huge_val) then
                  sss_prev(i,j,iblk) = sss_fallback
               endif
               if (ieee_is_nan(uocn_prev(i,j,iblk)) .or. abs(uocn_prev(i,j,iblk)) > huge_val) then
                  uocn_prev(i,j,iblk) = uv_fallback
               endif
               if (ieee_is_nan(vocn_prev(i,j,iblk)) .or. abs(vocn_prev(i,j,iblk)) > huge_val) then
                  vocn_prev(i,j,iblk) = uv_fallback
               endif
            enddo
         enddo
      enddo
      do iblk = 1, nblocks
         do j = 1, ny_block
            do i = 1, nx_block
               if (ieee_is_nan(sst(i,j,iblk)) .or. abs(sst(i,j,iblk)) > huge_val) then
                  sst(i,j,iblk) = c1intp * sst_prev(i,j,iblk) + c2intp * sst_curr(i,j,iblk)
               endif
            enddo
         enddo
      enddo
      if (debug_forcing .or. local_debug) then
         if (my_task == master_task) write(nu_diag,*) subname, ' global min/max after "NaN-fix in SST"'
         vmin = global_minval(sst, distrb_info, tmask)
         vmax = global_maxval(sst, distrb_info, tmask)
         if (my_task == master_task) write(nu_diag,*) subname, '  sst  ', vmin, vmax
      endif
      !-------------------------------------------------------------------
      ! 6) apply forcing using openmpi 
      !-------------------------------------------------------------------
      if (debug_forcing .or. local_debug) then
         if (my_task == master_task) write(nu_diag,*) subname, ' global min/max before "APPLY"'
         vmin = global_minval(sst, distrb_info, tmask)
         vmax = global_maxval(sst, distrb_info, tmask)
         if (my_task == master_task) write(nu_diag,*) subname, '  sst  ', vmin, vmax
      endif
      !$OMP PARALLEL DO PRIVATE(iblk, j, i, sst_interp, sss_interp, u_interp, v_interp)
      do iblk = 1, nblocks
         do j = 1, ny_block
            do i = 1, nx_block
               ! Interpolate forcing data
               sst_interp = c1intp * sst_prev(i,j,iblk) + c2intp * sst_curr(i,j,iblk)
               sss_interp = c1intp * sss_prev(i,j,iblk) + c2intp * sss_curr(i,j,iblk)
               u_interp   = c1intp * uocn_prev(i,j,iblk) + c2intp * uocn_curr(i,j,iblk)
               v_interp   = c1intp * vocn_prev(i,j,iblk) + c2intp * vocn_curr(i,j,iblk)
               ! CORRECTED MASK LOGIC: hm > p5 = OCEAN, else = LAND
               if (hm(i,j,iblk) > p5) then
                  ! OCEAN point
                  ! SST: Use direct assignment on first call, restoring afterwards
                  if (first_call) then
                     ! First call: directly assign forcing SST
                     sst(i,j,iblk) = sst_interp
                  elseif (restore_ocn) then
                     ! Subsequent calls with restoring: gradual adjustment
                     sst(i,j,iblk) = sst(i,j,iblk) + (sst_interp - sst(i,j,iblk)) * dt / trest
                  else
                     ! No restoring: direct assignment
                     sst(i,j,iblk) = sst_interp
                  endif
                  ! SSS, uocn, vocn: always direct assignment
                  sss(i,j,iblk) = sss_interp
                  uocn(i,j,iblk) = u_interp
                  vocn(i,j,iblk) = v_interp
                  ! Upper bound for SST
                  if (sst(i,j,iblk) > sst_thresh) sst(i,j,iblk) = sst_thresh
                  ! SSS bounds
                  if (sss(i,j,iblk) < sss_low_thresh)  sss(i,j,iblk) = sss_low_thresh
                  if (sss(i,j,iblk) > sss_high_thresh) sss(i,j,iblk) = sss_high_thresh
                  ! Velocity bounds
                  if (abs(uocn(i,j,iblk)) > ispd_thresh) uocn(i,j,iblk) = c0
                  if (abs(vocn(i,j,iblk)) > ispd_thresh) vocn(i,j,iblk) = c0
               else
                  ! LAND point
                  sst (i,j,iblk) = c0
                  sss (i,j,iblk) = c0
                  uocn(i,j,iblk) = c0
                  vocn(i,j,iblk) = c0
               endif
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
      ! After first call, disable direct assignment
      first_call = .false.
      !-------------------------------------------------------------------
      ! 7) update freezing temperature and write out globals to diagnostic
      !-------------------------------------------------------------------
      call ocn_freezing_temperature
      if (debug_forcing .or. local_debug) then
         if (my_task == master_task) write(nu_diag,*) subname, ' global min/max at "END" of subroutine'
         vmin = global_minval(sst, distrb_info, tmask)
         vmax = global_maxval(sst, distrb_info, tmask)
         if (my_task == master_task) write(nu_diag,*) subname, '  sst  ', vmin, vmax
         vmin = global_minval(sss, distrb_info, tmask)
         vmax = global_maxval(sss, distrb_info, tmask)
         if (my_task == master_task) write(nu_diag,*) subname, '  sss  ', vmin, vmax
         vmin = global_minval(uocn, distrb_info, umask)
         vmax = global_maxval(uocn, distrb_info, umask)
         if (my_task == master_task) write(nu_diag,*) subname, '  uocn ', vmin, vmax
         vmin = global_minval(vocn, distrb_info, umask)
         vmax = global_maxval(vocn, distrb_info, umask)
         if (my_task == master_task) write(nu_diag,*) subname, '  vocn ', vmin, vmax
      endif
   end subroutine AFIM_data

  !=======================================================================
  subroutine compute_shortwave(nx_block, ny_block, ilo, ihi, jlo, jhi, TLON, TLAT, hm, Qa, cldf, fsw)
    ! AOMIP shortwave forcing
    ! standard calculation using solar declination angle
    ! then shortwave is reduced using a function of cloud fraction
    integer (kind=int_kind), intent(in) :: nx_block, ny_block, & ! block dimensions
         ilo,ihi,jlo,jhi       ! beginning and end of physical domain
    real (kind=dbl_kind), dimension(nx_block,ny_block), intent(in) :: TLON, TLAT     , & ! longitude, latitude
         Qa             , & ! specific humidity
         cldf           , & ! cloud fraction
         hm                 ! land mask
    real (kind=dbl_kind), dimension(nx_block,ny_block), intent(inout) :: fsw             ! shortwave
    real (kind=dbl_kind) :: hour_angle, &
         solar_time, &
         declin    , &
         cosZ      , &
         e, d      , &
         sw0       , &
         secday    , &
         pi        , &
         lontmp    , &
         deg2rad
    integer (kind=int_kind) :: i, j
    character(len=*), parameter :: subname = '(compute_shortwave)'
    if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'
    call icepack_query_parameters(secday_out=secday, pi_out=pi)
    call icepack_warnings_flush(nu_diag)
    if (icepack_warnings_aborted()) call abort_ice(error_message=subname, file=__FILE__, line=__LINE__)
    do j=jlo,jhi
       do i=ilo,ihi
          deg2rad = pi/c180
          ! Convert longitude to range of -180 to 180 for LST calculation
          lontmp = mod(TLON(i,j)/deg2rad,c360)
          if (lontmp .gt. c180) lontmp = lontmp - c360
          if (lontmp .lt. -c180) lontmp = lontmp + c360
          solar_time = mod(real(msec,kind=dbl_kind),secday)/c3600 + lontmp/c15
          if (solar_time .ge. 24._dbl_kind) solar_time = solar_time - 24._dbl_kind
          hour_angle = (c12 - solar_time)*pi/c12
          ! use dayyr instead of c365???
          declin = 23.44_dbl_kind*cos((172._dbl_kind-yday) * c2*pi/c365)*deg2rad 
          cosZ   = sin(TLAT(i,j))*sin(declin) + cos(TLAT(i,j))*cos(declin)*cos(hour_angle)
          cosZ   = max(cosZ,c0)
          e      = 1.e5*Qa(i,j)/(0.622_dbl_kind + 0.378_dbl_kind*Qa(i,j))
          d      = (cosZ+2.7_dbl_kind)*e*1.e-5_dbl_kind+1.085_dbl_kind*cosZ+p1
          sw0    = 1353._dbl_kind*cosZ**2/d
          sw0    = max(sw0,c0)
          ! total downward shortwave for cice
          Fsw(i,j) = sw0*(c1-p6*cldf(i,j)**3)
          Fsw(i,j) = Fsw(i,j)*hm(i,j)
       enddo
    enddo
  end subroutine compute_shortwave

    !=======================================================================
    subroutine Qa_fixLY(nx_block, ny_block, Tair, Qa)
      ! prevents humidity from being super-saturated
      integer (kind=int_kind), intent(in) :: nx_block, ny_block ! block dimensions
      ! air temperature
      real (kind=dbl_kind), dimension(nx_block,ny_block), intent(in) :: Tair
      ! specific humidity
      real (kind=dbl_kind), dimension(nx_block,ny_block), intent(inout) :: Qa
      real (kind=dbl_kind), dimension (nx_block,ny_block) :: worka
      real (kind=dbl_kind) :: Tffresh, puny
      character(len=*), parameter :: subname = '(Qa_fixLY)'
      if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'
      call icepack_query_parameters(Tffresh_out=Tffresh, puny_out=puny)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, file=__FILE__, line=__LINE__)
      worka = Tair - Tffresh
      worka = c2 + (0.7859_dbl_kind + 0.03477_dbl_kind*worka)/(c1 + 0.00412_dbl_kind*worka) + 0.00422_dbl_kind*worka
      ! vapor pressure
      worka = (c10**worka)      ! saturated
      worka = max(worka,puny)   ! puny over land to prevent division by zero
      ! specific humidity
      worka = 0.622_dbl_kind*worka/(1.e5_dbl_kind-0.378_dbl_kind*worka)
      Qa = min(Qa, worka)
    end subroutine Qa_fixLY

!=======================================================================
    subroutine hadgem_files (yr)
      ! HadGEM or HadGAM atmospheric forcing
      ! Construct filenames based on selected model options
      !! Note: The year number in these filenames does not matter, because
      !       subroutine file_year will insert the correct year.
      ! author: Alison McLaren, Met Office
      ! current forcing year
      integer (kind=int_kind), intent(in) :: yr
      ! thickness category index
      integer (kind=int_kind) :: n
      logical (kind=log_kind) :: calc_strair, calc_Tsfc
      character(len=*), parameter :: subname = '(hadgem_files)'
      if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'
      call icepack_query_parameters(calc_strair_out=calc_strair, calc_Tsfc_out=calc_Tsfc)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, file=__FILE__, line=__LINE__)
      ! -----------------------------------------------------------
      ! Rainfall and snowfall
      ! -----------------------------------------------------------
      F_snow = trim(atm_data_dir)//'/MONTHLY/snowfall.1996.nc'
      call file_year_atm(F_snow,yr)
      F_rain = trim(atm_data_dir)//'/MONTHLY/rainfall.1996.nc'
      call file_year_atm(F_rain,yr)
      if (my_task == master_task) then
         write (nu_diag,*) ' '
         write (nu_diag,*) 'Atmospheric data files:'
         write (nu_diag,*) trim(F_rain)
         write (nu_diag,*) trim(F_snow)
      endif
      if (calc_strair) then
         ! --------------------------------------------------------
         ! Wind velocity
         ! --------------------------------------------------------
         F_uwind = trim(atm_data_dir)//'/MONTHLY/u_10.1996.nc'
         call file_year_atm(F_uwind,yr)
         F_vwind = trim(atm_data_dir)//'/MONTHLY/v_10.1996.nc'
         call file_year_atm(F_vwind,yr)
         if (my_task == master_task) then
            write (nu_diag,*) trim(F_uwind)
            write (nu_diag,*) trim(F_vwind)
         endif
      else
         ! --------------------------------------------------------
         ! Wind stress
         ! --------------------------------------------------------
         F_strax = trim(atm_data_dir)//'/MONTHLY/taux.1996.nc'
         call file_year_atm(F_strax,yr)
         F_stray = trim(atm_data_dir)//'/MONTHLY/tauy.1996.nc'
         call file_year_atm(F_stray,yr)
         if (my_task == master_task) then
            write (nu_diag,*) trim(F_strax)
            write (nu_diag,*) trim(F_stray)
         endif
         if (calc_Tsfc .or. oceanmixed_ice) then
            ! --------------------------------------------------
            ! Wind speed
            ! --------------------------------------------------
            F_wind = trim(atm_data_dir)//'/MONTHLY/wind_10.1996.nc'
            call file_year_atm(F_wind,yr)
            if (my_task == master_task) then
               write (nu_diag,*) trim(F_wind)
            endif
         endif   ! calc_Tsfc or oceanmixed_ice
      endif  ! calc_strair
      ! --------------------------------------------------------------
      ! Atmosphere properties.  Even if these fields are not
      ! being used to force the ice (i.e. calc_Tsfc=.false.), they
      ! are still needed to generate forcing for mixed layer model or
      ! to calculate wind stress
      ! --------------------------------------------------------------
       if (calc_Tsfc .or. oceanmixed_ice .or. calc_strair) then
          F_fsw = trim(atm_data_dir)//'/MONTHLY/SW_incoming.1996.nc'
          call file_year_atm(F_fsw,yr)
          F_flw = trim(atm_data_dir)//'/MONTHLY/LW_incoming.1996.nc'
          call file_year_atm(F_flw,yr)
          F_tair = trim(atm_data_dir)//'/MONTHLY/t_10.1996.nc'
          call file_year_atm(F_tair,yr)
          F_humid = trim(atm_data_dir)//'/MONTHLY/q_10.1996.nc'
          call file_year_atm(F_humid,yr)
          F_rhoa = trim(atm_data_dir)//'/MONTHLY/rho_10.1996.nc'
          call file_year_atm(F_rhoa,yr)
          if (my_task == master_task) then
             write (nu_diag,*) trim(F_fsw)
             write (nu_diag,*) trim(F_flw)
             write (nu_diag,*) trim(F_tair)
             write (nu_diag,*) trim(F_humid)
             write (nu_diag,*) trim(F_rhoa)
          endif ! master_task
       endif ! calc_Tsfc or oceanmixed_ice  or calc_strair
       if (.not. calc_Tsfc) then
          ! ------------------------------------------------------
          ! Sublimation, topmelt and botmelt
          ! ------------------------------------------------------
          do n = 1, ncat
             ! 'topmelt' = fsurf - fcondtop.
             write(topmelt_file(n), '(a,i1,a)') trim(atm_data_dir)//'/MONTHLY/topmeltn',n,'.1996.nc'
             call file_year_atm(topmelt_file(n),yr)
             ! 'botmelt' = fcondtop.
             write(botmelt_file(n), '(a,i1,a)') trim(atm_data_dir)//'/MONTHLY/botmeltn',n,'.1996.nc'
             call file_year_atm(botmelt_file(n),yr)
          enddo
          ! 'sublim' = - flat / Lsub.
          F_sublim = trim(atm_data_dir)//'/MONTHLY/sublim.1996.nc'
          call file_year_atm(F_sublim,yr)
          if (my_task == master_task) then
             do n = 1, ncat
                write (nu_diag,*) trim(topmelt_file(n))
                write (nu_diag,*) trim(botmelt_file(n))
             enddo
             write (nu_diag,*) trim(F_sublim)
          endif
       endif  ! .not. calc_Tsfc
     end subroutine hadgem_files

     !=======================================================================
     subroutine hadgem_data
       ! read HadGEM or HadGAM atmospheric data
       ! authors: Alison McLaren, Met Office
       use ice_domain, only: nblocks
       use ice_flux, only: fsnow, frain, uatm, vatm, strax, stray, wind, fsw, flw, Tair, rhoa, Qa, fcondtopn_f, fsurfn_f, flatn_f
       integer (kind=int_kind) :: i, j        , & ! horizontal indices
                                  n           , & ! thickness category index
                                  iblk        , & ! block index
                                  ixm,ixp     , & ! record numbers for neighboring months
                                  maxrec      , & ! maximum record number
                                  recslot     , & ! spline slot for current record
                                  midmonth        ! middle day of month

      logical (kind=log_kind) :: readm
      real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks) :: topmelt, & ! temporary fields
                                                                       botmelt, &
                                                                       sublim
      ! field name in netcdf file
      character (char_len) :: fieldname
      real (kind=dbl_kind) :: Lsub
      logical (kind=log_kind) :: calc_strair, calc_Tsfc
      character(len=*), parameter :: subname = '(hadgem_data)'
      if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'
      call icepack_query_parameters(Lsub_out=Lsub)
      call icepack_query_parameters(calc_strair_out=calc_strair, calc_Tsfc_out=calc_Tsfc)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, file=__FILE__, line=__LINE__)
      !-------------------------------------------------------------------
      ! monthly data
      ! Assume that monthly data values are located in the middle of the
      ! month.
      !-------------------------------------------------------------------
      midmonth = 15  ! data is given on 15th of every month
      !      midmonth = fix(p5 * real(daymo(mmonth)))  ! exact middle
      ! Compute record numbers for surrounding months
      maxrec = 12
      ixm  = mod(mmonth+maxrec-2,maxrec) + 1
      ixp  = mod(mmonth,         maxrec) + 1
      if (mday >= midmonth) ixm = -99  ! other two points will be used
      if (mday <  midmonth) ixp = -99
      ! Determine whether interpolation will use values 1:2 or 2:3
      ! recslot = 2 means we use values 1:2, with the current value (2)
      !  in the second slot
      ! recslot = 1 means we use values 2:3, with the current value (2)
      !  in the first slot
      recslot = 1                             ! latter half of month
      if (mday < midmonth) recslot = 2        ! first half of month
      ! Find interpolation coefficients
      call interp_coeff_monthly (recslot)
      ! Read 2 monthly values
      readm = .false.
      if (istep==1 .or. (mday==midmonth .and. msec==0)) readm = .true.
      ! -----------------------------------------------------------
      ! Rainfall and snowfall
      ! -----------------------------------------------------------
      fieldname='rainfall'
      call read_data_nc (readm, 0, fyear, ixm, mmonth, ixp, maxrec, F_rain, fieldname, frain_data, field_loc_center, field_type_scalar)
      fieldname='snowfall'
      call read_data_nc (readm, 0, fyear, ixm, mmonth, ixp, maxrec, F_snow, fieldname, fsnow_data, field_loc_center, field_type_scalar)
      ! Interpolate to current time step
      call interpolate_data (fsnow_data, fsnow)
      call interpolate_data (frain_data, frain)
      if (calc_strair) then
         ! --------------------------------------------------------
         ! Wind velocity
         ! --------------------------------------------------------
         fieldname='u_10'
         call read_data_nc (readm, 0, fyear, ixm, mmonth, ixp, maxrec, F_uwind, fieldname, uatm_data, field_loc_center, field_type_vector)
         fieldname='v_10'
         call read_data_nc (readm, 0, fyear, ixm, mmonth, ixp, maxrec, F_vwind, fieldname, vatm_data, field_loc_center, field_type_vector)
         ! Interpolate to current time step
         call interpolate_data (uatm_data, uatm)
         call interpolate_data (vatm_data, vatm)
      else
         ! --------------------------------------------------------
         ! Wind stress
         ! --------------------------------------------------------
         fieldname='taux'
         call read_data_nc (readm, 0, fyear, ixm, mmonth, ixp, maxrec, F_strax, fieldname, strax_data, field_loc_center, field_type_vector)
         fieldname='tauy'
         call read_data_nc (readm, 0, fyear, ixm, mmonth, ixp, maxrec, F_stray, fieldname, stray_data, field_loc_center, field_type_vector)
         ! Interpolate to current time step
         call interpolate_data (strax_data, strax)
         call interpolate_data (stray_data, stray)
         if (calc_Tsfc .or. oceanmixed_ice) then
            ! --------------------------------------------------
            ! Wind speed
            ! --------------------------------------------------
            fieldname='wind_10'
            call read_data_nc (readm, 0, fyear, ixm, mmonth, ixp, maxrec, F_wind, fieldname, wind_data, field_loc_center, field_type_scalar)
            ! Interpolate to current time step
            call interpolate_data (wind_data, wind)
         endif   ! calc_Tsfc or oceanmixed_ice
      endif      ! calc_strair
      ! -----------------------------------------------------------
      ! SW incoming, LW incoming, air temperature, density and humidity at 10m.
      ! Even if these fields are not being used to force the ice
      ! (i.e. calc_Tsfc=.false.), they are still needed to generate
      ! forcing for mixed layer model or to calculate wind stress
      ! -----------------------------------------------------------
      if (calc_Tsfc .or. oceanmixed_ice .or. calc_strair) then
         fieldname='SW_incoming'
         call read_data_nc (readm, 0, fyear, ixm, mmonth, ixp, maxrec, F_fsw, fieldname, fsw_data, field_loc_center, field_type_scalar)
         fieldname='LW_incoming'
         call read_data_nc (readm, 0, fyear, ixm, mmonth, ixp, maxrec, F_flw, fieldname, flw_data, field_loc_center, field_type_scalar)
         fieldname='t_10'
         call read_data_nc (readm, 0, fyear, ixm, mmonth, ixp, maxrec, F_tair, fieldname, Tair_data, field_loc_center, field_type_scalar)
         fieldname='rho_10'
         call read_data_nc (readm, 0, fyear, ixm, mmonth, ixp, maxrec, F_rhoa, fieldname, rhoa_data, field_loc_center, field_type_scalar)
         fieldname='q_10'
         call read_data_nc (readm, 0, fyear, ixm, mmonth, ixp, maxrec, F_humid, fieldname, Qa_data, field_loc_center, field_type_scalar)
         ! Interpolate onto current timestep
         call interpolate_data (fsw_data,   fsw)
         call interpolate_data (flw_data,  flw)
         call interpolate_data (Tair_data, Tair)
         call interpolate_data (rhoa_data, rhoa)
         call interpolate_data (Qa_data,   Qa)
      endif       ! calc_Tsfc or oceanmixed_ice or calc_strair
      if (.not. calc_Tsfc) then
         ! ------------------------------------------------------
         ! Sublimation, topmelt and botmelt
         ! ------------------------------------------------------
         fieldname='sublim'
         call read_data_nc (readm, 0, fyear, ixm, mmonth, ixp, maxrec, F_sublim, fieldname, sublim_data, field_loc_center, field_type_scalar)
         ! Interpolate to current time step
         call interpolate_data (sublim_data, sublim)
         do n = 1, ncat
            write(fieldname, '(a,i1)') 'topmeltn',n
            call read_data_nc (readm, 0, fyear, ixm, mmonth, ixp, maxrec, topmelt_file(n), fieldname, topmelt_data(:,:,:,:,n), field_loc_center, field_type_scalar)
            write(fieldname, '(a,i1)') 'botmeltn',n
            call read_data_nc (readm, 0, fyear, ixm, mmonth, ixp, maxrec, botmelt_file(n), fieldname, botmelt_data(:,:,:,:,n), field_loc_center, field_type_scalar)
            call interpolate_data (topmelt_data(:,:,:,:,n), topmelt)
            call interpolate_data (botmelt_data(:,:,:,:,n), botmelt)
            !--------------------------------------------------------
            ! Convert from UM variables to CICE variables
            !  topmelt = fsurf - fcondtop
            !  botmelt = fcondtop  (as zero layer)
            ! Convert UM sublimation data into CICE LH flux
            ! (sublim = - flatn / Lsub) and have same value for all
            ! categories
            !--------------------------------------------------------
            !$OMP PARALLEL DO PRIVATE(iblk,i,j)
            do iblk = 1, nblocks
               do j = 1, ny_block
                  do i = 1, nx_block
                     fcondtopn_f(i,j,n,iblk) = botmelt(i,j,iblk)
                     fsurfn_f(i,j,n,iblk)    = topmelt(i,j,iblk) + botmelt(i,j,iblk)
                     flatn_f(i,j,n,iblk)    = - sublim(i,j,iblk)*Lsub
                  enddo
               enddo
            enddo !$OMP END PARALLEL DO
         enddo  ! ncat
      endif   ! .not. calc_Tsfc
    end subroutine hadgem_data

    !=======================================================================
    subroutine monthly_files (yr)
      ! Construct filenames based on the LANL naming conventions for NCAR data.
      ! Edit for other directory structures or filenames.
      ! Note: The year number in these filenames does not matter, because
      !       subroutine file_year will insert the correct year.
      ! author: Elizabeth C. Hunke, LANL
      integer (kind=int_kind), intent(in) :: yr ! current forcing year
      character(len=*), parameter :: subname = '(monthly_files)'
      if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'
      F_flw   = trim(atm_data_dir)//'/MONTHLY/cldf.omip.dat'
      F_rain  = trim(atm_data_dir)//'/MONTHLY/prec.nmyr.dat'
      F_tair  = trim(atm_data_dir)//'/MONTHLY/t_10.1996.dat'
      call file_year_atm(F_tair,yr)
      F_humid = trim(atm_data_dir)//'/MONTHLY/q_10.1996.dat'
      call file_year_atm(F_humid,yr)
      ! stress/speed is used instead of wind components
      F_strax = trim(atm_data_dir)//'/MONTHLY/strx.1996.dat'
      call file_year_atm(F_strax,yr)
      F_stray = trim(atm_data_dir)//'/MONTHLY/stry.1996.dat'
      call file_year_atm(F_stray,yr)
      F_wind  = trim(atm_data_dir)//'/MONTHLY/wind.1996.dat'
      call file_year_atm(F_wind,yr)
      if (my_task == master_task) then
         write (nu_diag,*) ' '
         write (nu_diag,*) 'Forcing data year = ', fyear
         write (nu_diag,*) 'Atmospheric data files:'
         write (nu_diag,*) trim(F_flw)
         write (nu_diag,*) trim(F_rain)
         write (nu_diag,*) trim(F_tair)
         write (nu_diag,*) trim(F_humid)
         write (nu_diag,*) trim(F_uwind)
         write (nu_diag,*) trim(F_vwind)
      endif                     ! master_task
    end subroutine monthly_files

    !=======================================================================
    subroutine monthly_data
      ! read monthly atmospheric data
      use ice_blocks, only: block, get_block
      use ice_global_reductions, only: global_minval, global_maxval
      use ice_domain, only: nblocks, distrb_info, blocks_ice
      use ice_flux, only: fsnow, Tair, Qa, wind, strax, stray, fsw
      use ice_grid, only: hm, tlon, tlat, tmask, umask
      integer (kind=int_kind) :: i, j        , &
                                 ixm,ixp     , & ! record numbers for neighboring months
                                 maxrec      , & ! maximum record number
                                 recslot     , & ! spline slot for current record
                                 midmonth    , & ! middle day of month
                                 iblk        , & ! block index
                                 ilo,ihi,jlo,jhi ! beginning and end of physical domain
      real (kind=dbl_kind) :: vmin, vmax
      logical (kind=log_kind) :: readm
      type (block) :: this_block           ! block information for current block
      character(len=*), parameter :: subname = '(monthly_data)'
      if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'
      !-------------------------------------------------------------------
      ! Assume that monthly data values are located in the middle of the month.
      !-------------------------------------------------------------------
      midmonth = 15  ! data is given on 15th of every month
      !      midmonth = fix(p5 * real(daymo(mmonth)))  ! exact middle
      ! Compute record numbers for surrounding months
      maxrec = 12
      ixm  = mod(mmonth+maxrec-2,maxrec) + 1
      ixp  = mod(mmonth,         maxrec) + 1
      if (mday >= midmonth) ixm = -99  ! other two points will be used
      if (mday <  midmonth) ixp = -99
      ! Determine whether interpolation will use values 1:2 or 2:3
      ! recslot = 2 means we use values 1:2, with the current value (2)
      !  in the second slot
      ! recslot = 1 means we use values 2:3, with the current value (2)
      !  in the first slot
      recslot = 1                             ! latter half of month
      if (mday < midmonth) recslot = 2        ! first half of month
      ! Find interpolation coefficients
      call interp_coeff_monthly (recslot)
      ! Read 2 monthly values
      readm = .false.
      if (istep==1 .or. (mday==midmonth .and. msec==0)) readm = .true.
      call read_clim_data (readm, 0, ixm, mmonth, ixp, F_flw, cldf_data, field_loc_center, field_type_scalar)
      call read_clim_data (readm, 0, ixm, mmonth, ixp, F_rain, fsnow_data, field_loc_center, field_type_scalar)
      call read_clim_data (readm, 0, ixm, mmonth, ixp, F_tair, Tair_data, field_loc_center, field_type_scalar)
      call read_clim_data (readm, 0, ixm, mmonth, ixp, F_humid, Qa_data, field_loc_center, field_type_scalar)
      call read_clim_data (readm, 0, ixm, mmonth, ixp, F_wind, wind_data, field_loc_center, field_type_scalar)
      call read_clim_data (readm, 0, ixm, mmonth, ixp, F_strax, strax_data, field_loc_center, field_type_vector)
      call read_clim_data (readm, 0, ixm, mmonth, ixp, F_stray, stray_data, field_loc_center, field_type_vector)
      call interpolate_data (cldf_data, cldf)
      call interpolate_data (fsnow_data, fsnow)  ! units mm/s = kg/m^2/s
      call interpolate_data (Tair_data, Tair)
      call interpolate_data (Qa_data, Qa)
      call interpolate_data (wind_data, wind)
      call interpolate_data (strax_data, strax)
      call interpolate_data (stray_data, stray)
      !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block)
      do iblk = 1, nblocks
         call Qa_fixLY(nx_block, ny_block, Tair(:,:,iblk), Qa(:,:,iblk))
         do j = 1, ny_block
            do i = 1, nx_block
               Qa   (i,j,iblk) = Qa   (i,j,iblk) * hm(i,j,iblk)
               Tair (i,j,iblk) = Tair (i,j,iblk) * hm(i,j,iblk)
               wind (i,j,iblk) = wind (i,j,iblk) * hm(i,j,iblk)
               strax(i,j,iblk) = strax(i,j,iblk) * hm(i,j,iblk)
               stray(i,j,iblk) = stray(i,j,iblk) * hm(i,j,iblk)
            enddo
         enddo
         ! AOMIP
         this_block = get_block(blocks_ice(iblk),iblk)
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi
         call compute_shortwave(nx_block, ny_block, ilo, ihi, jlo, jhi, TLON (:,:,iblk), &
                                                                        TLAT (:,:,iblk), &
                                                                        hm   (:,:,iblk), &
                                                                        Qa   (:,:,iblk), &
                                                                        cldf (:,:,iblk), &
                                                                        fsw  (:,:,iblk))
      enddo  ! iblk !$OMP END PARALLEL DO
      if (debug_forcing) then
         if (my_task == master_task) write (nu_diag,*) 'LY_bulk_data'
         vmin = global_minval(fsw,distrb_info,tmask)
         vmax = global_maxval(fsw,distrb_info,tmask)
         if (my_task.eq.master_task) write (nu_diag,*) 'fsw',vmin,vmax
         vmin = global_minval(cldf,distrb_info,tmask)
         vmax = global_maxval(cldf,distrb_info,tmask)
         if (my_task.eq.master_task) write (nu_diag,*) 'cldf',vmin,vmax
         vmin =global_minval(fsnow,distrb_info,tmask)
         vmax =global_maxval(fsnow,distrb_info,tmask)
         if (my_task.eq.master_task) write (nu_diag,*) 'fsnow',vmin,vmax
         vmin = global_minval(Tair,distrb_info,tmask)
         vmax = global_maxval(Tair,distrb_info,tmask)
         if (my_task.eq.master_task) write (nu_diag,*) 'Tair',vmin,vmax
         vmin = global_minval(wind,distrb_info,umask)
         vmax = global_maxval(wind,distrb_info,umask)
         if (my_task.eq.master_task) write (nu_diag,*) 'wind',vmin,vmax
         vmin = global_minval(strax,distrb_info,umask)
         vmax = global_maxval(strax,distrb_info,umask)
         if (my_task.eq.master_task) write (nu_diag,*) 'strax',vmin,vmax
         vmin = global_minval(stray,distrb_info,umask)
         vmax = global_maxval(stray,distrb_info,umask)
         if (my_task.eq.master_task) write (nu_diag,*) 'stray',vmin,vmax
         vmin = global_minval(Qa,distrb_info,tmask)
         vmax = global_maxval(Qa,distrb_info,tmask)
         if (my_task.eq.master_task) write (nu_diag,*) 'Qa',vmin,vmax
      endif                   ! debug_forcing
    end subroutine monthly_data

    !=======================================================================
    subroutine oned_data
      ! Oned atmospheric data
      use ice_flux, only: uatm, vatm, Tair, fsw, fsnow, Qa, rhoa, frain
      ! local parameters
      character (char_len_long) :: met_file, &      ! netcdf filename
                                   fieldname        ! field name in netcdf file
      integer (kind=int_kind)   :: fid              ! file id for netCDF file
      real (kind=dbl_kind)      :: work             ! temporary variable
      logical (kind=log_kind)   :: diag
      integer (kind=int_kind)   :: status           ! status flag
      real (kind=dbl_kind)      :: Temp, &          ! air temperature (K), used to determine specific humidity
                                   rh, &            ! relative humidity (%)
                                   Psat, &          ! saturation vapour pressure (hPa)
                                   ws               ! saturation mixing ratio
      ! coefficients for Hyland-Wexler Qa
      real (kind=dbl_kind), parameter :: ps1 = 0.58002206e4_dbl_kind,    & ! (K)
                                         ps2 = 1.3914993_dbl_kind,       & !
                                         ps3 = 0.48640239e-1_dbl_kind,   & ! (K^-1)
                                         ps4 = 0.41764768e-4_dbl_kind,   & ! (K^-2)
                                         ps5 = 0.14452093e-7_dbl_kind,   & ! (K^-3)
                                         ps6 = 6.5459673_dbl_kind,       & !
                                         ws1 = 621.97_dbl_kind,          & ! for saturation mixing ratio
                                         Pair = 1020._dbl_kind             ! Sea level pressure (hPa)
      character(len=*), parameter :: subname = '(oned_data)'
      if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'
      diag = .false.   ! write diagnostic information
      if (trim(atm_data_format) == 'nc') then     ! read nc file
        ! hourly data beginning Jan 1, 1989, 01:00
        ! HARDWIRED for dt = 1 hour!
        met_file = F_uwind
        call ice_open_nc(met_file,fid)
        fieldname='Uatm'
        call ice_read_nc(fid,istep1,fieldname,work,diag)
        uatm(:,:,:) = work
        fieldname='Vatm'
        call ice_read_nc(fid,istep1,fieldname,work,diag)
        vatm(:,:,:) = work
        fieldname='Tair'
        call ice_read_nc(fid,istep1,fieldname,work,diag)
        Temp = work
        Tair(:,:,:) = Temp
        call ice_close_nc(fid)
        ! hourly solar data beginning Jan 1, 1989, 01:00
        met_file = F_fsw
        call ice_open_nc(met_file,fid)
        fieldname='fsw'
        call ice_read_nc(fid,istep1,fieldname,work,diag)
        fsw(:,:,:) = work
        call ice_close_nc(fid)
        ! hourly interpolated monthly  data beginning Jan 1, 1989, 01:00
        met_file = F_humid
        call ice_open_nc(met_file,fid)
        fieldname='rh'
        call ice_read_nc(fid,istep1,fieldname,work,diag)
        rh = work
        fieldname='fsnow'
        call ice_read_nc(fid,istep1,fieldname,work,diag)
        fsnow(:,:,:) = work
        call ice_close_nc(fid)
        !-------------------------------------------------------------------
        ! Find specific humidity using Hyland-Wexler formulation
        ! Hyland, R.W. and A. Wexler, Formulations for the Thermodynamic
        ! Properties of the saturated phases of H20 from 173.15K to 473.15K,
        ! ASHRAE Trans, 89(2A), 500-519, 1983
        !-------------------------------------------------------------------
        ! saturation vapour pressure
        Psat = exp(-ps1/Temp + ps2 - ps3*Temp + ps4*Temp**2 - ps5 * Temp**3 + ps6 * log(Temp))*p01
        ! saturation mixing ratio
        ws = ws1 * Psat/(Pair - Psat)
        ! specific humidity (kg/kg)
        Qa(:,:,:) = rh * ws * p01/(c1 + rh * ws * p01) * p001
     endif ! atm_data_format
     ! flw calculated in prepare_forcing
     rhoa (:,:,:) = 1.3_dbl_kind ! air density (kg/m^3)
     cldf (:,:,:) = p25          ! cloud fraction
     frain(:,:,:) = c0           ! this is available in hourlymet_rh file
   end subroutine oned_data

   !=======================================================================
   subroutine oned_files
     character(len=*), parameter :: subname = '(oned_files)'
     if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'
     F_fsw   = trim(atm_data_dir)//'/hourlysolar_brw1989_5yr.nc'
     F_rain  = trim(atm_data_dir)//'/hourlymet_rh_5yr.nc'
     F_uwind = trim(atm_data_dir)//'/hourlymet_brw1989_5yr.nc'
     F_vwind = trim(atm_data_dir)//'/hourlymet_brw1989_5yr.nc'
     F_tair  = trim(atm_data_dir)//'/hourlymet_brw1989_5yr.nc'
     F_humid = trim(atm_data_dir)//'/hourlymet_rh_5yr.nc'
     if (my_task == master_task) then
        write (nu_diag,*) ' '
        write (nu_diag,*) 'Atmospheric data files:'
        write (nu_diag,*) trim(F_fsw)
        write (nu_diag,*) trim(F_rain)
        write (nu_diag,*) trim(F_uwind)
        write (nu_diag,*) trim(F_vwind)
        write (nu_diag,*) trim(F_tair)
        write (nu_diag,*) trim(F_humid)
     endif                     ! master_task
   end subroutine oned_files

   !=======================================================================
   subroutine ocn_data_clim (dt)
     ! Climatological ocean forcing
     ! Interpolate monthly sss, sst data to timestep.
     ! Restore prognostic sst to data.
     ! Interpolate fields from U grid to T grid if necessary.
     ! author: Elizabeth C. Hunke and William H. Lipscomb, LANL
     use ice_domain, only: nblocks
     use ice_flux, only: sss, sst
     real (kind=dbl_kind), intent(in) :: dt      ! time step
     integer (kind=int_kind) :: i, j, iblk  , & ! horizontal indices
                                ixm,ixp     , & ! record numbers for neighboring months
                                maxrec      , & ! maximum record number
                                recslot     , & ! spline slot for current record
                                midmonth        ! middle day of month
     real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks) ::  sstdat ! data value toward which SST is restored
     logical (kind=log_kind) :: readm
     character(len=*), parameter :: subname = '(ocn_data_clim)'
     if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'
     if (my_task == master_task .and. istep == 1) then
        if (trim(ocn_data_type)=='clim') then
           write (nu_diag,*) ' '
           write (nu_diag,*) 'SSS data interpolated to timestep:'
           write (nu_diag,*) trim(F_sss)
           write (nu_diag,*) ' '
           write (nu_diag,*) 'SST data interpolated to timestep:'
           write (nu_diag,*) trim(F_sst)
           if (restore_ocn) write (nu_diag,*) 'SST restoring timescale (days) =', trestore
        endif
     endif ! my_task, istep
     !-------------------------------------------------------------------
     ! monthly data
     ! Assume that monthly data values are located in the middle of the
     ! month.
     !-------------------------------------------------------------------
     if (trim(ocn_data_type)=='clim') then
        midmonth = 15          ! data is given on 15th of every month
        !      midmonth = fix(p5 * real(daymo(mmonth)))  ! exact middle
        ! Compute record numbers for surrounding months
        maxrec = 12
        ixm  = mod(mmonth+maxrec-2,maxrec) + 1
        ixp  = mod(mmonth,         maxrec) + 1
        if (mday >= midmonth) ixm = -99 ! other two points will be used
        if (mday <  midmonth) ixp = -99
        ! Determine whether interpolation will use values 1:2 or 2:3
        ! recslot = 2 means we use values 1:2, with the current value (2)
        !  in the second slot
        ! recslot = 1 means we use values 2:3, with the current value (2)
        !  in the first slot
        recslot = 1            ! latter half of month
        if (mday < midmonth) recslot = 2 ! first half of month
        ! Find interpolation coefficients
        call interp_coeff_monthly (recslot)
        readm = .false.
        if (istep==1 .or. (mday==midmonth .and. msec==0)) readm = .true.
        !-------------------------------------------------------------------
        ! Read two monthly SSS values and interpolate.
        ! Note: SSS is restored instantaneously to data.
        !-------------------------------------------------------------------
        call read_clim_data (readm, 0, ixm, mmonth, ixp, F_sss, sss_data, field_loc_center, field_type_scalar)
        call interpolate_data (sss_data, sss)
        !$OMP PARALLEL DO PRIVATE(iblk,i,j)
        do iblk = 1, nblocks
           do j = 1, ny_block
              do i = 1, nx_block
                 sss(i,j,iblk) = max(sss(i,j,iblk), c0)
              enddo
           enddo
        enddo !$OMP END PARALLEL DO
        call ocn_freezing_temperature
     endif
     !-------------------------------------------------------------------
     ! Read two monthly SST values and interpolate.
     ! Restore toward interpolated value.
     !-------------------------------------------------------------------
     if (trim(ocn_data_type)=='clim') then
        call read_clim_data (readm, 0, ixm, mmonth, ixp, F_sst, sst_data, field_loc_center, field_type_scalar)
        call interpolate_data (sst_data, sstdat)
        if (restore_ocn) then
           !$OMP PARALLEL DO PRIVATE(iblk,i,j)
           do iblk = 1, nblocks
              do j = 1, ny_block
                 do i = 1, nx_block
                    sst(i,j,iblk) = sst(i,j,iblk) + (sstdat(i,j,iblk)-sst(i,j,iblk))*dt/trest
                 enddo
              enddo
           enddo
           !$OMP END PARALLEL DO
        endif
     endif
   end subroutine ocn_data_clim

   !=======================================================================
   subroutine ocn_data_ncar_init
     ! NCAR CESM M-configuration (AIO) ocean forcing
     ! Reads NCAR pop ocean forcing data set 'pop_frc_gx1v3_010815.nc'
     !
     ! List of ocean forcing fields: Note that order is important!
     ! (order is determined by field list in vname).
     !
     ! For ocean mixed layer-----------------------------units
     !
     ! 1  sst------temperature---------------------------(C)
     ! 2  sss------salinity------------------------------(ppt)
     ! 3  hbl------depth---------------------------------(m)
     ! 4  u--------surface u current---------------------(m/s)
     ! 5  v--------surface v current---------------------(m/s)
     ! 6  dhdx-----surface tilt x direction--------------(m/m)
     ! 7  dhdy-----surface tilt y direction--------------(m/m)
     ! 8  qdp------ocean sub-mixed layer heat flux-------(W/m2)
     !
     ! Fields 4, 5, 6, 7 are on the U-grid; 1, 2, 3, and 8 are
     ! on the T-grid.
     ! authors: Bruce Briegleb, NCAR
     !          Elizabeth Hunke, LANL
     use ice_blocks, only: nx_block, ny_block
     use ice_domain_size, only: max_blocks
#ifdef USE_NETCDF
     use netcdf
#endif
     integer (kind=int_kind) :: n   , & ! field index
                                m   , & ! month index
                                nrec, & ! record number for direct access
                                nbits
     character(char_len) :: vname(nfld) ! variable names to search for in file
     data vname /  'T', 'S', 'hblt', 'U', 'V', 'dhdx', 'dhdy', 'qdp' /
     integer (kind=int_kind) :: status  , & ! status flag
                                fid     , & ! file id
                                dimid   , & ! dimension id
                                nlat    , & ! number of longitudes of data
                                nlon        ! number of latitudes  of data
     real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: work1
     character(len=*), parameter :: subname = '(ocn_data_ncar_init)'
     if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'
     if (my_task == master_task) then
        write (nu_diag,*) 'WARNING: evp_prep calculates surface tilt'
        write (nu_diag,*) 'WARNING: stress from geostrophic currents,'
        write (nu_diag,*) 'WARNING: not data from ocean forcing file.'
        write (nu_diag,*) 'WARNING: Alter ice_dyn_evp.F90 if desired.'
        if (restore_ocn) write (nu_diag,*) 'SST restoring timescale = ',trestore,' days'
        F_sst = trim(ocn_data_dir)//'/'//trim(oceanmixed_file) ! not just sst
        !---------------------------------------------------------------
        ! Read in ocean forcing data from an existing file
        !---------------------------------------------------------------
        write (nu_diag,*) 'ocean mixed layer forcing data file = ', trim(F_sst)
     endif ! master_task
     if (trim(ocn_data_format) == 'nc') then
#ifdef USE_NETCDF
        if (my_task == master_task) then
           call ice_open_nc(F_sst, fid)
           status = nf90_inq_dimid(fid,'ni',dimid)
           status = nf90_inquire_dimension(fid,dimid,len=nlon)
           status = nf90_inq_dimid(fid,'nj',dimid)
           status = nf90_inquire_dimension(fid,dimid,len=nlat)
           if( nlon .ne. nx_global ) then
              call abort_ice (error_message=subname//'ice: ocn frc file nlon ne nx_global', file=__FILE__, line=__LINE__)
           endif
           if( nlat .ne. ny_global ) then
            call abort_ice (error_message=subname//'ice: ocn frc file nlat ne ny_global', file=__FILE__, line=__LINE__)
         endif
      endif ! master_task
      ! Read in ocean forcing data for all 12 months
      do n=1,nfld
         do m=1,12
            call ice_read_nc(fid, m, vname(n), work1, debug_forcing, field_loc_center, field_type_scalar)
            ocn_frc_m(:,:,:,n,m) = work1(:,:,:)
         enddo               ! month loop
      enddo               ! field loop
      if (my_task == master_task) call ice_close_nc(fid)
#else
      call abort_ice(subname//'ERROR: USE_NETCDF cpp not defined for '//trim(F_sst), file=__FILE__, line=__LINE__)
#endif
   else ! binary format
      nbits = 64
      call ice_open (nu_forcing, F_sst, nbits)
      nrec = 0
      do n=1,nfld
         do m=1,12
            nrec = nrec + 1
            if (n >= 4 .and. n <= 7) then
               call ice_read (nu_forcing, nrec, work1, 'rda8', debug_forcing, field_loc_NEcorner, field_type_vector)
            else
               call ice_read (nu_forcing, nrec, work1, 'rda8', debug_forcing, field_loc_center, field_type_scalar)
            endif
            ocn_frc_m(:,:,:,n,m) = work1(:,:,:)
         enddo               ! month loop
      enddo               ! field loop
      close (nu_forcing)
   endif
 end subroutine ocn_data_ncar_init

 !=======================================================================
 subroutine ocn_data_ncar_init_3D
   ! Reads NCAR pop ocean forcing data set 'oceanmixed_ice_depth.nc'
   !
   ! List of ocean forcing fields: Note that order is important!
   ! (order is determined by field list in vname).
   !
   ! For ocean mixed layer-----------------------------units
   !
   ! 1  sst------temperature---------------------------(C)
   ! 2  sss------salinity------------------------------(ppt)
   ! 3  hbl------depth---------------------------------(m)
   ! 4  u--------surface u current---------------------(m/s)
   ! 5  v--------surface v current---------------------(m/s)
   ! 6  dhdx-----surface tilt x direction--------------(m/m)
   ! 7  dhdy-----surface tilt y direction--------------(m/m)
   ! 8  qdp------ocean sub-mixed layer heat flux-------(W/m2)
   !
   ! All fields are on the T-grid.
   !
   ! authors: Bruce Briegleb, NCAR
   !          Elizabeth Hunke, LANL
   use ice_blocks, only: nx_block, ny_block
   use ice_domain_size, only: max_blocks
   use ice_grid, only: grid_average_X2Y, ANGLET
   use ice_read_write, only: ice_read_nc_uv
#ifdef USE_NETCDF
   use netcdf
#endif
#ifdef USE_NETCDF
   integer (kind=int_kind) :: n   , & ! field index
                              m   , & ! month index
                              nzlev   ! z level of currents
   character(char_len) :: vname(nfld) ! variable names to search for in file
   data vname / 'T', 'S', 'hblt', 'U', 'V', 'dhdx', 'dhdy', 'qdp' /
   integer (kind=int_kind) :: fid        , & ! file id
                              dimid          ! dimension id
   integer (kind=int_kind) :: status  , & ! status flag
                              nlat    , & ! number of longitudes of data
                              nlon        ! number of latitudes  of data
   real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: work1, work2
#endif
   character(len=*), parameter :: subname = '(ocn_data_ncar_init_3D)'
   if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'
   if (my_task == master_task) then
      write (nu_diag,*) 'WARNING: evp_prep calculates surface tilt'
      write (nu_diag,*) 'WARNING: stress from geostrophic currents,'
      write (nu_diag,*) 'WARNING: not data from ocean forcing file.'
      write (nu_diag,*) 'WARNING: Alter ice_dyn_evp.F if desired.'
      if (restore_ocn) write (nu_diag,*) 'SST restoring timescale = ',trestore,' days'
      F_sst = trim(ocn_data_dir)//'/'//trim(oceanmixed_file) ! not just sst
      !---------------------------------------------------------------
      ! Read in ocean forcing data from an existing file
      !---------------------------------------------------------------
      write (nu_diag,*) 'ocean mixed layer forcing data file = ', trim(F_sst)
      write (nu_diag,*)
   endif ! master_task
   if (trim(ocn_data_format) == 'nc') then
#ifdef USE_NETCDF
      if (my_task == master_task) then
         call ice_open_nc(F_sst, fid)
         status = nf90_inq_dimid(fid,'ni',dimid)
         status = nf90_inquire_dimension(fid,dimid,len=nlon)
         status = nf90_inq_dimid(fid,'nj',dimid)
         status = nf90_inquire_dimension(fid,dimid,len=nlat)
         if( nlon .ne. nx_global ) then
            call abort_ice (error_message=subname//'ice: ocn frc file nlon ne nx_global', file=__FILE__, line=__LINE__)
         endif
         if( nlat .ne. ny_global ) then
            call abort_ice (error_message=subname//'ice: ocn frc file nlat ne ny_global', file=__FILE__, line=__LINE__)
         endif
      endif ! master_task
      ! Read in ocean forcing data for all 12 months
      do n=1,nfld
         do m=1,12
            ! Note: netCDF does single to double conversion if necessary
            if (n == 4 .or. n == 5) then ! 3D currents
               nzlev = 1                 ! surface currents
               call ice_read_nc_uv(fid, m, nzlev, vname(n), work1, debug_forcing, field_loc_center, field_type_scalar)
            else
               call ice_read_nc(fid, m, vname(n), work1, debug_forcing, field_loc_center, field_type_scalar)
            endif
            ! the land mask used in ocean_mixed_depth.nc does not
            ! match our gx1v3 mask (hm)
            where (work1(:,:,:) < -900.) work1(:,:,:) = c0
            ocn_frc_m(:,:,:,n,m) = work1(:,:,:)
         enddo               ! month loop
      enddo               ! field loop
      if (my_task == master_task) call ice_close_nc(fid)
      ! Rotate vector quantities and shift to U-grid
      do n=4,6,2
         do m=1,12
            work1(:,:,:) = ocn_frc_m(:,:,:,n  ,m)
            work2(:,:,:) = ocn_frc_m(:,:,:,n+1,m)
            ocn_frc_m(:,:,:,n  ,m) = work1(:,:,:)*cos(ANGLET(:,:,:)) + work2(:,:,:)*sin(ANGLET(:,:,:))
            ocn_frc_m(:,:,:,n+1,m) = work2(:,:,:)*cos(ANGLET(:,:,:)) - work1(:,:,:)*sin(ANGLET(:,:,:))
            work1(:,:,:) = ocn_frc_m(:,:,:,n  ,m)
            work2(:,:,:) = ocn_frc_m(:,:,:,n+1,m)
            call grid_average_X2Y('F',work1,'T',ocn_frc_m(:,:,:,n  ,m),'U')
            call grid_average_X2Y('F',work2,'T',ocn_frc_m(:,:,:,n+1,m),'U')
         enddo               ! month loop
      enddo               ! field loop
#else
      call abort_ice(subname//'ERROR: USE_NETCDF cpp not defined', file=__FILE__, line=__LINE__)
#endif
   else  ! binary format
      call abort_ice (error_message=subname//'new ocean forcing is netcdf only', file=__FILE__, line=__LINE__)
   endif
 end subroutine ocn_data_ncar_init_3D

 !=======================================================================
 subroutine ocn_data_ncar(dt)
   ! Interpolate monthly ocean data to timestep.
   ! Restore sst if desired. sst is updated with surface fluxes in ice_ocean.F.
   use ice_blocks, only: nx_block, ny_block
   use ice_global_reductions, only: global_minval, global_maxval
   use ice_domain, only: nblocks, distrb_info
   use ice_domain_size, only: max_blocks
   use ice_flux, only: sss, sst, Tf, uocn, vocn, ss_tltx, ss_tlty, qdp, hmix
   use ice_restart_shared, only: restart
   use ice_grid, only: hm, tmask, umask
   real (kind=dbl_kind), intent(in) :: dt      ! time step
   integer (kind=int_kind) :: i, j, n, iblk   , &
                              ixm,ixp         , & ! record numbers for neighboring months
                              maxrec          , & ! maximum record number
                              recslot         , & ! spline slot for current record
                              midmonth            ! middle day of month
   real (kind=dbl_kind) :: vmin, vmax
   real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: work1
   character(len=*), parameter :: subname = '(ocn_data_ncar)'
   if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'
    !-------------------------------------------------------------------
    ! monthly data
    ! Assume that monthly data values are located in the middle of the
    ! month.
    !-------------------------------------------------------------------
   midmonth = 15  ! data is given on 15th of every month
   ! Compute record numbers for surrounding months
   maxrec = 12
   ixm  = mod(mmonth+maxrec-2,maxrec) + 1
   ixp  = mod(mmonth,         maxrec) + 1
   if (mday >= midmonth) ixm = -99  ! other two points will be used
   if (mday <  midmonth) ixp = -99
   ! Determine whether interpolation will use values 1:2 or 2:3
   ! recslot = 2 means we use values 1:2, with the current value (2)
   !  in the second slot
   ! recslot = 1 means we use values 2:3, with the current value (2)
   !  in the first slot
   recslot = 1                             ! latter half of month
   if (mday < midmonth) recslot = 2        ! first half of month
   ! Find interpolation coefficients
   call interp_coeff_monthly (recslot)
   sst_data(:,:,:,:) = c0
   do n = nfld, 1, -1
      do iblk = 1, nblocks
         ! use sst_data arrays as temporary work space until n=1
         if (ixm /= -99) then  ! first half of month
            sst_data(:,:,1,iblk) = ocn_frc_m(:,:,iblk,n,ixm)
            sst_data(:,:,2,iblk) = ocn_frc_m(:,:,iblk,n,mmonth)
         else ! second half of month
            sst_data(:,:,1,iblk) = ocn_frc_m(:,:,iblk,n,mmonth)
            sst_data(:,:,2,iblk) = ocn_frc_m(:,:,iblk,n,ixp)
         endif
      enddo
      call interpolate_data (sst_data,work1)
      ! masking by hm is necessary due to NaNs in the data file
      do j = 1, ny_block
         do i = 1, nx_block
            if (n == 2) sss    (i,j,:) = c0
            if (n == 3) hmix   (i,j,:) = c0
            if (n == 4) uocn   (i,j,:) = c0
            if (n == 5) vocn   (i,j,:) = c0
            if (n == 6) ss_tltx(i,j,:) = c0
            if (n == 7) ss_tlty(i,j,:) = c0
            if (n == 8) qdp    (i,j,:) = c0
            do iblk = 1, nblocks
               if (hm(i,j,iblk) == c1) then
                  if (n == 2) sss    (i,j,iblk) = work1(i,j,iblk)
                  if (n == 3) hmix   (i,j,iblk) = max(mixed_layer_depth_default,work1(i,j,iblk))
                  if (n == 4) uocn   (i,j,iblk) = work1(i,j,iblk)
                  if (n == 5) vocn   (i,j,iblk) = work1(i,j,iblk)
                  if (n == 6) ss_tltx(i,j,iblk) = work1(i,j,iblk)
                  if (n == 7) ss_tlty(i,j,iblk) = work1(i,j,iblk)
                  if (n == 8) qdp    (i,j,iblk) = work1(i,j,iblk)
               endif
            enddo
         enddo
      enddo
   enddo
   do j = 1, ny_block
      do i = 1, nx_block
         sss (i,j,:) = max (sss(i,j,:), c0)
         hmix(i,j,:) = max(hmix(i,j,:), c0)
      enddo
   enddo
   call ocn_freezing_temperature
   if (restore_ocn) then
      do j = 1, ny_block
         do i = 1, nx_block
            sst(i,j,:) = sst(i,j,:) + (work1(i,j,:)-sst(i,j,:))*dt/trest
         enddo
      enddo
      !     else sst is only updated in ice_ocean.F
   endif
   ! initialize sst properly on first step
   if (istep1 <= 1 .and. .not. (restart)) then
      call interpolate_data (sst_data,sst)
      !$OMP PARALLEL DO PRIVATE(iblk,i,j)
      do iblk = 1, nblocks
         do j = 1, ny_block
            do i = 1, nx_block
               if (hm(i,j,iblk) == c1) then
                  sst(i,j,iblk) =  max (sst(i,j,iblk), Tf(i,j,iblk))
               else
                  sst(i,j,iblk) = c0
               endif
            enddo
         enddo
      enddo !$OMP END PARALLEL DO
   endif
   if (debug_forcing) then
      if (my_task == master_task) write (nu_diag,*) 'ocn_data_ncar'
      vmin = global_minval(Tf,distrb_info,tmask)
      vmax = global_maxval(Tf,distrb_info,tmask)
      if (my_task.eq.master_task) write (nu_diag,*) 'Tf',vmin,vmax
      vmin = global_minval(sst,distrb_info,tmask)
      vmax = global_maxval(sst,distrb_info,tmask)
      if (my_task.eq.master_task) write (nu_diag,*) 'sst',vmin,vmax
      vmin = global_minval(sss,distrb_info,tmask)
      vmax = global_maxval(sss,distrb_info,tmask)
      if (my_task.eq.master_task) write (nu_diag,*) 'sss',vmin,vmax
      vmin = global_minval(hmix,distrb_info,tmask)
      vmax = global_maxval(hmix,distrb_info,tmask)
      if (my_task.eq.master_task) write (nu_diag,*) 'hmix',vmin,vmax
      vmin = global_minval(uocn,distrb_info,umask)
      vmax = global_maxval(uocn,distrb_info,umask)
      if (my_task.eq.master_task) write (nu_diag,*) 'uocn',vmin,vmax
      vmin = global_minval(vocn,distrb_info,umask)
      vmax = global_maxval(vocn,distrb_info,umask)
      if (my_task.eq.master_task) write (nu_diag,*) 'vocn',vmin,vmax
      vmin = global_minval(ss_tltx,distrb_info,umask)
      vmax = global_maxval(ss_tltx,distrb_info,umask)
      if (my_task.eq.master_task) write (nu_diag,*) 'ss_tltx',vmin,vmax
      vmin = global_minval(ss_tlty,distrb_info,umask)
      vmax = global_maxval(ss_tlty,distrb_info,umask)
      if (my_task.eq.master_task) write (nu_diag,*) 'ss_tlty',vmin,vmax
      vmin = global_minval(qdp,distrb_info,tmask)
      vmax = global_maxval(qdp,distrb_info,tmask)
      if (my_task.eq.master_task) write (nu_diag,*) 'qdp',vmin,vmax
   endif
 end subroutine ocn_data_ncar

 !=======================================================================
 subroutine ocn_data_oned
   ! ocean data for oned configuration
   ! Current (released) values are the same as the defaults (ice_flux.F90)
   use ice_flux, only: sss, sst, Tf, uocn, vocn, ss_tltx, ss_tlty, qdp, hmix, frzmlt
   character(len=*), parameter :: subname = '(ocn_data_oned)'
   if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'
   sss    (:,:,:) = 34.0_dbl_kind   ! sea surface salinity (ppt)
   call ocn_freezing_temperature
   sst    (:,:,:) = Tf(:,:,:)       ! sea surface temp (C)
   uocn   (:,:,:) = c0              ! surface ocean currents (m/s)
   vocn   (:,:,:) = c0
   ss_tltx(:,:,:) = c0              ! sea surface tilt (m/m)
   ss_tlty(:,:,:) = c0
   frzmlt (:,:,:) = c0              ! freezing/melting potential (W/m^2)
   qdp    (:,:,:) = c0              ! deep ocean heat flux (W/m^2)
   hmix   (:,:,:) = mixed_layer_depth_default   ! ocean mixed layer depth
 end subroutine ocn_data_oned

 !=======================================================================
 subroutine ocn_data_hadgem(dt)
   !  Reads in HadGEM ocean forcing data as required from netCDF files
   !  Current options (selected by ocn_data_type)
   !  hadgem_sst: read in sst only
   !  hadgem_sst_uvocn: read in sst plus uocn and vocn
   ! authors: Ann Keen, Met Office
   use ice_domain, only: nblocks
   use ice_domain_size, only: max_blocks
   use ice_flux, only: sst, uocn, vocn
   use ice_grid, only: grid_average_X2Y, ANGLET
   real (kind=dbl_kind), intent(in) :: dt      ! time step
   integer (kind=int_kind) :: i, j, iblk  , &
                              ixm,ixp     , & ! record numbers for neighboring months
                              maxrec      , & ! maximum record number
                              recslot     , & ! spline slot for current record
                              midmonth        ! middle day of month
   real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks) :: sstdat              ! data value toward which SST is restored
   real (kind=dbl_kind) :: workx, worky
   logical (kind=log_kind) :: readm
   character (char_len) :: fieldname     ! field name in netcdf file
   character (char_len_long) :: filename      ! name of netCDF file
   character(len=*), parameter :: subname = '(ocn_data_hadgem)'
   if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'
    !-------------------------------------------------------------------
    ! monthly data
    ! Assume that monthly data values are located in the middle of the
    ! month.
    !-------------------------------------------------------------------
   midmonth = 15  ! data is given on 15th of every month
   ! Compute record numbers for surrounding months
   maxrec = 12
   ixm  = mod(mmonth+maxrec-2,maxrec) + 1
   ixp  = mod(mmonth,         maxrec) + 1
   if (mday >= midmonth) ixm = -99  ! other two points will be used
   if (mday <  midmonth) ixp = -99
   ! Determine whether interpolation will use values 1:2 or 2:3
   ! recslot = 2 means we use values 1:2, with the current value (2)
   !  in the second slot
   ! recslot = 1 means we use values 2:3, with the current value (2)
   !  in the first slot
   recslot = 1                             ! latter half of month
   if (mday < midmonth) recslot = 2        ! first half of month
   ! Find interpolation coefficients
   call interp_coeff_monthly (recslot)
   ! Read 2 monthly values
   readm = .false.
   if (istep==1 .or. (mday==midmonth .and. msec==0)) readm = .true.
   if (my_task == master_task .and. istep == 1) then
      write (nu_diag,*) ' '
      write (nu_diag,*) 'SST data interpolated to timestep:'
      write (nu_diag,*) trim(ocn_data_dir)//'/MONTHLY/sst.1997.nc'
      if (restore_ocn) write (nu_diag,*) 'SST restoring timescale (days) =', trestore
      if (trim(ocn_data_type)=='hadgem_sst_uvocn') then
         write (nu_diag,*) ' '
         write (nu_diag,*) 'uocn and vocn interpolated to timestep:'
         write (nu_diag,*) trim(ocn_data_dir)//'/MONTHLY/uocn.1997.nc'
         write (nu_diag,*) trim(ocn_data_dir)//'/MONTHLY/vocn.1997.nc'
      endif
   endif                     ! my_task, istep
   ! -----------------------------------------------------------
   ! SST
   ! -----------------------------------------------------------
   F_sst = trim(ocn_data_dir)//'/MONTHLY/sst.1997.nc'
   fieldname='sst'
   call read_data_nc (readm, 0, fyear, ixm, mmonth, ixp, maxrec, F_sst, fieldname, sst_data, field_loc_center, field_type_scalar)
   ! Interpolate to current time step
   call interpolate_data (sst_data, sstdat)
   ! Restore SSTs if required
   if (restore_ocn) then
      !$OMP PARALLEL DO PRIVATE(iblk,i,j)
      do iblk = 1, nblocks
         do j = 1, ny_block
            do i = 1, nx_block
               sst(i,j,iblk) = sst(i,j,iblk) + (sstdat(i,j,iblk)-sst(i,j,iblk))*dt/trest
            enddo
         enddo
      enddo !$OMP END PARALLEL DO
   endif
   ! -----------------------------------------------------------
   ! Ocean currents
   ! --------------
   ! Values read in are on T grid and oriented geographically, hence
   ! vectors need to be rotated to model grid and then interpolated
   ! to U grid.
   ! Also need to be converted from cm s-1 (UM) to m s-1 (CICE)
   ! -----------------------------------------------------------
   if (trim(ocn_data_type)=='hadgem_sst_uvocn') then
      filename = trim(ocn_data_dir)//'/MONTHLY/uocn.1997.nc'
      fieldname='uocn'
      call read_data_nc (readm, 0, fyear, ixm, mmonth, ixp, maxrec, filename, fieldname, uocn_data, field_loc_center, field_type_vector)
      ! Interpolate to current time step
      call interpolate_data (uocn_data, uocn)
      filename = trim(ocn_data_dir)//'/MONTHLY/vocn.1997.nc'
      fieldname='vocn'
      call read_data_nc (readm, 0, fyear, ixm, mmonth, ixp, maxrec, filename, fieldname, vocn_data, field_loc_center, field_type_vector)
      ! Interpolate to current time step
      call interpolate_data (vocn_data, vocn)
      !-----------------------------------------------------------------
      ! Rotate zonal/meridional vectors to local coordinates and change  units
      !-----------------------------------------------------------------
      !$OMP PARALLEL DO PRIVATE(iblk,i,j,workx,worky)
      do iblk = 1, nblocks
         do j = 1, ny_block
            do i = 1, nx_block
               workx      = uocn(i,j,iblk)
               worky      = vocn(i,j,iblk)
               uocn(i,j,iblk) = workx*cos(ANGLET(i,j,iblk)) + worky*sin(ANGLET(i,j,iblk))
               vocn(i,j,iblk) = worky*cos(ANGLET(i,j,iblk)) - workx*sin(ANGLET(i,j,iblk))
               uocn(i,j,iblk) = uocn(i,j,iblk) * cm_to_m
               vocn(i,j,iblk) = vocn(i,j,iblk) * cm_to_m
            enddo   ! i
         enddo   ! j
      enddo      ! nblocks !$OMP END PARALLEL DO
      !-----------------------------------------------------------------
      ! Interpolate to U grid
      !-----------------------------------------------------------------
      ! tcraig, this is now computed in dynamics for consistency
   endif    !   ocn_data_type = hadgem_sst_uvocn
 end subroutine ocn_data_hadgem

 !=======================================================================
 subroutine ocn_data_hycom_init
   ! Read SSS+SST from a HYCOM file converted to NetCDF format.
   ! HYCOM binary2NetCDF: hcdata2ncdf2d (or hcdata2ncdf3z)
   !   + rename/link file
   use ice_blocks, only: nx_block, ny_block
   use ice_domain, only: nblocks
   use ice_flux, only: sss, sst, Tf
   integer (kind=int_kind) :: &
        i, j, iblk       , & ! horizontal indices
        fid                  ! file id for netCDF file
   character (char_len) :: &
        fieldname            ! field name in netcdf file
   character(len=*), parameter :: subname = '(ocn_data_hycom_init)'
   if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'
   if (trim(ocn_data_type) == 'hycom') then
      F_sss = trim(ocn_data_dir)//'ice.restart.surf.nc'
      if (my_task == master_task) then
         write (nu_diag,*)' '
         write (nu_diag,*)'Initial SSS file: ',trim(F_sss)
      endif
      fieldname = 'sss'
      call ice_open_nc (F_sss, fid)
      call ice_read_nc (fid, 1 , fieldname, sss, debug_forcing, field_loc_center, field_type_scalar)
      call ice_close_nc(fid)
      call ocn_freezing_temperature
      F_sst = trim(ocn_data_dir)//'ice.restart.surf.nc'
      if (my_task == master_task) then
         write (nu_diag,*)' '
         write (nu_diag,*)'Initial SST file: ',trim(F_sst)
      endif
      fieldname = 'sst'
      call ice_open_nc (F_sst, fid)
      call ice_read_nc (fid, 1 , fieldname, sst, debug_forcing, field_loc_center, field_type_scalar)
      call ice_close_nc(fid)
      ! Make sure sst is not less than freezing temperature Tf
      !$OMP PARALLEL DO PRIVATE(iblk,i,j)
      do iblk = 1, nblocks
         do j = 1, ny_block
            do i = 1, nx_block
               sst(i,j,iblk) = max(sst(i,j,iblk),Tf(i,j,iblk))
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
   endif
 end subroutine ocn_data_hycom_init

 !=======================================================================
 subroutine hycom_atm_files
   use ice_broadcast, only: broadcast_array, broadcast_scalar
   integer (kind = int_kind) :: &
        fid          ! File id
   character (char_len) :: &
        varname      ! variable name in netcdf file
   character(len=*), parameter :: subname = '(hycom_atm_files)'
   if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'
   F_fsw   = trim(atm_data_dir)//'/forcing.shwflx.nc'
   F_flw   = trim(atm_data_dir)//'/forcing.radflx.nc'
   F_rain  = trim(atm_data_dir)//'/forcing.precip.nc'
   F_uwind = trim(atm_data_dir)//'/forcing.wndewd.nc'
   F_vwind = trim(atm_data_dir)//'/forcing.wndnwd.nc'
   F_tair  = trim(atm_data_dir)//'/forcing.airtmp.nc'
   F_humid = trim(atm_data_dir)//'/forcing.vapmix.nc'
   ! Read time vector from "F_tair"
   call ice_open_nc(F_tair, fid)
   varname='MT'
   call ice_get_ncvarsize(fid,varname,Njday_atm)
   call broadcast_scalar(Njday_atm,master_task)
   allocate(jday_atm(Njday_atm))
   call ice_read_vec_nc(fid,Njday_atm, varname,jday_atm, .true.)
   call ice_close_nc(fid)
   call broadcast_array(jday_atm ,master_task)
   ! Write diag info
   if (my_task == master_task) then
      write (nu_diag,*) ' '
      write (nu_diag,*) 'CICE: Atm. (hycomdate) Start = ',jday_atm(1)
      write (nu_diag,*) 'CICE: Atm. (hycomdate) End   = ',jday_atm(Njday_atm)
      write (nu_diag,*) 'CICE: Total Atm timesteps    = ',Njday_atm
      write (nu_diag,*) 'CICE: Atmospheric forcing files:'
      write (nu_diag,*) trim(F_fsw)
      write (nu_diag,*) trim(F_flw)
      write (nu_diag,*) trim(F_rain)
      write (nu_diag,*) trim(F_uwind)
      write (nu_diag,*) trim(F_vwind)
      write (nu_diag,*) trim(F_tair)
      write (nu_diag,*) trim(F_humid)
   endif                     ! master_task
 end subroutine hycom_atm_files

 !=======================================================================
 subroutine hycom_atm_data
   use ice_flux, only: fsw, fsnow, Tair, uatm, vatm, Qa, flw
   use ice_domain, only: nblocks
   integer (kind=int_kind) :: &
        recnum       ! record number
   real (kind=dbl_kind) :: &
        hcdate         ! current time in HYCOM jday units
   logical (kind=log_kind) :: read6
   character (char_len) :: &
        fieldname    ! field name in netcdf file
   integer (kind=int_kind) :: &
        i, j, iblk      ! horizontal indices
   real (kind=dbl_kind) :: Tffresh, secday
   character(len=*), parameter :: subname = '(hycom_atm_data)'
   if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'
   call icepack_query_parameters(Tffresh_out=Tffresh)
   call icepack_query_parameters(secday_out=secday)
   ! current time in HYCOM jday units (HYCOM ref year: 1900,12,31,000000)
   hcdate = real(compute_days_between(1900,12,31,myear,mmonth,mday)) + msec/secday
   ! Init recnum try
   recnum=min(max(oldrecnum,1),Njday_atm-1)
   ! Find correct time in ATM data ... assume cont. incr. time-axis
   do while ( recnum<Njday_atm )
      if ( hcdate>=jday_atm(recnum) .and. &
           hcdate<=jday_atm(recnum+1) ) exit
      if ( abs(hcdate-jday_atm(recnum))<p001 ) exit  ! Accept within tolerance = 0.001 days
      if ( abs(hcdate-jday_atm(recnum+1))<p001 ) exit  ! Accept within tolerance = 0.001 days
      recnum=recnum+1
   enddo
   ! Last Atm date might be the same as last CICE date.
   recnum=min(recnum,Njday_atm-1)
   ! Check if current time do not exceed last forcing time
   ! + check forcing is available before (or at) current forcing time
   if ( hcdate>jday_atm(recnum+1)+p001 .or. hcdate<jday_atm(recnum)-p001) then
      write (nu_diag,*) 'ERROR: CICE: Atm forcing not available at hcdate =',hcdate
      write (nu_diag,*) 'ERROR: CICE: myear, yday ,msec = ',myear, yday, msec
      call abort_ice ('ERROR: CICE stopped')
   endif
   ! Compute interpolation coefficients
   call interp_coeff2 (hcdate, jday_atm(recnum), jday_atm(recnum+1) )
   ! Read
   read6 = .false.
   if (istep==1 .or. oldrecnum /= recnum) read6 = .true.
   if (trim(atm_data_format) == 'nc') then
      if (read6 .and. my_task == master_task) write(nu_diag,*) 'CICE: Atm. read: = ',jday_atm(recnum), jday_atm(recnum+1)
      fieldname = 'airtmp'
      call read_data_nc_hycom (read6, recnum, F_tair, fieldname, Tair_data, field_loc_center, field_type_scalar)
      fieldname = 'wndewd'
      call read_data_nc_hycom (read6, recnum, F_uwind, fieldname, uatm_data, field_loc_center, field_type_vector)
      fieldname = 'wndnwd'
      call read_data_nc_hycom (read6, recnum, F_vwind, fieldname, vatm_data, field_loc_center, field_type_vector)
      fieldname = 'vapmix'
      call read_data_nc_hycom (read6, recnum, F_humid, fieldname, Qa_data, field_loc_center, field_type_scalar)
      fieldname = 'shwflx'
      call read_data_nc_hycom (read6, recnum, F_fsw, fieldname, fsw_data, field_loc_center, field_type_scalar)
      fieldname = 'radflx'
      call read_data_nc_hycom (read6, recnum, F_flw, fieldname, flw_data, field_loc_center, field_type_scalar)
      fieldname = 'precip'
      call read_data_nc_hycom (read6, recnum, F_rain, fieldname, fsnow_data, field_loc_center, field_type_scalar)
   else
      call abort_ice(subname//'ERROR: atm_data_format unavailable for hycom')
   endif
   ! Interpolate
   if (debug_forcing) then
      if (my_task == master_task) then
         write(nu_diag,*)'CICE: Atm. interpolate: = ',hcdate,c1intp,c2intp
      endif
   endif
   call interpolate_data (Tair_data, Tair)
   call interpolate_data (uatm_data, uatm)
   call interpolate_data (vatm_data, vatm)
   call interpolate_data ( fsw_data, fsw)
   call interpolate_data ( flw_data, flw)
   call interpolate_data (  Qa_data, Qa)
   call interpolate_data (fsnow_data, fsnow)
   ! Adjust data forcing to CICE units
   !$OMP PARALLEL DO PRIVATE(iblk,i,j)
   do iblk = 1, nblocks
      do j = 1, ny_block
         do i = 1, nx_block
            ! Air temperature: Degrees --> Kelvin
            Tair(i,j,iblk) = Tair(i,j,iblk) + Tffresh
         enddo  ! i
      enddo  ! j
   enddo     ! nblocks  !$OMP END PARALLEL DO
   ! Save record number for next time step
   oldrecnum = recnum
 end subroutine hycom_atm_data

 !=======================================================================
 subroutine read_data_nc_point (flag, recd, yr, ixm, ixx, ixp, maxrec, data_file, fieldname, field_data, field_loc, field_type)
   !
   ! If data is at the beginning of a one-year record, get data from
   !  the previous year.
   ! If data is at the end of a one-year record, get data from the
   !  following year.
   ! If no earlier data exists (beginning of fyear_init), then
   !  (1) For monthly data, get data from the end of fyear_final.
   !  (2) For more frequent data, let the ixm value equal the
   !      first value of the year.
   ! If no later data exists (end of fyear_final), then
   !  (1) For monthly data, get data from the beginning of fyear_init.
   !  (2) For more frequent data, let the ixp value
   !      equal the last value of the year.
   ! In other words, we assume persistence when daily or 6-hourly
   !   data is missing, and we assume periodicity when monthly data
   !   is missing.
   !
   use ice_diagnostics, only: debug_model_step
   logical (kind=log_kind), intent(in) :: flag
   integer (kind=int_kind), intent(in) :: &
        recd                , & ! baseline record number
        yr                  , & ! year of forcing data
        ixm, ixx, ixp       , & ! record numbers of 3 data values
                                ! relative to recd
        maxrec                  ! maximum record value
   character (char_len_long), intent(in) :: &
        data_file               ! data file to be read
   character (char_len), intent(in) :: &
        fieldname               ! field name in netCDF file
   integer (kind=int_kind), intent(in) :: &
        field_loc, &      ! location of field on staggered grid
        field_type        ! type of field (scalar, vector, angle)
   real (kind=dbl_kind), dimension(2), intent(inout) :: &
        field_data              ! 2 values needed for interpolation
   integer (kind=int_kind) :: &
        nrec             , & ! record number to read
        n2, n4           , & ! like ixm and ixp, but
                                ! adjusted at beginning and end of data
        arg              , & ! value of time argument in field_data
        fid                  ! file id for netCDF routines
   character(len=*), parameter :: subname = '(read_data_nc_point)'
   if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'
   call ice_timer_start(timer_readwrite)  ! reading/writing
   field_data = c0 ! to satisfy intent(out) attribute
   if (istep1 > debug_model_step) debug_forcing = .true.  !! debugging
   if (my_task==master_task .and. (debug_forcing)) then
      write(nu_diag,*) '  ', trim(data_file)
   endif
   if (flag) then
      !-----------------------------------------------------------------
      ! Initialize record counters
      ! (n2, n4 will change only at the very beginning or end of
      !  a forcing cycle.)
      !-----------------------------------------------------------------
      n2 = ixm
      n4 = ixp
      arg = 0
      !-----------------------------------------------------------------
      ! read data
      !-----------------------------------------------------------------
      if (ixm /= -99) then
         ! currently in first half of data interval
         if (ixx <= 1) then
            if (yr > fyear_init) then ! get data from previous year
               !call file_year (data_file, yr-1)
            else             ! yr = fyear_init, no prior data exists
               if (maxrec > 12) then ! extrapolate from first record
                  if (ixx == 1) n2 = ixx
               else          ! go to end of fyear_final
                  ! call file_year (data_file, fyear_final)
               endif
            endif            ! yr > fyear_init
         endif               ! ixx <= 1
         ! write(nu_diag,*) '!! read_data_nc !!!', trim(data_file)
         ! write(nu_diag,*) 'istep  ', istep
         ! write(nu_diag,*) 'fyear_final  ', fyear_final
         ! write(nu_diag,*) 'fyear_init  ', fyear_init
         ! write(nu_diag,*) 'ixm, ixx, ixp  ', ixm, ixx, ixp
         ! write(nu_diag,*) 'maxrec ', maxrec
         ! write(nu_diag,*) 'fieldname  ', fieldname
         call ice_open_nc(data_file, fid)
         arg = 1
         nrec = recd + n2
         call ice_read_nc(fid, nrec, fieldname, field_data(arg), debug_forcing, field_loc, field_type)
         !if (ixx==1) call ice_close_nc(fid)
         call ice_close_nc(fid)
      endif                  ! ixm ne -99
      ! always read ixx data from data file for current year
      ! call file_year (data_file, yr)
      call ice_open_nc (data_file, fid)
      arg = arg + 1
      nrec = recd + ixx
      call ice_read_nc (fid, nrec, fieldname, field_data(arg), debug_forcing, field_loc, field_type)
      if (ixp /= -99) then
         ! currently in latter half of data interval
         if (ixx==maxrec) then
            if (yr < fyear_final) then ! get data from following year
               call ice_close_nc(fid)
               !call file_year (data_file, yr+1)
               call ice_open_nc (data_file, fid)
            else             ! yr = fyear_final, no more data exists
               if (maxrec > 12) then ! extrapolate from ixx
                  n4 = ixx
               else          ! go to beginning of fyear_init
                  call ice_close_nc(fid)
                  ! call file_year (data_file, fyear_init)
                  call ice_open_nc (data_file, fid)
               endif
            endif            ! yr < fyear_final
         endif               ! ixx = maxrec
         arg = arg + 1
         nrec = recd + n4
         call ice_read_nc (fid, nrec, fieldname, field_data(arg), debug_forcing, field_loc, field_type)
      endif                  ! ixp /= -99
      call ice_close_nc(fid)
   endif                     ! flag
   call ice_timer_stop(timer_readwrite)  ! reading/writing
 end subroutine read_data_nc_point

 !=======================================================================
 subroutine ISPOL_files
   character(len=*), parameter :: subname = '(ISPOL_files)'
   if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'
   F_fsw = trim(atm_data_dir)//'/fsw_sfc_4Xdaily.nc'
   F_flw = trim(atm_data_dir)//'/flw_sfc_4Xdaily.nc'
   F_rain = trim(atm_data_dir)//'/fsnow_sfc_daily_mod3.nc'
   F_uwind = trim(atm_data_dir)//'/uatm_10m_daily.nc'
   F_vwind = trim(atm_data_dir)//'/vatm_10m_daily.nc'
   F_tair = trim(atm_data_dir)//'/Tair_2m_daily.nc'
   F_humid = trim(atm_data_dir)//'/Qa_2m_daily.nc'
   if (my_task == master_task) then
      write (nu_diag,*) ' '
      write (nu_diag,*) 'Atmospheric data files:'
      write (nu_diag,*) trim(F_fsw)
      write (nu_diag,*) trim(F_flw)
      write (nu_diag,*) trim(F_rain)
      write (nu_diag,*) trim(F_uwind)
      write (nu_diag,*) trim(F_vwind)
      write (nu_diag,*) trim(F_tair)
      write (nu_diag,*) trim(F_humid)
   endif                     ! master_task
 end subroutine ISPOL_files

 !=======================================================================
 subroutine ISPOL_data
   ! Defines atmospheric data fields for Antarctic Weddell sea location
   ! authors: Nicole Jeffery, LANL
   !
   use ice_flux, only: uatm, vatm, Tair, fsw,  Qa, rhoa, &
        frain, fsnow, flw
   character (char_len_long) :: &
        met_file,   &    ! netcdf filename
        fieldname        ! field name in netcdf file
   real (kind=dbl_kind), dimension(2), save :: &
        Tair_data_p      , &      ! air temperature (K) for interpolation
        Qa_data_p,  fsnow_data_p, &
        fsw_data_p, flw_data_p, &
        uatm_data_p, vatm_data_p
   real (kind=dbl_kind), parameter :: & ! coefficients for Hyland-Wexler Qa
        ps1 = 0.58002206e4_dbl_kind,    & ! (K)
        ps2 = 1.3914993_dbl_kind,       & !
        ps3 = 0.48640239e-1_dbl_kind,   & ! (K^-1)
        ps4 = 0.41764768e-4_dbl_kind,   & ! (K^-2)
        ps5 = 0.14452093e-7_dbl_kind,   & ! (K^-3)
        ps6 = 6.5459673_dbl_kind,       & !
        ws1 = 621.97_dbl_kind,          & ! for saturation mixing ratio
        Pair = 1020._dbl_kind,          & ! Sea level pressure (hPa)
        lapse_rate = 0.0065_dbl_kind      ! (K/m) lapse rate over sea level
   ! for interpolation of hourly data
   integer (kind=int_kind) :: &
        ixm,ixx,ixp , &  ! record numbers for neighboring months
        maxrec      , &  ! maximum record number
        recslot     , &  ! spline slot for current record
        dataloc          ! = 1 for data located in middle of time interval
   ! = 2 for date located at end of time interval
   real (kind=dbl_kind) :: &
        secday    , &
        Qa_pnt
   real (kind=dbl_kind) :: &
        sec1hr           ! number of seconds in 1 hour
   logical (kind=log_kind) :: read1
   integer (kind=int_kind) :: &
        recnum      , & ! record number
        recnum4X        ! record number
   character(len=*), parameter :: subname = '(ISPOL_data)'
   if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'
   call icepack_query_parameters(secday_out=secday)
   call icepack_warnings_flush(nu_diag)
   if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
        file=__FILE__, line=__LINE__)
   if (trim(atm_data_format) == 'nc') then     ! read nc file
      !-------------------------------------------------------------------
      ! data from NCEP_DOE Reanalysis 2 and Bareiss et al 2008
      ! daily data located at the end of the 24-hour period.
      !-------------------------------------------------------------------
      dataloc = 2                          ! data located at end of interval
      sec1hr = secday                      ! seconds in day
      maxrec = 366                         !
      ! current record number
      recnum = int(yday)
      ! Compute record numbers for surrounding data (2 on each side)
      ixm = mod(recnum+maxrec-2,maxrec) + 1
      ixx = mod(recnum-1,       maxrec) + 1
      !     ixp = mod(recnum,         maxrec) + 1
      ! Compute interpolation coefficients
      ! If data is located at the end of the time interval, then the
      !  data value for the current record goes in slot 2
      recslot = 2
      ixp = -99
      call interp_coeff (recnum, recslot, sec1hr, dataloc)
      read1 = .false.
      if (istep==1 .or. oldrecnum .ne. recnum) read1 = .true.
      ! Daily 2m Air temperature 1991
      met_file = F_tair
      fieldname='Tair'
      call read_data_nc_point(read1, 0, fyear, ixm, ixx, ixp, &
           maxrec, met_file, fieldname, Tair_data_p, &
           field_loc_center, field_type_scalar)
      Tair(:,:,:) =  c1intp * Tair_data_p(1) &
           + c2intp * Tair_data_p(2) &
           - lapse_rate*8.0_dbl_kind
      met_file = F_humid
      fieldname='Qa'
      call read_data_nc_point(read1, 0, fyear, ixm, ixx, ixp, &
           maxrec, met_file, fieldname, Qa_data_p, &
           field_loc_center, field_type_scalar)
      Qa_pnt= c1intp * Qa_data_p(1) &
           + c2intp * Qa_data_p(2)
      Qa(:,:,:) = Qa_pnt
      met_file = F_uwind
      fieldname='uatm'
      call read_data_nc_point(read1, 0, fyear, ixm, ixx, ixp, &
           maxrec, met_file, fieldname, uatm_data_p, &
           field_loc_center, field_type_scalar)
      uatm(:,:,:) =  c1intp * uatm_data_p(1) &
           + c2intp * uatm_data_p(2)
      met_file = F_vwind
      fieldname='vatm'
      call read_data_nc_point(read1, 0, fyear, ixm, ixx, ixp, &
           maxrec, met_file, fieldname, vatm_data_p, &
           field_loc_center, field_type_scalar)
      vatm(:,:,:) =  c1intp * vatm_data_p(1) &
           + c2intp * vatm_data_p(2)
      met_file = F_rain
      fieldname='fsnow'
      call read_data_nc_point(read1, 0, fyear, ixm, ixx, ixp, &
           maxrec, met_file, fieldname, fsnow_data_p, &
           field_loc_center, field_type_scalar)
      fsnow(:,:,:) =  (c1intp * fsnow_data_p(1) + &
           c2intp * fsnow_data_p(2))
      !-----------------------------
      !fsw and flw are every 6 hours
      !------------------------------
      dataloc = 2                          ! data located at end of interval
      sec1hr = secday/c4                   ! seconds in 6 hours
      maxrec = 1460                        ! 366*4
      ! current record number
      recnum4X = 4*int(yday) - 3 + int(real(msec,kind=dbl_kind)/sec1hr)
      ! Compute record numbers for surrounding data (2 on each side)
      ixm = mod(recnum4X+maxrec-2,maxrec) + 1
      ixx = mod(recnum4X-1,       maxrec) + 1
      ! Compute interpolation coefficients
      ! If data is located at the end of the time interval, then the
      !  data value for the current record goes in slot 2
      recslot = 2
      ixp = -99
      call interp_coeff (recnum4X, recslot, sec1hr, dataloc)
      read1 = .false.
      if (istep==1 .or. oldrecnum4X .ne. recnum4X) read1 = .true.
      met_file = F_fsw
      fieldname='fsw'
      call read_data_nc_point(read1, 0, fyear, ixm, ixx, ixp, &
           maxrec, met_file, fieldname, fsw_data_p, &
           field_loc_center, field_type_scalar)
      fsw(:,:,:) =  c1intp * fsw_data_p(1) &
           + c2intp * fsw_data_p(2)
      met_file = F_flw
      fieldname='flw'
      call read_data_nc_point(read1, 0, fyear, ixm, ixx, ixp, &
           maxrec, met_file, fieldname, flw_data_p, &
           field_loc_center, field_type_scalar)
      flw(:,:,:) =  c1intp * flw_data_p(1) &
           + c2intp * flw_data_p(2)
   endif  !nc
   !flw   given cldf and Tair  calculated in prepare_atm_forcing
   !-----------------------------
   ! fixed data
   ! May not be needed
   !-----------------------------
   rhoa (:,:,:) = 1.3_dbl_kind ! air density (kg/m^3)
   cldf(:,:,:) =  c1  !0.25_dbl_kind ! cloud fraction
   frain(:,:,:) = c0            ! this is available in hourlymet_rh file
   ! Save record number for next time step
   oldrecnum = recnum
   oldrecnum4X = recnum4X
 end subroutine ISPOL_data

 !=======================================================================
 subroutine ocn_data_ispol_init
   ! Reads NCAR pop ocean forcing data set 'pop_frc_gx1v3_010815.nc'
   ! at the ISPOL location -67.4677N, 310.4375E
   !
   ! For ocean mixed layer-----------------------------units
   !
   ! 1  sst------temperature---------------------------(C)
   ! 2  sss------salinity------------------------------(ppt)
   ! 3  hbl------depth---------------------------------(m)
   ! 4  u--------surface u current---------------------(m/s)
   ! 5  v--------surface v current---------------------(m/s)
   ! 6  dhdx-----surface tilt x direction--------------(m/m)
   ! 7  dhdy-----surface tilt y direction--------------(m/m)
   ! 8  qdp------ocean sub-mixed layer heat flux-------(W/m2)
   !
   ! Fields 4, 5, 6, 7 are on the U-grid; 1, 2, 3, and 8 are
   ! on the T-grid.
   !
   ! authors: Nicole Jeffery, LANL
   !
   use ice_gather_scatter
   use ice_read_write
   integer (kind=int_kind) :: &
        n   , & ! field index
        m       ! month index
   character(char_len) :: &
        vname(nfld) ! variable names to search for in file
   data vname /  &
        'T',      'S',      'hblt',  'U',     'V', &
        'dhdx',   'dhdy',   'qdp' /
   real (kind=dbl_kind) :: &
        work
   integer (kind=int_kind) :: &
        fid         ! file id
   character(len=*), parameter :: subname = '(ocn_data_ispol_init)'
   if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'
   if (my_task == master_task) then
      if (restore_ocn) write (nu_diag,*)  &
           'SST restoring timescale = ',trestore,' days'
      F_sst = trim(ocn_data_dir)//'/'//trim(oceanmixed_file) ! not just sst
      !---------------------------------------------------------------
      ! Read in ocean forcing data from an existing file
      !---------------------------------------------------------------
      write (nu_diag,*) 'ocean mixed layer forcing data file = ', &
           F_sst
   endif ! master_task
   if (trim(ocn_data_format) == 'nc') then
      if (my_task == master_task) then
         call ice_open_nc(F_sst, fid)
      endif ! master_task
      ! Read in ocean forcing data for all 12 months
      do n=1,nfld
         do m=1,12
            ! Note: netCDF does single to double conversion if necessary
            if (n >= 4 .and. n <= 7) then
               call ice_read_nc(fid, m, vname(n), work, debug_forcing, &
                    field_loc_NEcorner, field_type_vector)
            else
               call ice_read_nc(fid, m, vname(n), work, debug_forcing, &
                    field_loc_center, field_type_scalar)
            endif
            ocn_frc_m(:,:,:,n,m) = work
         enddo               ! month loop
      enddo               ! field loop
      if (my_task == master_task) call ice_close_nc(fid)
   else  ! binary format
      call abort_ice (error_message=subname//'new ocean forcing is netcdf only', &
           file=__FILE__, line=__LINE__)
   endif
   !echmod - currents cause Fram outflow to be too large
   ocn_frc_m(:,:,:,4,:) = c0
   ocn_frc_m(:,:,:,5,:) = c0
   !echmod
 end subroutine ocn_data_ispol_init

 !=======================================================================
 subroutine box2001_data_atm
   ! wind fields as in Hunke, JCP 2001
   ! these are defined at the u point
   ! authors: Elizabeth Hunke, LANL
   use ice_domain, only: nblocks, blocks_ice
   use ice_calendar, only: timesecs
   use ice_blocks, only: block, get_block, nx_block, ny_block, nghost
   use ice_flux, only: uatm, vatm, wind, rhoa, strax, stray
   use ice_state, only: aice
   ! local parameters
   integer (kind=int_kind) :: &
        iblk, i,j           ! loop indices
   integer (kind=int_kind) :: &
        iglob(nx_block), & ! global indices
        jglob(ny_block)    ! global indices
   type (block) :: &
        this_block           ! block information for current block
   real (kind=dbl_kind) :: &
        secday, pi , puny, period, pi2, tau
   character(len=*), parameter :: subname = '(box2001_data_atm)'
   if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'
   call icepack_query_parameters(pi_out=pi, pi2_out=pi2, puny_out=puny)
   call icepack_query_parameters(secday_out=secday)
   period = c4*secday
   do iblk = 1, nblocks
      do j = 1, ny_block
         do i = 1, nx_block
            this_block = get_block(blocks_ice(iblk),iblk)
            iglob = this_block%i_glob
            jglob = this_block%j_glob
            uatm(i,j,iblk) = c5 + (sin(pi2*timesecs/period)-c3) &
                 * sin(pi2*real(iglob(i), kind=dbl_kind)  &
                 /real(nx_global,kind=dbl_kind)) &
                 * sin(pi *real(jglob(j), kind=dbl_kind)  &
                 /real(ny_global,kind=dbl_kind))
            vatm(i,j,iblk) = c5 + (sin(pi2*timesecs/period)-c3) &
                 * sin(pi *real(iglob(i), kind=dbl_kind)  &
                 /real(nx_global,kind=dbl_kind)) &
                 * sin(pi2*real(jglob(j), kind=dbl_kind)  &
                 /real(ny_global,kind=dbl_kind))
            ! wind stress
            wind(i,j,iblk) = sqrt(uatm(i,j,iblk)**2 + vatm(i,j,iblk)**2)
            tau = rhoa(i,j,iblk) * 0.0012_dbl_kind * wind(i,j,iblk)
            strax(i,j,iblk) = aice(i,j,iblk) * tau * uatm(i,j,iblk)
            stray(i,j,iblk) = aice(i,j,iblk) * tau * vatm(i,j,iblk)
         enddo
      enddo
   enddo ! nblocks
 end subroutine box2001_data_atm

 !=======================================================================
 subroutine box2001_data_ocn
   ! current fields as in Hunke, JCP 2001
   ! these are defined at the u point
   ! authors: Elizabeth Hunke, LANL
   use ice_domain, only: nblocks, blocks_ice
   use ice_blocks, only: block, get_block, nx_block, ny_block, nghost
   use ice_flux, only: uocn, vocn
   use ice_grid, only: uvm
   ! local parameters
   integer (kind=int_kind) :: &
        iblk, i,j           ! loop indices
   integer (kind=int_kind) :: &
        iglob(nx_block), & ! global indices
        jglob(ny_block)    ! global indices
   type (block) :: &
        this_block           ! block information for current block
   character(len=*), parameter :: subname = '(box2001_data_ocn)'
   if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'
   do iblk = 1, nblocks
      do j = 1, ny_block
         do i = 1, nx_block
            this_block = get_block(blocks_ice(iblk),iblk)
            iglob = this_block%i_glob
            jglob = this_block%j_glob
            ! ocean current
            ! constant in time, could be initialized in ice_flux.F90
            uocn(i,j,iblk) =  p2*real(jglob(j), kind=dbl_kind) &
                 / real(ny_global,kind=dbl_kind) - p1
            vocn(i,j,iblk) = -p2*real(iglob(i), kind=dbl_kind) &
                 / real(nx_global,kind=dbl_kind) + p1
            uocn(i,j,iblk) = uocn(i,j,iblk) * uvm(i,j,iblk)
            vocn(i,j,iblk) = vocn(i,j,iblk) * uvm(i,j,iblk)
         enddo
      enddo
   enddo ! nblocks
 end subroutine box2001_data_ocn

 !=======================================================================
 subroutine uniform_data_atm(dir,spd)
   !     uniform wind fields in some direction
   use ice_domain, only: nblocks
   use ice_blocks, only: nx_block, ny_block, nghost
   use ice_flux, only: uatm, vatm, wind, rhoa, strax, stray
   use ice_state, only: aice

   character(len=*), intent(in) :: dir
   real(kind=dbl_kind), intent(in), optional :: spd ! velocity

   ! local parameters

   integer (kind=int_kind) :: &
        iblk, i,j           ! loop indices

   real (kind=dbl_kind) :: &
        tau, &
        atm_val ! value to use for atm speed

   character(len=*), parameter :: subname = '(uniform_data_atm)'

   if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'

   ! check for optional spd
   if (present(spd)) then
      atm_val = spd
   else
      atm_val = c5 ! default
   endif

   ! wind components
   if (dir == 'NE') then
      uatm = atm_val
      vatm = atm_val
   elseif (dir == 'N') then
      uatm = c0
      vatm = atm_val
   elseif (dir == 'E') then
      uatm = atm_val
      vatm = c0
   elseif (dir == 'S') then
      uatm = c0
      vatm = -atm_val
   elseif (dir == 'W') then
      uatm = -atm_val
      vatm = c0
   else
      call abort_ice (subname//'ERROR: dir unknown, dir = '//trim(dir), &
           file=__FILE__, line=__LINE__)
   endif

   do iblk = 1, nblocks
      do j = 1, ny_block
         do i = 1, nx_block

            ! wind stress
            wind(i,j,iblk) = sqrt(uatm(i,j,iblk)**2 + vatm(i,j,iblk)**2)
            tau = rhoa(i,j,iblk) * 0.0012_dbl_kind * wind(i,j,iblk)
            strax(i,j,iblk) = aice(i,j,iblk) * tau * uatm(i,j,iblk)
            stray(i,j,iblk) = aice(i,j,iblk) * tau * vatm(i,j,iblk)

         enddo
      enddo
   enddo ! nblocks

 end subroutine uniform_data_atm

 !=======================================================================
 subroutine uniform_data_ocn(dir,spd)
   !     uniform current fields in some direction
   use ice_flux, only: uocn, vocn

   character(len=*), intent(in) :: dir

   real(kind=dbl_kind), intent(in), optional :: spd ! velocity

   ! local parameters

   real(kind=dbl_kind) :: &
        ocn_val ! value to use for ocean currents

   character(len=*), parameter :: subname = '(uniform_data_ocn)'

   if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'

   if (present(spd)) then
      ocn_val = spd
   else
      ocn_val = p1 ! default
   endif

   ! ocn components
   if (dir == 'NE') then
      uocn = ocn_val
      vocn = ocn_val
   elseif (dir == 'N') then
      uocn = c0
      vocn = ocn_val
   elseif (dir == 'E') then
      uocn = ocn_val
      vocn = c0
   else
      call abort_ice (subname//'ERROR: dir unknown, dir = '//trim(dir), &
           file=__FILE__, line=__LINE__)
   endif

 end subroutine uniform_data_ocn

 !=======================================================================
 subroutine get_wave_spec
   
   use ice_read_write, only: ice_read_nc_xyf
   use ice_arrays_column, only: wave_spectrum, &
        dwavefreq, wavefreq
   use ice_constants, only: c0
   use ice_domain_size, only: nfreq
   use ice_timers, only: ice_timer_start, ice_timer_stop, timer_fsd

   ! local variables
   integer (kind=int_kind) :: &
        fid                    ! file id for netCDF routines

   real(kind=dbl_kind), dimension(nfreq) :: &
        wave_spectrum_profile  ! wave spectrum

   character(char_len) :: wave_spec_type
   logical (kind=log_kind) :: wave_spec
   character(len=*), parameter :: subname = '(get_wave_spec)'

   if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'

   call ice_timer_start(timer_fsd)

   call icepack_query_parameters(wave_spec_out=wave_spec, &
        wave_spec_type_out=wave_spec_type)
   call icepack_warnings_flush(nu_diag)
   if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
        file=__FILE__, line=__LINE__)

   ! if no wave data is provided, wave_spectrum is zero everywhere
   wave_spectrum(:,:,:,:) = c0
   wave_spec_dir = ocn_data_dir
   debug_forcing = .false.

   ! wave spectrum and frequencies
   if (wave_spec) then
      ! get hardwired frequency bin info and a dummy wave spectrum profile
      ! the latter is used if wave_spec_type == profile
      call icepack_init_wave(nfreq,                 &
           wave_spectrum_profile, &
           wavefreq, dwavefreq)

      ! read more realistic data from a file
      if ((trim(wave_spec_type) == 'constant').OR.(trim(wave_spec_type) == 'random')) then
         if (trim(wave_spec_file(1:4)) == 'unkn') then
            call abort_ice (subname//'ERROR: wave_spec_file '//trim(wave_spec_file), &
                 file=__FILE__, line=__LINE__)
         else
#ifdef USE_NETCDF
            call ice_open_nc(wave_spec_file,fid)
            call ice_read_nc_xyf (fid, 1, 'efreq', wave_spectrum(:,:,:,:), debug_forcing, &
                 field_loc_center, field_type_scalar)
            call ice_close_nc(fid)
#else
            write (nu_diag,*) "wave spectrum file not available, requires cpp USE_NETCDF"
            write (nu_diag,*) "wave spectrum file not available, using default profile"
            call abort_ice (subname//'ERROR: wave_spec_file '//trim(wave_spec_file), &
                 file=__FILE__, line=__LINE__)
#endif
         endif
      endif
   endif

   call ice_timer_stop(timer_fsd)

 end subroutine get_wave_spec

 !=======================================================================
 subroutine init_snowtable
   ! initial snow aging lookup table
   !
   ! Dry snow metamorphism table
   ! snicar_drdt_bst_fit_60_c070416.nc
   ! Flanner (file metadata units mislabelled)
   ! drdsdt0 (10^-6 m/hr) tau (10^-6 m)
   !
   use ice_broadcast, only: broadcast_array, broadcast_scalar
   integer (kind=int_kind) :: &
        idx_T_max   , &  ! Table dimensions
        idx_rhos_max, &
        idx_Tgrd_max
   real (kind=dbl_kind), allocatable :: &
        snowage_rhos (:), &
        snowage_Tgrd (:), &
        snowage_T    (:), &
        snowage_tau  (:,:,:), &
        snowage_kappa(:,:,:), &
        snowage_drdt0(:,:,:)

   ! local variables

   logical (kind=log_kind) :: diag = .false.

   integer (kind=int_kind) :: &
        fid                  ! file id for netCDF file

   character (char_len) :: &
        snw_aging_table, &   ! aging table setting
        fieldname            ! field name in netcdf file

   character(len=*), parameter :: subname = '(init_snowtable)'

   !-----------------------------------------------------------------
   ! read table of snow aging parameters
   !-----------------------------------------------------------------

   call icepack_query_parameters(snw_aging_table_out=snw_aging_table, &
        isnw_rhos_out=idx_rhos_max, isnw_Tgrd_out=idx_Tgrd_max, isnw_T_out=idx_T_max)

   if (my_task == master_task) then
      write (nu_diag,*) ' '
      write (nu_diag,*) 'Snow aging file:', trim(snw_filename)
   endif

   if (snw_aging_table == 'snicar') then
      ! just read the 3d data and pass it in

      call ice_open_nc(snw_filename,fid)

      allocate(snowage_tau  (idx_rhos_max, idx_Tgrd_max, idx_T_max))
      allocate(snowage_kappa(idx_rhos_max, idx_Tgrd_max, idx_T_max))
      allocate(snowage_drdt0(idx_rhos_max, idx_Tgrd_max, idx_T_max))

      fieldname = trim(snw_tau_fname)
      call ice_read_nc(fid,fieldname,snowage_tau,  diag, &
           idx_rhos_max,idx_Tgrd_max,idx_T_max)
      fieldname = trim(snw_kappa_fname)
      call ice_read_nc(fid,fieldname,snowage_kappa,diag, &
           idx_rhos_max,idx_Tgrd_max,idx_T_max)
      fieldname = trim(snw_drdt0_fname)
      call ice_read_nc(fid,fieldname,snowage_drdt0,diag, &
           idx_rhos_max,idx_Tgrd_max,idx_T_max)

      call ice_close_nc(fid)

      call broadcast_array(snowage_tau  , master_task)
      call broadcast_array(snowage_kappa, master_task)
      call broadcast_array(snowage_drdt0, master_task)

      if (my_task == master_task) then
         write(nu_diag,*) subname,'  '
         write(nu_diag,*) subname,' Successfully read snow aging properties:'
         write(nu_diag,*) subname,' snw_aging_table = ',trim(snw_aging_table)
         write(nu_diag,*) subname,' idx_rhos_max = ',idx_rhos_max
         write(nu_diag,*) subname,' idx_Tgrd_max = ',idx_Tgrd_max
         write(nu_diag,*) subname,' idx_T_max    = ',idx_T_max
         write(nu_diag,*) subname,' Data at rhos, Tgrd, T at first index '
         write(nu_diag,*) subname,' snoage_tau (1,1,1)         = ',snowage_tau  (1,1,1)
         write(nu_diag,*) subname,' snoage_kappa (1,1,1)       = ',snowage_kappa(1,1,1)
         write(nu_diag,*) subname,' snoage_drdt0 (1,1,1)       = ',snowage_drdt0(1,1,1)
         write(nu_diag,*) subname,' Data at rhos, Tgrd, T at max index'
         write(nu_diag,*) subname,' snoage_tau (max,max,max)   = ',snowage_tau  (idx_rhos_max, idx_Tgrd_max, idx_T_max)
         write(nu_diag,*) subname,' snoage_kappa (max,max,max) = ',snowage_kappa(idx_rhos_max, idx_Tgrd_max, idx_T_max)
         write(nu_diag,*) subname,' snoage_drdt0 (max,max,max) = ',snowage_drdt0(idx_rhos_max, idx_Tgrd_max, idx_T_max)
      endif

      call icepack_init_parameters(        &
           snowage_tau_in   = snowage_tau,   &
           snowage_kappa_in = snowage_kappa, &
           snowage_drdt0_in = snowage_drdt0 )

      deallocate(snowage_tau)
      deallocate(snowage_kappa)
      deallocate(snowage_drdt0)

   else
      ! read everything and pass it in

      call ice_open_nc(snw_filename,fid)

      fieldname = trim(snw_rhos_fname)
      call ice_get_ncvarsize(fid,fieldname,idx_rhos_max)
      fieldname = trim(snw_Tgrd_fname)
      call ice_get_ncvarsize(fid,fieldname,idx_Tgrd_max)
      fieldname = trim(snw_T_fname)
      call ice_get_ncvarsize(fid,fieldname,idx_T_max)

      call broadcast_scalar(idx_rhos_max, master_task)
      call broadcast_scalar(idx_Tgrd_max, master_task)
      call broadcast_scalar(idx_T_max   , master_task)

      allocate(snowage_rhos (idx_rhos_max))
      allocate(snowage_Tgrd (idx_Tgrd_max))
      allocate(snowage_T    (idx_T_max))
      allocate(snowage_tau  (idx_rhos_max, idx_Tgrd_max, idx_T_max))
      allocate(snowage_kappa(idx_rhos_max, idx_Tgrd_max, idx_T_max))
      allocate(snowage_drdt0(idx_rhos_max, idx_Tgrd_max, idx_T_max))

      fieldname = trim(snw_rhos_fname)
      call ice_read_nc(fid,fieldname,snowage_rhos,  diag, &
           idx_rhos_max)
      fieldname = trim(snw_Tgrd_fname)
      call ice_read_nc(fid,fieldname,snowage_Tgrd,  diag, &
           idx_Tgrd_max)
      fieldname = trim(snw_T_fname)
      call ice_read_nc(fid,fieldname,snowage_T,  diag, &
           idx_T_max)

      fieldname = trim(snw_tau_fname)
      call ice_read_nc(fid,fieldname,snowage_tau,  diag, &
           idx_rhos_max,idx_Tgrd_max,idx_T_max)
      fieldname = trim(snw_kappa_fname)
      call ice_read_nc(fid,fieldname,snowage_kappa,diag, &
           idx_rhos_max,idx_Tgrd_max,idx_T_max)
      fieldname = trim(snw_drdt0_fname)
      call ice_read_nc(fid,fieldname,snowage_drdt0,diag, &
           idx_rhos_max,idx_Tgrd_max,idx_T_max)

      call ice_close_nc(fid)

      call broadcast_array(snowage_rhos , master_task)
      call broadcast_array(snowage_Tgrd , master_task)
      call broadcast_array(snowage_T    , master_task)
      call broadcast_array(snowage_tau  , master_task)
      call broadcast_array(snowage_kappa, master_task)
      call broadcast_array(snowage_drdt0, master_task)

      if (my_task == master_task) then
         write(nu_diag,*) subname,'  '
         write(nu_diag,*) subname,' Successfully read snow aging properties:'
         write(nu_diag,*) subname,' idx_rhos_max = ',idx_rhos_max
         write(nu_diag,*) subname,' idx_Tgrd_max = ',idx_Tgrd_max
         write(nu_diag,*) subname,' idx_T_max    = ',idx_T_max
         write(nu_diag,*) subname,' Data at rhos, Tgrd, T = ',snowage_rhos(1),snowage_Tgrd(1),snowage_T(1)
         write(nu_diag,*) subname,' snoage_tau (1,1,1)         = ',snowage_tau  (1,1,1)
         write(nu_diag,*) subname,' snoage_kappa (1,1,1)       = ',snowage_kappa(1,1,1)
         write(nu_diag,*) subname,' snoage_drdt0 (1,1,1)       = ',snowage_drdt0(1,1,1)
         write(nu_diag,*) subname,' Data at rhos, Tgrd, T = ', &
              snowage_rhos(idx_rhos_max),snowage_Tgrd(idx_Tgrd_max),snowage_T(idx_T_max)
         write(nu_diag,*) subname,' snoage_tau (max,max,max)   = ',snowage_tau  (idx_rhos_max, idx_Tgrd_max, idx_T_max)
         write(nu_diag,*) subname,' snoage_kappa (max,max,max) = ',snowage_kappa(idx_rhos_max, idx_Tgrd_max, idx_T_max)
         write(nu_diag,*) subname,' snoage_drdt0 (max,max,max) = ',snowage_drdt0(idx_rhos_max, idx_Tgrd_max, idx_T_max)
      endif

      call icepack_init_parameters(        &
           isnw_t_in        = idx_T_max,     &
           isnw_Tgrd_in     = idx_Tgrd_max,  &
           isnw_rhos_in     = idx_rhos_max,  &
           snowage_rhos_in  = snowage_rhos,  &
           snowage_Tgrd_in  = snowage_Tgrd,  &
           snowage_T_in     = snowage_T,     &
           snowage_tau_in   = snowage_tau,   &
           snowage_kappa_in = snowage_kappa, &
           snowage_drdt0_in = snowage_drdt0 )

      deallocate(snowage_rhos)
      deallocate(snowage_Tgrd)
      deallocate(snowage_T)
      deallocate(snowage_tau)
      deallocate(snowage_kappa)
      deallocate(snowage_drdt0)

   endif

 end subroutine init_snowtable

end module ice_forcing
