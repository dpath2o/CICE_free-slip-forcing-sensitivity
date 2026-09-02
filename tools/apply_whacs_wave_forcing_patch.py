#!/usr/bin/env python3
"""Patch ice_forcing.F90 for hourly monthly WHACS/CICE25 wave forcing.

The patch is deliberately narrow: it edits only the existing wave-forcing
hooks in cicecore/cicedyn/general/ice_forcing.F90 and leaves the ERA5, AFIM,
lateral-drag and harmonic-tide pathways untouched.

Reference implementation:
  Noah Day, CICE wave-forcing branch
  cicecore/cicedyn/general/ice_forcing.F90

Adaptations for this standalone Antarctic workflow:
  * monthly hourly files rather than one 6-hourly annual file;
  * true current-hour/next-hour interpolation;
  * forcing-cycle wrap at month/year boundaries;
  * dedicated wave record cache (does not reuse atmospheric oldrecnum);
  * explicit tmask guard during wave propagation;
  * static record-1 behaviour retained for ordinary wave_spec_file paths;
  * monthly mode selected by a YYYYMM token in wave_spec_file.

Example forcing_nml value:
  wave_spec_file = '/g/data/gv90/da1339/afim_input/CAWCR/CAWCR_efreq_for_CICE6_YYYYMM.nc'
"""

from __future__ import annotations

import re
from pathlib import Path

TARGET = Path("cicecore/cicedyn/general/ice_forcing.F90")


def replace_once(text: str, old: str, new: str, label: str) -> str:
    n = text.count(old)
    if n != 1:
        raise RuntimeError(f"{label}: expected exactly one anchor, found {n}")
    return text.replace(old, new, 1)


def main() -> None:
    text = TARGET.read_text()

    text = replace_once(
        text,
        "use ice_domain_size  , only: ncat, max_blocks, nx_global, ny_global",
        "use ice_domain_size  , only: ncat, max_blocks, nx_global, ny_global, nfreq",
        "ice_domain_size import",
    )

    text = replace_once(
        text,
        "  ! 4-d (categorical) field values at 2 temporal data points\n"
        "  real (kind=dbl_kind), dimension(:,:,:,:,:), allocatable, public :: topmelt_data, botmelt_data\n",
        "  ! 4-d (categorical) field values at 2 temporal data points\n"
        "  real (kind=dbl_kind), dimension(:,:,:,:,:), allocatable, public :: topmelt_data, botmelt_data\n"
        "  ! Hourly WHACS spectrum at two temporal records for interpolation.\n"
        "  ! Layout: (i,j,frequency,time_slot,block).\n"
        "  real (kind=dbl_kind), dimension(:,:,:,:,:), allocatable :: wave_spectrum_data\n"
        "  integer (kind=int_kind) :: &\n"
        "       wave_cache_year   = -9999, &\n"
        "       wave_cache_month  = -9999, &\n"
        "       wave_cache_recnum = -9999\n",
        "wave cache declaration",
    )

    text = replace_once(
        text,
        "         topmelt_file(ncat), &\n"
        "         botmelt_file(ncat), &\n"
        "         stat=ierr)",
        "         topmelt_file(ncat), &\n"
        "         botmelt_file(ncat), &\n"
        "         wave_spectrum_data(nx_block,ny_block,nfreq,2,max_blocks), &\n"
        "         stat=ierr)",
        "wave buffer allocation",
    )

    new_wave_block = r'''  !=======================================================================
  subroutine get_wave_spec
    ! Populate the Icepack/CICE wave spectrum.  For a wave_spec_file
    ! containing the token YYYYMM, read the shuga hourly monthly WHACS files;
    ! otherwise preserve the historical static record-1 behaviour.
    !
    ! Wave propagation and attenuation follow Noah Day's standalone CICE
    ! wave-forcing implementation, with an explicit ocean-mask guard.
    use ice_read_write, only: ice_read_nc_xyf
    use ice_arrays_column, only: wave_spectrum, dwavefreq, wavefreq
    use ice_timers, only: ice_timer_start, ice_timer_stop, timer_fsd

    integer (kind=int_kind) :: fid
    real(kind=dbl_kind), dimension(nfreq) :: wave_spectrum_profile
    character(char_len) :: wave_spec_type
    logical (kind=log_kind) :: wave_spec
    logical (kind=log_kind), parameter :: wave_propagation = .true.
    character(len=*), parameter :: subname = '(get_wave_spec)'

    call ice_timer_start(timer_fsd)

    call icepack_query_parameters(wave_spec_out=wave_spec, &
         wave_spec_type_out=wave_spec_type)
    call icepack_warnings_flush(nu_diag)
    if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

    wave_spectrum(:,:,:,:) = c0

    if (wave_spec) then
       call icepack_init_wave(nfreq, wave_spectrum_profile, wavefreq, dwavefreq)

       if ((trim(wave_spec_type) == 'constant') .or. &
           (trim(wave_spec_type) == 'random')) then

          if (index(trim(wave_spec_file),'YYYYMM') > 0) then
#ifdef USE_NETCDF
             call wave_spec_data_hourly
#else
             call abort_ice(subname//' ERROR monthly WHACS forcing requires USE_NETCDF', &
                  file=__FILE__, line=__LINE__)
#endif
          else
             ! Backwards-compatible static-file pathway used by earlier tests.
             if (len_trim(wave_spec_file) == 0 .or. &
                  trim(wave_spec_file(1:min(4,len_trim(wave_spec_file)))) == 'unkn') then
                call abort_ice(subname//' ERROR wave_spec_file '//trim(wave_spec_file), &
                     file=__FILE__, line=__LINE__)
             endif
#ifdef USE_NETCDF
             call ice_open_nc(trim(wave_spec_file),fid)
             call ice_read_nc_xyf(fid, 1, 'efreq', wave_spectrum(:,:,:,:), &
                  debug_forcing, field_loc=field_loc_center, field_type=field_type_scalar)
             call ice_close_nc(fid)
#else
             call abort_ice(subname//' ERROR static wave spectrum requires USE_NETCDF', &
                  file=__FILE__, line=__LINE__)
#endif
          endif
       endif

       if (wave_propagation) call propagate_waves
    endif

    call ice_timer_stop(timer_fsd)
  end subroutine get_wave_spec

  !=======================================================================
  subroutine wave_spec_data_hourly
    ! Read hourly shuga WHACS/CICE25 forcing and linearly interpolate between
    ! the current and following hourly records.
    !
    ! wave_spec_file is a template, for example:
    !   /g/data/gv90/da1339/afim_input/CAWCR/
    !   CAWCR_efreq_for_CICE6_YYYYMM.nc
    !
    ! Expected variable:
    !   efreq(time,nfreq,nj,ni)
    !
    ! Files are already on the active CICE T grid and exact CICE25 frequency
    ! grid, so no runtime horizontal or spectral remapping occurs here.
    use ice_read_write, only: ice_read_nc_xyf
    use ice_arrays_column, only: wave_spectrum

    integer (kind=int_kind) :: &
         fid, curr_year, curr_month, next_year, next_month, &
         recnum, next_recnum, maxrec, hour_index, modadj
    real (kind=dbl_kind) :: secday, sec1hr, frac
    character(char_len_long) :: curr_file, next_file
    logical (kind=log_kind) :: new_month
    character(len=*), parameter :: subname = '(wave_spec_data_hourly)'

    call icepack_query_parameters(secday_out=secday)
    call icepack_warnings_flush(nu_diag)
    if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

    sec1hr = secday / 24.0_dbl_kind

    ! Resolve model year onto the configured forcing cycle independently of
    ! atmospheric/ocean forcing record caches.
    if (ycycle == 0) then
       curr_year = myear
    else
       modadj = abs((min(0,myear-fyear_init)/ycycle+1)*ycycle)
       curr_year = fyear_init + mod(myear-fyear_init+modadj,ycycle)
    endif
    curr_month = mmonth

    maxrec = 24 * daymo(curr_month)
    hour_index = int(real(msec,kind=dbl_kind) / sec1hr)
    recnum = 24 * (mday - 1) + hour_index + 1
    recnum = max(1, min(maxrec, recnum))

    frac = (real(msec,kind=dbl_kind) - real(hour_index,kind=dbl_kind)*sec1hr) / sec1hr
    frac = max(c0, min(c1, frac))

    next_year   = curr_year
    next_month  = curr_month
    next_recnum = recnum + 1
    if (recnum >= maxrec) then
       next_recnum = 1
       next_month = curr_month + 1
       if (next_month > 12) then
          next_month = 1
          next_year = curr_year + 1
          if (ycycle > 0 .and. next_year > fyear_final) next_year = fyear_init
       endif
    endif

    call whacs_monthly_wave_file(curr_year, curr_month, curr_file)
    call whacs_monthly_wave_file(next_year, next_month, next_file)

    new_month = (wave_cache_year /= curr_year .or. wave_cache_month /= curr_month)

    ! Read only upon entering a new forcing hour.  Both slots are then reused
    ! by every sub-hourly CICE timestep within that hour.
    if (istep == 1 .or. wave_cache_year /= curr_year .or. &
         wave_cache_month /= curr_month .or. wave_cache_recnum /= recnum) then

       wave_spectrum_data = c0

       call ice_open_nc(trim(curr_file),fid)
       call ice_read_nc_xyf(fid, recnum, 'efreq', wave_spectrum_data(:,:,:,1,:), &
            debug_forcing, field_loc=field_loc_center, field_type=field_type_scalar)

       if (trim(next_file) == trim(curr_file)) then
          call ice_read_nc_xyf(fid, next_recnum, 'efreq', wave_spectrum_data(:,:,:,2,:), &
               debug_forcing, field_loc=field_loc_center, field_type=field_type_scalar)
          call ice_close_nc(fid)
       else
          call ice_close_nc(fid)
          call ice_open_nc(trim(next_file),fid)
          call ice_read_nc_xyf(fid, next_recnum, 'efreq', wave_spectrum_data(:,:,:,2,:), &
               debug_forcing, field_loc=field_loc_center, field_type=field_type_scalar)
          call ice_close_nc(fid)
       endif

       if (new_month .and. my_task == master_task) then
          write(nu_diag,*) subname//' current WHACS file = ', trim(curr_file)
          write(nu_diag,*) subname//' next WHACS file    = ', trim(next_file)
          write(nu_diag,*) subname//' nfreq              = ', nfreq
       endif

       wave_cache_year   = curr_year
       wave_cache_month  = curr_month
       wave_cache_recnum = recnum
    endif

    wave_spectrum(:,:,:,:) = (c1-frac) * wave_spectrum_data(:,:,:,1,:) + &
                                  frac  * wave_spectrum_data(:,:,:,2,:)
    where (wave_spectrum < c0) wave_spectrum = c0

  end subroutine wave_spec_data_hourly

  !=======================================================================
  subroutine whacs_monthly_wave_file(year, month, filename)
    ! Replace the first YYYYMM token in wave_spec_file with the forcing date.
    integer (kind=int_kind), intent(in) :: year, month
    character(char_len_long), intent(out) :: filename
    integer (kind=int_kind) :: ipos, nlen
    character(len=6) :: yyyymm
    character(len=*), parameter :: subname = '(whacs_monthly_wave_file)'

    nlen = len_trim(wave_spec_file)
    ipos = index(trim(wave_spec_file),'YYYYMM')
    if (ipos <= 0) then
       call abort_ice(subname//' ERROR wave_spec_file lacks YYYYMM token: '// &
            trim(wave_spec_file), file=__FILE__, line=__LINE__)
    endif

    write(yyyymm,'(i4.4,i2.2)') year, month
    filename = ' '

    if (ipos == 1 .and. ipos+5 == nlen) then
       filename = yyyymm
    elseif (ipos == 1) then
       filename = yyyymm//trim(wave_spec_file(ipos+6:nlen))
    elseif (ipos+5 == nlen) then
       filename = wave_spec_file(1:ipos-1)//yyyymm
    else
       filename = wave_spec_file(1:ipos-1)//yyyymm// &
            trim(wave_spec_file(ipos+6:nlen))
    endif
  end subroutine whacs_monthly_wave_file

  !=======================================================================
  subroutine propagate_waves
    ! Propagate open-water wave spectra into modeled ice using the empirical
    ! attenuation relation of Meylan, Bennetts & Kohout (2014).  The 10-pass
    ! neighbour structure follows Noah Day's 2025 standalone implementation.
    use ice_grid, only: tlat, tmask
    use ice_domain, only: nblocks, blocks_ice
    use ice_blocks, only: block, get_block
    use ice_arrays_column, only: wave_spectrum, dwavefreq, wavefreq
    use ice_state, only: aice

    type (block) :: this_block
    real(kind=dbl_kind), dimension(nfreq) :: attenuation_rate
    logical (kind=log_kind) :: new_cells_updated
    real(kind=dbl_kind) :: &
         delta_lat, max_delta_lat, local_sig_ht, neighbour_sig_ht, conc_obs
    integer (kind=int_kind) :: &
         i, j, idx_freq, pass, max_passes, iblk, ilo, ihi, jlo, jhi, &
         best_in, best_jn, best_dir, i_n, j_n, idx_d
    integer, dimension(4) :: di = [0, 0, -1, 1]
    integer, dimension(4) :: dj = [1, -1, 0, 0]
    integer, dimension(4) :: dir_code_list = [1, 3, 4, 2]

    conc_obs = 0.70_dbl_kind
    do idx_freq = 1, nfreq
       attenuation_rate(idx_freq) = fn_Attn_MBK(wavefreq(idx_freq)) / conc_obs
    enddo

    max_passes = 10

    ! Treat WHACS as an open-water source field.  Remove direct external wave
    ! energy under modeled ice before propagating attenuated spectra inward.
    do iblk = 1, nblocks
       this_block = get_block(blocks_ice(iblk),iblk)
       ilo = this_block%ilo
       ihi = this_block%ihi
       jlo = this_block%jlo
       jhi = this_block%jhi
       do j = jlo, jhi
          do i = ilo, ihi
             if (.not. tmask(i,j,iblk) .or. aice(i,j,iblk) >= 0.15_dbl_kind) &
                  wave_spectrum(i,j,:,iblk) = c0
          enddo
       enddo
    enddo

    do pass = 1, max_passes
       new_cells_updated = .false.

       do iblk = 1, nblocks
          this_block = get_block(blocks_ice(iblk),iblk)
          ilo = this_block%ilo
          ihi = this_block%ihi
          jlo = this_block%jlo
          jhi = this_block%jhi

          do j = jlo, jhi
             do i = ilo, ihi
                local_sig_ht = c4 * sqrt(sum(wave_spectrum(i,j,:,iblk)*dwavefreq(:)))

                if (tmask(i,j,iblk) .and. aice(i,j,iblk) >= 0.15_dbl_kind .and. &
                     local_sig_ht <= p1) then

                   max_delta_lat = -c1
                   best_in  = i
                   best_jn  = j
                   best_dir = -1

                   do idx_d = 1, 4
                      i_n = i + di(idx_d)
                      j_n = j + dj(idx_d)

                      if (i_n >= ilo .and. i_n <= ihi .and. &
                           j_n >= jlo .and. j_n <= jhi) then
                         neighbour_sig_ht = c4 * sqrt(sum( &
                              wave_spectrum(i_n,j_n,:,iblk)*dwavefreq(:)))

                         if (neighbour_sig_ht > p1) then
                            delta_lat = abs(tlat(i_n,j_n,iblk) - tlat(i,j,iblk))
                            if (delta_lat > max_delta_lat) then
                               max_delta_lat = delta_lat
                               best_in  = i_n
                               best_jn  = j_n
                               best_dir = dir_code_list(idx_d)
                            endif
                         endif
                      endif
                   enddo

                   if (max_delta_lat > c0) then
                      call increment_wave(i, j, iblk, best_in, best_jn, iblk, &
                           aice(best_in,best_jn,iblk), best_dir, attenuation_rate, &
                           wave_spectrum)
                      new_cells_updated = .true.
                   endif
                endif
             enddo
          enddo
       enddo

       if (.not. new_cells_updated) exit
    enddo

  end subroutine propagate_waves

  !=======================================================================
  subroutine increment_wave(i, j, iblk, src_i, src_j, src_iblk, conc, &
       dir_code, attenuation_rate, wave_spectrum)
    use ice_grid, only: HTE, HTN

    integer(kind=int_kind), intent(in) :: &
         i, j, iblk, src_i, src_j, src_iblk, dir_code
    real(kind=dbl_kind), intent(inout) :: wave_spectrum(:,:,:,:)
    real(kind=dbl_kind), intent(in) :: attenuation_rate(:)
    real(kind=dbl_kind), intent(in) :: conc

    integer (kind=int_kind) :: idx_freq
    real(kind=dbl_kind) :: exp_atten, dist

    if (i /= src_i .and. j /= src_j) then
       dist = sqrt(HTE(src_i,src_j,src_iblk)**2 + HTN(i,j,iblk)**2)
    elseif (i /= src_i) then
       dist = HTE(src_i,src_j,src_iblk)
    elseif (j /= src_j) then
       dist = HTN(i,j,iblk)
    else
       dist = c0
    endif

    do idx_freq = 1, nfreq
       exp_atten = exp(-attenuation_rate(idx_freq) * dist * conc)
       wave_spectrum(i,j,idx_freq,iblk) = &
            wave_spectrum(src_i,src_j,idx_freq,src_iblk) * exp_atten
    enddo

  end subroutine increment_wave

  !=======================================================================
  real(kind=dbl_kind) function fn_Attn_MBK(dum_freq)
    ! Meylan, Bennetts & Kohout (2014), GRL, DOI:10.1002/2014GL060809
    real(kind=dbl_kind), intent(in) :: dum_freq
    real(kind=dbl_kind), parameter :: &
         a = 2.12e-3_dbl_kind, &
         b = 4.59e-2_dbl_kind

    fn_Attn_MBK = a*dum_freq**2 + b*dum_freq**4
  end function fn_Attn_MBK
'''

    pattern = re.compile(
        r"  !=======================================================================\n"
        r"  subroutine get_wave_spec\b.*?"
        r"  end subroutine get_wave_spec\n",
        flags=re.DOTALL,
    )
    text, n = pattern.subn(new_wave_block, text, count=1)
    if n != 1:
        raise RuntimeError(f"get_wave_spec replacement: expected 1 match, found {n}")

    TARGET.write_text(text)
    print(f"patched {TARGET}")
    print("monthly mode: wave_spec_file must contain YYYYMM")
    print("next: git diff --check && git diff -- cicecore/cicedyn/general/ice_forcing.F90")


if __name__ == "__main__":
    main()
