program wrf2umwmgrid
  
  ! Reads longitude, latitude, and landmask from a WRF file 
  ! and writes them into a umwm.grid file.
  
  use netcdf
  use umwm_utils, only: nc_check
  use grid_utils, only: refineField, coarsenField
  
  implicit none

  integer :: idm_src, jdm_src
  integer :: idm_dst, jdm_dst

  integer :: stat, ncid, xdimid, ydimid, lonid, latid, maskid

  character(9999) :: wrfInputFile
  character(2) :: refinementFactorStr

  character(NF90_MAX_NAME) :: dimname

  real, allocatable :: lon_src(:,:), lat_src(:,:)
  real, allocatable :: lon_dst(:,:), lat_dst(:,:)

  integer, allocatable :: mask_src(:,:)
  real, allocatable :: mask_dst(:,:)

  integer :: refinementFactor = 1

  call get_command_argument(1, value=wrfInputFile, status=stat)

  ! Error checking:
  if (stat /=0 ) then
    write(*,*) 'wrf2umwmgrid: Error: Cannot retrieve argument 1'
    stop
  else
    write(*,*) 'wrf2umwmgrid: Processing file ' // trim(wrfInputFile)
  end if

  call get_command_argument(2, value=refinementFactorStr, status=stat)
  if (stat == 0) then
    read(refinementFactorStr,*) refinementFactor
  end if

  if (refinementFactor == 0)then
    write(*,*) 'wrf2umwmgrid: Error: refinementFactor must be a non-zero integer'
    stop
  end if

  if (refinementFactor >= 1) then
    write(*,'(a,I2)') 'wrf2umwmgrid: Refinement factor: ', refinementFactor
  else if(refinementFactor < 0)THEN
    write(*,'(a,I2)') 'wrf2umwmgrid: Refinement factor: 1 /', - refinementFactor
  end if

  write(*,*) 'wrf2umwmgrid: Reading input (WRF) grid.'

  call nc_check(nf90_open(trim(wrfInputFile), NF90_NOWRITE, ncid))
  call nc_check(nf90_inq_dimid(ncid, 'west_east', xdimid))
  call nc_check(nf90_inq_dimid(ncid, 'south_north', ydimid))
  call nc_check(nf90_inquire_dimension(ncid, xdimid, dimname, idm_src))
  call nc_check(nf90_inquire_dimension(ncid, ydimid, dimname, jdm_src))

  write(*,*) 'wrf2umwmgrid: WRF grid idm/jdm:', idm_src, jdm_src

  allocate(lon_src(idm_src, jdm_src))
  allocate(lat_src(idm_src, jdm_src))
  allocate(mask_src(idm_src, jdm_src))

  call nc_check(nf90_inq_varid(ncid, 'XLONG', lonid))
  call nc_check(nf90_inq_varid(ncid, 'XLAT', latid))
  call nc_check(nf90_inq_varid(ncid, 'XLAND', maskid))
  call nc_check(nf90_get_var(ncid, lonid, lon_src))
  call nc_check(nf90_get_var(ncid, latid, lat_src))
  call nc_check(nf90_get_var(ncid, maskid, mask_src))
  call nc_check(nf90_close(ncid))

  where (lon_src < 0) lon_src = lon_src + 360

  write(*,*)'wrf2umwmgrid: WRF grid longitude min/max:', &
            minval(lon_src),maxval(lon_src)
  write(*,*)'wrf2umwmgrid: WRF grid latitude min/max:', &
            minval(lat_src),maxval(lat_src)
  write(*,*)'wrf2umwmgrid: WRF grid number of sea points:', &
            count(mask_src == 0)
  write(*,*)'wrf2umwmgrid: WRF grid number of land points:', &
            count(mask_src == 1)

  ! Adjust land/sea values:
  mask_src = mask_src - 1

  if (refinementFactor > 1) then
    idm_dst = (idm_src - 1) * refinementFactor + 1
    jdm_dst = (jdm_src - 1) * refinementFactor + 1
  else if(refinementFactor < 0) then
    idm_dst = idm_src / abs(refinementFactor)
    jdm_dst = jdm_src / abs(refinementFactor)
  else
    idm_dst = idm_src
    jdm_dst = jdm_src
  end if

  write(*,*) 'wrf2umwmgrid: UMWM grid idm/jdm:', idm_dst, jdm_dst

  allocate(lon_dst(idm_dst, jdm_dst))
  allocate(lat_dst(idm_dst, jdm_dst))
  allocate(mask_dst(idm_dst, jdm_dst))

  lon_dst  = 0
  lat_dst  = 0
  mask_dst = 0

  if (refinementFactor > 1) then ! refine:
    call refineField(lon_src, lon_dst, refinementFactor)
    call refineField(lat_src, lat_dst, refinementFactor)
    call refineField(real(mask_src), mask_dst, refinementFactor)
  else if(refinementFactor == 1) then ! copy:
    lon_dst  = lon_src
    lat_dst  = lat_src
    mask_dst = mask_src
  else ! coarsen:
    call coarsenField(lon_src, lon_dst, - refinementFactor)
    call coarsenField(lat_src, lat_dst, - refinementFactor)
    call coarsenField(real(mask_src), mask_dst, - refinementFactor)
  end if

  write(*,*) 'wrf2umwmgrid: UMWM grid longitude min/max:', &
             minval(lon_dst),maxval(lon_dst)
  write(*,*) 'wrf2umwmgrid: UMWM grid latitude min/max:', &
             minval(lat_dst),maxval(lat_dst)
  write(*,*) 'wrf2umwmgrid: UMWM grid number of sea points:', &
             count(nint(mask_dst) == 1)
  write(*,*) 'wrf2umwmgrid: UMWM grid number of land points:', &
             count(nint(mask_dst) == 0)

  where (lon_dst > 180) lon_dst = lon_dst - 360

  write(*,*) 'wrf2umwmgrid: Writing umwm.grid file.'

  call nc_check(nf90_create('umwm.grid', NF90_CLOBBER, ncid))
  call nc_check(nf90_def_dim(ncid, 'x', idm_dst, xdimid))
  call nc_check(nf90_def_dim(ncid, 'y', jdm_dst, ydimid))
  call nc_check(nf90_def_var(ncid, 'lon', NF90_FLOAT, [xdimid,ydimid], lonid))
  call nc_check(nf90_def_var(ncid, 'lat', NF90_FLOAT, [xdimid,ydimid], latid))
  call nc_check(nf90_def_var(ncid, 'seamask', NF90_INT, [xdimid,ydimid], maskid))
  call nc_check(nf90_put_att(ncid, NF90_GLOBAL, name='history', values='Generated by wrf2umwmgrid'))
  call nc_check(nf90_enddef(ncid))
  call nc_check(nf90_put_var(ncid, lonid, lon_dst))
  call nc_check(nf90_put_var(ncid, latid, lat_dst))
  call nc_check(nf90_put_var(ncid, maskid, nint(mask_dst)))
  call nc_check(nf90_close(ncid))

  deallocate(lon_src, lat_src, mask_src)
  deallocate(lon_dst, lat_dst, mask_dst)

  write(*,*) 'wrf2umwmgrid: Success.'

end program wrf2umwmgrid
