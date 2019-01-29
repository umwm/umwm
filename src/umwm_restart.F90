module umwm_restart
!======================================================================!
!                                                                      !
! description: provides read/write routines for umwm restart files     !
!                                                                      !
! contains: restart_read                                               !
!           restart_write                                              !
!                                                                      !
!======================================================================!

use umwm_module
use netcdf
use umwm_util,only:raiseexception

#ifdef MPI
use mpi
use umwm_mpi
#endif 

!======================================================================!
implicit none
!======================================================================!
contains



subroutine restart_read(timestr)
!======================================================================!

! arguments
character(len=19),intent(in) :: timestr

character(len=19) :: timestrnew

character(len=9999) :: filename

integer :: stat
integer :: ncid,ustid,specid

!======================================================================!

if(nproc == 0)then
  write(*,'(a)')'umwm: restart_read: reading restart file for '//timestr
end if

timestrnew = timestr
timestrnew(11:11) = '_'

filename = 'restart/umwmrst_'//timestrnew//'.nc'

stat = nf90_open(trim(filename),nf90_share,ncid)

if(stat /= 0)then
  if(nproc == 0)then
    call raiseexception('abort','restart_read',nf90_strerror(stat))
    stop
  end if
end if

stat = nf90_inq_varid(ncid,'f',specid)
stat = nf90_inq_varid(ncid,'ust',ustid)
stat = nf90_get_var(ncid,specid,e(:,:,istart:iend),&
                    start=[1,1,istart],count=[om,pm,iend-istart+1])
stat = nf90_get_var(ncid,ustid,ustar(istart:iend),&
                    start=[istart],count=[iend-istart+1])
stat = nf90_close(ncid)

#ifdef MPI
call mpi_barrier(MPI_COMM_WORLD,ierr)
#endif

endsubroutine restart_read
!======================================================================!



subroutine restart_write(timestr)
!======================================================================!

! arguments
character(len=19),intent(in) :: timestr
  
integer :: i,o,nn,p
  
integer :: stat
integer :: ncid
integer :: xdimid,fdimid,thdimid
integer :: kid,xid,lonid,latid,freqid,thetaid,ustid,specid
  
real,dimension(im) :: lon_tmp,lat_tmp

!======================================================================!

if(nproc == 0)then

  stat = nf90_create('restart/umwmrst_'//timestr//'.nc',nf90_clobber,ncid)

  stat = nf90_def_dim(ncid,'x',im,xdimid)
  stat = nf90_def_dim(ncid,'f',om,fdimid)
  stat = nf90_def_dim(ncid,'th',pm,thdimid)

  stat = nf90_def_var(ncid,'lon',nf90_float,[xdimid],lonid)
  stat = nf90_put_att(ncid,lonid,name='description',values='longitude')
  stat = nf90_put_att(ncid,lonid,name='units',values='degrees east')

  stat = nf90_def_var(ncid,'lat',nf90_float,[xdimid],latid)
  stat = nf90_put_att(ncid,latid,name='description',values='latitude')
  stat = nf90_put_att(ncid,latid,name='units',values='degrees north')

  stat = nf90_def_var(ncid,'frequency',nf90_float,[fdimid],freqid)
  stat = nf90_put_att(ncid,freqid,name='description',values='frequency')
  stat = nf90_put_att(ncid,freqid,name='units',values='hz')

  stat = nf90_def_var(ncid,'theta',nf90_float,[thdimid],thetaid)
  stat = nf90_put_att(ncid,thetaid,name='description',values='directions')
  stat = nf90_put_att(ncid,thetaid,name='units',values='rad')

  stat = nf90_def_var(ncid,'ust',nf90_float,[xdimid],ustid)
  stat = nf90_put_att(ncid,ustid,name='description',values='friction velocity')
  stat = nf90_put_att(ncid,ustid,name='units',values='m s^-1')

  stat = nf90_def_var(ncid,'wavenumber',nf90_float,[fdimid,xdimid],kid)
  stat = nf90_put_att(ncid,kid,name='description',values='wavenumber')
  stat = nf90_put_att(ncid,kid,name='units',values='rad m^-1')

  stat = nf90_def_var(ncid,'f',nf90_float,[fdimid,thdimid,xdimid],specid)
  stat = nf90_put_att(ncid,specid,name='description',values='wave energy spectrum')
  stat = nf90_put_att(ncid,specid,name='units',values='m^4 rad^-1')

  stat = nf90_enddef(ncid)

  ! fill in lon and lat arrays:
  do i=1,im
    lon_tmp(i) = lon(mi(i),ni(i))
    lat_tmp(i) = lat(mi(i),ni(i))
  end do

  stat = nf90_put_var(ncid,lonid,lon_tmp)
  stat = nf90_put_var(ncid,latid,lat_tmp)
  stat = nf90_put_var(ncid,freqid,f)
  stat = nf90_put_var(ncid,thetaid,th)

  stat = nf90_close(ncid)

end if

! loop over processes in order
do nn = 0,mpisize-1

  ! write to file if it is my turn
  if(nproc == nn)then

    stat = nf90_open('restart/umwmrst_'//timestr//'.nc',nf90_write,ncid)
    stat = nf90_inq_dimid(ncid,'x',xdimid)
    stat = nf90_inq_dimid(ncid,'f',fdimid)
    stat = nf90_inq_dimid(ncid,'th',thdimid)
    stat = nf90_inq_varid(ncid,'f',specid)
    stat = nf90_inq_varid(ncid,'wavenumber',kid)
    stat = nf90_inq_varid(ncid,'ust',ustid)
    stat = nf90_put_var(ncid,specid,e(:,:,istart:iend),start=[1,1,istart],count=[om,pm,iend-istart+1])
    stat = nf90_put_var(ncid,kid,k(:,istart:iend),start=[1,istart],count=[om,iend-istart+1])
    stat = nf90_put_var(ncid,ustid,ustar(istart:iend),start=[istart],count=[iend-istart+1])
    stat = nf90_close(ncid)

  end if

  ! this call ensures that all processes wait for each other
#ifdef MPI
  call mpi_barrier(MPI_COMM_WORLD,ierr)
#endif

end do

if(nproc == 0)then
  write(*,'(a)')'umwm: restart_write: restart written to restart/umwmrst_'//timestr//'.nc'
end if

endsubroutine restart_write
!======================================================================!
end module umwm_restart
