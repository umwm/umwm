MODULE umwm_utils
!=======================================================================
USE netcdf

IMPLICIT NONE
!=======================================================================
CONTAINS



SUBROUTINE nc_check(stat)
!=======================================================================
! A convenience routine to check if stat/=NF90_NOERR and abort
! with a helpful message on stdout.
!=======================================================================
IMPLICIT NONE
INTEGER,INTENT(IN) :: stat
IF(stat/=NF90_NOERR)THEN
  WRITE(*,*)'nc_check: Error in NetCDF I/O: '//TRIM(NF90_STRERROR(stat))
  STOP
ENDIF
ENDSUBROUTINE nc_check

!=======================================================================
ENDMODULE umwm_utils
