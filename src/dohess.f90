PROGRAM DOHESS

  USE nr, ONLY : locate
  USE pixcmd_vars
  IMPLICIT NONE

  INTEGER :: i,i1,i2,stat,ndat
  CHARACTER(50) :: infile
  REAL(SP), DIMENSION(ndat_max) :: xdat=0.,ydat=0.

  !------------------------------------------------------------!

  IF (IARGC().LT.1) THEN
     infile='m31_bulge'
  ELSE
     CALL GETARG(1,infile)
  ENDIF

  CALL SETUP_MODELS()

  !read in the data
  OPEN(12,FILE=TRIM(PIXCMD_HOME)//'/data/'//TRIM(infile)//'.dat',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  DO i=1,ndat_max
     READ(12,*,IOSTAT=stat) xdat(i),ydat(i)
     IF (stat.NE.0) GOTO 20
  ENDDO
  WRITE(*,*) 'HESS_DATA ERROR: did not finish reading in data file'
  STOP
20 CONTINUE
  ndat = i-1
  CLOSE(12)

  !create a Hess diagram for the data
  hess_data=0.0
  DO i=1,ndat
     IF (xdat(i).LT.xhess(1).OR.xdat(i).GT.xhess(nx).OR.&
          ydat(i).LT.yhess(1).OR.ydat(i).GT.yhess(ny)) CYCLE
     i1 = locate(xhess,xdat(i))
     i2 = locate(yhess,ydat(i))
     hess_data(i1,i2) = hess_data(i1,i2)+1.
  ENDDO

  !save the Hess diagram to file
  OPEN(1,FILE=TRIM(PIXCMD_HOME)//'data/'//TRIM(infile)//'.hess',&
       FORM='UNFORMATTED',STATUS='REPLACE',access='direct',&
       recl=nx*ny*4)
  WRITE(1,rec=1) hess_data
  CLOSE(1)


END PROGRAM DOHESS
