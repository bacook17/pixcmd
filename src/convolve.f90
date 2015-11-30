FUNCTION CONVOLVE(arr,psf,na,nb,np)

  ! Routine to convolve a 2D image with a PSF
  
  USE nrtype
  IMPLICIT NONE

  INTEGER, INTENT(in) :: na,np,nb
  REAL(SP), DIMENSION(na,na,nb), INTENT(in) :: arr
  REAL(SP), DIMENSION(na+2*np,na+2*np,nb) :: padarr
  REAL(SP), DIMENSION(np,np), INTENT(in)  :: psf
  REAL(SP), DIMENSION(na,na,nb) :: convolve
  INTEGER :: i,j,k
  
  !------------------------------------------------------------!

  padarr(np+1:np+na,np+1:np+na,:) = arr

  !first pad Y
  padarr(np+1:np+na,1:np,:) = arr(:,na-np+1:na,:)
  padarr(np+1:np+na,na+np+1:na+2*np,:) = arr(:,1:np,:)

  !then pad X
  padarr(1:np,:,:) = padarr(na+np+1:na+2*np,:,:)
  padarr(na+np+1:na+2*np,:,:) = padarr(np+1:2*np,:,:)

  DO j=1,na
     DO i=1,na
        DO k=1,nb
           convolve(i,j,k) = &
                SUM(padarr(i+np-np/2:i+np+np/2,j+np-np/2:j+np+np/2,k)*psf)
        ENDDO
     ENDDO
  ENDDO
  

END FUNCTION CONVOLVE
