SUBROUTINE INTERP_ISO(zmet,ziso,zniso)

  USE pixcmd_vars; USE nrtype
  USE pixcmd_utils, ONLY : convolve,hist_2d,add_obs_err,myran,mypoidev,drawn
  USE nr, ONLY : locate,ran1
  IMPLICIT NONE

  REAL(SP), INTENT(in) :: zmet
  INTEGER, INTENT(inout) :: zniso
  TYPE(TISO), DIMENSION(niso_max) :: ziso
  INTEGER  :: f,ilo,i,khi1,klo1,khi2,klo2,kloi,khii,kmax,j
  REAL(SP) :: dz,age,age1

  !------------------------------------------------------------!

  ilo = MIN(MAX(locate(zmetarr,zmet),1),nz-1)
  dz  = (zmet-zmetarr(ilo))/(zmetarr(ilo+1)-zmetarr(ilo))

  khi1=0
  klo1=1
  khi2=0
  klo2=1
  kloi=1

  DO i=1,niso_age

     IF (i.GT.1) klo1=khi1+1
     IF (i.GT.1) klo2=khi2+1

     age  = age0+REAL(i-1)*0.2

     !find the relevant isochrone for zmet1
     age1 = age
     DO WHILE(ABS(age1-age).LT.0.01)
        khi1=khi1+1
        age1 = iso(ilo,khi1)%age
     ENDDO
     khi1=khi1-1

     !find the relevant isochrone for zmet2
     age1 = age
     DO WHILE(ABS(age1-age).LT.0.01)
        khi2=khi2+1
        age1 = iso(ilo+1,khi2)%age
     ENDDO
     khi2=khi2-1

     khii = MIN(khi1-klo1+1,khi2-klo2+1)
     kmax = NINT( ((khi1-klo1+1)+(khi2-klo2+1))/2.)
     kmax = MIN(kmax,khi2-klo2)

     ziso(kloi:kloi+kmax-1)%age = age
     IF (khi2.GT.khi1) THEN
        ziso(kloi:kloi+kmax-1)%aind = iso(ilo+1,klo2)%aind
     ELSE
        ziso(kloi:kloi+kmax-1)%aind = iso(ilo,klo1)%aind
     ENDIF

     !interpolate IMF
     ziso(kloi:kloi+khii-1)%imf = &
          10**( (1-dz)*LOG10(iso(ilo,klo1:klo1+khii-1)%imf) + &
          dz*LOG10(iso(ilo+1,klo2:klo2+khii-1)%imf) )
     IF (khi2.GT.khi1) THEN
        ziso(kloi+khii:kloi+kmax-1)%imf = &
             10**(LOG10(iso(ilo+1,klo2+khii:klo2+kmax-1)%imf) - &
             LOG10(iso(ilo+1,klo2+khii-1)%imf) + LOG10(ziso(kloi+khii-1)%imf) )
     ELSE IF (khi1.GT.khi2) THEN
        ziso(kloi+khii:kloi+kmax-1)%imf = &
             10**(LOG10(iso(ilo,klo1+khii:klo1+kmax-1)%imf) - &
             LOG10(iso(ilo,klo1+khii-1)%imf) + LOG10(ziso(kloi+khii-1)%imf) )
     ENDIF

     !interpolate mass
     ziso(kloi:kloi+khii-1)%mass = &
          (1-dz)*iso(ilo,klo1:klo1+khii-1)%mass + &
          dz*iso(ilo+1,klo2:klo2+khii-1)%mass 

     IF (khi2.GT.khi1) THEN
         ziso(kloi+khii:kloi+kmax-1)%mass = &
             iso(ilo+1,klo2+khii:klo2+kmax-1)%mass - &
             iso(ilo+1,klo2+khii-1)%mass + ziso(kloi+khii-1)%mass
     ELSE IF (khi1.GT.khi2) THEN
        ziso(kloi+khii:kloi+kmax-1)%mass = &
             iso(ilo,klo1+khii:klo1+kmax-1)%mass - &
             iso(ilo,klo1+khii-1)%mass + ziso(kloi+khii-1)%mass
     ENDIF

     !interpolate fluxes
     DO f=1,2
        ziso(kloi:kloi+khii-1)%bands(f) = &
             10**( (1-dz)*LOG10(iso(ilo,klo1:klo1+khii-1)%bands(f)) + &
             dz*LOG10(iso(ilo+1,klo2:klo2+khii-1)%bands(f)) )
        IF (khi2.GT.khi1) THEN
           ziso(kloi+khii:kloi+kmax-1)%bands(f) = &
                10**(LOG10(iso(ilo+1,klo2+khii:klo2+kmax-1)%bands(f)) - &
                LOG10(iso(ilo+1,klo2+khii-1)%bands(f)) + &
                LOG10(ziso(kloi+khii-1)%bands(f)) )
        ELSE IF (khi1.GT.khi2) THEN
           ziso(kloi+khii:kloi+kmax-1)%bands(f) = &
                10**(LOG10(iso(ilo,klo1+khii:klo1+kmax-1)%bands(f)) - &
                LOG10(iso(ilo,klo1+khii-1)%bands(f)) + &
                LOG10(ziso(kloi+khii-1)%bands(f)) )
        ENDIF
     ENDDO

     kloi=kloi+kmax

  ENDDO

  zniso = kloi-1


END SUBROUTINE INTERP_ISO
