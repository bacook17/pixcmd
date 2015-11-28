PRO RUN_PIXCMD

  ;grid of metallicity and mass
  FOR z=-1,1 DO BEGIN
     FOR m=0,5 DO BEGIN
        print,m,z
        pixcmd,mbin=10^float(m),ssp=10.,zh=z
     ENDFOR
  ENDFOR

stop

  ;run a giant image to be used for illustration purposes
  pixcmd,mbin=1e5,ssp=10.,nbin=1E3,/nosample,/more_eep
  pixcmd,mbin=1e4,ssp=10.,nbin=1E3,/nosample,/more_eep
  pixcmd,mbin=1e3,ssp=10.,nbin=1E3,/nosample,/more_eep
  pixcmd,mbin=1e2,ssp=10.,nbin=1E3,/nosample,/more_eep
  pixcmd,mbin=1e1,ssp=10.,nbin=1E3,/nosample,/more_eep
  pixcmd,mbin=1e0,ssp=10.,nbin=1E3,/nosample,/more_eep
  pixcmd,mbin=1e-1,ssp=10.,nbin=1E3,/nosample,/more_eep

  stop

  ;run models with a complex SFH
  FOR m=0,5 DO BEGIN
     pixcmd,mbin=10^float(m),ssp=10.,nbin=1E3
     pixcmd,mbin=10^float(m),sfh=1.,nbin=1E3
     pixcmd,mbin=10^float(m),sfh=2.,nbin=1E3
     pixcmd,mbin=10^float(m),sfh=5.,nbin=1E3
     pixcmd,mbin=10^float(m),sfh=10.,nbin=1E3
  ENDFOR

  ;grid of age and mass
  FOR t=0,3 DO BEGIN
     FOR m=0,5 DO BEGIN
        print,m,t
        pixcmd,mbin=10^float(m),ssp=7.+t,/more_eep
        pixcmd,mbin=10^float(m+0.3),ssp=7.+t,/more_eep
        pixcmd,mbin=10^float(m+0.5),ssp=7.+t,/more_eep
     ENDFOR
  ENDFOR



END
