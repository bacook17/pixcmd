PRO PIXCMD, mbin=mbin, zh=zh, ssp=ssp, sfh_tau=sfh_tau, $
            tag=tag,nbin=nbin,nosample=nosample,more_eep=more_eep

  ; Routine to compute simulated CMDs at the pixel level.  
  ; Inputs include the SFH (tau model or SSP), metallicity,
  ; a string to add to the output file, a flag indicating
  ; that the output array is *not* to be downsampled to nsample,
  ; and a flag specifying the type of isochrone file to use
  ; (either the "default" or an enhanced number of EEP points)

  ;a=read_binary('../src/hess.dat',data_dims=[44,51,101],data_type=4)

  nfix    = 512  ;default array size
  nsample = 5E4  ;default number of output pixels
  ;set the "background" flux level.  The faintest M dwarf
  ;has a flux of ~1E-7, so this is very faint.  Used to avoid NaNs
  bkgd = 1E-10

  IF NOT(keyword_set(ssp)) AND NOT(keyword_set(sfh_tau)) THEN ssp=10.
  IF keyword_set(ssp) AND keyword_set(sfh_tau) THEN BEGIN
     print,'make up your damn mind!  cant have SSP & SFH!'
     RETURN
  END

  IF NOT(keyword_set(mbin)) THEN mbin = 1E2
  IF NOT(keyword_set(zh))   THEN zh   = 5
  IF keyword_set(ssp) THEN BEGIN
     age  = float(ssp)
     IF age EQ 1.0 THEN age=10.
  ENDIF
  IF NOT(keyword_set(nbin)) THEN nbin = nfix
  IF NOT(keyword_set(tag))  THEN tag  = ''
  IF keyword_set(more_eep)  THEN isoc_tag = '_moreEEP' $
  ELSE isoc_tag=''
  
  zstr = ['m2.15','m1.13','m0.73','m0.52','m0.22','p0.00','p0.30',$
          'p0.50']

  ;read in the model isochrone and the ACS PSF
  a   = mrdfits('~/pixcmd/isoc/MIST_v29_Z'+zstr[zh]+$
                isoc_tag+'.fits',1,/sil)
  psf = mrdfits('~/DATA/HST/psf/psf_f814w_unbinned.fits',0,/sil)

  ;set up the isochrones
  IF keyword_set(ssp) THEN BEGIN

     ;SSP
     tt = a[where(a.logage EQ age,cta)]

  ENDIF ELSE BEGIN

     ;create weights for a SFH
     aa  = a[uniq(a.logage)].logage
     wha = where(aa LE 10.0,ctw)
     IF sfh_tau GE 10.0 THEN BEGIN
        ;constant SFH
        sfh = fltarr(ctw)*0. + 1/1E10
     ENDIF ELSE BEGIN
        ;tau model SFH
        sfh = exp(-(1E10-10^aa)/(sfh_tau*1E9))
        ;sfh = sfh / int_tabulated(reverse(1E10-10^aa),reverse(sfh))
        wgt = 0.0
        FOR i=1,ctw-1 DO BEGIN
           dt  = (10^(aa[i+1]<10.)-10^aa[i-1])/2.
           wgt = wgt + sfh[i]*dt
        ENDFOR
        sfh = sfh / wgt
     ENDELSE

     wgt  = 0.
     wgti = fltarr(ctw)
     FOR i=1,ctw-1 DO BEGIN
        wh  = where(a.logage EQ aa[wha[i]])
        dt  = (10^(aa[i+1]<10.)-10^aa[i-1])/2.
        wgti[i] = sfh[i]*dt
        wgt = wgt+sfh[i]*dt
        tmp = a[wh]
        tmp.logimfweight = tmp.logimfweight+alog10(sfh[i]*dt)
        tt  = (i EQ 1) ? tmp : [tt,tmp]
     ENDFOR

     ;mass-weighted average age
     ;agrees with analytic results
     mage = total(10^aa[wha]*wgti)

  ENDELSE

  ;remove Infs
  tt = tt[where(finite(tt.logimfweight) EQ 1)]

  ;compactify things into a structure with fewer tags
  sub = {b:0.0,v:0.0,i:0.0,j:0.0,h:0.0,logimfweight:0.0}
  sub = replicate(sub,n_elements(tt))
  sub.b = tt.acs_f475w
  sub.v = tt.acs_f555w
  sub.i = tt.acs_f814w
  sub.j = tt.wfc3_f110w
  sub.h = tt.wfc3_f160w
  sub.logimfweight = tt.logimfweight
  tt    = sub

  ;convert from mags to flux (last tag is imfweight)
  FOR i=0,n_tags(tt)-2 DO tt.(i) = 10^(-2./5*tt.(i)) 
  
  ;results structure
  all = replicate({b:0.0,v:0.0,i:0.0,j:0.0,h:0.0},nbin,nbin)

  ;spawn,'date'

  ;loop over each pixel
  FOR i=0,nbin-1 DO BEGIN
   ;  IF nbin GE 1E3 THEN IF i MOD 1E2 EQ 0 then print,i,nbin-i
     FOR j=0,nbin-1 DO BEGIN
        ;draw from a Poisson distribution
        nn = poidev(10^tt.logimfweight*mbin)
        ;nn = drawn(10^tt.logimfweight*mbin)
        ;add all of the stars to the pixel
        FOR k=0,n_tags(all)-1 DO $
           all[i,j].(k) = total(nn*tt.(k)) + bkgd
     ENDFOR 
  ENDFOR

  ;spawn,'date'

  ;convolve with ACS F814W PSF (convolve in flux space)
  ;NB: should be convolving with each PSF separately,
  ;but the F814W and F475W PSFs are not that different
  pall = all

  ;below, accounting for the fact that the stars are not always
  ;centered within each pixel
  nn = 4.
  FOR j=0,nn-1 DO BEGIN
     FOR k=0,nn-1 DO BEGIN
        FOR i=0,n_tags(pall)-1 DO BEGIN
           im = all[j*nbin/nn:(j+1)*nbin/nn-1,$
                    k*nbin/nn:(k+1)*nbin/nn-1].(i)
           pall[j*nbin/nn:(j+1)*nbin/nn-1,$
                k*nbin/nn:(k+1)*nbin/nn-1].(i) = $
              convol(im,fshift(psf,float(j)/nn,float(k)/nn),$
                  /center,/nan,/edge_wrap)
        ENDFOR
     ENDFOR
  ENDFOR

  ;spawn,'date'


  IF keyword_set(nosample) THEN BEGIN
     tpall = all
     FOR i=0,n_tags(tpall)-1 DO $
        tpall.(i) = convol(all.(i),psf,/center,/nan,/edge_wrap)
     FOR i=0,n_tags(tpall)-1 DO $
        tpall.(i) = -2.5*alog10(tpall.(i))
  ENDIF

  ;convert back to mags
  FOR i=0,n_tags(all)-1 DO $
     all.(i) = -2.5*alog10(all.(i))
  FOR i=0,n_tags(pall)-1 DO $
     pall.(i) = -2.5*alog10(pall.(i))


  ;-----------------------save results-------------------------;

  IF nbin NE nfix THEN nstr='_N'+strtrim(fix(nbin),2) ELSE nstr=''

  IF keyword_set(sfh_tau) THEN BEGIN
     file = '~/pixcmd/results/pixcmd_tau'+$
            strmid(strtrim(sfh_tau,2),0,4)+'_Z'+zstr[zh]+'_Mbin'+$
            strmid(strtrim(alog10(mbin),2),0,4)+nstr+tag+'.fits'
  ENDIF ELSE BEGIN
     file = '~/pixcmd/results/pixcmd_t'+$
            strmid(strtrim(age,2),0,4)+'_Z'+zstr[zh]+'_Mbin'+$
            strmid(strtrim(alog10(mbin),2),0,4)+nstr+tag+'.fits'
  ENDELSE

  ;sub-sample the raw model (we dont need so many independent samples)
  IF NOT(keyword_set(nosample)) THEN BEGIN
     IF n_elements(all) GE nsample  THEN all  = all[0:nsample-1]
     IF n_elements(pall) GE nsample THEN pall = pall[0:nsample-1]
  ENDIF

  mwrfits,all,file,/create,/silent
  mwrfits,pall,file,/silent
  IF keyword_set(nosample) THEN BEGIN
     mwrfits,all.b,file,/silent
     mwrfits,all.i,file,/silent
     mwrfits,pall.b,file,/silent
     mwrfits,pall.i,file,/silent
     mwrfits,tpall.i,file,/silent
  ENDIF

  
END
