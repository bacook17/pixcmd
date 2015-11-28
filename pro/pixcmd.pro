PRO PIXCMD, mbin=mbin, zh=zh, ssp=ssp, sfh_tau=sfh_tau, $
            tag=tag,nbin=nbin,nosample=nosample,more_eep=more_eep

  ; Routine to compute simulated CMDs at the pixel level.  
  ; Inputs include the SFH (tau model or SSP), metallicity,
  ; a string to add to the output file, a flag indicating
  ; that the output array is *not* to be downsampled to nsample,
  ; and a flag specifying the type of isochrone file to use
  ; (either the "default" or an enhanced number of EEP points)

  ;a=read_binary('../src/hess.dat',data_dims=[44,51,101],data_type=4)

  nfix    = 200  ;default array size
  nsample = 2E4  ;default number of output pixels
  ;set the "background" flux level.  The faintest M dwarf
  ;has a flux of ~1E-7, so this is very faint.  Used to avoid NaNs
  bkgd = 1E-10

  IF NOT(keyword_set(ssp)) AND NOT(keyword_set(sfh_tau)) THEN ssp=10.
  IF keyword_set(ssp) AND keyword_set(sfh_tau) THEN BEGIN
     print,'make up your damn mind!  cant have SSP & SFH!'
     RETURN
  END

  IF NOT(keyword_set(mbin)) THEN mbin = 1E4
  IF NOT(keyword_set(zh))   THEN zh   = 0.0
  IF keyword_set(ssp) THEN BEGIN
     age  = float(ssp)
     IF age EQ 1.0 THEN age=10.
  ENDIF
  IF NOT(keyword_set(nbin)) THEN nbin = nfix
  IF NOT(keyword_set(tag))  THEN tag  = ''
  IF keyword_set(more_eep)  THEN isoc_tag = 'moreEEP' $
  ELSE isoc_tag='default'
  

  IF zh EQ 0.0 THEN zstr='Z0.0190'
  IF zh GT 0.0 THEN zstr='Z0.0290'
  IF zh LT 0.0 THEN zstr='Z0.0040'

  ;read in the model isochrone and the ACS PSF
  cmdfile = 'tdsp/SSP_MISTv29_BaSeL_Salpeter_'+zstr+'_'+$
            isoc_tag+'.out.cmd'
  a       = cmd_add_variables(cmdfile,kband=kband)
  psf     = mrdfits('~/DATA/HST/psf_f814w_unbinned.fits',0,/sil)
  ;smooth the PSF in angle
  psf     = sqrt( psf * transpose(psf) )
  psf     = psf / total(psf)

  ;set up the isochrones
  IF keyword_set(ssp) THEN BEGIN

     ;SSP
     tt = a[where(a.logage EQ age,cta)]

  ENDIF ELSE BEGIN

     ;create weights for a SFH
     aa = a[uniq(a.logage)].logage
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

     wgt=0.
     FOR i=1,ctw-1 DO BEGIN
        wh  = where(a.logage EQ aa[wha[i]])
        dt  = (10^(aa[i+1]<10.)-10^aa[i-1])/2.
        wgt = wgt+sfh[i]*dt
        tmp = a[wh]
        tmp.logimfweight = tmp.logimfweight+alog10(sfh[i]*dt)
        tt  = (i EQ 1) ? tmp : [tt,tmp]
     ENDFOR

  ENDELSE

  ;remove Infs
  tt = tt[where(finite(tt.logimfweight) EQ 1)]

  ;compactify things into a structure with fewer tags
  sub = {b:0.0,v:0.0,i:0.0,j:0.0,h:0.0,logimfweight:0.0}
  sub = replicate(sub,n_elements(tt))
  struct_assign,tt,sub
  sub.b = tt.acs_f475w
  sub.i = tt.acs_f814w
  sub.j = tt.wfc3_f110w
  sub.h = tt.wfc3_f160w
  tt    = sub

  ;convert from mags to flux (last tag is imfweight)
  FOR i=0,n_tags(tt)-2 DO tt.(i) = 10^(-2./5*tt.(i)) 
  
  ;results structure
  all = replicate({b:0.0,v:0.0,i:0.0,j:0.0,h:0.0},nbin,nbin)

  spawn,'date'
  ;loop over each pixel
  FOR i=0,nbin-1 DO BEGIN
     IF nbin GE 1E3 THEN IF i MOD 1E2 EQ 0 then print,i,nbin-i
     FOR j=0,nbin-1 DO BEGIN
        ;draw from a Poisson distribution
        nn = poidev(10^tt.logimfweight*mbin)
        ;add all of the stars to the pixel
        FOR k=0,n_tags(all)-1 DO $
           all[i,j].(k) = total(nn*tt.(k)) + bkgd
     ENDFOR 
  ENDFOR
  spawn,'date'

  ;convolve with ACS F814W PSF (convolve in flux space)
  ;NB: need to convolve with each PSF separately
  pall = all
  FOR i=0,n_tags(pall)-1 DO $
     pall.(i) = convol(all.(i),psf,/center,/nan,/edge_wrap)

  ;convert back to mags
  FOR i=0,n_tags(all)-1 DO $
     all.(i) = -2.5*alog10(all.(i))
  FOR i=0,n_tags(pall)-1 DO $
     pall.(i) = -2.5*alog10(pall.(i))

  ;-----------------------save results-------------------------;

  IF nbin NE nfix THEN nstr='_N'+strtrim(fix(nbin),2) ELSE nstr=''

  IF keyword_set(sfh_tau) THEN BEGIN
     file = '~/sps/timedomain/pixcmd/results/pixcmd_tau'+$
            strmid(strtrim(sfh_tau,2),0,4)+'_'+zstr+'_Mbin'+$
            strmid(strtrim(alog10(mbin),2),0,4)+nstr+tag+'.fits'
  ENDIF ELSE BEGIN
     file = '~/sps/timedomain/pixcmd/results/pixcmd_t'+$
            strmid(strtrim(age,2),0,4)+'_'+zstr+'_Mbin'+$
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
  ENDIF

  
END
