PRO PLOT_47TUC

  ;[Z/H]=-0.7
  cmdfile = 'tdsp/SSP_MISTv29_BaSeL_Salpeter_Z0.0040_default.out.cmd'
  a       = cmd_add_variables(cmdfile,kband=kband)
  t1      = a[where(a.logage EQ 10.0,cta)]

  ;------------------------------Fornax------------------------------;

  im1  = mrdfits('~/DATA/HST/Fornax/hst_05917_03_wfpc2_f555w_wf_drz.fits',1,h1,/sil)
  tim2 = mrdfits('~/DATA/HST/Fornax/hst_05917_03_wfpc2_f814w_wf_drz.fits',1,h2,/sil)
  ;astrometrically align the two images
  hastrom,tim2,h2,im2,newhdr,h1,interp=2,cubic=-0.5,/silent
  
  zpt1 = 22.545 
  zpt2 = 21.639 +0.41 ;Vega to AB converstion
  dm   = 20.84

  exp1 = 5518.
  exp2 = 7720.
  sn1  = sqrt(im1*exp1)
  sn2  = sqrt(im2*exp2)
  im1  = -2.5*alog10(im1) + zpt1 - dm
  im2  = -2.5*alog10(im2) + zpt2 - dm

  ;generate mag err vs. mag relation in both bands
  wh = randomu(seed,1E5)*n_elements(im1)
  x  = im1[wh]
  y  = 2.5/sn1[wh] ;is this eqn correct?  (OK in small err limit)
  wh = where(finite(x) EQ 1 AND finite(y) EQ 1)
  p1 = poly_fit(x[wh],alog10(y[wh]),1)
  wh = randomu(seed,1E5)*n_elements(im2)
  x  = im2[wh]
  y  = 2.5/sn2[wh]
  wh = where(finite(x) EQ 1 AND finite(y) EQ 1)
  p2 = poly_fit(x[wh],alog10(y[wh]),1) 

  ;mass per pixel
  a    = read_fsps('SSP/SSP_Padova_BaSeL_Salpeter_Z0.0012.out.mags')
  m2l  = 10^a[92].logmass / 10^(2./5*(4.52-a[92].wfpc2_f814w))
  mpix = alog10( m2l * 10^(2./5*(4.52-im2)) )

  xc = 1034.
  yc = 639.
  nx = n_elements(im1[*,0])
  ny = n_elements(im1[0,*])
  xx = fltarr(nx,ny)
  yy = fltarr(nx,ny)
  FOR i=0,ny-1 DO xx[*,i] = findgen(nx)
  FOR i=0,nx-1 DO yy[i,*] = findgen(ny)
  ;distance from center in pixels
  rr = sqrt( (xx-xc)^2 + (yy-yc)^2 )
  wh = where(finite(mpix) EQ 0)
  rr[wh] = 1E5

  dir  = '~/sps/timedomain/pixcmd/results/'
  sage = '10.0'
  m1   = mrdfits(dir+'pixcmd_t'+sage+'_Z0.0040_Mbin1.00.fits',2,/sil)
  m2   = mrdfits(dir+'pixcmd_t'+sage+'_Z0.0040_Mbin2.00.fits',2,/sil)
  whm1 = where(tag_names(m1) EQ 'V')
  whm2 = where(tag_names(m1) EQ 'I')
  ;add photometric uncertainties
  m1 = add_phot_err(m1,whm1,whm2,p1,p2)
  m2 = add_phot_err(m2,whm1,whm2,p1,p2)

  !p.multi=[0,2,2]
  wh = where(rr LE 40 AND rr GE 10,ct2)
  plot,im1[wh]-im2[wh],im2[wh],ps=3,xr=[-1,3],yr=[6,-6]
  oplot,t1.wfpc2_f555w-t1.wfpc2_f814w,t1.wfpc2_f814w

  wh = where(rr LE 100 AND rr GE 40,ct1)
  plot,im1[wh]-im2[wh],im2[wh],ps=3,xr=[-1,3],yr=[6,-6]
  oplot,t1.wfpc2_f555w-t1.wfpc2_f814w,t1.wfpc2_f814w

  ii = ct2 < n_elements(m2)-1
  plot,m2[0:ii].v-m2[0:ii].i,m2[0:ii].i,ps=3,xr=[-1,3],yr=[6,-6]
  oplot,t1.wfpc2_f555w-t1.wfpc2_f814w,t1.wfpc2_f814w
 
  ii = ct1 < n_elements(m1)-1
  plot,m1[0:ii].v-m1[0:ii].i,m1[0:ii].i,ps=3,xr=[-1,3],yr=[6,-6]
  oplot,t1.wfpc2_f555w-t1.wfpc2_f814w,t1.wfpc2_f814w
  
  stop

  ;------------------------------47 Tuc-------------------------------;

  im1  = mrdfits('~/DATA/HST/N104/HST_11677_01_ACS_WFC_F606W_drz.fits',1,h1,/sil)
  tim2 = mrdfits('~/DATA/HST/N104/HST_11677_01_ACS_WFC_F814W_drz.fits',1,h2,/sil)
  ;astrometrically align the two images
  hastrom,tim2,h2,im2,newhdr,h1,interp=2,cubic=-0.5,/silent
  
  zpt1 = 26.493  ; F606W
  zpt2 = 25.9433 ; F814W
  dm   = 13.30   ; Richer et al. 2013

  im1 = -2.5*alog10(im1) + zpt1 - dm
  im2 = -2.5*alog10(im2) + zpt2 - dm

  xx = (im1-im2)[1300:1700,3000:3200]
  yy = im2[1300:1700,3000:3200]
  wh = where(finite(xx) EQ 1)
  xx = xx[wh]
  yy = yy[wh]



END
