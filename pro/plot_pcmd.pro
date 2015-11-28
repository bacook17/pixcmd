PRO PLOT_PCMD, psf=psf, age=age, mbin=mbin

  IF NOT(keyword_set(age)) THEN age = 10.0
  age  = float(age)
  sage = strmid(strtrim(age,2),0,4)
 
  IF keyword_set(psf) THEN BEGIN
     spsf='_psf' 
     pp=2
  ENDIF ELSE BEGIN
     spsf=''
     pp=1
  ENDELSE
  
  pdir = '~/sps/timedomain/pixcmd/plots/'
  rdir = '~/sps/timedomain/pixcmd/results/'

  cmdfile = 'tdsp/SSP_MISTv29_BaSeL_Salpeter_Z0.0190_default.out.cmd'
  a       = cmd_add_variables(cmdfile,kband=kband)
  t1      = a[where(a.logage EQ age,cta)]
 
  m0 = mrdfits(rdir+'pixcmd_t'+sage+'_Z0.0190_Mbin0.00.fits',pp,/sil)
  m1 = mrdfits(rdir+'pixcmd_t'+sage+'_Z0.0190_Mbin1.00.fits',pp,/sil)
  m2 = mrdfits(rdir+'pixcmd_t'+sage+'_Z0.0190_Mbin2.00.fits',pp,/sil)
  m3 = mrdfits(rdir+'pixcmd_t'+sage+'_Z0.0190_Mbin3.00.fits',pp,/sil)
  m4 = mrdfits(rdir+'pixcmd_t'+sage+'_Z0.0190_Mbin4.00.fits',pp,/sil)
  m5 = mrdfits(rdir+'pixcmd_t'+sage+'_Z0.0190_Mbin5.00.fits',pp,/sil)

  ;B-I vs. I

  begplot,name=pdir+'pixcmd_t'+sage+spsf+'_BI.eps',/col,$
          xsize=7,ysize=5,/quiet,/encap

    !p.multi=[0,3,2]
    !p.charsize=1.4

    plot,[0],yr=[10,-8],ys=1,xr=[0,4],xtit='B-I',ytit='I'
    oplot,t1.acs_f475w-t1.acs_f814w,t1.acs_f814w,col=!red,thick=4
    oplot,m0.b-m0.i,m0.i,ps=8,symsize=0.2
    legend,['log(M!Dpix!N)=0'],box=0,/right,charsize=0.7
    legend,['log(age)='+sage],box=0,charsize=1.1

    plot,[0],yr=[10,-8],ys=1,xr=[0,4],xtit='B-I',ytit='I'
    oplot,t1.acs_f475w-t1.acs_f814w,t1.acs_f814w,col=!red,thick=4
    oplot,m1.b-m1.i,m1.i,ps=8,symsize=0.2
    legend,['log(M!Dpix!N)=1'],box=0,/right,charsize=0.7

    plot,[0],yr=[10,-8],ys=1,xr=[0,4],xtit='B-I',ytit='I'
    oplot,t1.acs_f475w-t1.acs_f814w,t1.acs_f814w,col=!red,thick=4
    oplot,m2.b-m2.i,m2.i,ps=8,symsize=0.2
    legend,['log(M!Dpix!N)=2'],box=0,/right,charsize=0.7

    plot,[0],yr=[10,-8],ys=1,xr=[0,4],xtit='B-I',ytit='I'
    oplot,t1.acs_f475w-t1.acs_f814w,t1.acs_f814w,col=!red,thick=4
    oplot,m3.b-m3.i,m3.i,ps=8,symsize=0.2
    legend,['log(M!Dpix!N)=3'],box=0,/right,charsize=0.7

    plot,[0],yr=[10,-8],ys=1,xr=[0,4],xtit='B-I',ytit='I'
    oplot,t1.acs_f475w-t1.acs_f814w,t1.acs_f814w,col=!red,thick=4
    oplot,m4.b-m4.i,m4.i,ps=8,symsize=0.2
    legend,['log(M!Dpix!N)=4'],box=0,/right,charsize=0.7

    plot,[0],yr=[10,-8],ys=1,xr=[0,4],xtit='B-I',ytit='I'
    oplot,t1.acs_f475w-t1.acs_f814w,t1.acs_f814w,col=!red,thick=4
    oplot,m5.b-m5.i,m5.i,ps=8,symsize=0.2
    legend,['log(M!Dpix!N)=5'],box=0,/right,charsize=0.7

  endplot,/quiet

  ;J-H vs. H

  begplot,name=pdir+'pixcmd_t'+sage+spsf+'_JH.eps',/col,$
          xsize=7,ysize=5,/quiet,/encap

    !p.multi=[0,3,2]
    !p.charsize=1.4

    plot,[0],yr=[10,-8],ys=1,xr=[-0.5,1],xtit='J-H',ytit='H',xs=1
    oplot,t1.wfc3_f110w-t1.wfc3_f160w,t1.wfc3_f160w,col=!red,thick=2
    oplot,m0.j-m0.h,m0.j,ps=8,symsize=0.2
    legend,['log(M!Dpix!N)=0'],box=0,/right,charsize=0.7,/bottom
    legend,['log(age)='+sage],box=0,charsize=0.9

    plot,[0],yr=[10,-8],ys=1,xr=[-0.5,1],xtit='J-H',ytit='H',xs=1
    oplot,t1.wfc3_f110w-t1.wfc3_f160w,t1.wfc3_f160w,col=!red,thick=2
    oplot,m1.j-m1.h,m1.j,ps=8,symsize=0.2
    legend,['log(M!Dpix!N)=1'],box=0,/right,charsize=0.7,/bottom

    plot,[0],yr=[10,-8],ys=1,xr=[-0.5,1],xtit='J-H',ytit='H',xs=1
    oplot,t1.wfc3_f110w-t1.wfc3_f160w,t1.wfc3_f160w,col=!red,thick=2
    oplot,m2.j-m2.h,m2.j,ps=8,symsize=0.2
    legend,['log(M!Dpix!N)=2'],box=0,/right,charsize=0.7,/bottom

    plot,[0],yr=[10,-8],ys=1,xr=[-0.5,1],xtit='J-H',ytit='H',xs=1
    oplot,t1.wfc3_f110w-t1.wfc3_f160w,t1.wfc3_f160w,col=!red,thick=2
    oplot,m3.j-m3.h,m3.j,ps=8,symsize=0.2
    legend,['log(M!Dpix!N)=3'],box=0,/right,charsize=0.7,/bottom

    plot,[0],yr=[10,-8],ys=1,xr=[-0.5,1],xtit='J-H',ytit='H',xs=1
    oplot,t1.wfc3_f110w-t1.wfc3_f160w,t1.wfc3_f160w,col=!red,thick=2
    oplot,m4.j-m4.h,m4.j,ps=8,symsize=0.2
    legend,['log(M!Dpix!N)=4'],box=0,/right,charsize=0.7,/bottom

    plot,[0],yr=[10,-8],ys=1,xr=[-0.5,1],xtit='J-H',ytit='H',xs=1
    oplot,t1.wfc3_f110w-t1.wfc3_f160w,t1.wfc3_f160w,col=!red,thick=2
    oplot,m5.j-m5.h,m5.j,ps=8,symsize=0.2
    legend,['log(M!Dpix!N)=5'],box=0,/right,charsize=0.7,/bottom
 
  endplot,/quiet

  ;--------------------------vary SFH-------------------------------;

  FOR i=0,5 DO BEGIN

     mb = float(i)

     mbs = strmid(strtrim(mb,2),0,4)
     IF mb EQ -1. THEN mbs2='0.1'
     IF mb EQ 0. THEN mbs2='1'
     IF mb EQ 1. THEN mbs2='10'
     IF mb GT 1. THEN mbs2='10!U'+strtrim(fix(mb),2)+'!N'

     mt0  = mrdfits(rdir+'pixcmd_t10.0_Z0.0190_Mbin'+mbs+'.fits',pp,/sil)
     mt2  = mrdfits(rdir+'pixcmd_tau2.00_Z0.0190_Mbin'+mbs+'.fits',pp,/sil)
     mt5  = mrdfits(rdir+'pixcmd_tau5.00_Z0.0190_Mbin'+mbs+'.fits',pp,/sil)
     mt10 = mrdfits(rdir+'pixcmd_tau10.0_Z0.0190_Mbin'+mbs+'.fits',pp,/sil)

     flt0  = -2.5*alog10(median(10^(-2./5*mt0.i))*30*10)
     flt2  = -2.5*alog10(median(10^(-2./5*mt2.i))*30*10)
     flt5  = -2.5*alog10(median(10^(-2./5*mt5.i))*30*10)
     flt10 = -2.5*alog10(median(10^(-2./5*mt10.i))*30*10)
     
     t10 = a[where(a.logage EQ 10.,cta)]
     t8  = a[where(a.logage EQ 8.,cta)]
     
     file = pdir+'pixcmd_varySFH_Mbin'+strtrim(fix(mb),2)+spsf+'_BI.eps'
     begplot,name=file,/col,xsize=7,ysize=6,/quiet,/encap
     
       !p.multi=[0,2,2]
       !p.charsize=1.0

       plot,[0],yr=[10,-8],ys=1,xr=[-1,4],xtit='B-I',ytit='I'
       oplot,mt0.b-mt0.i,mt0.i,ps=8,symsize=0.2
       oplot,t10.acs_f475w-t10.acs_f814w,t10.acs_f814w,col=!red,thick=2
       oplot,t8.acs_f475w-t8.acs_f814w,t8.acs_f814w,col=!red,thick=2
       legend,['M!Dpix!N='+mbs2],box=0,charsize=1.0
       oplot,[-10,10],[1,1]*flt0,line=2,col=!dodgerblue
       legend,['10 Gyr (SSP)'],box=0,charsize=0.8,/right
       
       plot,[0],yr=[10,-8],ys=1,xr=[-1,4],xtit='B-I',ytit='I'
       oplot,mt2.b-mt2.i,mt2.i,ps=8,symsize=0.2
       oplot,t10.acs_f475w-t10.acs_f814w,t10.acs_f814w,col=!red,thick=2
       oplot,t8.acs_f475w-t8.acs_f814w,t8.acs_f814w,col=!red,thick=2
       oplot,[-10,10],[1,1]*flt2,line=2,col=!dodgerblue
       legend,[textoidl('\tau_{SF}=2 Gyr')],box=0,charsize=0.8,/right
       
       plot,[0],yr=[10,-8],ys=1,xr=[-1,4],xtit='B-I',ytit='I'
       oplot,mt5.b-mt5.i,mt5.i,ps=8,symsize=0.2
       oplot,t10.acs_f475w-t10.acs_f814w,t10.acs_f814w,col=!red,thick=2
       oplot,t8.acs_f475w-t8.acs_f814w,t8.acs_f814w,col=!red,thick=2
       oplot,[-10,10],[1,1]*flt5,line=2,col=!dodgerblue
       legend,[textoidl('\tau_{SF}=5 Gyr')],box=0,charsize=0.8,/right
       
       plot,[0],yr=[10,-8],ys=1,xr=[-1,4],xtit='B-I',ytit='I'
       oplot,mt10.b-mt10.i,mt10.i,ps=8,symsize=0.2
       oplot,t10.acs_f475w-t10.acs_f814w,t10.acs_f814w,col=!red,thick=2
       oplot,t8.acs_f475w-t8.acs_f814w,t8.acs_f814w,col=!red,thick=2
       oplot,[-10,10],[1,1]*flt10,line=2,col=!dodgerblue
       legend,[textoidl('SFH=constant')],box=0,charsize=0.8,/right
       
     endplot,/quiet

     spawn,'convert -density 500 '+file+' '+str_replace(file,'.eps','.png')

  ENDFOR 

  stop


  ;------------------vary age and metallicity-----------------------;


  m0 = mrdfits(rdir+'pixcmd_t'+sage+'_Z0.0190_Mbin0.00.fits',pp,/sil)
  m1 = mrdfits(rdir+'pixcmd_t'+sage+'_Z0.0190_Mbin1.00.fits',pp,/sil)
  m2 = mrdfits(rdir+'pixcmd_t'+sage+'_Z0.0190_Mbin2.00.fits',pp,/sil)
  m3 = mrdfits(rdir+'pixcmd_t'+sage+'_Z0.0190_Mbin3.00.fits',pp,/sil)
  m4 = mrdfits(rdir+'pixcmd_t'+sage+'_Z0.0190_Mbin4.00.fits',pp,/sil)
  m5 = mrdfits(rdir+'pixcmd_t'+sage+'_Z0.0190_Mbin5.00.fits',pp,/sil)

  m1zp = mrdfits(rdir+'pixcmd_t'+sage+'_Z0.0290_Mbin1.00.fits',pp,/sil)
  m2zp = mrdfits(rdir+'pixcmd_t'+sage+'_Z0.0290_Mbin2.00.fits',pp,/sil)
  m3zp = mrdfits(rdir+'pixcmd_t'+sage+'_Z0.0290_Mbin3.00.fits',pp,/sil)
  m5zp = mrdfits(rdir+'pixcmd_t'+sage+'_Z0.0290_Mbin5.00.fits',pp,/sil)
  m1zm = mrdfits(rdir+'pixcmd_t'+sage+'_Z0.0080_Mbin1.00.fits',pp,/sil)
  m2zm = mrdfits(rdir+'pixcmd_t'+sage+'_Z0.0080_Mbin2.00.fits',pp,/sil)
  m3zm = mrdfits(rdir+'pixcmd_t'+sage+'_Z0.0080_Mbin3.00.fits',pp,/sil)
  m5zm = mrdfits(rdir+'pixcmd_t'+sage+'_Z0.0080_Mbin5.00.fits',pp,/sil)
 
  m190 = mrdfits(rdir+'pixcmd_t9.00_Z0.0190_Mbin1.00.fits',pp,/sil)
  m180 = mrdfits(rdir+'pixcmd_t8.00_Z0.0190_Mbin1.00.fits',pp,/sil)
  m290 = mrdfits(rdir+'pixcmd_t9.00_Z0.0190_Mbin2.00.fits',pp,/sil)
  m280 = mrdfits(rdir+'pixcmd_t8.00_Z0.0190_Mbin2.00.fits',pp,/sil)
  m390 = mrdfits(rdir+'pixcmd_t9.00_Z0.0190_Mbin3.00.fits',pp,/sil)
  m590 = mrdfits(rdir+'pixcmd_t9.00_Z0.0190_Mbin5.00.fits',pp,/sil)
  m580 = mrdfits(rdir+'pixcmd_t8.00_Z0.0190_Mbin5.00.fits',pp,/sil)

  begplot,name=pdir+'pixcmd_t'+sage+spsf+'_vary1.eps',/col,$
          xsize=6,ysize=5,/quiet,/encap

    !p.multi=[0,2,2]
    !p.charsize=1.

    plot,[0],yr=[10,-8],ys=1,xr=[-1,4],xtit='B-I',ytit='I'
    oplot,m1.b-m1.i,m1.i,ps=8,symsize=0.2
    oplot,m1zp.b-m1zp.i,m1zp.i,ps=8,symsize=0.2,col=!red
    oplot,m1zm.b-m1zm.i,m1zm.i,ps=8,symsize=0.2,col=!blue
    legend,['log(M!Dpix!N)=1'],box=0,charsize=0.7
   
    plot,[0],yr=[10,-8],ys=1,xr=[-1,4],xtit='B-I',ytit='I'
    oplot,m1.b-m1.i,m1.i,ps=8,symsize=0.2
    oplot,m190.b-m190.i,m190.i,ps=8,symsize=0.2,col=!red
    oplot,m180.b-m180.i,m180.i,ps=8,symsize=0.2,col=!green
 
    plot,[0],yr=[10,-8],ys=1,xr=[-0.5,1],xtit='J-H',ytit='H'
    oplot,m1.j-m1.h,m1.h,ps=8,symsize=0.2
    oplot,m1zp.j-m1zp.h,m1zp.h,ps=8,symsize=0.2,col=!red
    oplot,m1zm.j-m1zm.h,m1zm.h,ps=8,symsize=0.2,col=!blue

    plot,[0],yr=[10,-8],ys=1,xr=[-0.5,1],xtit='J-H',ytit='H'
    oplot,m1.j-m1.h,m1.h,ps=8,symsize=0.2
    oplot,m190.j-m190.h,m190.h,ps=8,symsize=0.2,col=!red
    oplot,m180.j-m180.h,m180.h,ps=8,symsize=0.2,col=!green

  endplot,/quiet

  begplot,name=pdir+'pixcmd_t'+sage+spsf+'_vary2.eps',/col,$
          xsize=6,ysize=5,/quiet,/encap

    !p.multi=[0,2,2]
    !p.charsize=1.

    plot,[0],yr=[10,-8],ys=1,xr=[-1,4],xtit='B-I',ytit='I'
    oplot,m2.b-m2.i,m2.i,ps=8,symsize=0.2
    oplot,m2zp.b-m2zp.i,m2zp.i,ps=8,symsize=0.2,col=!red
    oplot,m2zm.b-m2zm.i,m2zm.i,ps=8,symsize=0.2,col=!blue
    legend,['log(M!Dpix!N)=1'],box=0,charsize=0.7
   
    plot,[0],yr=[10,-8],ys=1,xr=[-1,4],xtit='B-I',ytit='I'
    oplot,m2.b-m2.i,m2.i,ps=8,symsize=0.2
    oplot,m290.b-m290.i,m290.i,ps=8,symsize=0.2,col=!red
    oplot,m280.b-m280.i,m280.i,ps=8,symsize=0.2,col=!green
 
    plot,[0],yr=[10,-8],ys=1,xr=[-0.5,1],xtit='J-H',ytit='H'
    oplot,m2.j-m2.h,m2.h,ps=8,symsize=0.2
    oplot,m2zp.j-m2zp.h,m2zp.h,ps=8,symsize=0.2,col=!red
    oplot,m2zm.j-m2zm.h,m2zm.h,ps=8,symsize=0.2,col=!blue

    plot,[0],yr=[10,-8],ys=1,xr=[-0.5,1],xtit='J-H',ytit='H'
    oplot,m2.j-m2.h,m2.h,ps=8,symsize=0.2
    oplot,m290.j-m290.h,m290.h,ps=8,symsize=0.2,col=!red
    oplot,m280.j-m280.h,m280.h,ps=8,symsize=0.2,col=!green
  
  endplot,/quiet


  !P.multi=0
  stop

END

;-----------------------------------------------------------------;
;-----------------------------------------------------------------;

PRO WRITE_PIXCMD_GIF, inarr, file, ndim=ndim

  arr = bytscl(inarr-min(inarr))
  IF keyword_set(ndim) THEN arr = arr[0:ndim-1,0:ndim-1]

  writetifs,arr,file,/gif

  plotimage,arr,xs=4,ys=4,/PIXEL_ASPECT_RATIO

END

;-----------------------------------------------------------------;
;-----------------------------------------------------------------;

PRO PLOT_PIXVARY

  ;plot images and CMDs as a function of Mpix
 
  dir = '~/sps/timedomain/pixcmd/results/'
  pdir = '~/sps/timedomain/pixcmd/plots/'
  fm1 = 'pixcmd_t10.0_Z0.0190_Mbin-1.0_N5000.fits'
  f0  = 'pixcmd_t10.0_Z0.0190_Mbin0.00_N1000.fits'
  f1  = 'pixcmd_t10.0_Z0.0190_Mbin1.00_N1000.fits'
  f2  = 'pixcmd_t10.0_Z0.0190_Mbin2.00_N1000.fits'
  f3  = 'pixcmd_t10.0_Z0.0190_Mbin3.00_N1000.fits'
  f4  = 'pixcmd_t10.0_Z0.0190_Mbin4.00_N1000.fits'
  f5  = 'pixcmd_t10.0_Z0.0190_Mbin5.00_N1000.fits'
  bm1i = mrdfits(dir+fm1,3,/sil)
  im1i = mrdfits(dir+fm1,4,/sil)
  bm1p = mrdfits(dir+fm1,5,/sil)
  im1p = mrdfits(dir+fm1,6,/sil)
  b0i = mrdfits(dir+f0,3,/sil)
  i0i = mrdfits(dir+f0,4,/sil)
  b0p = mrdfits(dir+f0,5,/sil)
  i0p = mrdfits(dir+f0,6,/sil)
  b1i = mrdfits(dir+f1,3,/sil)
  i1i = mrdfits(dir+f1,4,/sil)
  b1p = mrdfits(dir+f1,5,/sil)
  i1p = mrdfits(dir+f1,6,/sil)
  b2i = mrdfits(dir+f2,3,/sil)
  i2i = mrdfits(dir+f2,4,/sil)
  b2p = mrdfits(dir+f2,5,/sil)
  i2p = mrdfits(dir+f2,6,/sil)
  b3i = mrdfits(dir+f3,3,/sil)
  i3i = mrdfits(dir+f3,4,/sil)
  b3p = mrdfits(dir+f3,5,/sil)
  i3p = mrdfits(dir+f3,6,/sil)
  b4i = mrdfits(dir+f4,3,/sil)
  i4i = mrdfits(dir+f4,4,/sil)
  b4p = mrdfits(dir+f4,5,/sil)
  i4p = mrdfits(dir+f4,6,/sil)
  b5i = mrdfits(dir+f5,3,/sil)
  i5i = mrdfits(dir+f5,4,/sil)
  b5p = mrdfits(dir+f5,5,/sil)
  i5p = mrdfits(dir+f5,6,/sil)

  write_pixcmd_gif,-2.5*im1p,pdir+'Mbin-1.',ndim=400
  write_pixcmd_gif,-2.5*i0p,pdir+'Mbin0.',ndim=400
  write_pixcmd_gif,-2.5*i1p,pdir+'Mbin1.',ndim=400
  write_pixcmd_gif,-2.5*i2p,pdir+'Mbin2.',ndim=400
  write_pixcmd_gif,-2.5*i3p,pdir+'Mbin3.',ndim=400
  write_pixcmd_gif,-2.5*i4p,pdir+'Mbin4.',ndim=400
  write_pixcmd_gif,-2.5*i5p,pdir+'Mbin5.',ndim=400
 
  ;source spread out over 30 resolution elements (Hogg 2001)
  ;HST ACS PSF contains ~60% of the power within 10 pixels
  flm1 = -2.5*alog10(median(10^(-2./5*im1p))*30*10)
  fl0  = -2.5*alog10(median(10^(-2./5*i0p))*30*10)
  fl1  = -2.5*alog10(median(10^(-2./5*i1p))*30*10)
  fl2  = -2.5*alog10(median(10^(-2./5*i2p))*30*10)
  fl3  = -2.5*alog10(median(10^(-2./5*i3p))*30*10)
  fl4  = -2.5*alog10(median(10^(-2./5*i4p))*30*10)
  
  ;SB_I mag/arcsec^2 at DM=18.5
  sbi = median(i0p)+18.5+2.5*alog10(0.05^2)
  vmi=1.0
  sbv = sbi +vmi
  limv = fl0+vmi

  
  ;read in model isochrone
  cmdfile = 'tdsp/SSP_MISTv29_BaSeL_Salpeter_Z0.0190_default.out.cmd'
  a       = cmd_add_variables(cmdfile,kband=kband)
  t1      = a[where(a.logage EQ 10.0,cta)]

  begplot,name=pdir+'pixvar_M-1.eps',/encap,xsize=7,ysize=3,/col,/quiet
    !p.multi=[0,2,1]
    !p.charsize=1.1
    plot,bm1i-im1i,im1i,ps=3,xr=[0,4],yr=[10,-8],xs=1,ys=1,xtit='B-I',ytit='I'
    legend,['M!Dpix!N=0.1'],box=0,charsize=1.,pos=[0.0,-7]
    legend,['no PSF'],box=0,/right,charsize=0.9
    oplot,t1.acs_f475w-t1.acs_f814w,t1.acs_f814w,col=!red,thick=2
    plot,bm1p-im1p,im1p,ps=3,xr=[0,4],yr=[10,-8],xs=1,ys=1,xtit='B-I',ytit='I'
    oplot,t1.acs_f475w-t1.acs_f814w,t1.acs_f814w,col=!red,thick=2
    oplot,[0,4],[1,1]*flm1,line=2,thick=3,col=!dodgerblue
    legend,['HST PSF'],box=0,/right,charsize=0.9
  endplot,/quiet

  begplot,name=pdir+'pixvar_M0.eps',/encap,xsize=9,ysize=3,/col,/quiet
    !p.multi=[0,3,1]
    !p.charsize=1.7
    inarr = -2.5*i0p
    arr = bytscl(inarr-min(inarr))
    arr = arr[0:400-1,0:400-1]
    loadct,0
    plotimage,arr,xs=4,ys=4,/preserve_aspect
    simpctable
    plot,b0i-i0i,i0i,ps=3,xr=[0,4],yr=[10,-8],xs=1,ys=1,xtit='B-I',ytit='I'
    oplot,t1.acs_f475w-t1.acs_f814w,t1.acs_f814w,col=!red,thick=2
    legend,['M!Dpix!N=1'],box=0,charsize=1.,pos=[0.0,-7]
    legend,['no PSF'],box=0,/right,charsize=0.9
    plot,b0p-i0p,i0p,ps=3,xr=[0,4],yr=[10,-8],xs=1,ys=1,xtit='B-I',ytit='I'
    oplot,t1.acs_f475w-t1.acs_f814w,t1.acs_f814w,col=!red,thick=2
    oplot,[0,4],[1,1]*fl0,line=2,thick=3,col=!dodgerblue
    legend,['HST PSF'],box=0,/right,charsize=0.9

  endplot,/quiet

  begplot,name=pdir+'pixvar_M1.eps',/encap,xsize=9,ysize=3,/col,/quiet
    !p.multi=[0,3,1]
    !p.charsize=1.7
    inarr = -2.5*i1p
    arr = bytscl(inarr-min(inarr))
    arr = arr[0:400-1,0:400-1]
    loadct,0
    plotimage,arr,xs=4,ys=4,/preserve_aspect
    simpctable
    plot,b1i-i1i,i1i,ps=3,xr=[0,4],yr=[10,-8],xs=1,ys=1,xtit='B-I',ytit='I'
    oplot,t1.acs_f475w-t1.acs_f814w,t1.acs_f814w,col=!red,thick=2
    legend,['M!Dpix!N=10'],box=0,charsize=1.,pos=[0.0,-7]
    legend,['no PSF'],box=0,/right,charsize=0.9
    plot,b1p-i1p,i1p,ps=3,xr=[0,4],yr=[10,-8],xs=1,ys=1,xtit='B-I',ytit='I'
    oplot,t1.acs_f475w-t1.acs_f814w,t1.acs_f814w,col=!red,thick=2
    oplot,[0,4],[1,1]*fl1,line=2,thick=3,col=!dodgerblue
    legend,['HST PSF'],box=0,/right,charsize=0.9
  endplot,/quiet

  begplot,name=pdir+'pixvar_M2.eps',/encap,xsize=9,ysize=3,/col,/quiet
    !p.multi=[0,3,1]
    !p.charsize=1.7
    inarr = -2.5*i2p
    arr = bytscl(inarr-min(inarr))
    arr = arr[0:400-1,0:400-1]
    loadct,0
    plotimage,arr,xs=4,ys=4,/preserve_aspect
    simpctable
    plot,b2i-i2i,i2i,ps=3,xr=[0,4],yr=[10,-8],xs=1,ys=1,xtit='B-I',ytit='I'
    oplot,t1.acs_f475w-t1.acs_f814w,t1.acs_f814w,col=!red,thick=2
    legend,['M!Dpix!N=10!U2!N'],box=0,charsize=1.,pos=[0.0,-7]
    legend,['no PSF'],box=0,/right,charsize=0.9
    plot,b2p-i2p,i2p,ps=3,xr=[0,4],yr=[10,-8],xs=1,ys=1,xtit='B-I',ytit='I'
    oplot,t1.acs_f475w-t1.acs_f814w,t1.acs_f814w,col=!red,thick=2
    oplot,[0,4],[1,1]*fl2,line=2,thick=3,col=!dodgerblue
    legend,['HST PSF'],box=0,/right,charsize=0.9
  endplot,/quiet

  begplot,name=pdir+'pixvar_M3.eps',/encap,xsize=9,ysize=3,/col,/quiet
    !p.multi=[0,3,1]
    !p.charsize=1.7
    inarr = -2.5*i3p
    arr = bytscl(inarr-min(inarr))
    arr = arr[0:400-1,0:400-1]
    loadct,0
    plotimage,arr,xs=4,ys=4,/preserve_aspect
    simpctable
    plot,b3i-i3i,i3i,ps=3,xr=[0,4],yr=[10,-8],xs=1,ys=1,xtit='B-I',ytit='I'
    oplot,t1.acs_f475w-t1.acs_f814w,t1.acs_f814w,col=!red,thick=2
    legend,['M!Dpix!N=10!U3!N'],box=0,charsize=1.,pos=[0.0,-7]
    legend,['no PSF'],box=0,/right,charsize=0.9
    plot,b3p-i3p,i3p,ps=3,xr=[0,4],yr=[10,-8],xs=1,ys=1,xtit='B-I',ytit='I'
    oplot,t1.acs_f475w-t1.acs_f814w,t1.acs_f814w,col=!red,thick=2
    oplot,[0,4],[1,1]*fl3,line=2,thick=3,col=!dodgerblue
    legend,['HST PSF'],box=0,/right,charsize=0.9
  endplot,/quiet

  begplot,name=pdir+'pixvar_M4.eps',/encap,xsize=9,ysize=3,/col,/quiet
    !p.multi=[0,3,1]
    !p.charsize=1.7
    inarr = -2.5*i4p
    arr = bytscl(inarr-min(inarr))
    arr = arr[0:400-1,0:400-1]
    loadct,0
    plotimage,arr,xs=4,ys=4,/preserve_aspect
    simpctable
    plot,b4i-i4i,i4i,ps=3,xr=[0,4],yr=[10,-8],xs=1,ys=1,xtit='B-I',ytit='I'
    oplot,t1.acs_f475w-t1.acs_f814w,t1.acs_f814w,col=!red,thick=2
    legend,['M!Dpix!N=10!U4!N'],box=0,charsize=1.,pos=[0.0,-7]
    legend,['no PSF'],box=0,/right,charsize=0.9
    plot,b4p-i4p,i4p,ps=3,xr=[0,4],yr=[10,-8],xs=1,ys=1,xtit='B-I',ytit='I'
    oplot,t1.acs_f475w-t1.acs_f814w,t1.acs_f814w,col=!red,thick=2
    oplot,[0,4],[1,1]*fl4,line=2,thick=3,col=!dodgerblue
    legend,['HST PSF'],box=0,/right,charsize=0.9
  endplot,/quiet

  begplot,name=pdir+'pixvar_M5.eps',/encap,xsize=9,ysize=3,/col,/quiet
    !p.multi=[0,3,1]
    !p.charsize=1.7
    inarr = -2.5*i5p
    arr = bytscl(inarr-min(inarr))
    arr = arr[0:400-1,0:400-1]
    loadct,0
    plotimage,arr,xs=4,ys=4,/preserve_aspect
    simpctable
    plot,b5i-i5i,i5i,ps=3,xr=[0,4],yr=[10,-8],xs=1,ys=1,xtit='B-I',ytit='I'
    legend,['M!Dpix!N=10!U5!N'],box=0,charsize=1.,pos=[0.0,-7]
    legend,['no PSF'],box=0,/right,charsize=0.9
    oplot,t1.acs_f475w-t1.acs_f814w,t1.acs_f814w,col=!red,thick=2
    plot,b5p-i5p,i5p,ps=3,xr=[0,4],yr=[10,-8],xs=1,ys=1,xtit='B-I',ytit='I'
    oplot,t1.acs_f475w-t1.acs_f814w,t1.acs_f814w,col=!red,thick=2
    legend,['HST PSF'],box=0,/right,charsize=0.9
  endplot,/quiet

  convert -density 500 pixvar_M0.eps pixvar_M0.png
  convert -density 500 pixvar_M1.eps pixvar_M1.png
  convert -density 500 pixvar_M2.eps pixvar_M2.png
  convert -density 500 pixvar_M3.eps pixvar_M3.png
  convert -density 500 pixvar_M4.eps pixvar_M4.png
  convert -density 500 pixvar_M5.eps pixvar_M5.png

  stop

END



