PRO PLOT_PCMD, psf=psf, age=age

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
  
  pdir = '~/pixcmd/plots/'
  rdir = '~/pixcmd/results/'

  a  = mrdfits('~/pixcmd/isoc/MIST_v29_Zp0.00.fits',1,/sil)
  t1 = a[where(a.logage EQ age,cta)]
 
  m0 = mrdfits(rdir+'pixcmd_t'+sage+'_Zp0.00_Mbin0.00.fits',pp,/sil)
  m1 = mrdfits(rdir+'pixcmd_t'+sage+'_Zp0.00_Mbin1.00.fits',pp,/sil)
  m2 = mrdfits(rdir+'pixcmd_t'+sage+'_Zp0.00_Mbin2.00.fits',pp,/sil)
  m3 = mrdfits(rdir+'pixcmd_t'+sage+'_Zp0.00_Mbin3.00.fits',pp,/sil)
  m4 = mrdfits(rdir+'pixcmd_t'+sage+'_Zp0.00_Mbin4.00.fits',pp,/sil)
  m5 = mrdfits(rdir+'pixcmd_t'+sage+'_Zp0.00_Mbin5.00.fits',pp,/sil)

  ;B-I vs. I

  file = pdir+'pixcmd_t'+sage+spsf+'_BI.eps'
  begplot,name=file,/col,xsize=7,ysize=5,/quiet,/encap

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

  spawn,'convert -density 500 '+file+' '+str_replace(file,'.eps','.png')

  ;J-H vs. H
  
  file = pdir+'pixcmd_t'+sage+spsf+'_JH.eps'
  begplot,name=file,/col,xsize=7,ysize=5,/quiet,/encap

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

  spawn,'convert -density 500 '+file+' '+str_replace(file,'.eps','.png')


  ;--------------------------------vary SFH------------------------------------;

  FOR i=1,4 DO BEGIN

     mb = float(i)

     mbs = strmid(strtrim(mb,2),0,4)
     IF mb EQ -1. THEN mbs2='0.1'
     IF mb EQ 0. THEN mbs2='1'
     IF mb EQ 1. THEN mbs2='10'
     IF mb GT 1. THEN mbs2='10!U'+strtrim(fix(mb),2)+'!N'

     mt0  = mrdfits(rdir+'pixcmd_t10.0_Zp0.00_Mbin'+mbs+'.fits',pp,/sil)
     mt2  = mrdfits(rdir+'pixcmd_tau2.00_Zp0.00_Mbin'+mbs+'.fits',pp,/sil)
     mt5  = mrdfits(rdir+'pixcmd_tau5.00_Zp0.00_Mbin'+mbs+'.fits',pp,/sil)
     mt10 = mrdfits(rdir+'pixcmd_tau10.0_Zp0.00_Mbin'+mbs+'.fits',pp,/sil)

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


  ;-----------------------vary metallicity----------------------------;


  m1   = mrdfits(rdir+'pixcmd_t'+sage+'_Zp0.00_Mbin1.00.fits',pp,/sil)
  m2   = mrdfits(rdir+'pixcmd_t'+sage+'_Zp0.00_Mbin2.00.fits',pp,/sil)
  m3   = mrdfits(rdir+'pixcmd_t'+sage+'_Zp0.00_Mbin3.00.fits',pp,/sil)
  m4   = mrdfits(rdir+'pixcmd_t'+sage+'_Zp0.00_Mbin4.00.fits',pp,/sil)
  m1zp = mrdfits(rdir+'pixcmd_t'+sage+'_Zp0.50_Mbin1.00.fits',pp,/sil)
  m2zp = mrdfits(rdir+'pixcmd_t'+sage+'_Zp0.50_Mbin2.00.fits',pp,/sil)
  m3zp = mrdfits(rdir+'pixcmd_t'+sage+'_Zp0.50_Mbin3.00.fits',pp,/sil)
  m4zp = mrdfits(rdir+'pixcmd_t'+sage+'_Zp0.50_Mbin4.00.fits',pp,/sil)
  m1zm = mrdfits(rdir+'pixcmd_t'+sage+'_Zm0.52_Mbin1.00.fits',pp,/sil)
  m2zm = mrdfits(rdir+'pixcmd_t'+sage+'_Zm0.52_Mbin2.00.fits',pp,/sil)
  m3zm = mrdfits(rdir+'pixcmd_t'+sage+'_Zm0.52_Mbin3.00.fits',pp,/sil)
  m4zm = mrdfits(rdir+'pixcmd_t'+sage+'_Zm0.52_Mbin4.00.fits',pp,/sil)

  ;Rv=3.1, below is A_X / E(B-V)
  eb = 3.268
  ei = 1.526
  ej = 0.881
  eh = 0.613
  ev = 2.742

  ebv = 0.3

  begplot,name=pdir+'mpix_vs_z.eps',/encap,/col,xsize=8,ysize=7,/quiet

    !p.multi=[0,3,3]
    !p.charsize=1.8
    loadct,2

    xr = [-0.5,4.]
    yr = [3.,-2.]
    xb = (max(xr)-min(xr))/150.
    yb = (max(yr)-min(yr))/150.

    dd = hist_2d(m2.i-m2.h,m2.i,bin1=xb,bin2=yb,$
                 min1=min(xr),max1=max(xr),min2=min(yr),max2=max(yr))
    ddzm = hist_2d(m2zm.i-m2zm.h,m2zm.i,bin1=xb,bin2=yb,$
                   min1=min(xr),max1=max(xr),min2=min(yr),max2=max(yr))
    ddzp = hist_2d(m2zp.i-m2zp.h,m2zp.i,bin1=xb,bin2=yb,$
                   min1=min(xr),max1=max(xr),min2=min(yr),max2=max(yr))

    plotimage,alog10(reverse(ddzm,2)),range=[2,-3],$
              imgxrange=xr,imgyrange=yr,xr=xr,yr=yr,xtit='I-H',ytit='H'
    arrow,0.2,2,0.2+(ei-eh)*ebv,2+eh*ebv,/data,thick=4,hsize=150
    legend,['M!Dpix!N=10!U2!N'],box=0,charsize=1.1,/right,pos=[4,-1.8]
    legend,['[Z/H]=-0.5'],/bottom,/right,box=0,charsize=0.8,pos=[4,2.8]
    plotimage,alog10(reverse(dd,2)),range=[2,-3],$
              imgxrange=xr,imgyrange=yr,xr=xr,yr=yr,xtit='I-H',ytit='H'
    arrow,0.2,2,0.2+(ei-eh)*ebv,2+eh*ebv,/data,thick=4,hsize=150
    legend,['[Z/H]=+0.0'],/bottom,/right,box=0,charsize=0.8,pos=[4,2.8]
    plotimage,alog10(reverse(ddzp,2)),range=[2,-3],$
              imgxrange=xr,imgyrange=yr,xr=xr,yr=yr,xtit='I-H',ytit='H'
    arrow,0.2,2,0.2+(ei-eh)*ebv,2+eh*ebv,/data,thick=4,hsize=150
    legend,['[Z/H]=+0.5'],/bottom,/right,box=0,charsize=0.8,pos=[4,2.8]
  
    xr = [-0.5,3.]
    yr = [0.,-3.]
    xb = (max(xr)-min(xr))/150.
    yb = (max(yr)-min(yr))/150.
   
    dd = hist_2d(m3.i-m3.h,m3.i,bin1=xb,bin2=yb,$
                 min1=min(xr),max1=max(xr),min2=min(yr),max2=max(yr))
    ddzm = hist_2d(m3zm.i-m3zm.h,m3zm.i,bin1=xb,bin2=yb,$
                 min1=min(xr),max1=max(xr),min2=min(yr),max2=max(yr))
    ddzp = hist_2d(m3zp.i-m3zp.h,m3zp.i,bin1=xb,bin2=yb,$
                   min1=min(xr),max1=max(xr),min2=min(yr),max2=max(yr))
    
    plotimage,alog10(reverse(ddzm,2)),range=[2,-3],/preserve,$
              imgxrange=xr,imgyrange=yr,xr=xr,yr=yr,xtit='I-H',ytit='H'
    ;oplot,[0.3,0.4,0.6],[-0.9,-0.65,-0.4],ps=-8,line=2
    arrow,0.2,-0.7,0.2+(ei-eh)*ebv,-0.7+eh*ebv,/data,thick=4,hsize=150
    legend,['M!Dpix!N=10!U3!N'],box=0,charsize=1.1,/right,pos=[3,-2.8]
    legend,['[Z/H]=-0.5'],/bottom,/right,box=0,charsize=0.8,pos=[3,-0.1]
    plotimage,alog10(reverse(dd,2)),range=[2,-3],/preserve,$
              imgxrange=xr,imgyrange=yr,xr=xr,yr=yr,xtit='I-H',ytit='H'
    ;oplot,[0.3,0.4,0.6],[-0.9,-0.65,-0.4],ps=-8,line=2
    arrow,0.2,-0.7,0.2+(ei-eh)*ebv,-0.7+eh*ebv,/data,thick=4,hsize=150
    legend,['[Z/H]=0.0'],/bottom,/right,box=0,charsize=0.8,pos=[3,-0.1]
    plotimage,alog10(reverse(ddzp,2)),range=[2,-3],/preserve,$
              imgxrange=xr,imgyrange=yr,xr=xr,yr=yr,xtit='I-H',ytit='H'
    ;oplot,[0.3,0.4,0.6],[-0.9,-0.65,-0.4],ps=-8,line=2
    arrow,0.2,-0.7,0.2+(ei-eh)*ebv,-0.7+eh*ebv,/data,thick=4,hsize=150
    legend,['[Z/H]=+0.5'],/bottom,/right,box=0,charsize=0.8,pos=[3,-0.1]
    
    xr = [-0.5,2.]
    yr = [-2.5,-5.]
    xb = (max(xr)-min(xr))/150.
    yb = (max(yr)-min(yr))/150.
   
    dd = hist_2d(m4.i-m4.h,m4.i,bin1=xb,bin2=yb,$
                 min1=min(xr),max1=max(xr),min2=min(yr),max2=max(yr))
    ddzm = hist_2d(m4zm.i-m4zm.h,m4zm.i,bin1=xb,bin2=yb,$
                 min1=min(xr),max1=max(xr),min2=min(yr),max2=max(yr))
    ddzp = hist_2d(m4zp.i-m4zp.h,m4zp.i,bin1=xb,bin2=yb,$
                   min1=min(xr),max1=max(xr),min2=min(yr),max2=max(yr))
    
    plotimage,alog10(reverse(ddzm,2)),range=[2,-3],/preserve,$
              imgxrange=xr,imgyrange=yr,xr=xr,yr=yr,xtit='I-H',ytit='H'
    arrow,0.2,-3.5,0.2+(ei-eh)*ebv,-3.5+eh*ebv,/data,thick=4,hsize=150
    legend,['M!Dpix!N=10!U4!N'],box=0,charsize=1.1,/right,pos=[2,-4.9]
    legend,['[Z/H]=-0.5'],/bottom,/right,box=0,charsize=0.8,pos=[2,-2.6]
    plotimage,alog10(reverse(dd,2)),range=[2,-3],/preserve,$
              imgxrange=xr,imgyrange=yr,xr=xr,yr=yr,xtit='I-H',ytit='H'
    arrow,0.2,-3.5,0.2+(ei-eh)*ebv,-3.5+eh*ebv,/data,thick=4,hsize=150
    legend,['[Z/H]=0.0'],/bottom,/right,box=0,charsize=0.8,pos=[2,-2.6]
    plotimage,alog10(reverse(ddzp,2)),range=[2,-3],/preserve,$
              imgxrange=xr,imgyrange=yr,xr=xr,yr=yr,xtit='I-H',ytit='H'
    arrow,0.2,-3.5,0.2+(ei-eh)*ebv,-3.5+eh*ebv,/data,thick=4,hsize=150
    legend,['[Z/H]=+0.5'],/bottom,/right,box=0,charsize=0.8,pos=[2,-2.6]
    
  endplot,/quiet
  simpctable


 
  !P.multi=0
  stop

END

;-----------------------------------------------------------------;
;-----------------------------------------------------------------;

PRO WRITE_PIXCMD_GIF, inarr, file, ndim=ndim

  arr = bytscl(inarr-min(inarr))
  IF keyword_set(ndim) THEN arr = arr[0:ndim-1,0:ndim-1]

  writetifs,arr,file,/gif
  ;plotimage,arr,xs=4,ys=4,/PIXEL_ASPECT_RATIO

END

;-----------------------------------------------------------------;
;-----------------------------------------------------------------;

FUNCTION READ_PCMDIM, file

  bi = mrdfits(file,3,/sil)
  nn = n_elements(bi[*,0])
  im = replicate({bi:0.0,ii:0.0,bp:0.0,ip:0.0,ip0:0.0},nn,nn)
  im.bi = bi
  im.ii = mrdfits(file,4,/sil)
  im.bp = mrdfits(file,5,/sil)
  im.ip = mrdfits(file,6,/sil)
  im.ip0 = mrdfits(file,7,/sil)

  RETURN,im

END

;-----------------------------------------------------------------;
;-----------------------------------------------------------------;

PRO PLOT_PCMD_OVERVIEW

  ;plot images and CMDs as a function of Mpix

  ;sub-space of the image to plot
  nc = 400

  dir  = '~/pixcmd/results/'
  pdir = '~/pixcmd/plots/'
  imm1 = read_pcmdim(dir+'pixcmd_t10.0_Zp0.00_Mbin-1.0_N1000.fits')
  im0  = read_pcmdim(dir+'pixcmd_t10.0_Zp0.00_Mbin0.00_N1000.fits')
  im1  = read_pcmdim(dir+'pixcmd_t10.0_Zp0.00_Mbin1.00_N1000.fits')
  im2  = read_pcmdim(dir+'pixcmd_t10.0_Zp0.00_Mbin2.00_N1000.fits')
  im3  = read_pcmdim(dir+'pixcmd_t10.0_Zp0.00_Mbin3.00_N1000.fits')
  im4  = read_pcmdim(dir+'pixcmd_t10.0_Zp0.00_Mbin4.00_N1000.fits')
  im5  = read_pcmdim(dir+'pixcmd_t10.0_Zp0.00_Mbin5.00_N1000.fits')
  im6  = read_pcmdim(dir+'pixcmd_t10.0_Zp0.00_Mbin6.00_N1000.fits')
 
  write_pixcmd_gif,-2.5*imm1.ip,pdir+'Mbin-1.',ndim=nc
  write_pixcmd_gif,-2.5*im0.ip,pdir+'Mbin0.',ndim=nc
  write_pixcmd_gif,-2.5*im1.ip,pdir+'Mbin1.',ndim=nc
  write_pixcmd_gif,-2.5*im2.ip,pdir+'Mbin2.',ndim=nc
  write_pixcmd_gif,-2.5*im3.ip,pdir+'Mbin3.',ndim=nc
  write_pixcmd_gif,-2.5*im4.ip,pdir+'Mbin4.',ndim=nc
  write_pixcmd_gif,-2.5*im5.ip,pdir+'Mbin5.',ndim=nc
 
  ;source spread out over 30 resolution elements (Hogg 2001)
  ;HST ACS PSF contains ~60% of the power within 10 pixels
  flm1 = -2.5*alog10(median(10^(-2./5*imm1.ip))*30*10)
  fl0  = -2.5*alog10(median(10^(-2./5*im0.ip))*30*10)
  fl1  = -2.5*alog10(median(10^(-2./5*im1.ip))*30*10)
  fl2  = -2.5*alog10(median(10^(-2./5*im2.ip))*30*10)
  fl3  = -2.5*alog10(median(10^(-2./5*im3.ip))*30*10)
  fl4  = -2.5*alog10(median(10^(-2./5*im4.ip))*30*10)
  
  ;SB_I mag/arcsec^2 at DM=18.5
  sbi = median(im0.ip)+18.5+2.5*alog10(0.05^2)
  vmi = 1.0
  sbv  = sbi+vmi
  limv = fl0+vmi

  
  ;read in model isochrone
  a  = mrdfits('~/pixcmd/isoc/MIST_v29_Zp0.00.fits',1,/sil)
  t1 = a[where(a.logage EQ 10.0,cta)]

  begplot,name=pdir+'pixvar_M-1.eps',/encap,xsize=9,ysize=2.8,/col,/quiet
    !p.multi=[0,3,1]
    !p.charsize=1.7
    inarr = -2.5*imm1.ip0
    arr = bytscl(inarr-min(inarr))
    arr = arr[0:nc-1,0:nc-1]
    loadct,0
    plotimage,arr,xs=4,ys=4,/preserve_aspect
    simpctable
    plot,imm1.bi-imm1.ii,imm1.ii,ps=3,xr=[0,4],yr=[10,-8],xs=1,ys=1,$
         xtit='B-I',ytit='I'
    legend,['M!Dpix!N=0.1'],box=0,charsize=1.,pos=[0.0,-7]
    legend,['no PSF'],box=0,/right,charsize=0.9
    oplot,t1.acs_f475w-t1.acs_f814w,t1.acs_f814w,col=!red,thick=2
    plot,imm1.bp-imm1.ip,imm1.ip,ps=3,xr=[0,4],yr=[10,-8],xs=1,ys=1,$
         xtit='B-I',ytit='I'
    oplot,t1.acs_f475w-t1.acs_f814w,t1.acs_f814w,col=!red,thick=2
    oplot,[0,4],[1,1]*flm1,line=2,thick=3,col=!dodgerblue
    legend,['HST PSF'],box=0,/right,charsize=0.9
  endplot,/quiet

  begplot,name=pdir+'pixvar_M0.eps',/encap,xsize=9,ysize=2.8,/col,/quiet
    !p.multi=[0,3,1]
    !p.charsize=1.7
    inarr = -2.5*im0.ip0
    arr = bytscl(inarr-min(inarr))
    arr = arr[0:nc-1,0:nc-1]
    loadct,0
    plotimage,arr,xs=4,ys=4,/preserve_aspect
    simpctable
    plot,im0.bi-im0.ii,im0.ii,ps=3,xr=[0,4],yr=[10,-8],xs=1,ys=1,$
         xtit='B-I',ytit='I'
    oplot,t1.acs_f475w-t1.acs_f814w,t1.acs_f814w,col=!red,thick=2
    legend,['M!Dpix!N=1'],box=0,charsize=1.,pos=[0.0,-7]
    legend,['no PSF'],box=0,/right,charsize=0.9
    plot,im0.bp-im0.ip,im0.ip,ps=3,xr=[0,4],yr=[10,-8],xs=1,ys=1,$
         xtit='B-I',ytit='I'
    oplot,t1.acs_f475w-t1.acs_f814w,t1.acs_f814w,col=!red,thick=2
    oplot,[0,4],[1,1]*fl0,line=2,thick=3,col=!dodgerblue
    legend,['HST PSF'],box=0,/right,charsize=0.9

  endplot,/quiet

  begplot,name=pdir+'pixvar_M1.eps',/encap,xsize=9,ysize=2.8,/col,/quiet
    !p.multi=[0,3,1]
    !p.charsize=1.7
    inarr = -2.5*im1.ip0
    arr = bytscl(inarr-min(inarr))
    arr = arr[0:nc-1,0:nc-1]
    loadct,0
    plotimage,arr,xs=4,ys=4,/preserve_aspect
    simpctable
    plot,im1.bi-im1.ii,im1.ii,ps=3,xr=[0,4],yr=[10,-8],xs=1,ys=1,$
         xtit='B-I',ytit='I'
    oplot,t1.acs_f475w-t1.acs_f814w,t1.acs_f814w,col=!red,thick=2
    legend,['M!Dpix!N=10'],box=0,charsize=1.,pos=[0.0,-7]
    legend,['no PSF'],box=0,/right,charsize=0.9
    plot,im1.bp-im1.ip,im1.ip,ps=3,xr=[0,4],yr=[10,-8],xs=1,ys=1,$
         xtit='B-I',ytit='I'
    oplot,t1.acs_f475w-t1.acs_f814w,t1.acs_f814w,col=!red,thick=2
    oplot,[0,4],[1,1]*fl1,line=2,thick=3,col=!dodgerblue
    legend,['HST PSF'],box=0,/right,charsize=0.9
  endplot,/quiet

  begplot,name=pdir+'pixvar_M2.eps',/encap,xsize=9,ysize=2.8,/col,/quiet
    !p.multi=[0,3,1]
    !p.charsize=1.7
    inarr = -2.5*im2.ip0
    arr = bytscl(inarr-min(inarr))
    arr = arr[0:nc-1,0:nc-1]
    loadct,0
    plotimage,arr,xs=4,ys=4,/preserve_aspect
    simpctable
    plot,im2.bi-im2.ii,im2.ii,ps=3,xr=[0,4],yr=[10,-8],xs=1,ys=1,$
         xtit='B-I',ytit='I'
    oplot,t1.acs_f475w-t1.acs_f814w,t1.acs_f814w,col=!red,thick=2
    legend,['M!Dpix!N=10!U2!N'],box=0,charsize=1.,pos=[0.0,-7]
    legend,['no PSF'],box=0,/right,charsize=0.9
    plot,im2.bp-im2.ip,im2.ip,ps=3,xr=[0,4],yr=[10,-8],xs=1,ys=1,$
         xtit='B-I',ytit='I'
    oplot,t1.acs_f475w-t1.acs_f814w,t1.acs_f814w,col=!red,thick=2
    oplot,[0,4],[1,1]*fl2,line=2,thick=3,col=!dodgerblue
    legend,['HST PSF'],box=0,/right,charsize=0.9
  endplot,/quiet

  begplot,name=pdir+'pixvar_M3.eps',/encap,xsize=9,ysize=2.8,/col,/quiet
    !p.multi=[0,3,1]
    !p.charsize=1.7
    inarr = -2.5*im3.ip0
    arr = bytscl(inarr-min(inarr))
    arr = arr[0:nc-1,0:nc-1]
    loadct,0
    plotimage,arr,xs=4,ys=4,/preserve_aspect
    simpctable
    plot,im3.bi-im3.ii,im3.ii,ps=3,xr=[0,4],yr=[10,-8],xs=1,ys=1,$
         xtit='B-I',ytit='I'
    oplot,t1.acs_f475w-t1.acs_f814w,t1.acs_f814w,col=!red,thick=2
    legend,['M!Dpix!N=10!U3!N'],box=0,charsize=1.,pos=[0.0,-7]
    legend,['no PSF'],box=0,/right,charsize=0.9
    plot,im3.bp-im3.ip,im3.ip,ps=3,xr=[0,4],yr=[10,-8],xs=1,ys=1,$
         xtit='B-I',ytit='I'
    oplot,t1.acs_f475w-t1.acs_f814w,t1.acs_f814w,col=!red,thick=2
    oplot,[0,4],[1,1]*fl3,line=2,thick=3,col=!dodgerblue
    legend,['HST PSF'],box=0,/right,charsize=0.9
  endplot,/quiet

  begplot,name=pdir+'pixvar_M4.eps',/encap,xsize=9,ysize=2.8,/col,/quiet
    !p.multi=[0,3,1]
    !p.charsize=1.7
    inarr = -2.5*im4.ip0
    arr = bytscl(inarr-min(inarr))
    arr = arr[0:nc-1,0:nc-1]
    loadct,0
    plotimage,arr,xs=4,ys=4,/preserve_aspect
    simpctable
    plot,im4.bi-im4.ii,im4.ii,ps=3,xr=[0,4],yr=[10,-8],xs=1,ys=1,$
         xtit='B-I',ytit='I'
    oplot,t1.acs_f475w-t1.acs_f814w,t1.acs_f814w,col=!red,thick=2
    legend,['M!Dpix!N=10!U4!N'],box=0,charsize=1.,pos=[0.0,-7]
    legend,['no PSF'],box=0,/right,charsize=0.9
    plot,im4.bp-im4.ip,im4.ip,ps=3,xr=[0,4],yr=[10,-8],xs=1,ys=1,$
         xtit='B-I',ytit='I'
    oplot,t1.acs_f475w-t1.acs_f814w,t1.acs_f814w,col=!red,thick=2
    oplot,[0,4],[1,1]*fl4,line=2,thick=3,col=!dodgerblue
    legend,['HST PSF'],box=0,/right,charsize=0.9
  endplot,/quiet

  begplot,name=pdir+'pixvar_M5.eps',/encap,xsize=9,ysize=2.8,/col,/quiet
    !p.multi=[0,3,1]
    !p.charsize=1.7
    inarr = -2.5*im5.ip0
    arr = bytscl(inarr-min(inarr))
    arr = arr[0:nc-1,0:nc-1]
    loadct,0
    plotimage,arr,xs=4,ys=4,/preserve_aspect
    simpctable
    plot,im5.bi-im5.ii,im5.ii,ps=3,xr=[0,4],yr=[10,-8],xs=1,ys=1,$
         xtit='B-I',ytit='I'
    legend,['M!Dpix!N=10!U5!N'],box=0,charsize=1.,pos=[0.0,-7]
    legend,['no PSF'],box=0,/right,charsize=0.9
    oplot,t1.acs_f475w-t1.acs_f814w,t1.acs_f814w,col=!red,thick=2
    plot,im5.bp-im5.ip,im5.ip,ps=3,xr=[0,4],yr=[10,-8],xs=1,ys=1,$
         xtit='B-I',ytit='I'
    oplot,t1.acs_f475w-t1.acs_f814w,t1.acs_f814w,col=!red,thick=2
    legend,['HST PSF'],box=0,/right,charsize=0.9
  endplot,/quiet

  begplot,name=pdir+'pixvar_M6.eps',/encap,xsize=9,ysize=2.8,/col,/quiet
    !p.multi=[0,3,1]
    !p.charsize=1.7
    inarr = -2.5*im6.ip0
    arr = bytscl(inarr-min(inarr))
    arr = arr[0:nc-1,0:nc-1]
    loadct,0
    plotimage,arr,xs=4,ys=4,/preserve_aspect
    simpctable
    plot,im6.bi-im6.ii,im6.ii,ps=3,xr=[0,4],yr=[10,-8],xs=1,ys=1,$
         xtit='B-I',ytit='I'
    legend,['M!Dpix!N=10!U5!N'],box=0,charsize=1.,pos=[0.0,-7]
    legend,['no PSF'],box=0,/right,charsize=0.9
    oplot,t1.acs_f475w-t1.acs_f814w,t1.acs_f814w,col=!red,thick=2
    plot,im6.bp-im6.ip,im6.ip,ps=3,xr=[0,4],yr=[10,-8],xs=1,ys=1,$
         xtit='B-I',ytit='I'
    oplot,t1.acs_f475w-t1.acs_f814w,t1.acs_f814w,col=!red,thick=2
    legend,['HST PSF'],box=0,/right,charsize=0.9
  endplot,/quiet

  spawn,'convert -density 500 '+pdir+'pixvar_M-1.eps '+pdir+'pixvar_M-1.png'
  spawn,'convert -density 500 '+pdir+'pixvar_M0.eps '+pdir+'pixvar_M0.png'
  spawn,'convert -density 500 '+pdir+'pixvar_M1.eps '+pdir+'pixvar_M1.png'
  spawn,'convert -density 500 '+pdir+'pixvar_M2.eps '+pdir+'pixvar_M2.png'
  spawn,'convert -density 500 '+pdir+'pixvar_M3.eps '+pdir+'pixvar_M3.png'
  spawn,'convert -density 500 '+pdir+'pixvar_M4.eps '+pdir+'pixvar_M4.png'
  spawn,'convert -density 500 '+pdir+'pixvar_M5.eps '+pdir+'pixvar_M5.png'
  spawn,'convert -density 500 '+pdir+'pixvar_M6.eps '+pdir+'pixvar_M6.png'

  stop

END



