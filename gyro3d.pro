PRO gyro3d ,result $
            ,NCube=NCube,r0=r0,dt=dt $
            ,lat=lat,lon=lon,incl=incl,azim=azim,length=length $
            ,n_Sections=n_Sections,date=date $
            ,delta=delta,energy=energy, Nel=Nel, map=map $
            ,bmag=bmag,np=np,temperature=temperature,pitchangle=m $
            ,freq=freq, ID=ID,geo=geo,plot=plot,show=show,data=data $
            ,idl=idl,coord=coord,cylinder=cylinder,f2tb=f2tb
;+
; gyro3d ,result $
;            ,NCube=NCube,r0=r0,dt=dt $
;            ,lat=lat,lon=lon,incl=incl,azim=azim,length=length $
;            ,n_Sections=n_Sections $
;            ,delta=delta,energy=energy, Nel=Nel, map=map $
;            ,bmag=bmag,np=np,temperature=temperature,pitchangle=m $
;            ,freq=freq, ID=ID,geo=geo,plot=plot,show=show,data=data,idl=idl
; IDL keyword: use IDL version of gyro, not the C version.
;
; ;; Model parameters:
; ;; ncube, n_sections
;
; ;; Loop geometry parameters:
; ;; LAT, LON, AZIM, INCL, LENGTH, R0, dt;
;
; ;; Source model
; ;; B(s), np(s), temperature(s)
;
; ;; Electron parameters:
; ;; delta, energy, pitchangle M, Nel
;
; ;; MW parameters:
; ;; freq[Hz]
;-

 if not keyword_set(length) then length=2e9
 if not keyword_set(ncube) then ncube=64
 if not keyword_set(r0) then r0=4
 if not keyword_set(delta) then delta=3.0
 if not keyword_set(energy) then energy=[10,5e3]
 if not keyword_set(nel) then nel=1e7
 if not keyword_set(m) then m=0.0
 if not keyword_set(bmag) then bmag=500.
 if not keyword_set(np) then np=1e9
 if not keyword_set(temperature) then temperature=10e6
 if not keyword_set(n_sections) then n_sections=10 
 IF NOT keyword_set(freq) THEN freq=10.^interpol(alog10([1,35.0]),15)*1e9 

 n_freq = n_elements(freq) 

if keyword_set(cylinder) then BEGIN
;; overrides the loop geometry for testing

build_cylinder,coord,loopmap,NCube=ncube,r0=r0 $
            ,length=length,angle=cylinder,bmag=bmag $
            ,geo=geo,n_sections=n_sections,date=date,map_struct=1 $
            ,_EXTRA=_EXTRA

n_sections=2

endif 

 jo = dblarr(n_freq,n_sections)
 jx = dblarr(n_freq,n_sections)
 ko = dblarr(n_freq,n_sections)
 kx = dblarr(n_freq,n_sections)

 IF n_elements(nel) ne n_sections THEN BEGIN 
    IF n_elements(Nel) EQ 1 THEN Nel1=replicate(Nel,n_Sections) ELSE Nel1=interpol(Nel,n_sections)
 ENDIF ELSE nel1=nel

 IF n_elements(np) NE n_sections THEN BEGIN  
    IF n_elements(np) EQ 1 THEN np1=replicate(np,n_Sections) ELSE np1=interpol(np,n_sections)
 ENDIF ELSE np1=np

 IF n_elements(bmag) NE n_sections THEN BEGIN  
    IF n_elements(bmag) EQ 1 THEN bmag1=replicate(bmag,n_Sections) ELSE bmag1=interpol(bmag,n_sections)
 ENDIF ELSE bmag1=bmag

 IF n_elements(temperature) NE n_sections THEN BEGIN  
    IF n_elements(temperature) EQ 1 THEN $
     temp=replicate(temperature,n_Sections) ELSE temp=interpol(temperature,n_sections)
 ENDIF ELSE temp=temperature

 if ~keyword_set(cylinder) then  build_loop,coord,loopmap,NCube=NCube,r0=r0,dt=dt $
            ,lat=lat,lon=lon,incl=incl,azim=azim,length=length $
            ,geo=geo,n_Sections=n_Sections                 $
            ,map_struct=1,box=box,quiet=quiet,bmag=bmag1 $
            ,date=date,radiuspx=radiuspx,px2cm=px2cm $
            ,_EXTRA=_EXTRA

 FOR i=0,n_sections-1 DO BEGIN 

    IF keyword_set(idl) THEN BEGIN  
       ;; change to gyro C++ later
       gs,freq,flux $
          ,delta=delta $
          ,energy=energy $
          ,nel=nel1[i]/max(nel1) $
          ,m=m $
          ,bmag=bmag1[i] $
          ,np=np1[i] $
          ,angle=geo.angle[i] $
          ,e1d=e1d,e2d=e2d,a1d=a1d,a2d=a2d $
          ,/bessel,/quiet
;       ,anor=anor $
;       ,ntot=ntot 
    ENDIF ELSE BEGIN 
       gyro,freq,flux $
            ,delta=delta $
            ,energy=energy $
            ,nel=nel1[i]/max(nel1) $
            ,m=m $
            ,bmag=bmag1[i] $
            ,np=np1[i] $
            ,angle=geo.angle[i] $
            ,e1d=e1d,e2d=e2d,a1d=a1d,a2d=a2d,quiet=0
    ENDELSE 
    
    jo[0:n_freq-1,i] = e1d
    jx[0:n_freq-1,i] = e2d
    ko[0:n_freq-1,i] = a1d
    kx[0:n_freq-1,i] = a2d  

    ;; plot,freq,e1d/a1d*(1-exp(-a1d*nel[i]*geo.cross_section_radius*2)),/xl,/yl
    ;; oplot,freq,e2d/a2d*(1-exp(-a2d*nel[i]*geo.cross_section_radius*2)),col='ff'x
    ;;  print,nel[i],nel[i]*geo.cross_section_radius*2

 ENDFOR 

 tr3d,coord,freq,jo,jx,ko,kx,nel1,Ncube $
      ,temperature=temp,np=np1,geo=geo $ 
      ,map=mapdata,flux=flux,mo=mo,mx=mx,f2tb=f2tb,pd=pd
 
 map=replicate(loopmap,n_freq)
 mapo=replicate(loopmap,n_freq)
 mapx=replicate(loopmap,n_freq)
 pol=replicate(loopmap,n_freq)
 map.data=mapdata
 mapo.data=mo
 mapx.data=mx
 pp=where(map.data gt 0.)
 pol.data=pd;/map.data

 ;; total number of electrons
 ntot=total(geo.sector_vol*nel1)
 vol=total(geo.sector_vol)

 IF ~keyword_set(ID) THEN ID='ID_string'
 result = {freq:freq,flux:flux,map:map,ID:ID,timestamp:systime(),ntot:ntot,vol:vol,mo:mapo,mx:mapx,pol:pol}

;; maps to Tb
 for i=0,n_elements(result.freq)-1 do BEGIN
 result.map[i].data = result.map[i].data * f2tb[i]
 result.mx[i].data = result.mx[i].data * f2tb[i]
 result.mo[i].data = result.mo[i].data * f2tb[i]
 result.pol[i].data = result.pol[i].data * f2tb[i]
 endfor

 IF keyword_set(plot) THEN BEGIN
    
    window,0,tit='GS results',xsize=800,ysize=600
    ps_start,file='gyro3d.eps',/encap,charsiz=1

    !p.multi=[0,2,2]
    plot,result.freq/1e9,result.flux,/xl,/yl,ps=-4,symsiz=0.25,xtit='Frequency [GHz]',ytit='Flux density [sfu]'

    IF keyword_set(data) THEN oplot,data.freq/1e9,data.flux,ps=4,col=cgcolor('red')

    loadct,3
    white=255
    
    sbox=[15,9,5]
    map1=map[[0,n_elements(freq)/2,n_elements(freq)-1]]
    ;; map1=replicate(loopmap,3)
    ;; map1[0].data=smooth(map[*,*,0],sbox[0])
    ;; map1[1].data=smooth(map[*,*,n_elements(freq)/2],sbox[1])
    ;; map1[2].data=smooth(map[*,*,n_elements(freq)-1],sbox[2])
    plot_map,map1[0],tit=string(freq[0]/1e9,format='(i0)')+'GHz',grid=5,/smooth,gcol=white
    plot_map,map1[0],/smooth,/over,/perc,lev=[30,50,70]
    plot_map,map1[1],tit=string(freq[n_elements(freq)/2]/1e9,format='(i0)')+'GHz',grid=5,/smooth,gcol=white
    plot_map,map1[1],/smooth,/over,/perc,lev=[30,50,70]
    plot_map,map1[2],tit=string(freq[n_elements(freq)-1]/1e9,format='(i0)')+'GHz',grid=5,/smooth,gcol=white
    plot_map,map1[2],/smooth,/over,/perc,lev=[30,50,70]
    ps_end
    !p.multi=0    
    wdelete,0
 ENDIF

 IF keyword_set(show) THEN BEGIN 

    window,2,tit='loop/plasma parameters',xsize=800,ysize=600
    ps_start,file='loop-plasma-params.eps',charsiz=1.5
    !p.multi=[0,3,2]
    xtit='loop sections'
    x=indgen(n_sections)
    plot,x,temp,ytit='T [K]',thic=3,xtit=xtit
    plot,x,np1,ytit='Np [cm^-3]',thic=3,xtit=xtit
    plot,x,nel1*geo.sector_vol,ytit='Ntot',thic=3,xtit=xtit
    plot,x,nel1,ytit='Nel [cm^-3]',thic=3,xtit=xtit
    plot,x,bmag1,ytit='Bmag [G]',thic=3,xtit=xtit
    plot,x,geo.angle,ytit='viewing angle [degrees]',thic=3,xtit=xtit
    !p.multi=0    
    ps_end
 ENDIF 

 

END 
