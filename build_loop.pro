pro build_loop,p,v,NCube=NCube,r0=r0,dt=dt $
                ,lat=lat,lon=lon,incl=incl,azim=azim,length=length $
                ,geo=geo,n_Sections=n_Sections                 $
                ,map_struct=map_struct,box=box,quiet=quiet $
                ,date=date,radiuspx=radiuspx,bmag=bmag,d_arr=d_arr $
                ,px2cm=px2cm,ratiob=ratiob,show=show ;;,_EXTRA=_EXTRA

;+
; build_loop
;
; USE:
;
;  build_loop, p, v
;
;primer: build_loop,a,b,azim=80,lon=80,ncube=256,r0=10,bmag=[9.,1.,9.]
;tvcon,total(b,3)
; Creates a filled pixelized semi-circular loop/tube model.
; The output is used with resuts from gyrosync,
; to perform 3D radiative transfer.
;
; OUTPUT:
;
; p: data structure
;    x,y,z: voxel coordinates (voxels, sections)
;    px2cm: resolution cm/px
;
; v: datacube (ncube x ncube x ncube) or map_struct (for SSW plot_map) if MAP_STRUCT=1
;    NOTE: map_struct is just for easy visualization of the loop
;    model; solar coordinates are not fully tested.
;
; INPUT KEYWORDS:
;
; LAT: latitudinal rotation [degrees]
; LON: longitudinal rotation [degrees]
; INCL: inclination angle (local) [degrees]
; AZIM: azimuthal angle (to the "equator") [degrees]
; LENGTH: full length of the tube
; N_SECTIONS: number of divisions of the tube
; R0: cross-section of the tube [voxels] (def: r0=R/10)
; nCube: size of the datacube (nCube x nCube x nCube) Def=128
; DT: with dt=0, the loop will be a semi-circunference (pi);
;     setting dt to a degree value the loop increases: pi+2*dt
;
; OPTIONAL KEYWORDS:
; BOX: smooth boxcar for map_struct image
; quiet: don't print info on screen
;
; TESTING:
; radiuspx: loop radius in voxels (def: R=Ncube/4)
; px2cm: resolution cm/px. Override internal px2cm, use it with
; radiuspx and length values.
;
; OUTPUT KEYWORDS:
; geo: structure with several info and data (loop and voxel volume, length,
; height, viewing angles, position angles, projected area, angular
; diameter (arcsec) for circular source of the same area, etc)
;
; PJSimoes
;-

show=keyword_set(show) 

 AU = 1.49597870d13             ; Astronomic Unit
 arc2cm = !dtor/3600d0 * AU     ; arcsec to cm in Sun

; cube size (pixels)
 if n_elements(nCube) eq 0 then n = 128 else n=nCube
 nCube=n

; main radius
 IF ~keyword_set(radiuspx) THEN R = N/4 ELSE R = radiuspx

; cross-section radius
 if n_elements(r0) eq 0 then r0 = round(R/10.0)>1.0

 if ~keyword_set(bmag) then bmag=[1.0,1.0]

; center
 c=[n,n,n]/2

; volume
 v = fltarr(n,n,n)
 l = fltarr(n,n,n)

 M = n*10
 if n_elements(dt) eq 0 then dt = 0.0
 dt1=dt
 ;; IF n_elements(dt2) EQ 0 THEN dt2=dt

 ;; from here, dt is in RAD
 dt = dt*!dtor
 ;; dt2 = dt2*!dtor

 ;; ratio of local cross-section radius and r0 (pixels)
 ratioB=(min(bmag)/bmag)^0.5
 ratioB=(interpol(ratioB,M))
 ratio= (ratiob-min(ratiob))*!pi/2.
 ;;ratio= (ratiob)*!pi/2.

                                ;di = (r0*ratioB)
 foot1=(min(r0*ratiob[0:m/2-1]))
 foot2=(min(r0*ratiob[m/2:*]))
 di1=interpol([foot1,r0],m/2)
 di2=interpol([foot2,r0],m/2)
 ;;di1=ratiob[0:m/2-1]*r0
 ;;di2=reverse(ratiob[m/2:*])*r0
 

 ;;stop

 print,'cross-section radius:'
 print,'r0 main: ',r0
 print,'foot 1: ',foot1
 print,'foot 2: ',foot2
                                ;t = interpol([0.0-dt,!pi+dt],M)

; d_arr=0

; vou de xc-r0 a xc+r0, e encontro onde existe arco
 for xi=c[0]-r0,c[0]+r0 do begin ; subst por vetor de xc-r0,xc+r0 com step=1

    img = fltarr(n,n)
    len = fltarr(n,n)

    for ri=R-r0,R+r0 do begin

; coordenadas do pont no apex do loop, em relacao ao ponto central do apex:
;apex = [c[0],0.0,c[2]+r]
       dx = xi-c[0]          ;; posicao X (ou seja, xi) - Xc (ou c[0])
       dy = 0.0              ;; posicao Y (ou ri*cos(90)+c[1] - Yc (ou c[1])
       dz = ri-R             ;; posicao Z (ou ri*sin(90)+c[2]) - (Zc (ou c[2]+R) + Raio R)
; distancia entre o ponto e o centro do apex
       d = (sqrt((dx)^2.+(dy)^2.+(dz)^2.))

                                ;  d_arr = [d_arr,d]
;d=round(d)
       IF d LE r0 THEN BEGIN ;; draw loop if inside r0

          ;;  d=d<r0
          if d le foot1 then wid1=0.0 else $
           wid1=interpol(ratio[0:M/2-1],di1,d)
          if d le foot2 then wid2=0.0 else $
           wid2=interpol(ratio[M/2:*],reverse(di2),d)
          ;; if d gt max([foot1,foot2]) and wid2/!dtor gt 0 and wid1/!dtor gt 0 $
          ;; then print,foot1,d,foot2,'   ',wid1/!dtor,wid2/!dtor else continue
;if d ge 9 then continue

          ;;IF d ge r0 THEN
          ;; print,d,wid1,wid2

          t = interpol([0.0-dt+wid1,!pi+dt-wid2],M)
          y=ri*cos(t)+c[1]
          z=ri*sin(t)+c[2]
          y=round(y)
          z=round(z)
          value = 1.0

          img[y,z] = value
          v[*,xi,*]=img         ; nao estao erradas as posicoes [*,xi,*]; para a rotacao
          ;; da geometria funcionar, esta deve ser a orientacao do arco.
          ;;v[c[0],c[1],c[2]]=max(v)

                                ; tv,bytscl(rebin(total(v,1),256,256)),0
                                ; tv,bytscl(congrid(total(v,2),256,256)),1
                                ; tv,bytscl(rebin(total(v,3),256,256)),2
                                ; wait,.005

          len[y,z]= atan(z-c[2],y-c[0])
          l[*,xi,*] = len
          ;; apex = [mean(y),max(z)]
          ;; linearDistance = sqrt((y-apex[0])^2 + (z-apex[1])^2)
          ;; alpha = 2.0*asin(linearDistance/2.0/ri) * (y-apex[0])/abs(y-apex[0])
          ;; NaN = where(finite(alpha) eq 0)
          ;; if NaN[0] NE -1 then alpha[NaN]=0.0
          ;; comp = alpha          ;;* ri ; lenght=angle*radius (L=2!PI*R)
          ;; len[y,z] = comp
          ;; l[*,xi,*] = len
       endif
    endfor
 endfor

 p = where(v NE 0.0)
 if p[0] eq -1 then begin
    print,'NO VALID POINTS.'
    RETURN
 endif
 p = array_indices(v,p)
 ;; loop built. p[*,3] array is x,y,z coordinates, before rotation

; this is important to identify loop sections:
; sorting by polar angle:
 ll=L[p[0,*],p[1,*],p[2,*]]
 sections=sort(ll)
 pll = where(ll LE -1.0*!pi/2.,pllcount)
 IF pllcount NE 0 THEN ll[pll] = (ll[pll])+2.0*!pi
 ;; ll = (ll+!pi/2.0) ;; polar angle from (0-dt) to (!pi+dt)

 x=reform(p[0,*])
 y=reform(p[1,*])
 z=reform(p[2,*])

; no. de secoes do arco (definido no fkrplk)
 if n_elements(n_Sections) eq 0 then nSections = 15 else nSections=n_Sections

 ;; FINDING THE VOXELS IN EACH SECTION: THE SECTIONS ARE DEFINED BY
 ;; CONSTANT LENGHT (OR CONSTANT POLAR ANGLE STEP)
 ;; values in voxels
 
; h=histogram(reform(ll),loc=loc,rev=rr,nbins=nSections+1,min=-dt*0.99,max=(!pi+dt)*1.01)
 ;;h=histogram(reform(ll),loc=loc,rev=rr,nbins=nSections+1,min=-dt,max=pi_aprox+dt)

;; binsize=(max(ll)-min(ll))/(nsections)
;; h=histogram(reform(ll),loc=loc,rev=rr,binsize=binsize)

 nloc=nsections
h=hist(reform(ll),nsections,rr,binsize,loc)
;;col=cgcolor(['blue','magenta'])
;;plot,[1,1],[0,0],xr=[0,n-1],yr=[0,n-1],/xs,/ys,/nodata,/iso                          
;;for i=0,nsections-1 do oplot,x[pos[ind[i]:ind[i+1]-1]],z[pos[ind[i]:ind[i+1]-1]],ps=3,col=col[i mod 2]

 IF show THEN BEGIN
    col=cgcolor(['blue','magenta'])
    plot,x,z,ps=3,xr=[0,n-1],yr=[0,n-1],/xs,/ys,/iso
    for i=0,nloc-1 do BEGIN
      ; IF rr[rr[i]] NE rr[rr[i+1]-1] THEN $
        oplot,x[rr[rr[i]:rr[i+1]-1]],z[rr[rr[i]:rr[i+1]-1]],col=col[i MOD 2],ps=3
    ENDFOR
    for i=0,nloc-1 do plots,[c[0],(r+r0)*cos(loc[i]+binsize/2)+c[0]],[c[2],(r+r0)*sin(loc[i]+binsize/2)+c[2]],lines=1
    for i=0,nloc-1 do xyouts,mean(x[rr[rr[i]:rr[i+1]-1]]),mean(z[rr[rr[i]:rr[i+1]-1]]),string(i,format='(i2)'),align=0.5
 ENDIF

;;stop
 ;; FORGET THIS
 ;; centering the loop base at the center of the cube
 ;; z = z-min(z)+c[2]

 ;; ROTATION
 if n_elements(lat) eq 0 then lat=0.
 if n_elements(lon) eq 0 then lon=0.
 if n_elements(incl) eq 0 then incl=0.
 if n_elements(azim) eq 0 then azim=0.

; GIRA AS COORDENADAS

 p = rot_voxel(x,y,z,lat,lon,incl=incl,azim=azim,size=size,b0=b0,plot=plot)

 xp=reform(p[*,0])
 yp=reform(p[*,1])
 zp=reform(p[*,2])

 IF (0) THEN BEGIN
    col=cgcolor(['blue','magenta'])
    plot,xp,yp,ps=3,xr=[0,n-1],yr=[0,n-1],/xs,/ys,/iso
    for i=0,nSections-1 do BEGIN
  ;;     IF rr[rr[i]] NE rr[rr[i+1]-1] THEN $
        oplot,xp[rr[rr[i]:rr[i+1]-1]],yp[rr[rr[i]:rr[i+1]-1]],col=col[i MOD 2],ps=3
    ENDFOR
;for i=0,nsections do plots,[c[0],(r+r0)*cos(loc[i])+c[0]],[c[2],(r+r0)*sin(loc[i])+c[2]],lines=1
    for i=0,nsections-1 do xyouts,mean(xp[rr[rr[i]:rr[i+1]-1]]),mean(yp[rr[rr[i]:rr[i+1]-1]]),string(i,format='(i2)'),align=0.5
 ENDIF

 ;; /* OLD VERSION
 ;;setting outputs:
 ;; x = p[sections,0]
 ;; y = p[sections,1]
 ;; z = p[sections,2]
 ;; nVoxels = n_elements(x)        not right. pnts[x,y,z] overlap....
 ;; i = indgen(nVoxels/nSections)
 ;; xp = fltarr(n_elements(i),nSections)
 ;; yp = fltarr(n_elements(i),nSections)
 ;; zp = fltarr(n_elements(i),nSections)
 ;; for j=0,nSections-1 do begin
 ;;    pos = i+nVoxels/nSections*j
 ;;    xp[*,(nSections-1)-j]=x[pos]
 ;;    yp[*,(nSections-1)-j]=y[pos]
 ;;    zp[*,(nSections-1)-j]=z[pos]
 ;; endfor
 ;; OLD VERSION */

; ajuste para volume discretizado:
 ;; v = fltarr(n,n,n)
 ;; v[xp,yp,zp] = 1.0
 ;; nVoxels = total(v)

 ;;[rr[rr[i]:rr[1+i]-1]]
 nvoxels=lonarr(nsections)
 FOR i=0,nloc-1 DO BEGIN
    v = fltarr(n,n,n)
    v[xp[rr[rr[i]:rr[1+i]-1]],yp[rr[rr[i]:rr[1+i]-1]],zp[rr[rr[i]:rr[1+i]-1]]] = 1.0
    nVoxels[i] = total(v)
 ENDFOR

 ;; SCALING TO REAL WORLD
 ;; if n_elements(dt) eq 0 then dt = 0.
 ;; LENGTH is FULL LOOP distance.
 if n_elements(LENGTH) eq 0 then LENGTH = 4.0E9 ; standard length.
 RADIUS = LENGTH / (!pi+(dt+dt))
 IF ~keyword_set(px2cm) THEN px2cm = RADIUS/R

 dV = px2cm^3.
 Volume = total(nVoxels) * dV
 Sect_Vol = nVoxels * dV
 ;; number of voxels in each section (not in voxels yet)
 ;;nvox = rr[1:nsections] - rr[0:nsections-1]
;;;;; nvol = nvox * dV

;stop

 ;; old
 ;; Sect_Vol = Volume / nSections

 ;; arco: (input) dt em graus
 arco,LAT,LON,AZIM,INCL,LENGTH,cs,coord,npts_output=nSections,dt=dt/!dtor
 angle = acos(cs)/!dtor

; os valores de volume, raio e raio cross-section

 v = fltarr(n,n,n)
 v[xp,yp,zp] = 1.0

 area=total(v,3,/nan)
 pp=where(area GE 1,nV)
 area=nV*px2cm^2
 omega=2.0*sqrt(area/!pi)/arc2cm
 px2arc = px2cm/arc2cm

 IF ~keyword_set(quiet) THEN BEGIN
    print,'------ LOOP GEOMETRIC PARAMETERS ------'
    print,'px2cm',px2cm,' cm/px'
    print,'px2arc',px2arc,' arc/px'
    print,''
    print,'volume',volume,' cm^3'
    print,'radius',px2cm*r,' cm'
    print,'cross-section',px2cm*r0,' cm'
    print,'full length',px2cm*r*(!pi+dt+dt),' cm'
    print,'projected area',area,' cm^2'
    print,'equivalent angular diameter',omega,' arcsec'
    print,'voxel volume',px2cm^3,' cm^3'
    print,'---------------------------------------'
 ENDIF
; voxel information
 p = {$
     x:xp,y:yp,z:zp,n:n,r:r,r0:r0,center:c,coord:coord $
     ,nVoxels:nVoxels,px2cm:px2cm,nSections:nSections,ri:rr}
 ;; ri is reverse index from histogram, which defines the sections.

 geo = { $
       vol:volume ,px2arc:px2arc, radius:px2cm*p.r,cross_section_radius:px2cm*p.r0 $
       ,height:(p.r+p.r0)*px2cm,length:px2cm*p.r*(!pi+dt+dt) $
       ,eq_angular_diameter:omega, angle:angle, dt:dt/!dtor $
       ,proj_area:area,Sector_Vol:Sect_Vol,lat:lat,lon:lon,azim:azim,incl:incl}

 if keyword_set(map_struct) then begin
    if not keyword_set(box) then box=1
    data = smooth(total(v,3),box)
    IF n_elements(date) EQ 0 THEN b0=0.0
    pos = hel2arcmin(LAT,LON,B0=b0,date=date) * 60.0
    ;; xc = pos[0] - (px2cm*p.r/ 72527093.)*sin(dt)*sin(lon*!dtor)
    ;; yc = pos[1] - (px2cm*p.r/ 72527093.)*sin(dt)*sin(lat*!dtor)
    xc = pos[0]                 ;- (px2cm*p.r/ 72527093.)*sin(dt)*sin(lon*!dtor)
    yc = pos[1]                 ;- (px2cm*p.r/ 72527093.)*sin(dt)*sin(lat*!dtor)
    dx = px2cm / 72527093.
    dy = dx
    get_utc,utc,/ccs            ;'2002-08-24T01:00:41.480'
    IF keyword_set(date) THEN time=date ELSE time = utc
    dur = 1.0
    units = 'arcsecs'
    ID = 'Loop Structure'
    roll_center = 0.0
    roll_angle = 0.0
    v = {data:data,xc:xc,yc:yc,dx:dx,dy:dy,time:time,dur:dur,units:units $
         ,ID:ID,ROLL_CENTER:roll_center,ROLL_ANGLE:roll_angle}

 endif

 dt=dt1

 ;; wid1=interpol(ratio[0:M/2-1],di1,d_arr)/2.
 ;; wid2=interpol(ratio[M/2:*],reverse(di2),d_arr)/2.
 ;; pp=where(d_arr le foot1)
 ;; wid1[pp]=0
 ;; pp=where(d_arr le foot2)
 ;; wid2[pp]=0

 ;; plot,d_arr,wid1,ps=2
 ;; oplot,d_arr,wid2,ps=4,col=180
 ;; plots,foot1,!y.crange
 ;; plots,foot2,!y.crange
 ;; plots,r0,!y.crange
;plot,scale_vector(g.sector_vol),/xs & oplot,scale_vector(1./interpol(b,15)),ps=-2
;stop
end

