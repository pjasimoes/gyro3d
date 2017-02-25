pro arco,LAT,LON,AZIM,INCL,LENGTH,cs,coord,npts_output=npts_output,plot=plot,dt=dt

; input: valores em graus
; LAT - latitude heliografica do arco
; LON - longitude  heliografica do arco
; AZIM - giro em azimute
; INCL - inclinacao do arco

; output:
; cs - cos dos angulos de visada
; coord - coordenadas x,y,z dos pontos do arco

  if n_elements(dt) eq 0 then dt = 0.
  if not keyword_set(npts_output) then npts_output = 51

  arc2cm = 1d0 * !dtor/3600d0 * 1.49597870d13

; input de LENGTH e' FULL LENGTH
  Lc = LENGTH
; raio do arco:
  Rc = Lc / (!pi+dt*!dtor)

;cria arco hi-res p/ determinar ângulo de visada:
  pts = 110;;npts_output*10.0        ; no. de pontos no arco hi-res
  r = rc                        ;r=110
  t = interpol([0.-dt,180.+dt]*!dtor,pts)

; eixo central
  xp = r*cos(t)
  zp = r*sin(t)
  yp = xp*0.
;;  xp = xp - xp[0]

  B0 = 0                        ; inclinação do equador (compõe com 'incl')
  c = rot_voxel(xp,yp,zp,LAT,LON,azim=AZIM,b0=B0,incl=incl)

; calcula angulo de visada de cada trecho do arco.

  xpm = c[*,0]
  ypm = c[*,1]
  zpm = c[*,2]

  ii = indgen(pts-1)+1
;pos = p+i
  xx = xpm[ii] - xpm[ii-1]
  yy = ypm[ii] - ypm[ii-1]
  zz = zpm[ii] - zpm[ii-1]

  OBS = [0,0,-1.]
  angs = fltarr(pts)

  for i=0,pts-2 do angs[i] = ang_vect([xx[i],yy[i],zz[i]],OBS)
  angs[109] = angs[108]

  cs = fltarr(11)
  ii = indgen(10)
  jj=0
  for i=0,100,10 do begin
     cs[jj] = cos(mean(angs[indgen(10)+i]))
;print,mean(angs[indgen(10)+i])
     jj++
  endfor

;stop


  x = interpol(xpm,npts_output)
  y = interpol(ypm,npts_output)
  z = interpol(zpm,npts_output)
  cs = interpol(cs,npts_output)

  if keyword_set(PLOT) then begin

     xr=[-15,15]
     yr=[-15,15]

     window,1,xsiz=350,ysiz=700
     !p.multi=[0,1,2]
     plot,x/arc2cm+9,y/arc2cm,xtit='X [arcsec]',/ys,/iso $
          ,symsiz=2,thic=5,/xs,ytit='Y [arcsec]',ps=0,yr=yr,xr=xr,charsiz=1.2
     oplot,x/arc2cm+9,y/arc2cm,ps=4,col=200,symsiz=2
     plots,circle(x[0]/arc2cm+9,y[0]/arc2cm,1)

     plot,x/arc2cm+9,acos(cs)/!dtor,xtit='X [arcsec]'$
          ,symsiz=1.5,/xs,ytit='View Angle [degree]',ps=-4,xr=xr,charsiz=1.2,/yno
     !p.multi=0

     plots,circle(x[0]/arc2cm+9,acos(cs[0])/!dtor,1)

  endif

  coord = [[x],[y],[z]]

end
