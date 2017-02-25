FUNCTION rot_voxel,x,y,z,lat,lon,incl=incl,azim=azim,size=size,b0=b0,plot=plot

;+
; rot_voxel,x,y,z,LAT,LON,INCL=INCL,AZIM=AZIM
;
; rotates voxel coordinates.
;
; Written by PJ Simoes. 2004-2008
;-

  x0 = (max(x)+min(x))/2.       ; centro de giro
  y0 = (max(y)+min(y))/2.
  z0 = min(z)

; Calculo da matriz de tranformacoes tmp
  T3D, /Reset, Translate = -[x0,y0,z0],matrix=tmp                  ; centro de giro no pe' do loop
  IF KEYWORD_SET(incl) THEN T3D, tmp,Rotate=[-incl,0.,0.],matrix=tmp ; inclinacao do loop
  IF KEYWORD_SET(b0) THEN T3D, tmp,Rotate=[-b0,0.,0.],matrix=tmp     ; inclinacao do Equador
  T3D, tmp,Rotate=[0.,0.,azim],matrix=tmp                            ; giro em azimute do arco
  T3D, tmp,Rotate=[-lat,0,0.],matrix=tmp                             ; giro ate' a posicao do arco
  T3D, tmp,Rotate=[0,lon,0.],matrix=tmp                              ; giro ate' a posicao do arco
  T3D, tmp, Translate=[x0,y0,z0],matrix=tmp ; origem das coord. em 0,0,0 de novo

; Coordenadas dos arcos numa matriz a ser transformada por tmp
  c_l = [[x],[y],[z],[replicate(1.,n_elements(x))]]
  c_l = c_l # tmp               ; Tranformacao (transladar, girar etc)

  coord = c_l[*,0:2]

  return, coord

end
