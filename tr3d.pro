pro tr3d,p,freq,jo,jx,ko,kx,nel,N,geo=geo $
         ,temperature=temperature,np=np $
         ,map=map,flux=flux,mo=mo,mx=mx,pd=pd,f2tb=f2tb $
         ,opa_o=opa_o,opa_x=opa_x,jff_i=jff_i,kff_i=kff_i ;,joV=joV,jxV=jxV,koV=koV,kxV=kxV

; nel: non-thermal electron number density (real, comes from solflare.pro)

 AU = 1.49597870d13
 arc2cm = !dtor/3600d0 * AU
 cm2arc = 1.0/arc2cm
 kb=1.38e-16
 c=3e10

 nf = n_elements(freq)

 map = fltarr(N,N,nf)
 mo = map
 mx = map
 pd = map
 opa_o = map
 opa_x = map

 ;if n_elements(temperature) ne 0 then mt = map else mt = 0.0
 flux = fltarr(nf)
 joV = fltarr(n,n,n)
 jxV = fltarr(n,n,n)
 koV = fltarr(n,n,n)
 kxV = fltarr(n,n,n)

 L = p.px2cm
 omega = L^2.0 / au^2
 cgs2sfu = 1e19
 flux2tb = 1./cgs2sfu * c^2 /kb/(freq)^2/(L^2/au^2)
 f2tb=flux2tb

 ;; free-free emission and absorption (only used if temperature is set)
 jff = 0.0d
 kff = 0.0d
 jff_i=dblarr(nf,p.nsections)
 kff_i=dblarr(nf,p.nsections)

 IF keyword_set(temperature) THEN print,'FREE-FREE'
 for f=0,nf-1 do BEGIN
    for i=0,p.nsections-1 do BEGIN
       
;; REMOVING FREE-FREE FOR NOW TO DEBUG (2013-FEB)
       IF keyword_set(temperature) THEN freefree,temperature[i],np[i],freq[f],jff,kff

; summation of free-free and GS:
       ;; OLD VERSION:
       ;; joV[p.x[*,i],p.y[*,i],p.z[*,i]]=jo[f,i]*nel[i]+jff
       ;; 		jff_i(f,i)=jff
       ;; koV[p.x[*,i],p.y[*,i],p.z[*,i]]=ko[f,i]*nel[i]+kff
       ;; 		kff_i(f,i)=kff
       ;; jxV[p.x[*,i],p.y[*,i],p.z[*,i]]=jx[f,i]*nel[i]+jff
       ;; kxV[p.x[*,i],p.y[*,i],p.z[*,i]]=kx[f,i]*nel[i]+kff

       joV[p.x[p.ri[p.ri[i]:p.ri[i+1]-1]],p.y[p.ri[p.ri[i]:p.ri[i+1]-1]],p.z[p.ri[p.ri[i]:p.ri[i+1]-1]]] $
        =jo[f,i]*nel[i]+jff
       koV[p.x[p.ri[p.ri[i]:p.ri[i+1]-1]],p.y[p.ri[p.ri[i]:p.ri[i+1]-1]],p.z[p.ri[p.ri[i]:p.ri[i+1]-1]]] $
        =ko[f,i]*nel[i]+kff
       jxV[p.x[p.ri[p.ri[i]:p.ri[i+1]-1]],p.y[p.ri[p.ri[i]:p.ri[i+1]-1]],p.z[p.ri[p.ri[i]:p.ri[i+1]-1]]] $
        =jx[f,i]*nel[i]+jff
       kxV[p.x[p.ri[p.ri[i]:p.ri[i+1]-1]],p.y[p.ri[p.ri[i]:p.ri[i+1]-1]],p.z[p.ri[p.ri[i]:p.ri[i+1]-1]]] $
        =kx[f,i]*nel[i]+kff

       jff_i[f,i]=jff
       kff_i[f,i]=kff

    endfor


v=intarr(p.n,p.n,p.n)
for i=0,p.nsections-1 do $
       v[p.x[p.ri[p.ri[i]:p.ri[i+1]-1]],p.y[p.ri[p.ri[i]:p.ri[i+1]-1]],p.z[p.ri[p.ri[i]:p.ri[i+1]-1]]] $
        =cos(geo.angle[i]*!dtor)/abs(cos(geo.angle[i]*!dtor))

sgn=float(total(v,3))
ltz=where(sgn lt 0)
gtz=where(sgn gt 0)
sgn[ltz]=-1
sgn[gtz]=+1


;VM: jo,jx and ko, kx were calculated in the programs gscore.pro and gyrosync.pro for the specific normalization
; n1=n(i)/nt(i) number density=1.
; So, here we multiply them by nel that is a real electron number density (obtained using the normalization expression)

    ;;  if keyword_set(temperature) then begin
    ;;        ;; j_ff = eta_ff(freq[f],np,temperature)
    ;;        ;; k_ff = kappa_ff(freq[f],np,temperature)
    ;;        j_ff = jff(freq[f],np,temperature)
    ;;        k_ff = kff(freq[f],np,temperature)
    ;;
    ;;        joV[where(jo ne 0)>0] += j_ff[0]
    ;;        jxV[where(jx ne 0)>0] += j_ff[0]
    ;;        koV[where(ko ne 0)>0] += k_ff[0]
    ;;        kxV[where(kx ne 0)>0] += k_ff[0]
    ;;
    ;;        jvff = fltarr(n,n,n)
    ;;        kvff = fltarr(n,n,n)
    ;;        jvff[p.x,p.y,p.z]=j_ff[0]
    ;;        kvff[p.x,p.y,p.z]=k_ff[0]
    ;;
    ;;        imgT = transfrad(jvff,kvff,L)*omega*cgs2sfu*flux2tb[f]
    ;;        mt[*,*,f] = imgT
    ;;     endif

                                ;imgO = transfrad(joV,koV,L)*omega*cgs2sfu
                                ;imgX = transfrad(jxV,kxV,L)*omega*cgs2sfu

    imgO = transfrad(joV,koV,L)*omega*cgs2sfu
    imgX = transfrad(jxV,kxV,L)*omega*cgs2sfu
    img = imgO + imgX
    pol = sgn*(imgo - imgx)
    flux[f] = total(img,/nan)
    map[*,*,f] = img            ;* flux2tb[f]
    mo[*,*,f] = imgO            ;* flux2tb[f]
    mx[*,*,f] = imgX            ;* flux2tb[f]
    pd[*,*,f] = pol
;    opa_o[*,*,f] = opac_o
;    opa_x[*,*,f] = opac_x

 endfor

end
