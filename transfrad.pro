FUNCTION transfrad,jv,kv,L,Bsign=Bsign,pol=pol
;+
; map = transfrad(jv,kv,L)
; jv: array de 3 dimensões com coeficientes de emissão
; kv:                                          absorção
; L: tamanho de cada elemento do volume (em cm)
;-

 t0=systime(/sec)
 s = size(jv)
 nlayer = s[3]
 map = fltarr(s[1],s[2])
 pol = fltarr(s[1],s[2])

 FOR i = 0, nlayer-1 DO BEGIN

    emission = jv[*,*,i]
    absorption = kv[*,*,i]
    opacity = absorption * L

    thn = where(opacity LT 0.001) ; if tau<<1 then e^-tau ~ 1-tau
    nor = where((opacity GE 0.001) and (opacity lt 10.))
    thk = where(opacity GE 10.) ; optc. thick

    map = map * exp(-opacity)

    IF thn[0] NE -1 THEN map[thn] = map[thn] + emission[thn] * L

    IF nor[0] NE -1 THEN map[nor] = map[nor] + emission[nor]/absorption[nor] $
                                    *(1.0-exp(-opacity[nor]))

    IF thk[0] NE -1 THEN map[thk] = emission[thk]/absorption[thk]

    IF keyword_set(bsign) THEN BEGIN 

       sign = bsign[*,*,i]
       pol = pol * exp(-opacity)
       IF thn[0] NE -1 THEN pol[thn] = pol[thn] + emission[thn] * L
       
       IF nor[0] NE -1 THEN pol[nor] = pol[nor] + emission[nor]/absorption[nor] $
                                       *(1.0-exp(-opacity[nor]))
       
       IF thk[0] NE -1 THEN pol[thk] = emission[thk]/absorption[thk]

       pp=where(sign NE 0.0)
       pol[pp] *= sign[pp]

    ENDIF 

 ENDFOR

 ;;print,'TRANSFRAD: elapsed time: ',systime(/sec)-t0
 RETURN, map

END

