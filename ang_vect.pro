FUNCTION ang_vect,a,b
;+
; obtencao do angulo entre dois vetores (produto escalar).
; angle_rad = ang_vect([a1,a2,a3],[b1,b2,b3])
; 
; INPUT: a,b : vetores.
; OUTPUT: theta : angulo entre os dois vetores em radianos.
;-
;theta = acos((total(double(a)*b)) / (sqrt(total(double(b)^2.)) * sqrt(total(double(a)^2.))))/!dtor
;

  RETURN, acos(total((a / sqrt(total(a^2.))) * (b / sqrt(total(b^2.)))))

END
