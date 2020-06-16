function asymptotic_relation_nonradial,enn,numax,ell,param
; param(0) = Dnu
; param(1) = epsilon
; param(2) = d0ell
; param(3) = alpha (Dnu curvature)
; param(4) = beta  (small spacing curvature)
freq_asymp = param(0)*(enn + param(1) + ell/2. + param(3)/2.*(enn - numax/param(0))^2 - param(4)*(enn - numax/param(0))) - param(2)
 
return,freq_asymp
end