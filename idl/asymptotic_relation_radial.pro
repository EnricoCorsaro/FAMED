function asymptotic_relation_radial,enn,numax,param
; param(0) = Dnu
; param(1) = epsilon
; param(2) = alpha
freq_asymp = param(0)*(enn + param(1) + param(2)/2.*(enn - numax/param(0))^2)
  
return,freq_asymp
end