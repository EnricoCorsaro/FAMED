function get_modeid,freq,dnu,epsi_input,d01,numax,alpha
; -------------------------------------------------------------------------------------------------------
; Author:  Enrico Corsaro
; e-mail:  enrico.corsaro@inaf.it
; Date:    August 2019
; Place:   Catania, Italy
; Purpose: This function obtains the order number and angular degree of an input frequency value.
;          It assumes that the input epsilon parameter is reliable and incorporates the curvature
;          term, which plays no effect if set to zero.
; -------------------------------------------------------------------------------------------------------
radial_order_array = findgen(60)

epsi_radial_array = (freq - dnu*radial_order_array - dnu*alpha/2.*(radial_order_array - numax/dnu)^2)/dnu
diff = epsi_radial_array - epsi_input

tmp = where(diff ge -0.25)
diff = diff(tmp)
min_diff = min(diff,j)

good_order = radial_order_array(j)
epsi_local_radial = (freq - dnu*good_order - dnu*alpha/2.*(good_order - numax/dnu)^2)/dnu
epsi_local_dipole = (freq - dnu*(good_order+0.5) - dnu*alpha/2.*(good_order - numax/dnu)^2 + d01)/dnu
diff_radial = abs(epsi_local_radial - epsi_input)
diff_dipole = abs(epsi_local_dipole - epsi_input)
order_number = round(good_order)

if abs(diff_radial) lt abs(diff_dipole) then begin
    ; In this case it is a radial mode

    angular_degree = 0
endif else begin
    ; In this case it is a dipole mode

    angular_degree = 1
endelse

mode_id = {    order:      order_number,   $
               degree:     angular_degree  $
          }

return,mode_id
end
