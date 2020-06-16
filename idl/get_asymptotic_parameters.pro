function get_asymptotic_parameters,numax,dnu,teff
; -------------------------------------------------------------------------------------------------------
; Author:  Enrico Corsaro
; e-mail:  enrico.corsaro@inaf.it
; Date:    August 2019
; Place:   Catania, Italy
; Purpose: This function compute some parameters from asymptotic values of nuMax and DeltaNu, including
;          an estimate for Teff. It provides estimates of stellar mass, radius, as well as the small
;          frequency spacings d02, d01, d03.
;          NOTE: The asymptotic parameters d02, d01, d03, are reliable only for RG stars. When using
;          this routine for SG and MS, make sure that proper corrections are applied.
; -------------------------------------------------------------------------------------------------------
COMMON CONFIG,cp

; Raw mass estimate from scaling relation

mass = (numax/cp.numax_sun)^3 * (dnu/cp.dnu_sun)^(-4) * (teff/cp.teff_sun)^1.5

; Raw radius estimate from scaling relation

radius = (numax/cp.numax_sun) * (dnu/cp.dnu_sun)^(-2) * (teff/cp.teff_sun)^0.5

; Scaling relation Dnu - d02 (Corsaro et al. 2012b)

slope = cp.d02_mass_offset + cp.d02_mass_slope*mass
d02 = slope*dnu + cp.d02_offset
               
; Scaling relation Dnu - d01 (Corsaro et al. 2012b)

slope = cp.d01_mass_offset + cp.d01_mass_slope*mass
d01 = slope*dnu + cp.d01_offset

; Scaling relation Dnu - d03 (Huber et al. 2010)     [1,20] microHz

d03 = cp.d03_slope*dnu + cp.d03_offset

asymptotic_parameters = { d02:     d02,    $
                          d01:     d01,    $
                          d03:     d03,    $
                          mass:    mass,   $
                          radius:  radius  $ 
                        }

return, asymptotic_parameters
end
