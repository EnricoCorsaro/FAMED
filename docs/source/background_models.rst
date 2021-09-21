Background Models
=================
The current version of FAMED relies on a preliminary step that involves the fitting of a background model to the stellar power spectral density. In this section the models available are all introduced in terms of their analytical expressions and of the list and order of free parameters estimates that are required as input by FAMED. 

The names of the background models that can be used within FAMED are:

1. ``Flat``

2. ``Original``

3. ``OneHarvey``

4. ``OneHarveyColor``

5. ``TwoHarvey``

6. ``TwoHarveyColor``

7. ``ThreeHarvey``

8. ``ThreeHarveyColor``

Except for the components related to the instrumental noise, each component appearing in the background model is re-modulated by the apodization signal (or response function) :math:`R(\nu)` of the given observations. This is expressed as

:math:`R(\nu) = \mbox{sinc}^2 \left(\frac{\pi \nu}{2 \nu_\mathrm{Nyq}}\right)`

where :math:`\nu_\mathrm{Nyq}` is the Nyquist frequency of the dataset. 

All the background models include by default a Gaussian envelope to reproduce the region of the solar-like oscillations. This Gaussian envelope is represented by the Gaussian function

:math:`G(\nu) = H_\mathrm{osc} \exp \left[- \frac{(\nu - \nu_\mathrm{max})^2}{2 \sigma_\mathrm{env}^2} \right]`

with :math:`H_\mathrm{osc}` the height of the oscillation envelope (in :math:`\mbox{ppm}^2/\mu\mbox{Hz}`), :math:`\nu_\mathrm{max}` the frequency of maximum oscillation power (in :math:`\mu\mbox{Hz}`), and :math:`\sigma_\mathrm{env}` the standard deviation of the Gaussian (in :math:`\mu\mbox{Hz}`), proportional to the width of the oscillation region.

Flat
^^^^
This model can be identified with the label ``Flat``. This is the simplest of the background models available because it comprises only the white noise level in addition to the Gaussian envelope parameters. It is expressed by the equation:

:math:`P_\mathrm{bkg}(\nu) = W + R(\nu)G(\nu)`

There are 4 free parameters that can be fitted for this model. Their estimates must be supplied in the following order:

1. :math:`W`

2. :math:`H_\mathrm{osc}`

3. :math:`\nu_\mathrm{max}`

4. :math:`\sigma_\mathrm{env}`

Original
^^^^^^^^
This model can be identified with the label ``Original``. It comprises one Harvey profile in its original form, namely using an exponent set to 2. It is an outdated version of more recent background models. It can in some cases be adopted to test the effect of a change in the Harvey-profile slope. It is expressed by the equation:

:math:`P_\mathrm{bkg}(\nu) = W + R(\nu) \left[ B(\nu) + G(\nu) \right]`

where

:math:`B(\nu) = \frac{4 a^2/b}{1 + (2 \pi \nu/b)^2}`

This model comprises 6 free parameters, whose estimates must be supplied in the following order:

1. :math:`W`

2. :math:`a`

3. :math:`b`

4. :math:`H_\mathrm{osc}`

5. :math:`\nu_\mathrm{max}`

6. :math:`\sigma_\mathrm{env}`

One Harvey
^^^^^^^^^^
This model can be identified with the label ``OneHarvey``. It incorporates the white noise and a single Harvey-like component. The Harvey-like profile is a modified version of the classic Harvey profile (the one implemented in the model ``Original``) that has an exponent set to 4 (instead of 2). It has been shown by Kallinger et al. (2014) that a modified Harvey profile provides a better fit to the observed granulation-related signal, and therefore more evidence to the quality of the model in the light of the data. This model is suited for stars having a low oscillation frequency and/or stars with a high level of white noise (e.g. faint stars or with a low signal-to-noise ratio oscillations). It is expressed by the equation:

:math:`P_\mathrm{bkg}(\nu) = W + R(\nu) \left[ B(\nu) + G(\nu) \right]`

where 

:math:`B(\nu) = \frac{2 \sqrt{2}}{\pi} \frac{a^2/b}{1 + (\nu/b)^4}`

is the Harvey-like component that typically reproduces a granulation signal for stars with very low oscillation frequencies (e.g. below 10 :math:`\mu\mbox{Hz}`), defined by its rms amplitude :math:`a` (in ppm) and its characteristic frequency :math:`b` (in :math:`\mu\mbox{Hz}).
There are 6 free parameters that can be fitted for this model. Their estimates must be supplied in the following order:

1. :math:`W`

2. :math:`a`

3. :math:`b`

4. :math:`H_\mathrm{osc}`

5. :math:`\nu_\mathrm{max}`

6. :math:`\sigma_\mathrm{env}`

One Harvey Color
^^^^^^^^^^^^^^^^
This model can be identified with the label ``OneHarveyColor``. It is a more complex version of the ``OneHarvey`` model, including an additional Harvey-like profile to reproduce the colored noise at low frequency. This model is only suited for stars with very low oscillation frequencies (e.g. below 10 :math:`\mu\mbox{Hz}`) and a high frequency resolution and signal-to-noise ratio oscillations. This is mostly the case of stars observed for a long period (more than one year), where the signatures of a colored noise component can become significant. It is expressed by the equation:

:math:`P_\mathrm{bkg}(\nu) = N(\nu) + R(\nu) \left[ B(\nu) + G(\nu) \right]`

where again

:math:`B(\nu) = \frac{2 \sqrt{2}}{\pi} \frac{a^2/b}{1 + (\nu/b)^4}`

and where

:math:`N(\nu) = W + 2 \pi \frac{a_n^2/b_n}{1 + (\nu/b_n)^2}`

is a component that models the instrumental noise, now comprising both a white and a colored noise, the latter one represented by a Harvey profile with a rms amplitude :math:`a_n` (in ppm) and a characteristic frequency :math:`b_n` (in :math:`\mu\mbox{Hz}`). There are 8 free parameters that can be fitted for this model. Their estimates must be supplied in the following order:

1. :math:`W`

2. :math:`a_n`

3. :math:`b_n`

4. :math:`a`

5. :math:`b`

6. :math:`H_\mathrm{osc}`

7. :math:`\nu_\mathrm{max}`

8. :math:`\sigma_\mathrm{env}`

Two Harvey
^^^^^^^^^^
This model can be identified with the label ``TwoHarvey``. It comprises two Harvey-like components, typically used to model signal of medium quality for asteroseismology (e.g. TESS or K2 data). The two components are related to a low-frequency signal (e.g. activity, rotational modulation, super-granulation) and to a granulation activity (mostly referring to the meso-granulation signal). It is expressed by the equation:

:math:`P_\mathrm{bkg}(\nu) = W + R(\nu) \left[ B(\nu) + G(\nu) \right]`

where

:math:`B(\nu) = \frac{2 \sqrt{2}}{\pi} \left[ \frac{a_1^2/b_1}{1 + (\nu/b_1)^4} + \frac{a_2^2/b_2}{1 + (\nu/b_2)^4} \right]`

is the two-component term of the Harvey-like profiles. There are 8 free parameters that can be fitted for this model. Their estimates must be supplied in the following order:

1. :math:`W`

2. :math:`a_1`

3. :math:`b_1`

4. :math:`a_2`

5. :math:`b_2`

6. :math:`H_\mathrm{osc}`

7. :math:`\nu_\mathrm{max}`

8. :math:`\sigma_\mathrm{env}`

Two Harvey Color
^^^^^^^^^^^^^^^^
This model can be identified with the label ``TwoHarveyColor``. It is a more complex version of the ``TwoHarvey`` model. Possible applications comprise stars with low oscillation frequencies (e.g. red clump stars, or below 30 :math:`\mu\mbox{Hz}`) that exhibit two clear background components. The colored-noise component will allow a more reliable estimation of the granulation component, provided that the quality of the data is sufficiently high to justify the adoption of this model. It is expressed by the equation:

:math:`P_\mathrm{bkg}(\nu) = N(\nu) + R(\nu) \left[ B(\nu) + G(\nu) \right]`

where

:math:`B(\nu) = \frac{2 \sqrt{2}}{\pi} \left[ \frac{a_1^2/b_1}{1 + (\nu/b_1)^4} + \frac{a_2^2/b_2}{1 + (\nu/b_2)^4} \right]`

It therefore comprises 10 free parameters, whose estimates must be supplied in the following order:

1. :math:`W`

2. :math:`a_n`

3. :math:`b_n`

4. :math:`a_1`

5. :math:`b_1`

6. :math:`a_2`

7. :math:`b_2`

8. :math:`H_\mathrm{osc}`

9. :math:`\nu_\mathrm{max}`

10. :math:`\sigma_\mathrm{env}`

Three Harvey
^^^^^^^^^^^^
This model can be identified with the label ``ThreeHarvey``. It is one of the most commonly used when referring to datasets from NASA Kepler (long observations, more than one year), and especially for main sequence, subgiant, and RGB stars. It incorporates three different Harvey-like profiles, where the low-frequency one refers to a signal of potential stellar activity, rotational modulation, and super-granulation, while the two additional profiles are aimed to model the meso-granulation and granulation signal (see Corsaro et al. 2017b for more details). It is expressed by the equation:

:math:`P_\mathrm{bkg}(\nu) = W + R(\nu) \left[ B(\nu) + G(\nu) \right]`

where in this case

:math:`B(\nu) = \frac{2 \sqrt{2}}{\pi} \left[ \frac{a_1^2/b_1}{1 + (\nu/b_1)^4} + \frac{a_2^2/b_2}{1 + (\nu/b_2)^4} + \frac{a_3^2/b_3}{1 + (\nu/b_3)^4} \right]`

This model comprises 10 free parameters. Their estimates must be supplied in the following order:

1. :math:`W`

2. :math:`a_1`

3. :math:`b_1`

4. :math:`a_2`

5. :math:`b_2`

6. :math:`a_3`

7. :math:`b_3`

8. :math:`H_\mathrm{osc}`

9. :math:`\nu_\mathrm{max}`

10. :math:`\sigma_\mathrm{env}`

Three Harvey Color
^^^^^^^^^^^^^^^^^^
This model can be identified with the label ``ThreeHarveyColor``. It is the most complete model among those implemented and a more complex version of the ``ThreeHarvey`` model. It is adopted only for stars that are observed for a long period of time (more than one year), which comprises stars observed by the nominal NASA Kepler mission and stars observed by NASA TESS for one year. It can be applied only in conditions of a good signal-to-noise ratio of the overall astrophysical signal. It is expressed by the equation:

:math:`P_\mathrm{bkg}(\nu) = N(\nu) + R(\nu) \left[ B(\nu) + G(\nu) \right]`

where again

:math:`B(\nu) = \frac{2 \sqrt{2}}{\pi} \left[ \frac{a_1^2/b_1}{1 + (\nu/b_1)^4} + \frac{a_2^2/b_2}{1 + (\nu/b_2)^4} + \frac{a_3^2/b_3}{1 + (\nu/b_3)^4} \right]`

It accounts for 12 free parameters, whose estimates must be supplied in the following order:

1. :math:`W`

2. :math:`a_n`

3. :math:`b_n`

4. :math:`a_1`

5. :math:`b_1`

6. :math:`a_2`

7. :math:`b_2`

8. :math:`a_3`

9. :math:`b_3`

10. :math:`H_\mathrm{osc}`

11. :math:`\nu_\mathrm{max}`

12. :math:`\sigma_\mathrm{env}`