.. _ipython_directive:

========================================================================================================
Tutorial microbial organism removal (mbo-removal)
========================================================================================================

Overview
========

Sutra includes the functionailty to calculate the advective removal of microbial organisms 
(also called 'pathogens') from source to end_point.

Main features:
 - Includes database of removal parameters for microbial organisms. 
 - Calculate the removal and concentration of the microbial organism over distance, and with time   

The code for microbial organism removal is based on the equations given in chapter 6.7 of 
BTO 2012.015 "Gezondheidsrisico's van fecale verontreiniging" [in Dutch], pp 71-74. The equations used are
described below.

Background information
======================

During transport in the subsurface, microbial organism removal takes place by both attachment to the soil matrix and by inactivation.
The virus concentration 'C' [m-3] through steady-state transport of microbial organisms along pathlines in the saturated
groundwater can be approximated by:

.. image:: _images/mrlp_20211018_equation1.PNG
  :width: 200
  :alt: Equation1: dC/dx + ((k_att + μ_1) / v) * C = 0

Where k_att + mu1 --> equal to removal rate 'lambda'
'k_att': attachment coefficient [day-1]
'mu1': inactivation coefficient [day-1] 
x: the distance traveled [m] 
v: the porewater velocity [m day-1] or 'darcyflux divided by the effective porosity'

Assuming that the background concentration of the relevant microbial organism equals 0,
the relative removal 'C_x/C0' can be calculated as follows.

.. image:: _images/mrlp_20211018_equation2.PNG
  :width: 250
  :alt: Equation2: log(C_x/C_0) = -(((k_att + μ_1))/ln⁡(10)) * (x/v)

The attachment coefficient 'k_att' depends on the effective porosity 'epsilon', the grain diameter of the sediment 'd_c',
'sticky coefficient' alpha [day-1], the porosity dependent Happel's parameter 'A_s', diffusion constant 'D_BM' [m2 day-1], and
the porewater velocity [m day-1].

.. image:: _images/mrlp_20211018_equation3.PNG
  :width: 300
  :alt: Equation3: k_att = 3/2 *(((1-ε))/d_c) * α*4*A_s^(1⁄3)*(D_BM/(d_c*ε*v))^(2⁄3) * v

The sticky coefficient alpha is determined by coefficient 'alpha0', which depends on both the soil type and the type of organism.
Alpha0 is being determined for a reference pH [pH0], e.g. pH=7.5.
Alpha relates to alpha0 as follows [corrected for different pH].

.. image:: _images/mrlp_20211018_equation4.PNG
  :width: 250
  :alt: Equation4: α = α_0 * 0.9^(((pH - pH_0) / 0.1))

The other parameters are calculated as follows:

Happel's porosity dependent parameter

.. image:: _images/mrlp_20211018_equation5.PNG
  :width: 250
  :alt: Equation5: A_s = 2 * ((1 - γ^5))/((2 - 3γ + 3γ^5 - 2γ^6))

Where:

.. image:: _images/mrlp_20211018_equation6.PNG
  :width: 250
  :alt: Equation6: γ = (1 - ε)^(1⁄3)

Boltzmann diffusion coefficient:

.. image:: _images/mrlp_20211018_equation7.PNG
  :width: 250
  :alt: Equation7: D_BM = (K_B * (T + 273))/(3π * d_p * μ) * 86400

with Boltzmann constant K_B [1,38 × 10-23 J K-1], organism diameter d_p [m], water temperature T [degr C], 
and conversion factor 86,400 [s day-1].

The dynamic viscosity 'mu' [kg m-1 s-1] depends on the groundwater density 'rho'.
The water density is assumed to be 999.7 [kg m-3], representative for fresh groundwater in the Netherlands under a reference
temperature of 12 degrees centigrade.

.. image:: _images/mrlp_20211018_equation8.PNG
  :width: 250
  :alt: Equation8: μ = (ρ * 497*10^(-6))/(T + 42.5)^(3⁄2) 