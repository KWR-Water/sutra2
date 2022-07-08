===============
Getting Started
===============
sutra2 (SUbsurface TRAnsport) is a python package to model the behavior of Organic
MicroPollutants (OMPs) and pathogens for 4 standard types of Public Supply Well
Fields (PSWFs).

Main features:
 - Calculate the travel time distribition for 4 PSWF types using an analytical or numerical (Modflow) solution
 - Calculate the travel time for different aquifer zones
 - Numerical solution using FloPy
 - Includes database of substance paramters
 - Calculate the removal and concentration of the substance or pathogen is interest in each aquifer zone and the well

..
    @Steven/MartinvdS anythign to add here? for modflow? (#AH @MartinK, somehting about QSAR here)

Install
-------
sutra2 requires Python 3.6 or later.

To get sutra2, use the following command::

  pip install git+https://github.com/KWR-Water/sutra2.git

..
  #AH @MartinK -> check how to do this. @ALEX: do we want to have this on pypi?

Philosophy
----------

..
  #AH AH @MartinvdS @MartinK ...  what here?

sutra2 calculates the behavior of OMPs and pathogens for 4 standard types of PSWFs, which cover the most frequently occurring and most vulnerable
groundwater resources for drinking water supply in the Netherlands (and Flanders).
One of the aims of this approach is to forecast the behavior of new OMPs in
groundwater. Groundwater is often overlooked in current environmental risk
assessment methods, which are a priori or a posteriori applied when new organic
chemicals appear on the market.


Conventions
-----------
..
  #AH AH @MartinvdS @MartinK ...  what here?
  #@ALEX: MWK: dont know. ask MartinvdS to see if he wants to indicate we are using flopy convention here.

The naming conventions follow the naming conventions for... FloPy?