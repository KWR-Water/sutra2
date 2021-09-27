============
Known Issues
============

1. Phreatic schematisation with analytical solution: 
If the vadose zone and shallow aquifer have small thicknesses then it could happen
that the drawdown results in a water level at the well which is lower than the top 
of the target aquifer. In that case the travel time in the shallow aquifer close to the well is zero.
This can result in a travel time distribution which initially decreases before increasing again 
with radial distance from the well. The Modflow solution uses a fixed head, therefore
if this issue needs to be avoided use the Modflow solution instead. 