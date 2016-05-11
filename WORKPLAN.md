The idea is to be able to handle mergers and sinking galactic satellites by 
having each subsystem centered on one (massless?) particle, and have the 
expansion calculated separately in each center of mass system (only for 
particles belonging to it). After the internal forces are calculated for each 
subsystem, the mutual forces are calculated as well, and all the particles are 
integrated in a global coordinate system.

Technically, most of the SCF module has to become a class. Some of the arrays 
like RadCoeff and AngCoeff (anything else?) should remain global and would be 
initialized by some global initialization routine. Then each subsystem would be 
a different instance of the SCF class, with its own list of coefficients and 
particle cache. It seems that all subsystems can (should?) share the integrator 
object. Maybe (probably?) the integrator should be agnostic to how particles are 
divided into subsystems, and the call to CalculateGravity() would take care of 
everything. In this case we are talking about a class for the current SCF 
routines but also some wrapper that knows how to CalculateGravity() by calling 
the appropriate methods for each subsystem's object, and also work out the 
mutual gravity. Need a good name for this intermediate level between the 
integrator and the SCF core.

This intermediate level or even a different class would deal how to associate 
particles to the different subsystems and how to figure out the initial 
conditions. Also, the central particles might need some special treatment 
(possibly the integrator considers them normal particles but this other level 
knows not to use them for force calculation; I should think whether those should 
be omitted from output, but it might be helpful to keep them).
