This code is setup as library so that it can be easily integrated with other codes that provide rotor dynamics, acoustics, or even things like optimizers for rotor optimization. As such there is no end executable to run. An example scripts are provided to get feel for how to setup and run a simulation.

Normalizations
==============

Internally OpenCOPTER uses non-dimensional values. Currently it expects all inputs to be normalized appropriately and all state must be de-normalized if required. The normalizations used by OpenCOPTER are as follows:

Forces
------
Forces are normalized by the ambient density, rotor tip speed, and rotor area. For example the thrust coefficient for the rotor :math:`i`:

.. math::
    C_{T_i} = \frac{T_i}{\rho \pi R_i^2 (\Omega_i R_i)^2}

Similarly for the spanwise sectional thrust coefficient:

.. math::
    dC_{T_i} = \frac{dT_i}{\rho \pi R_i^2 (\Omega_i R_i)^2}

Velocities
----------
Velocities are normalized by the angular velocity and rotor radius of the rotor that produced them. For example the advance ratio for rotor :math:`i` with a radius :math:`R_i` and angular velocity :math:`\Omega_i`:

.. math::
    \mu = \frac{V_{\infty}}{\Omega_i R_i}

Distances
---------
Distances are non-dimensionalized by the rotor radius of the current rotor. For example, on rotor :math:`i` with radius :math:`R_i` the blade spanwise coordinate is normalized as:

.. math::
    r = \frac{y}{R_i}

Internal Data Structure
=======================

Internally OpenCOPTER represents all aircraft geometry and state in chunks of 8. Chunks of data that are used close to each other are then grouped together. These two things allow OpenCOPTER to make very effective use of a CPU's vector extensions (avx, avx2, avx512) without the need for inline assembly while maintaining good cache locality. However this means that interacting with the internal data directly can be a bit cumbersome. To alieviate this, a number of functions are provided that allow for getting and setting internal data using simple linear arrays. These functions can be seen in the :ref:`api documentation <API Documentation>`

State vs Input State
====================

As this code is just an aerodynamics code that is intended to be integrated with other codes to take care of things like vehicle/rotor dynamics we have separate structures for the models current aerodynamic state and the state that might be input from a different code or just prescribed. As such the :class:`libopencopter.AircraftState` / :class:`libopencopter.RotorState` / :class:`libopencopter.BladeState` structures are both outputs from the current iteration as well as inputs to the next and the :class:`libopencopter.AircraftInputState` and :class:`libopencopter.RotorInputState` structures are purely inputs to the model. Separate from these structures is the :class:`libopencopter.Wake` structure which can be seen as a state structure as well.

Units
=====

While many of the parameters are non-dimensional there are ones that are not. In these instances the specific units generally do not matter as long as they are applied consistently everywhere. The few units that are specific are angles which are always in radians and time which is always in seconds.
