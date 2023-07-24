Simulation
=============

Classes
-------

.. autoclass:: libopencopter.Atmosphere
   :members:
   :undoc-members:

Functions
---------

.. autofunction:: libopencopter.step(ac_state, aircraft, ac_input_state, inflows, wake_history, atmo, iteration, dt)

.. autofunction:: libopencopter.basic_aircraft_rotor_dynamics(input_state, dt)

.. autofunction:: libopencopter.basic_single_rotor_dynamics(input_state, dt)

.. autofunction:: libopencopter.chunk_size

Geometry
========

Classes
-------

.. autoclass:: libopencopter.Aircraft
   :members:
   :undoc-members:

.. autoclass:: libopencopter.RotorGeometry
   :members:
   :undoc-members:

.. autoclass:: libopencopter.BladeGeometry
   :members:
   :undoc-members:

.. autoclass:: libopencopter.BladeGeometryChunk
   :members:
   :undoc-members:

Functions
---------

.. autofunction:: libopencopter.set_r(bg, data)

.. autofunction:: libopencopter.set_twist(bg, data)

.. autofunction:: libopencopter.set_chord(bg, data)

.. autofunction:: libopencopter.set_C_l_alpha(bg, data)

.. autofunction:: libopencopter.set_alpha_0(bg, data)

.. autofunction:: libopencopter.set_sweep(bg, data)

State
=====

Classes
-------

.. autoclass:: libopencopter.AircraftState
   :members:
   :undoc-members:

.. autoclass:: libopencopter.RotorState
   :members:
   :undoc-members:

.. autoclass:: libopencopter.BladeState
   :members:
   :undoc-members:

.. autoclass:: libopencopter.BladeStateChunk
   :members:
   :undoc-members:

Functions
---------

.. autofunction:: libopencopter.get_dC_T(blade_state)

.. autofunction:: libopencopter.get_dC_Q(blade_state)

.. autofunction:: libopencopter.get_aoa(blade_state)

.. autofunction:: libopencopter.get_gamma(blade_state)

.. autofunction:: libopencopter.get_d_gamma(blade_state)

.. autofunction:: libopencopter.get_dC_L(blade_state)

.. autofunction:: libopencopter.get_dC_D(blade_state)

Input State
===========

Classes
-------

.. autoclass:: libopencopter.AircraftInputState
   :members:
   :undoc-members:

.. autoclass:: libopencopter.RotorInputState
   :members:
   :undoc-members:

Functions
---------

Inflow Models
=============

Classes
-------
.. autoclass:: libopencopter.HuangPeters
   :members:
   :undoc-members:

Wake
====

Classes
-------

.. autoclass:: libopencopter.WakeHistory
   :members:
   :undoc-members:

.. autoclass:: libopencopter.Wake
   :members:
   :undoc-members:

.. autoclass:: libopencopter.RotorWake
   :members:
   :undoc-members:

.. autoclass:: libopencopter.ShedVortex
   :members:
   :undoc-members:

.. autoclass:: libopencopter.VortexFilament
   :members:
   :undoc-members:

.. autoclass:: libopencopter.FilamentChunk
   :members:
   :undoc-members:

Functions
---------

.. autofunction:: libopencopter.get_wake_x_component(votex_filament)

.. autofunction:: libopencopter.get_wake_y_component(votex_filament)

.. autofunction:: libopencopter.get_wake_z_component(votex_filament)

.. autofunction:: libopencopter.get_wake_gamma_component(votex_filament)

.. autofunction:: libopencopter.get_wake_r_c_component(votex_filament)

.. autofunction:: libopencopter.get_wake_v_z_component(votex_filament)
   

Vtk integration
===============

Classes
-------

.. autoclass:: libopencopter.VtkWake
   :members:
   :undoc-members:

.. autoclass:: libopencopter.VtkRotor
   :members:
   :undoc-members:

Functions
---------

.. autofunction:: libopencopter.build_base_vtu_rotor(rotor)

.. autofunction:: libopencopter.write_rotor_vtu(base_filename, iteration, rotor_idx, vtk_rotor, rotor_state, rotor_input_state)

.. autofunction:: libopencopter.build_base_vtu_wake(wake)

.. autofunction:: libopencopter.write_wake_vtu(base_filename, iteration, vtk_wake, wake)

