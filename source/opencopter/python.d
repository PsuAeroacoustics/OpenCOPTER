module opencopter.python;

static import opencopter.aircraft;
import opencopter.airfoilmodels;
import opencopter.atmosphere;
import opencopter.math;
import opencopter.memory;

static import opencopter.vtk;

static import opencopter.bladeelement;
static import opencopter.wake;
static import opencopter.inflow;


import pyd.pyd;

import std.conv : to;
import std.exception : enforce;
import std.math : abs, fmod, PI;
import std.traits : isBasicType;

alias BI = opencopter.inflow.BeddosInflow!(ArrayContainer.array);
alias HP = opencopter.inflow.HuangPetersInflowT!(ArrayContainer.array);

size_t chunk_size() {
	static import opencopter.config;
	return opencopter.config.chunk_size;
}

opencopter.inflow.Direction Direction_clockwise() {
	return opencopter.inflow.Direction.clockwise;
}

opencopter.inflow.Direction Direction_counter_clockwise() {
	return opencopter.inflow.Direction.counter_clockwise;
}

//alias Direction = opencopter.inflow.Direction;
struct Direction {
	static opencopter.inflow.Direction clockwise() {
		return opencopter.inflow.Direction.clockwise;
	}

	static opencopter.inflow.Direction counter_clockwise() {
		return opencopter.inflow.Direction.counter_clockwise;
	}
}

void basic_aircraft_rotor_dynamics(PyAircraftInputState* ac_input, double dt) {
	foreach(r_idx, ref rotor; ac_input.rotor_inputs) {
		//double angle = rotor.angular_velocity*dt + rotor.angular_accel*dt*dt;

		//rotor.frame.rotate(Vec3(0, 0, 1), angle);
		rotor.azimuth += rotor.angular_velocity*dt + rotor.angular_accel*dt*dt;

		// Keep the azimuth between 0 and 2*PI so we don't
		// lose fp precicion as the sim marches in time and
		// the azimuth grows unbounded.
		if(rotor.azimuth > 2.0*PI) {
			rotor.azimuth = fmod(abs(rotor.azimuth), 2.0*PI);
		}
	}
}

double basic_single_rotor_dynamics(PyRotorInputState* input_state, double dt) {
	double angle = input_state.angular_velocity*dt + input_state.angular_accel*dt*dt;

	//rotor.frame.rotate(Vec3(0, 0, 1), angle);
	input_state.azimuth += angle;//input_state.angular_velocity*dt + input_state.angular_accel*dt*dt;

	// Keep the azimuth between 0 and 2*PI so we don't
	// lose fp precicion as the sim marches in time and
	// the azimuth grows unbounded.
	if(input_state.azimuth > 2.0*PI) {
		input_state.azimuth = fmod(abs(input_state.azimuth), 2.0*PI);
	}

	return angle;
}

/++
 +	This is here purely to act as in interface, but python does not have those. abstract class doesn't work
 +	either as link errors are thrown when the module is loaded complaining about the abstract class not
 +	having an implementation it can link to. So here we are.
 +/
class Inflow {
	void update(double C_T, PyRotorInputState* rotor, PyRotorState* rotor_state, double advance_ratio, double axial_advance_ratio, double dt) { assert(0); }
	void update(double C_T, PyRotorInputState rotor, PyRotorState rotor_state, double advance_ratio, double axial_advance_ratio, double dt) { assert(0); }
	Chunk inflow_at(immutable Chunk x, immutable Chunk y, immutable Chunk z, immutable Chunk x_e, double angle_of_attack) { assert(0); }
	Chunk inflow_at(immutable Chunk r, immutable double cos_azimuth, immutable double sin_azimuth) { assert(0); }
	double wake_skew() { assert(0); }
	Frame* frame() { assert(0); }
}


class HuangPeters : Inflow {
	private HP huang_peters;

	this(long Mo, long Me, PyRotorGeometry* rotor, double dt) {
		huang_peters = new HP(Mo, Me, rotor, dt);
	}

	this(PyRotorGeometry* rotor, double dt) {
		huang_peters = new HP(4, 2, rotor, dt);
	}

	override void update(double C_T, PyRotorInputState* rotor, PyRotorState* rotor_state, double advance_ratio, double axial_advance_ratio, double dt) {
		huang_peters.update(C_T, rotor, rotor_state, advance_ratio, axial_advance_ratio, dt);
	}

	override void update(double C_T, PyRotorInputState rotor, PyRotorState rotor_state, double advance_ratio, double axial_advance_ratio, double dt) {
		huang_peters.update(C_T, rotor, rotor_state, advance_ratio, axial_advance_ratio, dt);
	}

	override Chunk inflow_at(immutable Chunk x, immutable Chunk y, immutable Chunk z, immutable Chunk x_e, double angle_of_attack) {
		return huang_peters.inflow_at(x, y, z, x_e, angle_of_attack);
	}

	override Chunk inflow_at(immutable Chunk r, immutable double cos_azimuth, immutable double sin_azimuth) {
		return huang_peters.inflow_at(r, cos_azimuth, sin_azimuth);
	}

	override double wake_skew() {
		return huang_peters.wake_skew();
	}

	override Frame* frame() {
		return huang_peters.frame();
	}
}

class Beddoes : Inflow {
	private BI beddoes;

	this() {
		beddoes = new BI();
	}

	this(PyRotorGeometry* rotor, double dt) {
		beddoes = new BI();
	}

	override void update(double C_T, PyRotorInputState* rotor, PyRotorState* rotor_state, double advance_ratio, double axial_advance_ratio, double dt) {
		beddoes.update(C_T, rotor, rotor_state, advance_ratio, axial_advance_ratio, dt);
	}

	override void update(double C_T, PyRotorInputState rotor, PyRotorState rotor_state, double advance_ratio, double axial_advance_ratio, double dt) {
		beddoes.update(C_T, rotor, rotor_state, advance_ratio, axial_advance_ratio, dt);
	}

	override Chunk inflow_at(immutable Chunk x, immutable Chunk y, immutable Chunk z, immutable Chunk x_e, double angle_of_attack) {
		return beddoes.inflow_at(x, y, z, x_e, angle_of_attack);
	}

	override Chunk inflow_at(immutable Chunk r, immutable double cos_azimuth, immutable double sin_azimuth) {
		return beddoes.inflow_at(r, cos_azimuth, sin_azimuth);
	}

	override double wake_skew() {
		return beddoes.wake_skew();
	}

	override Frame* frame() {
		return null;
	}
}

alias PyWakeHistory = opencopter.wake.WakeHistoryT!(ArrayContainer.array);
alias PyWake = opencopter.wake.WakeT!(ArrayContainer.array);
alias PyRotorWake = opencopter.wake.RotorWakeT!(ArrayContainer.array);
alias PyVortexFilament = opencopter.wake.VortexFilamentT!(ArrayContainer.array);
alias PyShedVortex = opencopter.wake.ShedVortexT!(ArrayContainer.array);

alias PyAircraftTimehistory = opencopter.aircraft.AircraftTimehistoryT!(ArrayContainer.array);
alias PyAircraftState = opencopter.aircraft.AircraftStateT!(ArrayContainer.array);
alias PyRotorState = opencopter.aircraft.RotorStateT!(ArrayContainer.array);
alias PyBladeState = opencopter.aircraft.BladeStateT!(ArrayContainer.array);

alias PyAircraft = opencopter.aircraft.AircraftT!(ArrayContainer.array);
alias PyAircraftInputState = opencopter.aircraft.AircraftInputStateT!(ArrayContainer.array);
alias PyRotorGeometry = opencopter.aircraft.RotorGeometryT!(ArrayContainer.array);
alias PyBladeGeometry = opencopter.aircraft.BladeGeometryT!(ArrayContainer.array);
//alias PyFrame = opencopter.aircraft.Frame!(ArrayContainer.array);
alias Frame = opencopter.aircraft.Frame;

alias PyRotorInputState = opencopter.aircraft.RotorInputStateT!(ArrayContainer.array);

alias set_geometry_array = opencopter.aircraft.set_geometry_array;
alias get_geometry_array = opencopter.aircraft.get_geometry_array;
alias get_state_array = opencopter.aircraft.get_state_array;

void set_twist(ref PyBladeGeometry bg, double[] data) {
	bg.set_geometry_array!"twist"(data);
}

double[] get_twist(ref PyBladeGeometry bg) {
	return bg.get_geometry_array!"twist";
}

void set_chord(ref PyBladeGeometry bg, double[] data) {
	bg.set_geometry_array!"chord"(data);
}

double[] get_chord(ref PyBladeGeometry bg) {
	return bg.get_geometry_array!"chord";
}

void set_r(ref PyBladeGeometry bg, double[] data) {
	bg.set_geometry_array!"r"(data);
}

double[] get_r(ref PyBladeGeometry bg) {
	return bg.get_geometry_array!"r";
}

void set_C_l_alpha(ref PyBladeGeometry bg, double[] data) {
	bg.set_geometry_array!"C_l_alpha"(data);
}

void set_alpha_0(ref PyBladeGeometry bg, double[] data) {
	bg.set_geometry_array!"alpha_0"(data);
}

void set_sweep(ref PyBladeGeometry bg, double[] data) {
	bg.set_geometry_array!"sweep"(data);
}

void set_xi(ref PyBladeGeometry bg, double[] data) {
	bg.set_geometry_array!"xi"(data);
}

void set_xi_p(ref PyBladeGeometry bg, double[] data) {
	bg.set_geometry_array!"xi_p"(data);
}

double[] get_dC_T(ref PyBladeState blade) {
	return blade.get_state_array!"dC_T";
}

double[] get_dC_N(ref PyBladeState blade) {
	return blade.get_state_array!"dC_N";
}

double[] get_dC_c(ref PyBladeState blade) {
	return blade.get_state_array!"dC_c";
}

void fill_dC_N(ref PyBladeState blade, double[] data) {
	return blade.get_state_array!"dC_N"(data);
}

void fill_dC_c(ref PyBladeState blade, double[] data) {
	return blade.get_state_array!"dC_c"(data);
}

void fill_dC_D(ref PyBladeState blade, double[] data) {
	return blade.get_state_array!"dC_D"(data);
}

void fill_dC_Nf(ref PyBladeState blade, float[] data) {
	return blade.get_state_array!"dC_N"(data);
}

void fill_dC_cf(ref PyBladeState blade, float[] data) {
	return blade.get_state_array!"dC_c"(data);
}

void fill_dC_Df(ref PyBladeState blade, float[] data) {
	return blade.get_state_array!"dC_D"(data);
}

void fill_dC_Tf(ref PyBladeState blade, float[] data) {
	return blade.get_state_array!"dC_T"(data);
}

void fill_dC_Td(ref PyBladeState blade, double[] data) {
	return blade.get_state_array!"dC_T"(data);
}

void fill_dC_Qf(ref PyBladeState blade, float[] data) {
	return blade.get_state_array!"dC_Q"(data);
}

void fill_dC_Qd(ref PyBladeState blade, double[] data) {
	return blade.get_state_array!"dC_Q"(data);
}

double[] get_dC_Q(ref PyBladeState blade) {
	return blade.get_state_array!"dC_Q";
}

double[] get_dC_L(ref PyBladeState blade) {
	return blade.get_state_array!"dC_L";
}

double[] get_dC_D(ref PyBladeState blade) {
	return blade.get_state_array!"dC_D";
}

double[] get_aoa(ref PyBladeState blade) {
	return blade.get_state_array!"aoa";
}

double[] get_theta(ref PyBladeState blade) {
	return blade.get_state_array!"theta";
}

void fill_aoad(ref PyBladeState blade, double[] data) {
	return blade.get_state_array!"aoa"(data);
}

void fill_aoaf(ref PyBladeState blade, float[] data) {
	return blade.get_state_array!"aoa"(data);
}

void fill_u_td(ref PyBladeState blade, double[] data) {
	return blade.get_state_array!"u_t"(data);
}

void fill_u_tf(ref PyBladeState blade, float[] data) {
	return blade.get_state_array!"u_t"(data);
}

void fill_u_pd(ref PyBladeState blade, double[] data) {
	return blade.get_state_array!"u_p"(data);
}

void fill_u_pf(ref PyBladeState blade, float[] data) {
	return blade.get_state_array!"u_p"(data);
}

double[] get_gamma(ref PyBladeState blade) {
	return blade.get_state_array!"gamma";
}

double[] get_u_p(ref PyBladeState blade) {
	return blade.get_state_array!"u_p";
}

double[] get_u_t(ref PyBladeState blade) {
	return blade.get_state_array!"u_t";
}

double[] get_inflow_angle(ref PyBladeState blade) {
	return blade.get_state_array!"inflow_angle";
}

double[] get_aoa_eff(ref PyBladeState blade) {
	return blade.get_state_array!"aoa_eff";
}

double[] get_d_gamma(ref PyBladeState blade) {
	return blade.get_state_array!"d_gamma";
}

double[] get_wake_x_component(ref PyVortexFilament filament) {
	return opencopter.wake.get_wake_component!"x"(filament);
}

double[] get_wake_y_component(ref PyVortexFilament filament) {
	return opencopter.wake.get_wake_component!"y"(filament);
}

double[] get_wake_z_component(ref PyVortexFilament filament) {
	return opencopter.wake.get_wake_component!"z"(filament);
}

double[] get_wake_gamma_component(ref PyVortexFilament filament) {
	return opencopter.wake.get_wake_component!"gamma"(filament);
}

double[] get_wake_r_c_component(ref PyVortexFilament filament) {
	return opencopter.wake.get_wake_component!"r_c"(filament);
}

double[] get_wake_v_z_component(ref PyVortexFilament filament) {
	return opencopter.wake.get_wake_component!"v_z"(filament);
}


void fill_wake_x_component(ref PyVortexFilament filament, double[] data) {
	return opencopter.wake.get_wake_component!"x"(filament, data);
}

void fill_wake_y_component(ref PyVortexFilament filament, double[] data) {
	return opencopter.wake.get_wake_component!"y"(filament, data);
}

void fill_wake_z_component(ref PyVortexFilament filament, double[] data) {
	return opencopter.wake.get_wake_component!"z"(filament, data);
}

void fill_wake_gamma_component(ref PyVortexFilament filament, double[] data) {
	return opencopter.wake.get_wake_component!"gamma"(filament, data);
}

void fill_wake_r_c_component(ref PyVortexFilament filament, double[] data) {
	return opencopter.wake.get_wake_component!"r_c"(filament, data);
}

void fill_wake_v_z_component(ref PyVortexFilament filament, double[] data) {
	return opencopter.wake.get_wake_component!"v_z"(filament, data);
}

string FrameType_aircraft() {
	return opencopter.aircraft.FrameType.aircraft.to!string;
}

string FrameType_connection() {
	return opencopter.aircraft.FrameType.connection.to!string;
}

string FrameType_rotor() {
	return opencopter.aircraft.FrameType.rotor.to!string;
}

string FrameType_blade() {
	return opencopter.aircraft.FrameType.blade.to!string;
}

string FrameType_wing() {
	return opencopter.aircraft.FrameType.wing.to!string;
}


void step(PyAircraftState* ac_state, PyAircraft* aircraft, PyAircraftInputState* ac_input_state, Inflow[] inflows, PyWakeHistory* wake_history, Atmosphere* atmo, size_t iteration, double dt) {
	opencopter.bladeelement.step(*ac_state, *aircraft, *ac_input_state, inflows, *wake_history, *atmo, iteration, dt);
}

opencopter.vtk.VtkRotor build_base_vtu_rotor(PyRotorGeometry* rotor) {
	return opencopter.vtk.build_base_vtu_rotor(rotor);
}

void write_rotor_vtu(string base_filename, size_t iteration, size_t rotor_idx, opencopter.vtk.VtkRotor rotor, PyRotorState* rotor_state, PyRotorInputState* rotor_input, PyRotorGeometry* rotor_geom) {
	opencopter.vtk.write_rotor_vtu(base_filename, iteration, rotor_idx, rotor, rotor_state, rotor_input, rotor_geom);
}

void write_rotors_vtu(string base_filename, size_t iteration, opencopter.vtk.VtkRotor[] rotors, PyAircraftState* ac_state, PyAircraftInputState* ac_input, PyAircraft* aircraft) {
	foreach(r_idx; 0..aircraft.rotors.length()) {
		opencopter.vtk.write_rotor_vtu(base_filename, iteration, r_idx, rotors[r_idx], ac_state.rotor_states[r_idx], ac_input.rotor_inputs[r_idx], aircraft.rotors[r_idx]);
	}
}

opencopter.vtk.VtkWake build_base_vtu_wake(PyWake* wake) {
	return opencopter.vtk.build_base_vtu_wake(wake);
}

void write_wake_vtu(string base_filename, size_t iteration, opencopter.vtk.VtkWake vtk_wake, PyWake* wake) {
	opencopter.vtk.write_wake_vtu(base_filename, iteration, vtk_wake, wake);
}

void write_wake_field_vtu(string filename, PyAircraftState ac_state, PyWake wake, Vec3 delta, Vec3 starts, size_t num_x, size_t num_y, size_t num_z) {
	opencopter.vtk.write_wake_field_vtu(filename, ac_state, wake, delta, starts, num_x, num_y, num_z);
}

alias write_inflow_vtu = opencopter.vtk.write_inflow_vtu!(Inflow, Array!PyRotorGeometry);

void wrap_array(T)() {
	static if(isBasicType!T) {
		wrap_struct!(
			Array!T,
			OpIndex!(),
			OpIndexAssign!(T, size_t),
			Def!(Array!(T).next),
			Def!(Array!(T).__iter__),
			Def!(Array!(T).__len__),
			Def!(Array!(T).len),
			Def!(Array!(T).length),
			Repr!(Array!(T).toString)
		);

		ex_python_to_d((PyObject* o) {
			auto size = PyList_Size(o);
			auto arr = Array!(T)(size.to!size_t);

			foreach(idx; 0..size) {
				auto item = PyList_GetItem(o, idx.to!Py_ssize_t);
				arr[idx] = python_to_d!(T)(cast(PyObject*)item);
			}
			return arr;
		});
	} else {
		wrap_struct!(
			Array!T,
			OpIndex!(),
			OpIndexAssign!(T*, size_t),
			Def!(Array!(T).next),
			Def!(Array!(T).__iter__),
			Def!(Array!(T).__len__),
			Def!(Array!(T).len),
			Def!(Array!(T).length),
			Repr!(Array!(T).toString)
		);

		ex_python_to_d((PyObject* o) {
			
			auto size = PyList_Size(o);
			auto arr = Array!(T)(size.to!size_t);

			foreach(idx; 0..size) {
				auto item = PyList_GetItem(o, idx.to!Py_ssize_t);
				arr[idx] = python_to_d!(T*)(cast(PyObject*)item);
			}
			return arr;
		});
	}
}

opencopter.wake.InducedVelocities compute_wake_induced_velocities(ref PyWake wake, immutable Chunk x, immutable Chunk y, immutable Chunk z, ref PyAircraftState ac_state, size_t rotor_idx, bool single_rotor = false) {
	return opencopter.wake.compute_wake_induced_velocities(wake, x, y, z, ac_state, rotor_idx, 0, single_rotor);
}

void compute_blade_vectors(ref PyBladeGeometry blade) {
	opencopter.aircraft.compute_blade_vectors(blade);
}

Mat4 Mat4_identity() {
	return Mat4.identity();
}

Mat3 Mat3_identity() {
	return Mat3.identity();
}

extern(C) void PydMain() {
	import std.stdio : writeln;

	def!FrameType_aircraft;
	def!FrameType_connection;
	def!FrameType_rotor;
	def!FrameType_blade;
	def!FrameType_wing;

	def!(Mat4_identity);
	def!(Mat3_identity);

	def!(compute_blade_vectors);

	def!(compute_wake_induced_velocities);

	def!(chunk_size, Docstring!q{
		Returns the chunk size used by the internal data structures
	});

	def!(build_base_vtu_rotor, Docstring!q{
		Construct and allocate all required data for a :class:`VtkRotor`.

		:param rotor: An instance of a :class:`RotorGeometry` object
		:return: An initialized :class:`VtkRotor` object
	});
	def!(write_rotor_vtu, Docstring!q{
		Takes the :class:`RotorState` and :class:`RotorInput` data, converts it a :class:`VtkRotor`,
		and saves it to a vtu file.

		The filaname that the data is saved to is constructed by the function such that::

			filename = base_filename + "_" + rotor_idx + "_" + iteration + ".vtu"

		This ensures that timeseries files can be iterpreted by paraview correctly

		.. caution::
			
			This function has a relatively high runtime cost, use sparingly

		:param base_filename: The filename to save the data to
		:param iteration: The iteration of the saved data
		:param rotor_idx: The index of the rotor being saved
		:param vtk_rotor: The :class:`VtkRotor` object to save out
		:param rotor_state: The :class:`RotorState` object containing the rotor state data to save
		:param rotor_input: the :class:`RotorInput` object containing the rotor position data
		:param rotor_geom: the :class:`RotorGeometry` object containing the rotor position data
	});

	def!(write_rotors_vtu, Docstring!q{
		Takes the :class:`AircraftState` and :class:`AircraftInput` data, converts it :class:`VtkRotor`s,
		and saves it to a vtu file.

		The filaname that the data is saved to is constructed by the function such that::

			filename = base_filename + "_" + rotor_idx + "_" + iteration + ".vtu"

		This ensures that timeseries files can be iterpreted by paraview correctly

		.. caution::
			
			This function has a relatively high runtime cost, use sparingly

		:param base_filename: The filename to save the data to
		:param iteration: The iteration of the saved data
		:param rotor_idx: The index of the rotor being saved
		:param vtk_rotor: The :class:`VtkRotor[]` objects to save out
		:param rotor_state: The :class:`AircraftState` object containing the rotor state data to save
		:param rotor_input: the :class:`AircraftInput` object containing the rotor position data
		:param rotor_geom: the :class:`Aircraft` object containing the rotor position data
	});


	def!(build_base_vtu_wake, Docstring!q{
		Construct and allocate all required data for a :class:`VtkWake`.

		:param wake: An instance of a :class:`Wake` object
		:return: An initialized :class:`VtkWake` object
	});
	def!(write_wake_vtu, Docstring!q{
		Takes the :class:`Wake` data, converts it a :class:`VtkWake`,
		and saves it to several vtu files. One file is saved per blade tip vortex filament,
		and one file is saved per blade shed wake.

		The filanames that the data is saved to is constructed by the function such that::

			tip_vortex_filename = base_filename + "_" + rotor_idx + "_" + blade_idx + iteration + ".vtu"
			
			shed_wake_filename = base_filename + "_shed_" + rotor_idx + "_" + blade_idx + iteration + ".vtu"

		This ensures that timeseries files can be iterpreted by paraview correctly

		.. caution::
			
			This function has a very high runtime cost, use sparingly

		:param base_filename: The filename to save the data to
		:param iteration: The iteration of the saved data
		:param vtk_wake: The :class:`VtkWake` object to save out
		:param wake: The :class:`Wake` object containing the wake data
	});

	def!(write_wake_field_vtu, Docstring!q{
		Takes the :class:`Wake` data and writes the induced velocity field to a vtu file.

		.. caution::
			
			This function has a very high runtime cost, use sparingly

		:param filename: The filename to save the data to
		:param ac_state: The aircraft state object to save
		:param wake: The :class:`Wake` object containing the wake data
	});

	def!(basic_aircraft_rotor_dynamics, Docstring!q{
		Update all rotor azimuthal positions using a basic dynamics equation with simple
		forward Euler integration:

		:math:`\frac{\mathrm{d}\psi}{\mathrm{d}t} = \Omega + \dot\Omega \Delta t`
		
		:param input_state: The :class:`AircraftInputState` to update
		:param dt: The current time step size
	});

	def!(basic_single_rotor_dynamics, Docstring!q{
		Updates a single rotor azimuthal positions using a basic dynamics equation with simple
		forward Euler integration:

		:math:`\frac{\mathrm{d}\psi}{\mathrm{d}t} = \Omega + \dot\Omega \Delta t`
		
		:param input_state: The :class:`RotorInputState` to update
		:param dt: The current time step size
	});

	def!(opencopter.aircraft.generate_radius_points, Docstring!q{
		Generate a list of spanwise radial stations with the appropriate spacing

		.. attention::

			As we use different internal representation the size of the returned array
			may be different than the requested size so that it is a multipe of :func:`chunk_size`.

		:param n_sections: The desired number of radial stations.
		:return: The list of radial stations. The length of this list may be different than n_sections
	});

	def!(set_twist, Docstring!q{
		Set the spanwise twist distribution of a :class:`BladeGeometry` from a linear array.

		:param bg: :class:`BladeGeometry` object to apply twist to.
		:param data: List of twist angles (in radians) for each radial station
	});

	def!(get_twist, Docstring!q{
		Get the spanwise twist distribution of a :class:`BladeGeometry`

		:param bg: :class:`BladeGeometry` object get twist from.
	});

	def!(set_chord, Docstring!q{
		Set the spanwise chord distribution (non-dimensional) of a :class:`BladeGeometry` from a linear array.

		:param bg: :class:`BladeGeometry` object to apply twist to.
		:param data: List of chord lengths (non-dimensional) for each radial station
	});

	def!(get_chord, Docstring!q{
		Get the spanwise chord distribution (non-dimensional) of a :class:`BladeGeometry`

		:param bg: :class:`BladeGeometry` object get chord from.
	});

	def!(get_r, Docstring!q{
		Get the spanwise radial station distribution (non-dimensional) of a :class:`BladeGeometry`

		:param bg: :class:`BladeGeometry` object get radial stations from.
	});

	def!(set_r, Docstring!q{
		Set the spanwise radial stations (non-dimensional) of a :class:`BladeGeometry` from a linear array.

		:param bg: :class:`BladeGeometry` object to apply twist to.
		:param data: List of radial stations
	});
	
	def!(set_C_l_alpha, Docstring!q{
		Set the spanwise lift curve slope (:math:`C_{L_{\alpha}}`)(per radian) distribution of a
		:class:`BladeGeometry` from a linear array.

		:param bg: :class:`BladeGeometry` object to apply twist to.
		:param data: List of lift curve slope (per radians) for each radial station
	});
	
	def!(set_alpha_0, Docstring!q{
		Set the spanwise 0 lift angle (:math:`\alpha_0`)(in radians) distribution of a
		:class:`BladeGeometry` from a linear array.

		:param bg: :class:`BladeGeometry` object to apply twist to.
		:param data: List of :math:`\alpha_0` (in radians) for each radial station
	});

	def!(set_sweep, Docstring!q{
		Set the spanwise sweep angle (in radians) distribution of a :class:`BladeGeometry` from a linear array.

		:param bg: :class:`BladeGeometry` object to apply twist to.
		:param data: List of sweep angles (in radians) for each radial station
	});

	def!(set_xi, Docstring!q{
		Set the spanwise quarter chord position distribution of a :class:`BladeGeometry` from a linear array.

		:param bg: :class:`BladeGeometry` object to apply twist to.
		:param data: List of sweep angles (in radians) for each radial station
	});

	def!(set_xi_p, Docstring!q{
		Set the spanwise quarter chord derivative position distribution of a :class:`BladeGeometry`
		from a linear array.

		:param bg: :class:`BladeGeometry` object to apply twist to.
		:param data: List of sweep angles (in radians) for each radial station
	});

	def!(opencopter.aircraft.sweep_from_quarter_chord, Docstring!q{
		Set the spanwise quarter chord derivative position distribution of a :class:`BladeGeometry`
		from a linear array.

		:param r: spanwise radial station
		:param x: x location of the quarter chord
	});

	def!(get_dC_N, double[] function(ref PyBladeState), Docstring!q{
		Extract blade spanwise normal force coefficient to a linear array.

		:param blade_state: the :class:`BladeState` to extract the spanwise :math:`dC_N` from
		:return: List of spanwise :math:`dC_N` values
	});

	def!(fill_dC_N, void function(ref PyBladeState, double[]), Docstring!q{
		Extract blade spanwise normal force coefficient to a linear array.

		:param blade_state: the :class:`BladeState` to extract the spanwise :math:`dC_N` from
		:return: List of spanwise :math:`dC_N` values
	});

	def!(fill_dC_Nf, void function(ref PyBladeState, float[]), Docstring!q{
		Extract blade spanwise normal force coefficient to a linear array.

		:param blade_state: the :class:`BladeState` to extract the spanwise :math:`dC_N` from
		:return: List of spanwise :math:`dC_N` values
	});

	def!(get_dC_c, double[] function(ref PyBladeState), Docstring!q{
		Extract blade spanwise chord wise force coefficient to a linear array.

		:param blade_state: the :class:`BladeState` to extract the spanwise :math:`dC_c` from
		:return: List of spanwise :math:`dC_c` values
	});

	def!(fill_dC_c, void function(ref PyBladeState, double[]), Docstring!q{
		Extract blade spanwise chord wise force coefficient to a linear array.

		:param blade_state: the :class:`BladeState` to extract the spanwise :math:`dC_c` from
		:return: List of spanwise :math:`dC_c` values
	});

	def!(fill_dC_D, void function(ref PyBladeState, double[]), Docstring!q{
		Extract blade spanwise chord wise force coefficient to a linear array.

		:param blade_state: the :class:`BladeState` to extract the spanwise :math:`dC_D` from
		:return: List of spanwise :math:`dC_D` values
	});

	def!(fill_dC_cf, void function(ref PyBladeState, float[]), Docstring!q{
		Extract blade spanwise chord wise force coefficient to a linear array.

		:param blade_state: the :class:`BladeState` to extract the spanwise :math:`dC_c` from
		:return: List of spanwise :math:`dC_c` values
	});

	def!(fill_dC_Df, void function(ref PyBladeState, float[]), Docstring!q{
		Extract blade spanwise chord wise force coefficient to a linear array.

		:param blade_state: the :class:`BladeState` to extract the spanwise :math:`dC_D` from
		:return: List of spanwise :math:`dC_D` values
	});

	def!(get_dC_T, double[] function(ref PyBladeState), Docstring!q{
		Extract blade spanwise thrust coefficient to a linear array.

		:param blade_state: the :class:`BladeState` to extract the spanwise :math:`dC_T` from
		:return: List of spanwise :math:`dC_T` values
	});

	def!(fill_dC_Td, void function(ref PyBladeState, double[]), Docstring!q{
		Extract blade spanwise thrust coefficient to a linear double precision floating point array.

		:param blade_state: the :class:`BladeState` to extract the spanwise :math:`dC_T` from
		:param data: numpy slice to fill with :math:`dC_T` values
	});

	def!(fill_dC_Tf, void function(ref PyBladeState, float[]), Docstring!q{
		Extract blade spanwise thrust coefficient to a linear single precision floating point array.

		:param blade_state: the :class:`BladeState` to extract the spanwise :math:`dC_T` from
		:param data: numpy slice to fill with :math:`dC_T` values
	});

	def!(fill_dC_Qd, void function(ref PyBladeState, double[]), Docstring!q{
		Extract blade spanwise torque coefficient to a linear double precision floating point array.

		:param blade_state: the :class:`BladeState` to extract the spanwise :math:`dC_Q` from
		:param data: numpy slice to fill with :math:`dC_Q` values
	});

	def!(fill_dC_Qf, void function(ref PyBladeState, float[]), Docstring!q{
		Extract blade spanwise torque coefficient to a linear single precision floating point array.

		:param blade_state: the :class:`BladeState` to extract the spanwise :math:`dC_Q` from
		:param data: numpy slice to fill with :math:`dC_Q` values
	});

	def!(get_dC_Q, double[] function(ref PyBladeState), Docstring!q{
		Extract blade spanwise torque coefficient to a linear array.

		:param blade_state: the :class:`BladeState` to extract the spanwise :math:`Cd_Q` from
		:return: List of spanwise :math:`dC_Q` values
	});
	
	def!(get_aoa, double[] function(ref PyBladeState), Docstring!q{
		Extract blade spanwise angle of attack (:math:`{\alpha}`) to a linear array.

		:param blade_state: the :class:`BladeState` to extract the spanwise :math:`{\alpha}` from
		:return: List of spanwise :math:`{\alpha}`
	});

	def!(get_theta, double[] function(ref PyBladeState), Docstring!q{
		Extract blade spanwise angle of attack (:math:`{\theta}`) to a linear array.

		:param blade_state: the :class:`BladeState` to extract the spanwise :math:`{\theta}` from
		:return: List of spanwise :math:`{\theta}`
	});

	def!(fill_aoad, void function(ref PyBladeState, double[]), Docstring!q{
		Extract blade spanwise angle of attack to a linear double precision floating point array.

		:param blade_state: the :class:`BladeState` to extract the spanwise :math:`dC_T` from
		:param data: numpy slice to fill with :math:`dC_T` values
	});

	def!(fill_aoaf, void function(ref PyBladeState, float[]), Docstring!q{
		Extract blade spanwise angle of attack to a linear single precision floating point array.

		:param blade_state: the :class:`BladeState` to extract the spanwise :math:`dC_T` from
		:param data: numpy slice to fill with :math:`dC_T` values
	});

	def!(fill_u_td, void function(ref PyBladeState, double[]), Docstring!q{
		Extract blade spanwise tangential velocity component to a linear double precision floating point array.

		:param blade_state: the :class:`BladeState` to extract the spanwise :math:`dC_T` from
		:param data: numpy slice to fill with :math:`dC_T` values
	});

	def!(fill_u_tf, void function(ref PyBladeState, float[]), Docstring!q{
		Extract blade spanwise tangential velocity component to a linear single precision floating point array.

		:param blade_state: the :class:`BladeState` to extract the spanwise :math:`dC_T` from
		:param data: numpy slice to fill with :math:`dC_T` values
	});

	def!(fill_u_pd, void function(ref PyBladeState, double[]), Docstring!q{
		Extract blade spanwise perpendicular velocity component to a linear double precision floating point array.

		:param blade_state: the :class:`BladeState` to extract the spanwise :math:`dC_T` from
		:param data: numpy slice to fill with :math:`dC_T` values
	});

	def!(fill_u_pf, void function(ref PyBladeState, float[]), Docstring!q{
		Extract blade spanwise perpendicular velocity component to a linear single precision floating point array.

		:param blade_state: the :class:`BladeState` to extract the spanwise :math:`dC_T` from
		:param data: numpy slice to fill with :math:`dC_T` values
	});

	def!(get_gamma, double[] function(ref PyBladeState), Docstring!q{
		Extract blade spanwise bound circulation (:math:`{\Gamma}`) to a linear array.

		:param blade_state: the :class:`BladeState` to extract the spanwise :math:`{\Gamma}` from
		:return: List of spanwise circulation values
	});

	def!(get_u_p, double[] function(ref PyBladeState), Docstring!q{
		Extract blade spanwise bound circulation (:math:`{u_p}`) to a linear array.

		:param blade_state: the :class:`BladeState` to extract the spanwise :math:`{u_p}` from
		:return: List of spanwise circulation values
	});

	def!(get_u_t, double[] function(ref PyBladeState), Docstring!q{
		Extract blade spanwise bound circulation (:math:`{u_t}`) to a linear array.

		:param blade_state: the :class:`BladeState` to extract the spanwise :math:`{u_p}` from
		:return: List of spanwise circulation values
	});

	def!(get_inflow_angle, double[] function(ref PyBladeState), Docstring!q{
		Extract blade spanwise bound circulation (:math:`{\Phi}`) to a linear array.

		:param blade_state: the :class:`BladeState` to extract the spanwise :math:`{\Phi}` from
		:return: List of spanwise circulation values
	});

	def!(get_aoa_eff, double[] function(ref PyBladeState), Docstring!q{
		Extract blade spanwise effective angle of attack (:math:`{\alpha_{eff}}`) to a linear array.

		:param blade_state: the :class:`BladeState` to extract the spanwise :math:`{\alpha_{eff}}` from
		:return: List of spanwise effective angle of attack values
	});

	def!(get_d_gamma, double[] function(ref PyBladeState), Docstring!q{
		Extract blade spanwise change in bound circulation (:math:`{\Delta}`:math:`{\Gamma}`) to a linear array.

		:param blade_state: the :class:`BladeState` to extract the spanwise :math:`{\Delta}`:math:`{\Gamma}` from
		:return: List of spanwise :math:`{\Delta}`:math:`{\Gamma}` values
	});

	def!(get_dC_L, double[] function(ref PyBladeState), Docstring!q{
		Extract blade spanwise lift coefficient (:math:`dC_L`) to a linear array.

		:param blade_state: the :class:`BladeState` to extract the spanwise :math:`C_L` from
		:return: List of spanwise :math:`dC_L` values
	});
	
	def!(get_dC_D, double[] function(ref PyBladeState), Docstring!q{
		Extract blade spanwise drag coefficient (:math:`dC_D`) to a linear array.

		:param blade_state: the :class:`BladeState` to extract the spanwise :math:`C_D` from
		:return: List of spanwise :math:`dC_D` values
	});
	
	def!(get_wake_x_component, double[] function(ref PyVortexFilament), Docstring!q{
		Extract filament x position component to a linear array.

		:param votex_filament: the :class:`VortexFilament` to extract the x position component from
		:return: List of filament x positions
	});
	
	def!(get_wake_y_component, double[] function(ref PyVortexFilament), Docstring!q{
		Extract filament y position component to a linear array.

		:param votex_filament: the :class:`VortexFilament` to extract the y position component from
		:return: List of filament y positions
	});
	
	def!(get_wake_z_component, double[] function(ref PyVortexFilament), Docstring!q{
		Extract filament z position component to a linear array.

		:param votex_filament: the :class:`VortexFilament` to extract the z position component from
		:return: List of filament z positions
	});
	
	def!(get_wake_gamma_component, double[] function(ref PyVortexFilament), Docstring!q{
		Extract filament cirulation strength (:math:`\Gamma`) to a linear array.

		:param votex_filament: the :class:`VortexFilament` to extract :math:`\Gamma` from
		:return: List of filament :math:`\Gamma`
	});

	def!(get_wake_r_c_component, double[] function(ref PyVortexFilament), Docstring!q{
		Extract filament vortex core size (:math:`r_c`) to a linear array.

		:param votex_filament: the :class:`VortexFilament` to extract :math:`r_c` from
		:return: List of filament core sizes
	});

	def!(get_wake_v_z_component, double[] function(ref PyVortexFilament), Docstring!q{
		Extract induced velocity (:math:`v_z`) acting upon the filament to a linear array.

		:param votex_filament: the :class:`VortexFilament` to extract :math:`v_z` from
		:return: List of induced velocities
	});
	


	def!(fill_wake_x_component, void function(ref PyVortexFilament, double[]), Docstring!q{
		Extract filament x position component to a linear array.

		:param votex_filament: the :class:`VortexFilament` to extract the x position component from
		:return: List of filament x positions
	});
	
	def!(fill_wake_y_component, void function(ref PyVortexFilament, double[]), Docstring!q{
		Extract filament y position component to a linear array.

		:param votex_filament: the :class:`VortexFilament` to extract the y position component from
		:return: List of filament y positions
	});
	
	def!(fill_wake_z_component, void function(ref PyVortexFilament, double[]), Docstring!q{
		Extract filament z position component to a linear array.

		:param votex_filament: the :class:`VortexFilament` to extract the z position component from
		:return: List of filament z positions
	});
	
	def!(fill_wake_gamma_component, void function(ref PyVortexFilament, double[]), Docstring!q{
		Extract filament cirulation strength (:math:`\Gamma`) to a linear array.

		:param votex_filament: the :class:`VortexFilament` to extract :math:`\Gamma` from
		:return: List of filament :math:`\Gamma`
	});

	def!(fill_wake_r_c_component, void function(ref PyVortexFilament, double[]), Docstring!q{
		Extract filament vortex core size (:math:`r_c`) to a linear array.

		:param votex_filament: the :class:`VortexFilament` to extract :math:`r_c` from
		:return: List of filament core sizes
	});

	def!(fill_wake_v_z_component, void function(ref PyVortexFilament, double[]), Docstring!q{
		Extract induced velocity (:math:`v_z`) acting upon the filament to a linear array.

		:param votex_filament: the :class:`VortexFilament` to extract :math:`v_z` from
		:return: List of induced velocities
	});


	def!(step, Docstring!(q{
		Step the simulation by one timestep

		:param ac_state: :class:`AircraftState` object holding the current state of the aircraft
		:param aircraft: :class:`Aircraft` object describing the aircraft geometry
		:param ac_input_state: :class:`AircraftInputState` describing the aircrafts inputs
		:param inflows: a list of :class:`Inflow` objects. One for each rotor
		:param wake_history: :class:`WakeHistory` object containing the current and past wake states
		:param atmo: :class:`Atmosphere` object describing ambient atmospheric conditions
		:param iteration: the current iteration
		:param dt: the current timestep size
	}));

	def!(create_aerodas_from_xfoil_polar);

	def!(load_c81_file);

	def!(write_inflow_vtu);

	module_init;

	wrap_class!(
		AirfoilModel
	);

	wrap_class!(
		BladeAirfoil,
		Init!(AirfoilModel[], size_t[2][])
	);

	wrap_class!(
		AeroDAS,
		Init!(double[], double[], double[], double, double),
		Def!(AeroDAS.get_Cl, double function(double, double)),
		Def!(AeroDAS.get_Cd, double function(double, double)),
		Member!("CL"),
		Member!("CD"),
		Member!("alpha")
	);

	wrap_class!(
		C81,
		Init!(string, double[], double[], double[][],
			  double[], double[], double[][],
			  double[], double[], double[][]
		),
		Def!(C81.get_Cl, double function(double, double)),
		Def!(C81.get_Cd, double function(double, double))
	);

	wrap_class!(
		ThinAirfoil,
		Init!(double),
		Def!(ThinAirfoil.get_Cl, double function(double, double)),
		Def!(ThinAirfoil.get_Cd, double function(double, double))
	);

	wrap_class!(
		opencopter.vtk.VtkRotor,
		Docstring!q{
			This is an opaque type that holds all the relevant data
			for writing out a vtu file that holds the rotor geometric data
			as well as various spanwise data.

			.. attention::
				
				Do not instantiate one of these yourself. Instead use :func:`build_base_vtu_rotor`
		}
	);

	wrap_class!(
		opencopter.vtk.VtkWake,
		Docstring!q{
			This is an opaque type that holds all the relevant data
			for writing out a vtu file that holds the wake geometric data
			as well as various data in a filament.

			.. attention::
			
				Do not instantiate one of these yourself. Instead use :func:`build_base_vtu_wake`
		}
	);

	
	wrap_struct!(
		opencopter.wake.InducedVelocities,
		Member!"v_x",
		Member!"v_y",
		Member!"v_z"
	);

	wrap_struct!(
		opencopter.wake.FilamentChunk,
		Member!("x", Docstring!q{A chunk of x positions along the filament}),
		Member!("y", Docstring!q{A chunk of y positions along the filament}),
		Member!("z", Docstring!q{A chunk of z positions along the filament}),
		Member!("gamma", Docstring!q{A chunk of circulation strengths along the filament}),
		//Member!"phi",
		Member!("r_c", Docstring!q{A chunk of vortex core size along the filament}),
		//Member!"r_0",
		//Member!"dx",
		//Member!"dy",
		//Member!"dz",
		//Member!"l_0",
		Member!("v_z", Docstring!q{A chunk of induced velocities along the filament}),
		//Member!"volume",
		//Member!"d_volume"
	);

	wrap_struct!(
		PyVortexFilament,
		PyName!"VortexFilament",
		Init!size_t,
		Member!("chunks", Docstring!q{An array of :class:`FilamentChunk`}),
		//Member!"length"
	);

	wrap_struct!(
		PyShedVortex,
		PyName!"ShedVortex",
		//Init!(size_t, size_t),
		Member!("shed_filaments", Docstring!q{An array of :class:`VortexFilament`}),
	);

	wrap_struct!(
		PyRotorWake,
		PyName!"RotorWake",
		Init!(size_t),
		Member!("tip_vortices", Docstring!q{An array of :class:`VortexFilament`}),
		Member!("shed_vortices", Docstring!q{An array of :class:`ShedVortex`})
	);

	wrap_struct!(
		PyWake,
		PyName!"Wake",
		Init!(size_t, size_t[], size_t[], size_t, size_t),
		Member!("rotor_wakes", Docstring!q{An array of :class:`RotorWake`})
	);

	wrap_struct!(
		PyWakeHistory,
		PyName!"WakeHistory",
		Init!(size_t, size_t[], size_t[], size_t, size_t, size_t, double),
		Docstring!q{
			Top level structure for holding the wake and its history.
			
			Contructor:

			:param num_rotors: Number of rotors the aircraft has.
			:param num_blades: A list of the number of blades per rotor. The length of the list must match the number of rotors.
			:param wake_history: The length of the wake history per rotor. The length of the list must match the number of rotors.
			:param time_history: The number of timesteps to store. 2 is the minimum number.
			:param radial_elements: The number of radial elements down the blade span.
			:param shed_history: The number of timesteps to store the shed wake for.
			:param a1: Wake core growth constant.
		},
		Member!("history", Docstring!q{An array of :class:`Wake` s, one for each timestep})
	);

	wrap_struct!(
		Atmosphere, 
		Docstring!q{
			Object describing basic atmosphere attributes.
			Make sure units are consistent with everything else.

			Constructor:
			
			:param density: Atmospheric density.
			:param dynamic_viscosity: Atmospheric dynamic viscosity.
		},
		Init!(double, double),
		Member!("density", Mode!"r", Docstring!"atmospheric density"),
		Member!("dynamic_viscosity", Mode!"r", Docstring!"atmospheric dynamic viscosity"),
		Member!("kinematic_viscosity", Mode!"r", Docstring!"atmospheric kinematic viscosity")
	);

	wrap_struct!(
		PyRotorInputState,
		PyName!("RotorInputState"),
		Member!("angular_velocity", Docstring!q{The angular velocity of the rotor in :math:`mathrm{rad}/s`}),
		Member!("angular_accel", Docstring!q{The angular acceleration of the rotor in :math:`mathrm{rad}/s^2`}),
		Member!("azimuth", Docstring!q{The current azimuthal position of the rotor in radians}),
		Member!("r_0", Docstring!q{The initial non-dimensional tip vortex core size}),
		Member!("blade_pitches", Docstring!q{An array of blade pitches taken at the 75% spanwise location in radians}),
		Member!("blade_flapping_rate", Docstring!q{An array of blade flapping rates in :math:`mathrm{rad}/s`}),
		Member!("blade_flapping", Docstring!q{An array of blade flapping angle in radians}),
	);

	wrap_struct!(
		PyAircraftInputState,
		PyName!"AircraftInputState",
		Init!(size_t, size_t[]),
		Docstring!q{
			This class represents the input parameters for the entire aircraft.
			
			Contructor:

			:param num_rotors: Number of rotors the aircraft has.
			:param num_blades: The number of blades each rotor has.
		},
		Member!("rotor_inputs", Docstring!q{An array of :class:`RotorInputState`, one for each rotor})
	);

	wrap_struct!(
		opencopter.aircraft.BladeGeometryChunk,
		Member!"twist",
		Member!"chord",
		Member!"r",
		Member!"C_l_alpha",
		Member!"alpha_0"
	);

	wrap_struct!(
		PyBladeGeometry,
		PyName!"BladeGeometry",
		Init!(size_t, double, double, BladeAirfoil, double),
		Docstring!q{
			This class allocates and holds the blade geomteric parameters.

			Constructor:

			:param num_elements: The number of spanwise blade elements.
			:param azimuth_offset: The azimuthal offset that this blade has relative to its rotor.
			:param average_chord: The average chord of the blade.
			:param r_c: The dimensional root cutout of the blade.
		},
		Member!"chunks",
		Member!"r_c",
		Member!"airfoil",
		Member!("azimuth_offset", Docstring!q{
			The azimuthal offset for this blade. This is added to the :class:`RotorInputState` azimuth.
		}),
		Member!("average_chord", Docstring!q{The dimensional average chord of the blade}),
		Member!("blade_length", Docstring!q{The dimensional actual length of the blade not accounting for root cutout}),
		Member!"frame"
	);

	wrap_struct!(
		PyRotorGeometry,
		PyName!"RotorGeometry",
		Init!(size_t, Vec3, double, double),
		Docstring!q{
			This class allocates and holds all the blades and other rotor geometric data.

			Constructor:

			:param num_blades: Number of blades the rotor has.
			:param origin: The global origin of the rotor. This is where the center of the hub is located.
			:param radius: The dimensional radius of the rotor.
			:param solidity: The solidity of the rotor.
		},
		Member!("blades", Docstring!q{An array of :class:`BladeGeomtery`, one for each blade of the rotor}),
		Member!("origin", Docstring!q{The global origin of the rotor. This is where the center of the hub is located.}),
		Member!("radius", Docstring!q{The dimensional radius of the rotor.}),
		Member!("solidity", Docstring!q{The solidity of the rotor.}),
		Member!"frame"
	);

	wrap_struct!(
		PyAircraft,
		PyName!"Aircraft",
		Init!size_t,
		Docstring!q{
			This class allocates and holds an array of rotors for the aircraft.

			Constructor:

			:param num_rotors: The number of rotors the aircraft has.
		},
		Member!("rotors", Docstring!q{An array of :class:`RotorGeometry`, one for each rotor on the aircraft}),
		Member!"root_frame"
	);

	wrap_struct!(
		opencopter.aircraft.BladeStateChunk,
		Member!("dC_T", Docstring!q{A chunk of spanwise sectional thrust coefficient}),
		Member!("dC_T_dot", Docstring!q{A chunk of spanwise sectional thrust coefficient time derivative}),
		Member!("dC_Q", Docstring!q{A chunk of spanwise sectional power/torque coefficient}),
		Member!("dC_L", Docstring!q{A chunk of spanwise sectional lift coefficient}),
		Member!("dC_L_dot", Docstring!q{A chunk of spanwise sectional lift coefficient time derivative}),
		Member!("dC_N", Docstring!q{A chunk of spanwise sectional normal force coefficient}),
		Member!("dC_c", Docstring!q{A chunk of spanwise sectional chordwise force coefficient}),
		Member!("dC_Mx", Docstring!q{A chunk of spanwise sectional x moment coefficient}),
		Member!("dC_My", Docstring!q{A chunk of spanwise sectional y moment coefficient}),
		Member!("u_p", Docstring!q{A chunk of spanwise sectional perpendicular inflow}),
		Member!("u_t", Docstring!q{A chunk of spanwise sectional tangential inflow}),
		Member!("aoa", Docstring!q{A chunk of spanwise sectional angle of attack}),
		Member!("aoa_eff", Docstring!q{A chunk of spanwise sectional effective angle of attack}),
		Member!("inflow_angle", Docstring!q{A chunk of spanwise sectional inflow angle}),
		Member!("gamma", Docstring!q{A chunk of spanwise sectional circulation}),
		Member!("d_gamma", Docstring!q{A chunk of spanwise sectional change in circulation}),
		Member!("x", Docstring!q{The current x position of the blade in world space}),
		Member!("y", Docstring!q{The current y position of the blade in world space}),
		Member!("z", Docstring!q{The current z position of the blade in world space}),
		Member!("theta", Docstring!q{The currect pitch of the spanwise element}),
		Member!("d_gamma", Docstring!q{Spanwise change in circulation}),
		Member!("gamma", Docstring!q{Spanwise circulation}),
	);

	wrap_struct!(
		PyBladeState,
		PyName!"BladeState",
		//Init!(size_t, PyBladeGeometry*),
		Docstring!q{
			This class holds the current aerodynamic state of the blade.

			.. attention::

				This never needs to be explicitly constructed as it will be constructed by :class:`AircraftState`
		},
		Member!"chunks",
		Member!("azimuth", Docstring!q{The current azimuthal angle of this blade}),
		Member!("C_T", Docstring!q{The current thrust coefficient contribution from this blade}),
		Member!("C_Q", Docstring!q{The current power/torque coefficient contribution from this blade}),
		Member!("C_L", Docstring!q{The current lift coefficient contribution from this blade}),
		Member!("C_D", Docstring!q{The current drag coefficient contribution from this blade}),
		Member!("C_Mx", Docstring!q{The current x moment coefficient contribution from this blade}),
		Member!("C_My", Docstring!q{The current y moment coefficient contribution from this blade})
	);

	wrap_struct!(
		PyRotorState,
		PyName!"RotorState",
		//Init!(size_t, size_t, PyRotorGeometry*),
		Docstring!q{
			This class holds the current aerodynamic state of the rotor.

			.. attention::

				This never needs to be explicitly constructed as it will be constructed by :class:`AircraftState`
		},
		Member!("blade_states", Docstring!q{An array of :class:`BladeState`, one for each blade of this rotor}),
		Member!("C_T", Docstring!q{The current thrust coefficient of this rotor}),
		Member!("C_Q", Docstring!q{The current power/torque coefficient of this rotor}),
		Member!("C_Mx", Docstring!q{The current x moment coefficient of this rotor}),
		Member!("C_My", Docstring!q{The current y moment coefficient of this rotor}),
		Member!("advance_ratio", Docstring!q{The current advance ratio of this rotor}),
		Member!("axial_advance_ratio", Docstring!q{The current axial advance ratio of this rotor})
	);

	wrap_struct!(
		PyAircraftState,
		PyName!"AircraftState",
		Init!(size_t, size_t[], size_t, PyAircraft*),
		Docstring!q{
			This is the top level class that holds the current aerodynamic state of an aircraft.

			Constructor:

			:param num_rotors: The number of rotors on this aircraft
			:param num_blades: The number of blades on each rotor
			:param num_elements: The requested number of spanwise elements. This will be rounded up to the nearest :func:`chunk_size`.
			:param ac: The :class:`Aircraft` object associated with this state
		},
		Member!("rotor_states", Docstring!q{An array of :class:`RotorState`, one for each rotor on the aircraft}),
		Member!"freestream"
	);

	wrap_class!(
		Inflow,
		Def!(Inflow.update),
		Def!(Inflow.inflow_at, Chunk function(immutable Chunk, immutable Chunk, immutable Chunk, immutable Chunk, double), PyName!"inflow_at_xyz"),
		Def!(Inflow.inflow_at, Chunk function(immutable Chunk, immutable double, immutable double), PyName!"inflow_at_r"),
		Def!(Inflow.wake_skew)
	)();

	wrap_class!(
		HuangPeters,
		Init!(long, long, PyRotorGeometry*, double),
		Docstring!q{
			This class instantiates a dynamic inflow model for a single rotor.
			One of these will be needed for each rotor in the aircraft.

			Constructor:

			:param Mo: The number of odd modes used in the model. 4 typically works well.
			:param Me: The number of even modes used in the model. 2 typically works well.
			:param rotor: The geometry of the rotor this inflow model is modeling.
			:param dt: The timestep of the simulation.
		},
		Def!(HuangPeters.update, Docstring!q{
			Updates the inflow model by one timestep

			.. attention

				You will likely never have to call this function yourself. This is called automaticall
				in the :func:`step` function.

			:param C_T: Current rotor thrust coefficient. This does nothing for this inflow model.
			:param rotor: The current input state for the rotor.
			:param rotor_state: The current rotor state.
			:param advance_ratio: The current advance ratio for the rotor.
			:param axial_advance_ratio: The current axial advance ratio for the rotor.
			:param dt: The current timestep size.
		}),
		Def!(HuangPeters.inflow_at,
			Chunk function(immutable Chunk, immutable Chunk, immutable Chunk, immutable Chunk, double),
			PyName!"inflow_at",
			Docstring!q{
				Computes the rotor induced flow at the requested location.

				:param x: A chunk of x positions to compute the induced velocity at.
				:param y: A chunk of y positions to compute the induced velocity at.
				:param z: A chunk of z positions to compute the induced velocity at.
				:param x_e: A chunk of x_e positions to compute the induced velocity at. Unused in this inflow model.
				:param angle_of_attack: Current angle of attack of the rotor. Unused in this inflow model.
				:return: A chunk of z induced velocities.
			}
		),
		Def!(HuangPeters.wake_skew, Docstring!q{
			:return: The current wake skew angle of the rotor in radians.
		})
	)();

	wrap_class!(
		Beddoes,
		Init!(),
		Docstring!q{
			This class instantiates a dynamic inflow model for a single rotor.
			One of these will be needed for each rotor in the aircraft.

			Constructor:

			:param Mo: The number of odd modes used in the model. 4 typically works well.
			:param Me: The number of even modes used in the model. 2 typically works well.
			:param rotor: The geometry of the rotor this inflow model is modeling.
			:param dt: The timestep of the simulation.
		},
		Def!(Beddoes.update, Docstring!q{
			Updates the inflow model by one timestep

			.. attention

				You will likely never have to call this function yourself. This is called automaticall
				in the :func:`step` function.

			:param C_T: Current rotor thrust coefficient. This does nothing for this inflow model.
			:param rotor: The current input state for the rotor.
			:param rotor_state: The current rotor state.
			:param advance_ratio: The current advance ratio for the rotor.
			:param axial_advance_ratio: The current axial advance ratio for the rotor.
			:param dt: The current timestep size.
		}),
		Def!(Beddoes.inflow_at,
			Chunk function(immutable Chunk, immutable Chunk, immutable Chunk, immutable Chunk, double),
			PyName!"inflow_at",
			Docstring!q{
				Computes the rotor induced flow at the requested location.

				:param x: A chunk of x positions to compute the induced velocity at.
				:param y: A chunk of y positions to compute the induced velocity at.
				:param z: A chunk of z positions to compute the induced velocity at.
				:param x_e: A chunk of x_e positions to compute the induced velocity at. Unused in this inflow model.
				:param angle_of_attack: Current angle of attack of the rotor. Unused in this inflow model.
				:return: A chunk of z induced velocities.
			}
		),
		//Def!(HuangPeters.inflow_at, Chunk function(immutable Chunk, immutable double, immutable double), PyName!"inflow_at_r"),
		Def!(Beddoes.wake_skew, Docstring!q{
			:return: The current wake skew angle of the rotor in radians.
		})
	)();

	import std.meta : AliasSeq, staticMap;
	import std.traits : FunctionTypeOf;

	wrap_struct!(
		Vec3,
		PyName!"Vec3",
		Init!(double[3]),
		OpIndex!(staticMap!(FunctionTypeOf, __traits(getOverloads, Vec3, "opIndex"))[1]),
	);

	wrap_struct!(
		Vec4,
		PyName!"Vec4",
		Init!(double[4]),
		OpIndex!(staticMap!(FunctionTypeOf, __traits(getOverloads, Vec4, "opIndex"))[1]),
	);

	wrap_struct!(
		Mat3,
		PyName!"Mat3",
		Init!(double[9]),
		OpIndex!(staticMap!(FunctionTypeOf, __traits(getOverloads, Mat3, "opIndex"))[1])
	);

	wrap_struct!(
		Mat4,
		PyName!"Mat4",
		Init!(double[16]),
		OpIndex!(staticMap!(FunctionTypeOf, __traits(getOverloads, Mat4, "opIndex"))[1])
	);

	wrap_struct!(
		Frame,
		Init!(Vec3, double, Vec3, Frame*, string, string),
		Member!"global_matrix",
		Member!"local_matrix",
		Member!"name",
		Member!"children",
		Member!"parent",
		Def!(Frame.update, void function(Mat4*), PyName!"update"),
		Def!(Frame.get_frame_type),
		Def!(Frame.set_frame_type),
		Def!(Frame.set_rotation),
		Def!(Frame.local_rotation_axis),
		Def!(Frame.local_rotation_angle),
		Def!(Frame.local_position),
		//Def!(Frame.set_parent),
		Def!(Frame.rotate),
		Def!(Frame.translate),
		Def!(Frame.global_position),
	);

	wrap_array!double;
	wrap_array!PyRotorGeometry;
	wrap_array!PyBladeGeometry;
	wrap_array!(opencopter.aircraft.BladeGeometryChunk);
	wrap_array!PyRotorInputState;
	wrap_array!PyAircraftState;
	wrap_array!PyRotorState;
	wrap_array!PyBladeState;
	wrap_array!(opencopter.aircraft.BladeStateChunk);
	wrap_array!PyWake;
	wrap_array!PyRotorWake;
	wrap_array!PyShedVortex;
	wrap_array!PyVortexFilament;
	wrap_array!(opencopter.wake.FilamentChunk);
	wrap_array!PyShedVortex;
}

// Boilerplate take from Pyd so I don't have to use distutils to build

import pyd.def;
import pyd.exception;
import pyd.thread;

//extern(C) void PydMain();

version(Python_3_0_Or_Later) {
    import deimos.python.Python;
    extern(C) export PyObject* PyInit_libopencopter() {
        return pyd.exception.exception_catcher(delegate PyObject*() {
                pyd.thread.ensureAttached();
                pyd.def.pyd_module_name = "opencopter";
                PydMain();
                return pyd.def.pyd_modules[""];
                });
    }
}else version(Python_2_4_Or_Later) {
    extern(C) export void initlibopencopter() {
        pyd.exception.exception_catcher(delegate void() {
                pyd.thread.ensureAttached();
                pyd.def.pyd_module_name = "opencopter";
                PydMain();
                });
    }
}else static assert(false);

extern(C) void _Dmain(){
    // make druntime happy
}

version(Windows) {
	import core.sys.windows.windows;
	import core.sys.windows.dll;

	__gshared HINSTANCE g_hInst;


	extern (Windows)
	BOOL DllMain(HINSTANCE hInstance, ULONG ulReason, LPVOID pvReserved)
	{
		switch (ulReason)
		{
			case DLL_PROCESS_ATTACH:
				g_hInst = hInstance;
				dll_process_attach( hInstance, true );
				break;

			case DLL_PROCESS_DETACH:
				dll_process_detach( hInstance, true );
				break;

			case DLL_THREAD_ATTACH:
				dll_thread_attach( true, true );
				break;

			case DLL_THREAD_DETACH:
				dll_thread_detach( true, true );
				break;

		default:
			assert(0);
		}

		return true;
	}
} else {
	// This file requires the .so be compiled with '-nostartfiles'.
	// Also note that this is inferior to the Windows version: it does not call the
	// static constructors or unit tests. As far as I can tell, this can't be done
	// until Phobos is updated to explicitly allow it.
	extern(C) shared bool _D2rt6dmain212_d_isHaltingOb;
	alias _D2rt6dmain212_d_isHaltingOb _d_isHalting;
	extern(C) {

		void rt_init();
		void rt_term();
		extern (C) void gc_init();
		extern (C) void gc_init_nothrow();
		extern (C) void _d_register_conservative_gc();
		extern (C) void  rt_moduleTlsDtor() nothrow;
		extern (C) void  rt_moduleTlsCtor() nothrow;

		version(LDC) {
			pragma(LDC_global_crt_ctor)
				void hacky_init() {
					rt_init();
				}

			pragma(LDC_global_crt_dtor)
				void hacky_fini() {
					if(!_d_isHalting){
						rt_term();
					}
				}
		}else{
			void hacky_init() {
				rt_init();
			}

			void hacky_fini() {
				if(!_d_isHalting){
					rt_term();
				}
			}
		}

	} /* extern(C) */
}
