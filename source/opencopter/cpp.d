module opencopter.cpp;

import opencopter.aircraft;
import opencopter.atmosphere;
static import opencopter.inflow;
import opencopter.math;
import opencopter.memory;
import opencopter.wake;

extern(C++) {

	RotorGeometry create_rotor_geometry(size_t num_blades, double x, double y, double z, double radius, double solidity) {
		return RotorGeometry(num_blades, Vec3(x, y, z), radius, solidity);
	}

	/++
	 +	This is really for fortran
	 +/
	Aircraft* create_aircraft(size_t num_rotors) {
		return new Aircraft(num_rotors);
	}

	/*void set_twist(ref BladeGeometry bg, double[] data) {
		bg.set_geometry_array!"twist"(data);
	}

	void set_chord(ref BladeGeometry bg, double[] data) {
		bg.set_geometry_array!"chord"(data);
	}

	void set_r(ref BladeGeometry bg, double[] data) {
		bg.set_geometry_array!"r"(data);
	}

	void set_C_l_alpha(ref BladeGeometry bg, double[] data) {
		bg.set_geometry_array!"C_l_alpha"(data);
	}

	void set_alpha_0(ref BladeGeometry bg, double[] data) {
		bg.set_geometry_array!"alpha_0"(data);
	}*/


	double[] get_dC_T(ref BladeState blade) {
		return blade.get_state_array!"dC_T";
	}

	void get_dC_T(ref BladeState blade, double* state_array, size_t length) {
		blade.get_state_array!"dC_T"(state_array[0..length]);
	}

	/+double[] get_dC_P(ref BladeState blade, size_t elements) {
		return blade.get_state_array!"dC_P"(elements);
	}

	void get_dC_P(ref BladeState blade, double* state_array, size_t length) {
		blade.get_state_array!"dC_P"(state_array[0..length]);
	}+/

	double[] get_dC_Q(ref BladeState blade) {
		return blade.get_state_array!"dC_Q";
	}

	void get_dC_Q(ref BladeState blade, double* state_array, size_t length) {
		blade.get_state_array!"dC_Q"(state_array[0..length]);
	}

	/+double[] get_inflow(ref BladeState blade, size_t elements) {
		return blade.get_state_array!"inflow"(elements);
	}

	void get_inflow(ref BladeState blade, double* state_array, size_t length) {
		blade.get_state_array!"inflow"(state_array[0..length]);
	}+/

	double[] get_dC_L(ref BladeState blade) {
		return blade.get_state_array!"dC_L";
	}

	void get_dC_L(ref BladeState blade, double* state_array, size_t length) {
		blade.get_state_array!"dC_L"(state_array[0..length]);
	}

	double[] get_dC_D(ref BladeState blade) {
		return blade.get_state_array!"dC_D";
	}

	void get_dC_D(ref BladeState blade, double* state_array, size_t length) {
		blade.get_state_array!"dC_D"(state_array[0..length]);
	}

	double[] get_wake_x_component(ref VortexFilament filament) {
		return opencopter.wake.get_wake_component!"x"(filament);
	}

	double[] get_wake_y_component(ref VortexFilament filament) {
		return opencopter.wake.get_wake_component!"y"(filament);
	}

	double[] get_wake_z_component(ref VortexFilament filament) {
		return opencopter.wake.get_wake_component!"z"(filament);
	}

	double[] get_wake_x_e_component(ref VortexFilament filament) {
		return opencopter.wake.get_wake_component!"x_e"(filament);
	}

	double[] get_wake_gamma_component(ref VortexFilament filament) {
		return opencopter.wake.get_wake_component!"gamma"(filament);
	}

	void get_wake_x_component(ref VortexFilament filament, double* d, size_t length) {
		opencopter.wake.get_wake_component!"x"(filament, d[0..length]);
	}

	void get_wake_y_component(ref VortexFilament filament, double* d, size_t length) {
		opencopter.wake.get_wake_component!"y"(filament, d[0..length]);
	}

	void get_wake_z_component(ref VortexFilament filament, double* d, size_t length) {
		opencopter.wake.get_wake_component!"z"(filament, d[0..length]);
	}

	void get_wake_x_e_component(ref VortexFilament filament, double* d, size_t length) {
		opencopter.wake.get_wake_component!"x_e"(filament, d[0..length]);
	}

	void get_wake_gamma_component(ref VortexFilament filament, double* d, size_t length) {
		opencopter.wake.get_wake_component!"gamma"(filament, d[0..length]);
	}


	struct PosChunk {
		Chunk x;
		Chunk y;
		Chunk z;
		Chunk x_e;
	}

	struct CChunk {
		Chunk c;
		alias c this;
	}

	struct Inflows {
		Inflow* inflows;
		size_t length;

		extern (D) Inflow[] get_d_array() {
			return inflows[0..length];
		}
	}

	interface Inflow {
		@nogc void update(double C_T, ref RotorInputState rotor, ref RotorState rotor_state, double advance_ratio, double axial_advance_ratio, double dt);
		@nogc void inflow_at(ref PosChunk pos, double angle_of_attack, ref CChunk output);
		extern (D) @nogc Chunk inflow_at(immutable Chunk x, immutable Chunk y, immutable Chunk z, immutable Chunk x_e, double angle_of_attack);
		@nogc double wake_skew();
	}

	class BeddosInflow : Inflow {
		static import opencopter.inflow.beddos;
		private opencopter.inflow.beddos.BeddosInflow inflow;

		this() {
			inflow = new opencopter.inflow.beddos.BeddosInflow();
		}

		@nogc void update(double C_T, ref RotorInputState rotor, ref RotorState rotor_state, double advance_ratio, double axial_advance_ratio, double dt) {
			inflow.update(C_T, rotor, rotor_state, advance_ratio, axial_advance_ratio, dt);
		}

		@nogc void inflow_at(ref PosChunk pos, double angle_of_attack, ref CChunk output) {
			output.c[] = inflow.inflow_at(pos.x, pos.y, pos.z, pos.x_e, angle_of_attack);
		}

		extern (D) @nogc Chunk inflow_at(immutable Chunk x, immutable Chunk y, immutable Chunk z, immutable Chunk x_e, double angle_of_attack) {
			return CChunk(inflow.inflow_at(x, y, z, x_e, angle_of_attack));
		}

		@nogc double wake_skew() {
			return inflow.wake_skew();
		}
	}

	alias UpdateCB = @nogc void function(double, ref RotorInputState, double, double, double);
	alias InflowAtCB = @nogc void function(ref PosChunk pos, double angle_of_attack, ref CChunk output);
	alias WakeSkewCB = @nogc double function();
	
	struct UpdateDelegate {
		UpdateCB func;
		void* data;

		alias func this;
	}

	struct InflowAtDelegate {
		InflowAtCB func;
		void* data;

		alias func this;
	}

	struct WakeSkewDelegate {
		WakeSkewCB func;
		void* data;

		alias func this;
	}

	class InflowInterface : Inflow {

		private UpdateDelegate update_cb;
		private InflowAtDelegate inflow_at_cb;
		private WakeSkewDelegate wake_skew_cb;

		this(UpdateDelegate _update_cb, InflowAtDelegate _inflow_at_cb, WakeSkewDelegate _wake_skew_cb) {
			update_cb = _update_cb;
			inflow_at_cb = _inflow_at_cb;
			wake_skew_cb = _wake_skew_cb;
		}

		@nogc void update(double C_T, ref RotorInputState rotor, ref RotorState rotor_state, double advance_ratio, double axial_advance_ratio, double dt) {
			update_cb(C_T, rotor, advance_ratio, axial_advance_ratio, dt);
		}

		@nogc void inflow_at(ref PosChunk pos, double angle_of_attack, ref CChunk output) {
			inflow_at_cb(pos, angle_of_attack, output);
		}

		extern (D) @nogc Chunk inflow_at(immutable Chunk x, immutable Chunk y, immutable Chunk z, immutable Chunk x_e, double angle_of_attack) {
			PosChunk pos = {
				x: x,
				y: y,
				z: z,
				x_e: x_e
			};
			CChunk output;
			inflow_at_cb(pos, angle_of_attack, output);
			return output.c;
		}

		@nogc double wake_skew() {
			return wake_skew_cb();
		}
	}


	static import opencopter.bladeelement;

	void compute_blade_properties(ref BladeGeometry blade, ref BladeState blade_state, ref RotorGeometry rotor, ref RotorInputState rotor_input, ref RotorState rotor_state, ref AircraftState ac_state, Inflow inflow, ref Wake wake, double time, double dt, size_t rotor_idx, size_t blade_idx) {
		import opencopter.liftmodels : steady_lift_model;
		opencopter.bladeelement.compute_blade_properties!steady_lift_model(blade, blade_state, rotor, rotor_input, rotor_state, ac_state, inflow, wake, time, dt, rotor_idx, blade_idx);
	}

	void compute_rotor_properties(ref RotorGeometry rotor, ref RotorState rotor_state, ref RotorInputState rotor_input, ref AircraftState ac_state, Inflow inflow, ref Wake wake, double C_Ti, double time, double dt, size_t rotor_idx) {
		import opencopter.liftmodels : steady_lift_model;
		opencopter.bladeelement.compute_rotor_properties!steady_lift_model(rotor, rotor_state, rotor_input, ac_state, inflow, wake, C_Ti, time, dt, rotor_idx);
	}

	void step(ref AircraftState ac_state, ref Aircraft aircraft, ref AircraftInputState ac_input_state, ref Inflows inflows, ref WakeHistory wake_history, immutable Atmosphere atmo, size_t iteration, double dt) {
		opencopter.bladeelement.step(ac_state, aircraft, ac_input_state, inflows.get_d_array, wake_history, atmo, iteration, dt);
	}

	void update_wake(ref Aircraft ac, ref AircraftState ac_state, ref AircraftInputState ac_input_state, ref WakeHistory wake_history, ref Inflows inflows, immutable Atmosphere atmo, size_t time_step, double dt) {
		opencopter.wake.update_wake(ac, ac_state, ac_input_state, wake_history, inflows.get_d_array, atmo, time_step, dt);
	}

	void compute_wake_induced_velocities(ref Wake wake, ref PosChunk pos, ref AircraftState ac_state, double angular_velocity, ref InducedVelocities output, size_t rotor_idx, size_t blade_idx) {
		output = opencopter.wake.compute_wake_induced_velocities(wake, pos.x, pos.y, pos.z, ac_state, angular_velocity, rotor_idx, blade_idx);
	}
}
