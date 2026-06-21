module opencopter.cpp;

import opencopter.aircraft;
import opencopter.atmosphere;
static import opencopter.inflow;
import opencopter.math;
import opencopter.memory;
import opencopter.wake;

//module opencopter.python;

static import opencopter.aircraft;
import opencopter.aircraft : get_wing_state_array, get_wing_state_matrix, WingPartGeometryChunk, WingPartCtrlPointChunk, WingPartStateChunk, WingPartCtrlPointStateChunk;

import opencopter.airfoilmodels;
import opencopter.atmosphere;
import opencopter.io;
import opencopter.math;
import opencopter.memory;
import opencopter.vortexlattice;

static import opencopter.vtk;

static import opencopter.aircraft.geometry;
static import opencopter.bladeelement;
static import opencopter.wake;
static import opencopter.inflow;
static import opencopter.bwi;

import std.algorithm;
import std.array;
import std.conv : to;
import std.exception : enforce;
import std.math : abs, fmod, PI, sgn;
import std.traits : isBasicType;

import numd.linearalgebra.matrix;

extern(C++) {

	import core.stdcpp.array : array;
	import core.stdcpp.vector;

	alias CppChunk = array!(double, chunk_size());

	alias BI = opencopter.inflow.BeddosInflow!(ArrayContainer.none);
	alias HP = opencopter.inflow.HuangPetersInflowT!(ArrayContainer.none);
	alias NI = opencopter.inflow.NullInflow!(ArrayContainer.none);
	alias SW = opencopter.inflow.SimpleWingT!(ArrayContainer.none);
	alias WI = opencopter.inflow.WingInflowT!(ArrayContainer.none);

	alias IV = opencopter.wake.InducedVelocities;

	alias Inflow_D = opencopter.inflow.InflowT!(ArrayContainer.none);

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

	void basic_aircraft_rotor_dynamics(AircraftInputState* ac_input, double dt) {
		foreach(r_idx, ref rotor; ac_input.rotor_inputs) {
			rotor.azimuth += rotor.angular_velocity*dt + rotor.angular_accel*dt*dt;
			auto sign = sgn(rotor.azimuth);
			// Keep the azimuth between 0 and 2*PI so we don't
			// lose fp precicion as the sim marches in time and
			// the azimuth grows unbounded.
			if(abs(rotor.azimuth) > 2.0*PI) {
				rotor.azimuth = sign * fmod(abs(rotor.azimuth), 2.0*PI);
			}
		}
	}

	double basic_single_rotor_dynamics(RotorInputState* input_state, double dt) {
		double angle = input_state.angular_velocity*dt + input_state.angular_accel*dt*dt;

		input_state.azimuth += angle;

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
	extern(C++, class) class Inflow {
		Inflow_D get_wrapped_inflow(){ assert (0);}
		void update(AircraftState* ac_state,Wake* wake, double dt) { assert(0); }
		CppChunk inflow_at(CppChunk x, CppChunk y, CppChunk z) { assert(0); }
		//Chunk inflow_at(Vector!(4, Chunk) xyz) { assert(0); }
		void update_wing_circulation(WingState* wing_state) { assert(0); }
		void update_wing_dC_L(WingState* wing_state) { assert(0); }
		IV compute_wing_induced_vel_on_blade(CppChunk x, CppChunk y, CppChunk z) { assert(0); }
		double wake_skew() { assert(0); }
		Frame* frame() { assert(0); }
		Mat4 inverse_global_frame() { assert(0); }
	}

	class HuangPeters : Inflow {
		private HP huang_peters;

		this(long Mo, long Me, RotorGeometry* rotor, RotorInputState* rotor_input, double dt) {
			huang_peters = new HP(Mo, Me, rotor, rotor_input, dt);
		}

		this(RotorGeometry* rotor, RotorInputState* rotor_input,  double dt) {
			huang_peters = new HP(4, 2, rotor, rotor_input, dt);
		}

		override void update(AircraftState* ac_state, Wake* wake, double dt) {
			huang_peters.update(*ac_state, *wake, dt);
		}

		override CppChunk inflow_at(CppChunk x, CppChunk y, CppChunk z) {
			auto xyz = Vector!(4, Chunk)(0);

			xyz[0][] = x[];
			xyz[1][] = y[];
			xyz[2][] = z[];

			return CppChunk(huang_peters.inflow_at(xyz));
		}

		// override Chunk inflow_at(immutable Vector!(4, Chunk) xyz) {
		// 	return huang_peters.inflow_at(xyz);
		// }

		override double wake_skew() {
			return huang_peters.wake_skew();
		}

		override Frame* frame() {
			return huang_peters.frame();
		}

		override Mat4 inverse_global_frame() {
			return huang_peters.inverse_global_frame;
		}

		override IV compute_wing_induced_vel_on_blade(CppChunk x, CppChunk y, CppChunk z){
			return huang_peters.compute_wing_induced_vel_on_blade(x, y, z);
		}

		override Inflow_D get_wrapped_inflow(){
			return huang_peters;
		}

	}

	class NullInflow : Inflow {
		private NI null_inflow;

		this(RotorGeometry* rotor, RotorInputState* rotor_input) {
			null_inflow = new NI(rotor, rotor_input);
		}

		override void update(AircraftState* ac_state, Wake* wake, double dt) {
			null_inflow.update(*ac_state, *wake, dt);
		}

		override CppChunk inflow_at(CppChunk x, CppChunk y, CppChunk z) {
			auto xyz = Vector!(4, Chunk)(0);

			xyz[0][] = x[];
			xyz[1][] = y[];
			xyz[2][] = z[];

			return CppChunk(null_inflow.inflow_at(xyz));
		}

		// override Chunk inflow_at(immutable Vector!(4, Chunk) xyz) {
		// 	return null_inflow.inflow_at(xyz);
		// }

		override double wake_skew() {
			return null_inflow.wake_skew();
		}

		override Frame* frame() {
			return null_inflow.frame();
		}

		override Mat4 inverse_global_frame() {
			return null_inflow.frame.inverse_global_matrix;
		}

		override IV compute_wing_induced_vel_on_blade(CppChunk x, CppChunk y, CppChunk z){
			return null_inflow.compute_wing_induced_vel_on_blade(x, y, z);
		}

		override Inflow_D get_wrapped_inflow(){
			return null_inflow;
		}

	}

	class WingInflow: Inflow{
		private WI wing_inflow;

		this(WingGeometry* _wing, WingInputState* _wing_inputs, WingLiftSurf* _wing_lift_surf) {
			wing_inflow = new WI(_wing, _wing_inputs, _wing_lift_surf);
		}

		override void update(AircraftState* ac_state, Wake* wake, double dt) {
			wing_inflow.update(*ac_state, *wake, dt);
		}

		override CppChunk inflow_at(CppChunk x, CppChunk y, CppChunk z) {
			auto xyz = Vector!(4, Chunk)(0);

			xyz[0][] = x[];
			xyz[1][] = y[];
			xyz[2][] = z[];

			return CppChunk(wing_inflow.inflow_at(xyz));
		}

		// override Chunk inflow_at(immutable Vector!(4, Chunk) xyz) {
		// 	return wing_inflow.inflow_at(xyz);
		// }

		override IV compute_wing_induced_vel_on_blade(CppChunk x, CppChunk y, CppChunk z){
			return wing_inflow.compute_wing_induced_vel_on_blade(x, y, z);
		}

		override void update_wing_circulation(WingState* wing_state){
			wing_inflow.update_wing_circulation(*wing_state);
		}

		override void update_wing_dC_L(WingState* wing_state){
			wing_inflow.update_wing_dC_L(*wing_state);
		}

		override Inflow_D get_wrapped_inflow(){
			return wing_inflow;
		}
	}

	alias WakeHistory = opencopter.wake.WakeHistoryT!(ArrayContainer.none);
	alias Wake = opencopter.wake.WakeT!(ArrayContainer.none);
	alias RotorWake = opencopter.wake.RotorWakeT!(ArrayContainer.none);
	alias VortexFilament = opencopter.wake.VortexFilamentT!(ArrayContainer.none);
	alias BWIinputs = opencopter.bwi.TipVortexInteractionT!(ArrayContainer.none);
	alias BladeVortexInteraction = opencopter.bwi.VortexInteractionT!(ArrayContainer.none);
	alias ShedVortex = opencopter.wake.ShedVortexT!(ArrayContainer.none);
	alias InteractionPerRotor = opencopter.bwi.VortexInteraction_multiRotorT!(ArrayContainer.none);

	alias VortexLattice = opencopter.vortexlattice.VortexLatticeT!(ArrayContainer.none);
	alias WingLiftSurf = opencopter.vortexlattice.WingLiftSurfT!(ArrayContainer.none);
	alias WingVortexFilament = opencopter.vortexlattice.WingVortexFilamentT!(ArrayContainer.none);
	alias WingPartLiftingSurf = opencopter.vortexlattice.WingPartLiftingSurfT!(ArrayContainer.none);

	alias AircraftTimehistory = opencopter.aircraft.AircraftTimehistoryT!(ArrayContainer.none);
	alias AircraftState = opencopter.aircraft.AircraftStateT!(ArrayContainer.none);
	alias RotorState = opencopter.aircraft.RotorStateT!(ArrayContainer.none);
	alias BladeState = opencopter.aircraft.BladeStateT!(ArrayContainer.none);

	alias Aircraft = opencopter.aircraft.AircraftT!(ArrayContainer.none);
	alias AircraftInputState = opencopter.aircraft.AircraftInputStateT!(ArrayContainer.none);
	alias RotorGeometry = opencopter.aircraft.RotorGeometryT!(ArrayContainer.none);
	alias BladeGeometry = opencopter.aircraft.BladeGeometryT!(ArrayContainer.none);
	//alias PyFrame = opencopter.aircraft.Frame!(ArrayContainer.none);
	alias Frame = opencopter.aircraft.Frame;

	alias RotorInputState = opencopter.aircraft.RotorInputStateT!(ArrayContainer.none);

	alias set_geometry_array = opencopter.aircraft.set_geometry_array;
	alias get_geometry_array = opencopter.aircraft.get_geometry_array;
	alias get_state_array = opencopter.aircraft.get_state_array;

	alias WingGeometry = opencopter.aircraft.WingGeometryT!(ArrayContainer.none);
	alias WingPartGeometry = opencopter.aircraft.WingPartGeometryT!(ArrayContainer.none);
	alias WingState = opencopter.aircraft.WingStateT!(ArrayContainer.none);
	alias WingPartState = opencopter.aircraft.WingPartStateT!(ArrayContainer.none);
	alias WingInputState = opencopter.aircraft.WingInputStateT!(ArrayContainer.none);

	alias DWingPartGeometry = opencopter.aircraft.WingPartGeometryT!(ArrayContainer.none);

	//void set_twist(ref BladeGeometry bg, core.stdcpp.vector.vector!double data) {
	void set_twist(ref BladeGeometry bg, double* data, size_t len) {
		double[] d_data = data[0..len];
		bg.set_geometry_array!"twist"(d_data);
	}

	void set_chord(ref BladeGeometry bg, double* data, size_t len) {
		double[] d_data = data[0..len];
		bg.set_geometry_array!"chord"(d_data);
	}

	void set_r(ref BladeGeometry bg, double* data, size_t len) {
		double[] d_data = data[0..len];
		bg.set_geometry_array!"r"(d_data);
	}

	void set_C_l_alpha(ref BladeGeometry bg, double* data, size_t len) {
		double[] d_data = data[0..len];
		bg.set_geometry_array!"C_l_alpha"(d_data);
	}

	void set_alpha_0(ref BladeGeometry bg, double* data, size_t len) {
		double[] d_data = data[0..len];
		bg.set_geometry_array!"alpha_0"(d_data);
	}

	void set_sweep(ref BladeGeometry bg, double* data, size_t len) {
		double[] d_data = data[0..len];
		bg.set_geometry_array!"sweep"(d_data);
	}

	void set_xi(ref BladeGeometry bg, double* data, size_t len) {
		double[] d_data = data[0..len];
		bg.set_geometry_array!"xi"(d_data);
	}

	void set_thickness(ref BladeGeometry bg, double* data, size_t len) {
		double[] d_data = data[0..len];
		bg.set_geometry_array!"thickness"(d_data);
	}

	void set_xi_p(ref BladeGeometry bg, double* data, size_t len) {
		double[] d_data = data[0..len];
		bg.set_geometry_array!"xi_p"(d_data);
	}

	void set_wing_chord(ref WingPartGeometry wg, double* data, size_t len) {
		double[] d_data = data[0..len];
		wg.set_geometry_array!"chord"(d_data);
	}

	void set_wing_twist(ref WingPartGeometry wg, double* data, size_t len) {
		double[] d_data = data[0..len];
		wg.set_geometry_array!"twist"(d_data);
	}

	void set_wing_sweep(ref WingPartGeometry wg, double* data, size_t len) {
		double[] d_data = data[0..len];
		wg.set_geometry_array!"sweep"(d_data);
	}

	void set_wing_y_span(ref WingPartGeometry wg, double* data, size_t len){
		double[] d_data = data[0..len];
		wg.set_geometry_array!"y_span"(d_data);
	}

	void fill_dC_N(ref BladeState blade, double* data, size_t len) {
		double[] d_data = data[0..len];
		blade.get_state_array!"dC_N"(d_data);
	}

	void fill_dC_c(ref BladeState blade, double* data, size_t len) {
		double[] d_data = data[0..len];
		blade.get_state_array!"dC_c"(d_data);
	}

	void fill_dC_D(ref BladeState blade, double* data, size_t len) {
		double[] d_data = data[0..len];
		blade.get_state_array!"dC_D"(d_data);
	}

	void fill_dC_Db(ref BladeState blade, double* data, size_t len) {
		double[] d_data = data[0..len];
		blade.get_state_array!"dC_Db"(d_data);
	}

	void fill_dynamic_dC_Db_profile(ref BladeState blade, double* data, size_t len) {
		double[] d_data = data[0..len];
		blade.get_state_array!"dynamic_dC_Db_profile"(d_data);
	}

	void fill_dynamic_dC_Db_induced(ref BladeState blade, double* data, size_t len) {
		double[] d_data = data[0..len];
		blade.get_state_array!"dynamic_dC_Db_induced"(d_data);
	}

	void fill_dC_Dbf(ref BladeState blade, float* data, size_t len) {
		float[] f_data = data[0..len];
		blade.get_state_array!"dC_Db"(f_data);
	}

	void fill_dC_Nf(ref BladeState blade, float* data, size_t len) {
		float[] f_data = data[0..len];
		blade.get_state_array!"dC_N"(f_data);
	}

	void fill_dC_cf(ref BladeState blade, float* data, size_t len) {
		float[] f_data = data[0..len];
		blade.get_state_array!"dC_c"(f_data);
	}

	void fill_dC_Df(ref BladeState blade, float* data, size_t len) {
		float[] f_data = data[0..len];
		blade.get_state_array!"dC_D"(f_data);
	}

	void fill_dC_Tf(ref BladeState blade, float* data, size_t len) {
		float[] f_data = data[0..len];
		blade.get_state_array!"dC_T"(f_data);
	}

	void fill_dC_Td(ref BladeState blade, double* data, size_t len) {
		double[] d_data = data[0..len];
		blade.get_state_array!"dC_T"(d_data);
	}

	void fill_dC_Qf(ref BladeState blade, float* data, size_t len) {
		float[] f_data = data[0..len];
		blade.get_state_array!"dC_Q"(f_data);
	}

	void fill_dC_Qd(ref BladeState blade, double* data, size_t len) {
		double[] d_data = data[0..len];
		blade.get_state_array!"dC_Q"(d_data);
	}

	void fill_aoad(ref BladeState blade, double* data, size_t len) {
		double[] d_data = data[0..len];
		blade.get_state_array!"aoa"(d_data);
	}

	void fill_aoa_effd(ref BladeState blade, double* data, size_t len) {
		double[] d_data = data[0..len];
		blade.get_state_array!"aoa_eff"(d_data);
	}

	void fill_aoaf(ref BladeState blade, float* data, size_t len) {
		float[] f_data = data[0..len];
		blade.get_state_array!"aoa"(f_data);
	}

	void fill_u_td(ref BladeState blade, double* data, size_t len) {
		double[] d_data = data[0..len];
		blade.get_state_array!"u_t"(d_data);
	}

	void fill_u_tf(ref BladeState blade, float* data, size_t len) {
		float[] f_data = data[0..len];
		blade.get_state_array!"u_t"(f_data);
	}

	void fill_u_pd(ref BladeState blade, double* data, size_t len) {
		double[] d_data = data[0..len];
		blade.get_state_array!"u_p"(d_data);
	}

	void fill_u_pf(ref BladeState blade, float* data, size_t len) {
		float[] f_data = data[0..len];
		blade.get_state_array!"u_p"(f_data);
	}

	void fill_dynamic_u_pd(ref BladeState blade, double* data, size_t len) {
		double[] d_data = data[0..len];
		blade.get_state_array!"dynamic_u_p"(d_data);
	}

	void fill_dynamic_u_pf(ref BladeState blade, float* data, size_t len) {
		float[] f_data = data[0..len];
		blade.get_state_array!"dynamic_u_p"(f_data);
	}

	void fill_interactionPt_wake_idx(ref BWIinputs VortexInteraction, size_t* outArray, size_t len) {
		size_t[] outArray_slice = outArray[0..len];
		opencopter.bwi.fill_indexArray!"wake_idx"(VortexInteraction, outArray_slice);
	}

	void fill_interactionPt_bladeSec_idx(ref BWIinputs VortexInteraction, size_t* outArray, size_t len) {
		size_t[] outArray_slice = outArray[0..len];
		opencopter.bwi.fill_indexArray!"bladeSec_idx"(VortexInteraction, outArray_slice);
	}

	void fill_interaction_point_r_blade(ref BWIinputs VortexInteraction, double* outArray, size_t len) {
		double[] outArray_slice = outArray[0..len];
		opencopter.bwi.fill_interaction_point_directionVec!("r_blade")(VortexInteraction, [outArray_slice]);
	}

	void fill_interaction_point_r_blade_v(ref BWIinputs VortexInteraction, double* outArray, size_t len) {
		double[] outArray_slice = outArray[0..len];
		opencopter.bwi.fill_interaction_point_directionVec!("r_blade_v")(VortexInteraction, [outArray_slice]);
	}

	void fill_interaction_point_blade_normal(ref BWIinputs VortexInteraction, double* outArray, size_t len) {
		double[] outArray_slice = outArray[0..len];
		opencopter.bwi.fill_interaction_point_directionVec!("normal")(VortexInteraction, [outArray_slice]);
	}

	void fill_interaction_point_r_vortex(ref BWIinputs VortexInteraction, double* outArray, size_t len) {
		double[] outArray_slice = outArray[0..len];
		opencopter.bwi.fill_interaction_point_directionVec!("r_vortex")(VortexInteraction, [outArray_slice]);
	}

	void fill_interaction_point_gamma_w(ref BWIinputs VortexInteraction, double* outArray, size_t len) {
		double[] outArray_slice = outArray[0..len];
		opencopter.bwi.fill_interaction_point_components!"gamma_w"(VortexInteraction, outArray_slice);
	}

	void fill_interaction_point_secLen(ref BWIinputs VortexInteraction, double* outArray, size_t len) {
		double[] outArray_slice = outArray[0..len];
		opencopter.bwi.fill_interaction_point_components!"secLen"(VortexInteraction, outArray_slice);
	}
	void fill_interaction_point_r_c(ref BWIinputs VortexInteraction, double* outArray, size_t len) {
		double[] outArray_slice = outArray[0..len];
		opencopter.bwi.fill_interaction_point_components!"r_c"(VortexInteraction, outArray_slice);
	}
	void fill_interaction_point_miss_dist(ref BWIinputs VortexInteraction, double* outArray, size_t len) {
		double[] outArray_slice = outArray[0..len];
		opencopter.bwi.fill_interaction_point_components!"miss_dist"(VortexInteraction, outArray_slice);
	}
	void fill_interaction_point_gamma_sec(ref BWIinputs VortexInteraction, double* outArray, size_t len) {
		double[] outArray_slice = outArray[0..len];
		opencopter.bwi.fill_interaction_point_components!"gamma_sec"(VortexInteraction, outArray_slice);
	}
	void fill_interaction_point_C_d(ref BWIinputs VortexInteraction, double* outArray, size_t len) {
		double[] outArray_slice = outArray[0..len];
		opencopter.bwi.fill_interaction_point_components!"C_d"(VortexInteraction, outArray_slice);
	}
	void fill_interaction_point_l(ref BWIinputs VortexInteraction, double* outArray, size_t len) {
		double[] outArray_slice = outArray[0..len];
		opencopter.bwi.fill_interaction_point_components!"l"(VortexInteraction, outArray_slice);
	}

	void fill_wake_x_component(ref VortexFilament filament, double* data, size_t len) {
		double[] data_slice = data[0..len];
		opencopter.wake.get_wake_component!"x"(filament, data_slice);
	}

	void fill_wake_y_component(ref VortexFilament filament, double* data, size_t len) {
		double[] data_slice = data[0..len];
		opencopter.wake.get_wake_component!"y"(filament, data_slice);
	}

	void fill_wake_z_component(ref VortexFilament filament, double* data, size_t len) {
		double[] data_slice = data[0..len];
		opencopter.wake.get_wake_component!"z"(filament, data_slice);
	}

	void fill_wake_gamma_component(ref VortexFilament filament, double* data, size_t len) {
		double[] data_slice = data[0..len];
		opencopter.wake.get_wake_component!"gamma"(filament, data_slice);
	}

	void fill_wake_r_c_component(ref VortexFilament filament, double* data, size_t len) {
		double[] data_slice = data[0..len];
		opencopter.wake.get_wake_component!"r_c"(filament, data_slice);
	}

	void fill_wake_v_z_component(ref VortexFilament filament, double* data, size_t len) {
		double[] data_slice = data[0..len];
		opencopter.wake.get_wake_component!"v_z"(filament, data_slice);
	}

	void fill_wake_xyz_rotor_frame(ref RotorGeometry rotor, ref VortexFilament filament, double[] x, double[] y, double[] z) {
		opencopter.wake.fill_wake_xyz_rotor_frame(rotor, filament, x, y, z);
	}

	void set_blade_pitch(ref AircraftInputState ac_input, size_t rotor_idx, size_t blade_idx, double pitch) {
		ac_input.rotor_inputs[rotor_idx].blade_pitches[blade_idx] = pitch;
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

	void step(AircraftState* ac_state, Aircraft* aircraft, AircraftInputState* ac_input_state, WakeHistory* wake_history, Atmosphere* atmo, size_t iteration, double dt, bool trackBWIevents, bool converged) {
		//void step(ArrayContainer AC = ArrayContainer.None)(ref AircraftStateT!AC ac_state, ref AircraftT!AC aircraft, ref AircraftInputStateT!AC ac_input_state, ref WakeHistoryT!AC wake_history, immutable Atmosphere atmo, size_t iteration, double dt, bool trackBWIevents, bool converged) {
		opencopter.bladeelement.step(*ac_state, *aircraft, *ac_input_state, *wake_history, *atmo, iteration, dt, trackBWIevents, converged);
	}

	void set_wing_vortex_geometry(WingLiftSurf* wing_lift_surf, WingGeometry* wing, size_t spanwise_chunks, size_t chordwise_nodes){
		opencopter.vortexlattice.set_wing_vortex_geometry(wing_lift_surf, wing, spanwise_chunks, chordwise_nodes);
	}

	void set_wing_ctrl_pt_geometry(WingGeometry* wing, size_t spanwise_nodes, size_t chordwise_nodes, double camber){
		opencopter.aircraft.set_wing_ctrl_pt_geometry(wing, spanwise_nodes, chordwise_nodes, camber);
	}

	opencopter.vtk.VtkRotor build_base_vtu_rotor(RotorGeometry* rotor) {
		return opencopter.vtk.build_base_vtu_rotor(rotor);
	}

	void write_rotor_vtu(string base_filename, size_t iteration, size_t rotor_idx, opencopter.vtk.VtkRotor rotor, RotorState* rotor_state, RotorGeometry* rotor_geom) {
		opencopter.vtk.write_rotor_vtu(base_filename, iteration, rotor_idx, rotor, rotor_state, rotor_geom);
	}

	void write_rotors_vtu(string base_filename, size_t iteration, opencopter.vtk.VtkRotor[] rotors, AircraftState* ac_state, Aircraft* aircraft) {
		foreach(r_idx; 0..aircraft.rotors.length()) {
			opencopter.vtk.write_rotor_vtu(base_filename, iteration, r_idx, rotors[r_idx], ac_state.rotor_states[r_idx], aircraft.rotors[r_idx]);
		}
	}

	opencopter.vtk.VtkWing build_base_vtu_wing(WingGeometry* wing){
		return opencopter.vtk.build_base_vtu_wing(wing);
	}

	void write_wing_vtu(string base_filename, size_t iteration, size_t wing_idx, opencopter.vtk.VtkWing wing, WingState* wing_state, WingGeometry* wing_geom){
		opencopter.vtk.write_wing_vtu(base_filename, iteration, wing_idx, wing, wing_state, wing_geom);
	}

	opencopter.vtk.VtkWake build_base_vtu_wake(Wake* wake) {
		return opencopter.vtk.build_base_vtu_wake(wake);
	}

	void write_wake_vtu(string base_filename, size_t iteration, opencopter.vtk.VtkWake vtk_wake, Wake* wake) {
		opencopter.vtk.write_wake_vtu(base_filename, iteration, vtk_wake, wake);
	}

	opencopter.vtk.VtkWingWake build_base_vtu_wing_wake(WingGeometry* wing, WingLiftSurf* wing_lift_surf){
		return opencopter.vtk.build_base_vtu_wing_wake(wing, wing_lift_surf);
	}

	void write_wing_wake_vtu(string base_filename, size_t iteration, size_t wing_idx,  opencopter.vtk.VtkWingWake wing_wake, WingGeometry* wing_geom, WingLiftSurf* wing_lift_surf, WingInputState wing_input){
		opencopter.vtk.write_wing_wake_vtu(base_filename, iteration, wing_idx, wing_wake, wing_geom, wing_lift_surf, wing_input);
	}

	void write_wake_field_vtu(string filename, AircraftState ac_state, Wake wake, Vec3 delta, Vec3 starts, size_t num_x, size_t num_y, size_t num_z) {
		opencopter.vtk.write_wake_field_vtu(filename, ac_state, wake, delta, starts, num_x, num_y, num_z);
	}

	alias write_inflow_vtu = opencopter.vtk.write_inflow_vtu!(Inflow, Array!RotorGeometry);

	IV compute_wake_induced_velocities(ref Wake wake, immutable Chunk x, immutable Chunk y, immutable Chunk z, ref AircraftState ac_state, size_t rotor_idx, bool single_rotor = false) {
		return opencopter.wake.compute_wake_induced_velocities(wake, x, y, z, ac_state, rotor_idx, 0.to!size_t, single_rotor);
	}

	AircraftState* CreateAircraftState(size_t num_rotors, size_t[] num_blades, size_t num_elements, size_t num_wings, size_t[] num_wing_parts, size_t num_span_nodes, size_t num_chord_nodes, Aircraft* ac, Inflow[] rotor_inflows, Inflow[] wing_inflows, double[] direction) {
		auto oc_rotor_inflows = rotor_inflows.map!(a => a.get_wrapped_inflow()).array;
		auto oc_wing_inflows = wing_inflows.map!(a => a.get_wrapped_inflow()).array;
		return new AircraftState(num_rotors, num_blades, num_elements, num_wings, num_wing_parts, num_span_nodes, num_chord_nodes, *ac, oc_rotor_inflows, oc_wing_inflows, direction);
	}

	void compute_blade_vectors(ref BladeGeometry blade) {
		opencopter.aircraft.compute_blade_vectors(blade);
	}

	Mat4 Mat4_identity() {
		return Mat4.identity();
	}

	Mat3 Mat3_identity() {
		return Mat3.identity();
	}

	/*Mat4 get_local_matrix(ref Frame comp_frame){
		return comp_frame.local_matrix;
	}

	Mat4 get_global_matrix(ref Frame comp_frame){
		return comp_frame.global_matrix;
	}*/

	alias WingLoc = opencopter.aircraft.geometry.Location;

	struct Location{
		private WingLoc loc;

		string toString() const {
			return loc;
		}
	}

	Location location_right_wing(){
		return Location(opencopter.aircraft.Location.right);
	}

	Location location_left_wing(){
		return Location(opencopter.aircraft.Location.left);
	}

	WingPartGeometry build_wing_part_geometry(size_t span_elements, size_t chordwise_nodes, Vec3 wing_root_origin, double average_chord, double wing_root_chord, double wing_tip_chord, double le_sweep_angle, double te_sweep_angle, double wing_span, Location Pyloc){
		auto loc = Pyloc.loc;
		auto wing_part = WingPartGeometry(span_elements, chordwise_nodes, wing_root_origin, average_chord, wing_root_chord, wing_tip_chord, le_sweep_angle, te_sweep_angle, wing_span, loc);

		return wing_part;
	}

}
