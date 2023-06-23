module opencopter.wake;

import opencopter.aircraft;
import opencopter.atmosphere;
import opencopter.config;
import opencopter.inflow;
import opencopter.math;
import opencopter.memory;

import std.algorithm;
import std.math;
import std.stdio;
import std.traits;
import std.typecons;

extern (C++) struct FilamentChunk {
	Chunk x;
	Chunk y;
	Chunk z;
	Chunk x_e;
	Chunk gamma; // Vortex circulation
	Chunk phi; // Time marched wake strain integral
	Chunk r_c; // vortex core radius.
	Chunk r_0; // vortex core radius.
	Chunk dx;
	Chunk dy;
	Chunk dz;
	Chunk l_0;
	Chunk v_z;
	Chunk volume;
	Chunk d_volume;
}

extern (C++) struct VortexFilamentT(ArrayContainer AC) {

	mixin ArrayDeclMixin!(AC, FilamentChunk, "chunks");

	size_t length;

	this(size_t wake_history) {

		assert(wake_history % chunk_size == 0);

		immutable num_chunks = wake_history/chunk_size;

		mixin(array_ctor_mixin!(AC, "FilamentChunk", "chunks", "num_chunks"));

		foreach(ref chunk; chunks) {
			chunk.gamma[] = 0;
			chunk.phi[] = 0;
			chunk.r_c = 0.00001;
		}

		length = 0;
	}

	this(FilamentChunk* slice, size_t slice_length) {
		mixin(array_ctor_mixin_slice!(AC, "FilamentChunk", "chunks", "cast(FilamentChunk[])slice[0..slice_length]"));

		foreach(ref chunk; chunks) {
			chunk.gamma[] = 0;
			chunk.phi[] = 0;
			chunk.r_c = 0.000001;
		}

		length = 0;
	}
}

extern(C++) struct ShedVortexT(ArrayContainer AC) {
	mixin ArrayDeclMixin!(AC, VortexFilamentT!(AC), "shed_filaments");

	private size_t shed_wake_length;
	this(size_t radial_elements, size_t _shed_wake_length) {
		shed_wake_length = _shed_wake_length;
		immutable radial_chunks = radial_elements/chunk_size;

		auto shed_filament_mem = new FilamentChunk[radial_chunks*shed_wake_length];
		mixin(array_ctor_mixin!(AC, "VortexFilamentT!(AC)", "shed_filaments", "shed_wake_length"));

		foreach(i, ref shed_filament; shed_filaments) {
			shed_filament = VortexFilamentT!AC(&shed_filament_mem[i*radial_chunks], radial_chunks);
		}
	}
}

alias VortexFilament = VortexFilamentT!(ArrayContainer.none);

alias RotorWake = RotorWakeT!(ArrayContainer.none);

extern(C++) struct RotorWakeT(ArrayContainer AC) {
	mixin ArrayDeclMixin!(AC, VortexFilamentT!(AC), "tip_vortices");
	mixin ArrayDeclMixin!(AC, ShedVortexT!(AC), "shed_vortices");

	size_t current_shed_idx;
	alias tip_vortices this;

	this(size_t num_blades) {
		mixin(array_ctor_mixin!(AC, "VortexFilamentT!(AC)", "tip_vortices", "num_blades"));
		mixin(array_ctor_mixin!(AC, "ShedVortexT!(AC)", "shed_vortices", "num_blades"));
		current_shed_idx = 0;
	}
}

template is_wake(A) {
	enum bool is_wake = {
		static if(isPointer!(A)) {
			return isInstanceOf!(WakeT, PointerTarget!A);
		} else {
			return isInstanceOf!(WakeT, A);
		}
	}();
}

alias Wake = WakeT!(ArrayContainer.none);

/++
 +	Holds the wake data for the specified number of rotors. This represents a single timestep
 +/
struct WakeT(ArrayContainer AC) {

	mixin ArrayDeclMixin!(AC, RotorWakeT!(AC), "rotor_wakes");

	this(size_t num_rotors, size_t num_blades, size_t wake_history, size_t radial_elements, size_t shed_history) {

		immutable actual_wake_history = wake_history%chunk_size == 0 ? wake_history : wake_history + (chunk_size - wake_history%chunk_size);
		immutable actual_radial_elements = radial_elements%chunk_size == 0 ? radial_elements : radial_elements + (chunk_size - radial_elements%chunk_size);

		mixin(array_ctor_mixin!(AC, "RotorWakeT!(AC)", "rotor_wakes", "num_rotors"));
		foreach(ref rotor_wake; rotor_wakes) {
			rotor_wake = RotorWakeT!AC(num_blades);
			foreach(ref filament; rotor_wake.tip_vortices) {
				filament = VortexFilamentT!AC(actual_wake_history);
			}

			foreach(ref shed_vortex; rotor_wake.shed_vortices) {
				shed_vortex = ShedVortexT!AC(actual_radial_elements, shed_history);
			}
		}
	}

	this(size_t num_rotors, size_t num_blades, size_t[] wake_history, size_t radial_elements, size_t shed_history) {

		//immutable actual_wake_history = wake_history%chunk_size == 0 ? wake_history : wake_history + (chunk_size - wake_history%chunk_size);

		auto actual_wake_history = wake_history.map!(w => w%chunk_size == 0 ? w : w + (chunk_size - w%chunk_size)).array;

		immutable actual_radial_elements = radial_elements%chunk_size == 0 ? radial_elements : radial_elements + (chunk_size - radial_elements%chunk_size);

		mixin(array_ctor_mixin!(AC, "RotorWakeT!(AC)", "rotor_wakes", "num_rotors"));
		foreach(r_idx, ref rotor_wake; rotor_wakes) {

			rotor_wake = RotorWakeT!AC(num_blades);
			foreach(ref filament; rotor_wake.tip_vortices) {
				filament = VortexFilamentT!AC(actual_wake_history[r_idx]);
			}

			foreach(ref shed_vortex; rotor_wake.shed_vortices) {
				shed_vortex = ShedVortexT!AC(actual_radial_elements, shed_history);
			}
		}
	}
}

template is_wake_history(A) {
	enum bool is_wake_history = {
		static if(isPointer!(A)) {
			return isInstanceOf!(WakeHistoryT, PointerTarget!A);
		} else {
			return isInstanceOf!(WakeHistoryT, A);
		}
	}();
}

alias WakeHistory = WakeHistoryT!(ArrayContainer.none);

/++
 +	Contains the history of the wake at each timestep. The number of timesteps
 +	to store is configurable, but at least 2 are required.
 +/
struct WakeHistoryT(ArrayContainer AC) {

	mixin ArrayDeclMixin!(AC, WakeT!(AC), "history");
	alias history this;

	this(size_t num_rotors, size_t num_blades, size_t wake_history, size_t time_history, size_t radial_elements, size_t shed_history = 20) {

		mixin(array_ctor_mixin!(AC, "WakeT!(AC)", "history", "time_history"));

		foreach(ref wake; history) {
			wake = WakeT!AC(num_rotors, num_blades, wake_history, radial_elements, shed_history);
		}
	}

	this(size_t num_rotors, size_t num_blades, size_t[] wake_history, size_t time_history, size_t radial_elements, size_t shed_history = 20) {

		mixin(array_ctor_mixin!(AC, "WakeT!(AC)", "history", "time_history"));

		foreach(ref wake; history) {
			wake = WakeT!AC(num_rotors, num_blades, wake_history, radial_elements, shed_history);
		}
	}

	void push_back_wake() {
		for(size_t i = history.length - 1; i > 0; i--) {
			foreach(rotor_idx, ref rotor_wake; history[i - 1].rotor_wakes) {
				//history[i].rotor_wakes[rotor_idx].shed_vortices[] = rotor_wake.shed_vortices[];
				foreach(blade_idx, ref filament; rotor_wake.tip_vortices) {
					//history[i].rotor_wakes[rotor_idx].shed_vortices[blade_idx].shed_filaments[] = rotor_wake.shed_vortices[blade_idx].shed_filaments[];
					foreach(f_idx, ref shed_filament; rotor_wake.shed_vortices[blade_idx].shed_filaments) {
						foreach(chunk_idx, ref chunk; shed_filament.chunks) {
							//auto rotor = &history[i].rotor_wakes[rotor_idx];
							//auto shed_wake = &history[i].rotor_wakes[rotor_idx].shed_vortices[blade_idx];
							auto fil = history[i].rotor_wakes[rotor_idx].shed_vortices[blade_idx].shed_filaments[f_idx];
							//import std.stdio : writeln;
							//writeln(chunk_idx, ": ", fil.chunks.length, ", ", filament.chunks.length);
							fil.chunks[chunk_idx] = chunk;
						}
					}
					foreach(chunk_idx, ref chunk; filament.chunks) {
						history[i].rotor_wakes[rotor_idx].tip_vortices[blade_idx].chunks[chunk_idx] = chunk;
					}
				}
			}
		}
	}
}

double[] get_wake_component(string value, VF)(auto ref VF filament) {
	
	immutable elements = filament.chunks.length * chunk_size;

	double[] d = new double[elements];

	foreach(c_idx, ref chunk; filament.chunks) {
		immutable out_start_idx = c_idx*chunk_size;

		immutable remaining = elements - out_start_idx;
		
		immutable out_end_idx = remaining > chunk_size ? (c_idx + 1)*chunk_size : out_start_idx + remaining;
		immutable in_end_idx = remaining > chunk_size ? chunk_size : remaining;

		mixin("d[out_start_idx..out_end_idx] = chunk."~value~"[0..in_end_idx];");
	}
	return d;
}

void get_wake_component(string value, VF)(auto ref VF filament, double[] d) {
	foreach(c_idx, ref chunk; filament.chunks) {
		immutable out_start_idx = c_idx*chunk_size;

		immutable remaining = d.length - out_start_idx;
		
		immutable out_end_idx = remaining > chunk_size ? (c_idx + 1)*chunk_size : out_start_idx + remaining;
		immutable in_end_idx = remaining > chunk_size ? chunk_size : remaining;

		mixin("d[out_start_idx..out_end_idx] = chunk."~value~"[0..in_end_idx];");
	}
}

import std.algorithm : map, fold, maxElement, minElement;
import std.array;
import std.range;

extern(C++) struct InducedVelocities {
	Chunk v_x;
	Chunk v_y;
	Chunk v_z;
}

import std.math : PI;
immutable Chunk one_over_four_pi = 1.0/(4.0*PI);

import std.parallelism;

InducedVelocities compute_filament_induced_velocities(FC)(auto ref FC chunks, immutable Chunk x, immutable Chunk y, immutable Chunk z, size_t chunk_offset) {
	InducedVelocities ret;
	Chunk[chunk_size] v_x, v_y, v_z;

	foreach(n_idx; 0..Chunk.length) {
		v_x[n_idx][] = 0;
		v_y[n_idx][] = 0;
		v_z[n_idx][] = 0;
	}

	foreach(i_c_idx, ref chunk_i; chunks[chunk_offset..$]) {

		i_c_idx += chunk_offset;

		Chunk x_b, y_b, z_b;
		Chunk x_a, y_a, z_a;

		Chunk x_bm1, y_bm1, z_bm1;
		Chunk x_am1, y_am1, z_am1;

		Chunk x_bp1, y_bp1, z_bp1;
		Chunk x_ap1, y_ap1, z_ap1;

		x_am1[1..$] = chunk_i.x[0..$-1];
		y_am1[1..$] = chunk_i.y[0..$-1];
		z_am1[1..$] = chunk_i.z[0..$-1];

		x_ap1[0..$-1] = chunk_i.x[1..$];
		y_ap1[0..$-1] = chunk_i.y[1..$];
		z_ap1[0..$-1] = chunk_i.z[1..$];

		x_bp1[0..$-2] = chunk_i.x[2..$];
		y_bp1[0..$-2] = chunk_i.y[2..$];
		z_bp1[0..$-2] = chunk_i.z[2..$];

		if(i_c_idx == 0) {
			x_am1[0] = x_am1[1];
			z_am1[0] = y_am1[1];
			y_am1[0] = z_am1[1];
		} else {
			x_am1[0] = chunks[i_c_idx - 1].x[$-1];
			y_am1[0] = chunks[i_c_idx - 1].y[$-1];
			z_am1[0] = chunks[i_c_idx - 1].z[$-1];
		}

		if(i_c_idx == chunks.length - 1) {
			x_ap1[$-1] = x_ap1[$-2];
			y_ap1[$-1] = y_ap1[$-2];
			z_ap1[$-1] = z_ap1[$-2];

			x_bp1[$-2..$] = x_bp1[$-2];
			y_bp1[$-2..$] = y_bp1[$-2];
			z_bp1[$-2..$] = z_bp1[$-2];
		} else {
			x_ap1[$-1] = chunks[i_c_idx + 1].x[0];
			y_ap1[$-1] = chunks[i_c_idx + 1].y[0];
			z_ap1[$-1] = chunks[i_c_idx + 1].z[0];

			x_bp1[$-2..$] = chunks[i_c_idx + 1].x[0..2];
			y_bp1[$-2..$] = chunks[i_c_idx + 1].y[0..2];
			z_bp1[$-2..$] = chunks[i_c_idx + 1].z[0..2];
		}

		x_a[] = chunk_i.x[];
		y_a[] = chunk_i.y[];
		z_a[] = chunk_i.z[];

		x_b[0..$-1] = chunk_i.x[1..$];
		y_b[0..$-1] = chunk_i.y[1..$];
		z_b[0..$-1] = chunk_i.z[1..$];

		if(i_c_idx != chunks.length - 1) {
			x_b[$-1] = chunks[i_c_idx + 1].x[0];
			y_b[$-1] = chunks[i_c_idx + 1].y[0];
			z_b[$-1] = chunks[i_c_idx + 1].z[0];
		} else {
			chunk_i.gamma[$-1] = 0;
			x_b[$-1] = 0;
			y_b[$-1] = 0;
			z_b[$-1] = 0;
		}

		bool all_nan = chunk_i.x[].map!(a => a.isNaN).fold!((res, a) => res && a)(true);

		if(all_nan) {
			break;
		}

		Chunk gamma = chunk_i.gamma[];
		Chunk r_c = chunk_i.r_c[];
		// put some dummy data in values that
		// haven't been populated so we don't
		// NaN infect the below calcs
		foreach(n_idx; 0..Chunk.length) {
			if(x_b[n_idx].isNaN || x_a[n_idx].isNaN) {
				r_c[n_idx] = 0.001;
				gamma[n_idx] = 0;
				x_a[n_idx] = 100;
				y_a[n_idx] = 100;
				z_a[n_idx] = 100;

				x_am1[n_idx] = 99;
				y_am1[n_idx] = 99;
				z_am1[n_idx] = 99;

				x_b[n_idx] = 101;
				y_b[n_idx] = 101;
				z_b[n_idx] = 101;
			}

			if(x_bp1[n_idx].isNaN || x_ap1[n_idx].isNaN) {
				x_ap1[n_idx] = 101;
				y_ap1[n_idx] = 101;
				z_ap1[n_idx] = 101;

				x_bp1[n_idx] = 102;
				y_bp1[n_idx] = 102;
				z_bp1[n_idx] = 102;
			}
		}

		x_bm1[] = x_a[];
		y_bm1[] = y_a[];
		z_bm1[] = z_a[];

		import std.stdio : writeln;
		//writeln(x_am1);
		//writeln(x_bm1);
		//writeln(x_ap1);
		//writeln(x_bp1);
		//writeln;
		// Calc vortex segment mid point
		immutable Chunk x_m = 0.5*(x_a[] + x_b[]);
		immutable Chunk y_m = 0.5*(y_a[] + y_b[]);
		immutable Chunk z_m = 0.5*(z_a[] + z_b[]);

		immutable Chunk x_mm1 = 0.5*(x_am1[] + x_bm1[]);
		immutable Chunk y_mm1 = 0.5*(y_am1[] + y_bm1[]);
		immutable Chunk z_mm1 = 0.5*(z_am1[] + z_bm1[]);

		immutable Chunk x_mp1 = 0.5*(x_ap1[] + x_bp1[]);
		immutable Chunk y_mp1 = 0.5*(y_ap1[] + y_bp1[]);
		immutable Chunk z_mp1 = 0.5*(z_ap1[] + z_bp1[]);

		// Calc vortex segment length vec
		//immutable Chunk dx = chunk_i.dx[];//(x_b[] - x_a[]);
		//immutable Chunk dy = chunk_i.dy[];//(y_b[] - y_a[]);
		//immutable Chunk dz = chunk_i.dz[];//(z_b[] - z_a[]);
		immutable Chunk dx = (x_b[] - x_a[]);
		immutable Chunk dy = (y_b[] - y_a[]);
		immutable Chunk dz = (z_b[] - z_a[]);

		immutable Chunk dx_m1 = (x_bm1[] - x_am1[]);
		immutable Chunk dy_m1 = (y_bm1[] - y_am1[]);
		immutable Chunk dz_m1 = (z_bm1[] - z_am1[]);

		immutable Chunk dx_p1 = (x_bp1[] - x_ap1[]);
		immutable Chunk dy_p1 = (y_bp1[] - y_ap1[]);
		immutable Chunk dz_p1 = (z_bp1[] - z_ap1[]);

		// Perform Biot-savart calc for each passed in point.
		foreach(n_idx; 0..Chunk.length) {
			immutable Chunk x_p = x[n_idx];
			immutable Chunk y_p = y[n_idx];
			immutable Chunk z_p = z[n_idx];

			immutable Chunk r_x1 = x_p[] - x_mm1[];
			immutable Chunk r_y1 = y_p[] - y_mm1[];
			immutable Chunk r_z1 = z_p[] - z_mm1[];

			immutable Chunk r_x2 = x_p[] - x_m[];
			immutable Chunk r_y2 = y_p[] - y_m[];
			immutable Chunk r_z2 = z_p[] - z_m[];

			immutable Chunk r_x3 = x_p[] - x_mp1[];
			immutable Chunk r_y3 = y_p[] - y_mp1[];
			immutable Chunk r_z3 = z_p[] - z_mp1[];

			immutable Chunk x_cross1 = r_y1[]*dz_m1[] - r_z1[]*dy_m1[];
			immutable Chunk y_cross1 = r_z1[]*dx_m1[] - r_x1[]*dz_m1[];
			immutable Chunk z_cross1 = r_x1[]*dy_m1[] - r_y1[]*dx_m1[];

			immutable Chunk x_cross2 = r_y2[]*dz[] - r_z2[]*dy[];
			immutable Chunk y_cross2 = r_z2[]*dx[] - r_x2[]*dz[];
			immutable Chunk z_cross2 = r_x2[]*dy[] - r_y2[]*dx[];

			immutable Chunk x_cross3 = r_y3[]*dz_p1[] - r_z3[]*dy_p1[];
			immutable Chunk y_cross3 = r_z3[]*dx_p1[] - r_x3[]*dz_p1[];
			immutable Chunk z_cross3 = r_x3[]*dy_p1[] - r_y3[]*dx_p1[];

			immutable Chunk x_cross = (1.0/3.0)*(x_cross1[] + x_cross2[] + x_cross3[]);
			immutable Chunk y_cross = (1.0/3.0)*(y_cross1[] + y_cross2[] + y_cross3[]);
			immutable Chunk z_cross = (1.0/3.0)*(z_cross1[] + z_cross2[] + z_cross3[]);

			immutable Chunk miss_dist1 = r_x1[]*r_x1[] + r_y1[]*r_y1[] + r_z1[]*r_z1[];
			immutable Chunk miss_dist2 = r_x2[]*r_x2[] + r_y2[]*r_y2[] + r_z2[]*r_z2[];
			immutable Chunk miss_dist3 = r_x3[]*r_x3[] + r_y3[]*r_y3[] + r_z3[]*r_z3[];

			immutable Chunk miss_dist = (1.0/3.0)*(miss_dist1[] + miss_dist2[] + miss_dist3[]);

			immutable Chunk unorm = 1.0/(miss_dist[] + r_c[]*r_c[]);

			import std.stdio : writeln;

			immutable Chunk norm = sqrt(unorm)[]*unorm[];

			immutable Chunk circulation = norm[]*gamma[]*one_over_four_pi[];

			// Well WTF. Keeping in this random 1.0 multiplication factor
			// makes things faster over just remove the multiply.............
			// I uh, i don't even know
			immutable Chunk random_mult_by_one = 1.0;
			
			immutable Chunk tmp_v_x = random_mult_by_one[]*circulation[]*x_cross[];
			immutable Chunk tmp_v_y = random_mult_by_one[]*circulation[]*y_cross[];
			immutable Chunk tmp_v_z = random_mult_by_one[]*circulation[]*z_cross[];

			v_x[n_idx][] += tmp_v_x[];
			v_y[n_idx][] += tmp_v_y[];
			v_z[n_idx][] += tmp_v_z[];
		}
	}

	foreach(n_idx; 0..Chunk.length) {
		ret.v_x[n_idx] = v_x[n_idx].sum;
		ret.v_y[n_idx] = v_y[n_idx].sum;
		ret.v_z[n_idx] = v_z[n_idx].sum;
	}

	return ret;
}

InducedVelocities compute_wake_induced_velocities(W, AS)(auto ref W wake, immutable Chunk x, immutable Chunk y, immutable Chunk z, auto ref AS ac_state, double angular_velocity, size_t rotor_idx, size_t blade_idx, bool single_rotor = false)
	if(is_wake!W && is_aircraft_state!AS)
{
	static import std.math;
	
	InducedVelocities ret;
	InducedVelocities ret_shed;
	ret.v_x[] = 0;
	ret.v_y[] = 0;
	ret.v_z[] = 0;

	ret_shed.v_x[] = 0;
	ret_shed.v_y[] = 0;
	ret_shed.v_z[] = 0;

	import std.math : isNaN;
	import std.stdio : writeln;

	immutable chunk_offset = 0;

	foreach(i_rotor_idx, ref i_rotor; ac_state.rotor_states) {
		foreach(i_blade_idx; 0..ac_state.rotor_states[i_rotor_idx].blade_states.length) {

			if((i_rotor_idx != rotor_idx) || ((i_rotor_idx == rotor_idx) && (i_blade_idx != blade_idx))) {
				auto ind_vel = compute_filament_induced_velocities(ac_state.rotor_states[i_rotor_idx].blade_states[i_blade_idx].chunks, x, y, z, 0);
				ret_shed.v_x[] += ind_vel.v_x[];
				ret_shed.v_y[] += ind_vel.v_y[];
				ret_shed.v_z[] += ind_vel.v_z[];
			}

			foreach(fil_idx, ref shed_filament; wake.rotor_wakes[i_rotor_idx].shed_vortices[i_blade_idx].shed_filaments[1..$]) {
				auto ind_vel = compute_filament_induced_velocities(shed_filament.chunks, x, y, z, 0);
				ret_shed.v_x[] += ind_vel.v_x[];
				ret_shed.v_y[] += ind_vel.v_y[];
				ret_shed.v_z[] += ind_vel.v_z[];
			}
		}
	}

	if(single_rotor) {
		foreach(i_blade_idx; 0..ac_state.rotor_states[rotor_idx].blade_states.length) {
			auto ind_vel = compute_filament_induced_velocities(wake.rotor_wakes[rotor_idx].tip_vortices[i_blade_idx].chunks, x, y, z, chunk_offset);
			ret.v_x[] += ind_vel.v_x[];
			ret.v_y[] += ind_vel.v_y[];
			ret.v_z[] += ind_vel.v_z[];
		}
	} else {
		foreach(i_rotor_idx, ref i_rotor; ac_state.rotor_states) {
			foreach(i_blade_idx; 0..i_rotor.blade_states.length) {
				auto ind_vel = compute_filament_induced_velocities(wake.rotor_wakes[i_rotor_idx].tip_vortices[i_blade_idx].chunks, x, y, z, chunk_offset);
				/+if(i_rotor_idx != rotor_idx) {
					ret.v_x[] -= ind_vel.v_x[];
					ret.v_y[] -= ind_vel.v_y[];
					ret.v_z[] -= ind_vel.v_z[];
				} else {+/
					ret.v_x[] += ind_vel.v_x[];
					ret.v_y[] += ind_vel.v_y[];
					ret.v_z[] += ind_vel.v_z[];
				//}
				
			}
		}
	}
	
	ret.v_x[] += ret_shed.v_x[];
	ret.v_y[] += ret_shed.v_y[];
	ret.v_z[] += ret_shed.v_z[];

	//ret.v_x[] /= (2.*PI);//(0.01*angular_velocity);
	//ret.v_y[] /= (2.*PI);//(0.01*angular_velocity);
	//ret.v_z[] /= (2.*PI);//(0.01*angular_velocity);

	//ret.v_x[] *= (2.0*PI);
	//ret.v_y[] *= (2.0*PI);
	//ret.v_z[] *= (2.0*PI);

	//ret.v_x[] /= (angular_velocity);
	//ret.v_y[] /= (angular_velocity);
	//ret.v_z[] /= (angular_velocity);

	//ret.v_x[] /= (10);
	//ret.v_y[] /= (10);
	//ret.v_z[] /= (10);
	
	return ret;
}

immutable double alpha_l = 1.25643;
immutable double a1 =  6.5e-5;

import core.thread;
import core.sync.barrier;

void update_wake(I, ArrayContainer AC = ArrayContainer.None)(ref AircraftT!AC ac, ref AircraftStateT!AC ac_state, ref AircraftInputStateT!AC ac_input_state, ref WakeHistoryT!AC wake_history, I[] inflows, immutable Atmosphere atmo, size_t time_step, double dt) {

	static import std.math;
	import std.stdio : writeln;
	import std.traits : ForeachType;
	
	wake_history.push_back_wake;

	//double vortex_offset = 35*PI/180.0;
	double vortex_offset = 0;
	
	foreach(rotor_idx, ref rotor; ac_state.rotor_states) {

		immutable omega_sgn = std.math.sgn(ac_input_state.rotor_inputs[rotor_idx].angular_velocity);

		immutable double inboard_factor = ac.rotors[rotor_idx].blades[0].chunks[$-1].r[$-1];

		immutable Chunk free_stream_x = -ac_input_state.rotor_inputs[rotor_idx].freestream_velocity/(ac.rotors[rotor_idx].radius);
		
		immutable sin_aoa =  ac_input_state.rotor_inputs[rotor_idx].sin_aoa;
		immutable cos_aoa =  ac_input_state.rotor_inputs[rotor_idx].cos_aoa;

		foreach(blade_idx, ref blade; rotor.blade_states) {

			immutable cos_beta = std.math.cos(ac_input_state.rotor_inputs[rotor_idx].blade_flapping[blade_idx]);
			immutable sin_beta = std.math.sin(ac_input_state.rotor_inputs[rotor_idx].blade_flapping[blade_idx]);

			immutable x_f = cos_beta*inboard_factor;
			immutable z_f = -sin_beta*inboard_factor;

			immutable x_tpp = ac.rotors[rotor_idx].origin[0] + x_f*std.math.cos(blade.azimuth - vortex_offset) + omega_sgn*ac.rotors[rotor_idx].blades[blade_idx].chunks[$-1].xi[$-1]*std.math.sin(blade.azimuth);
			immutable y_tpp = ac.rotors[rotor_idx].origin[1] + inboard_factor*std.math.sin(blade.azimuth - vortex_offset) - omega_sgn*ac.rotors[rotor_idx].blades[blade_idx].chunks[$-1].xi[$-1]*std.math.cos(ac_input_state.rotor_inputs[rotor_idx].blade_pitches[blade_idx])*std.math.cos(blade.azimuth);
			immutable z_tpp = ac.rotors[rotor_idx].origin[2] + z_f - ac.rotors[rotor_idx].blades[blade_idx].chunks[$-1].xi[$-1]*std.math.sin(ac_input_state.rotor_inputs[rotor_idx].blade_pitches[blade_idx]);
			immutable x = x_tpp*cos_aoa + z_tpp*sin_aoa;
			immutable y = y_tpp;
			immutable z = -x_tpp*sin_aoa + z_tpp*cos_aoa;

			static if(AC == ArrayContainer.none) {
				auto current_tip_filament = &wake_history.history[0].rotor_wakes[rotor_idx].tip_vortices[blade_idx];
				auto blade_geo = &ac.rotors[rotor_idx].blades[blade_idx];
				auto blade_state = &ac_state.rotor_states[rotor_idx].blade_states[blade_idx];
				auto current_shed_filament = &wake_history.history[0].rotor_wakes[rotor_idx].shed_vortices[blade_idx].shed_filaments[0];
			} else {
				auto current_tip_filament = wake_history.history[0].rotor_wakes[rotor_idx].tip_vortices[blade_idx];
				auto blade_geo = ac.rotors[rotor_idx].blades[blade_idx];
				auto blade_state = ac_state.rotor_states[rotor_idx].blade_states[blade_idx];
				auto current_shed_filament = wake_history.history[0].rotor_wakes[rotor_idx].shed_vortices[blade_idx].shed_filaments[0];
			}

			immutable omega = std.math.abs(ac_input_state.rotor_inputs[rotor_idx].angular_velocity);
			immutable rotation_dir = -1.0*std.math.sgn(ac_input_state.rotor_inputs[rotor_idx].angular_velocity);

			foreach(c_idx, ref chunk; current_shed_filament.chunks) {
				//immutable Chunk x_fc = cos_beta*blade_geo.chunks[c_idx].r[];
				//immutable Chunk z_fc = -sin_beta*blade_geo.chunks[c_idx].r[];

				//immutable Chunk x_tppc = ac.rotors[rotor_idx].origin[0] + x_fc[]*std.math.cos(blade.azimuth);
				//immutable Chunk y_tppc = ac.rotors[rotor_idx].origin[1] + blade_geo.chunks[c_idx].r[]*std.math.sin(blade.azimuth);
				//immutable Chunk z_tppc = ac.rotors[rotor_idx].origin[2] + z_fc[];
				//chunk.x[] = x_tppc[]*cos_aoa + z_tppc[]*sin_aoa;
				//chunk.y[] = y_tppc[];
				//chunk.z[] = -x_tppc[]*sin_aoa + z_tppc[]*cos_aoa;

				chunk.x[] = blade_state.chunks[c_idx].x[];
				chunk.y[] = blade_state.chunks[c_idx].y[];
				chunk.z[] = blade_state.chunks[c_idx].z[];

				chunk.gamma[] = /+rotation_dir*+/blade_state.chunks[c_idx].d_gamma[];
				chunk.l_0[] = 0;
				chunk.r_0[] = 0.0001;//ac_input_state.rotor_inputs[rotor_idx].r_0[blade_idx];
				chunk.r_c[] = 0.0001;
				chunk.x_e[] = 0;
				chunk.phi[] = 0;

				chunk.dx[0..$-1] = chunk.x[1..$] - chunk.x[0..$-1];
				chunk.dy[0..$-1] = chunk.y[1..$] - chunk.y[0..$-1];
				chunk.dz[0..$-1] = chunk.z[1..$] - chunk.z[0..$-1];

				chunk.dx[$ - 1] = 0;
				chunk.dy[$ - 1] = 0;
				chunk.dz[$ - 1] = 0;

				if(c_idx > 0) {
					current_shed_filament.chunks[c_idx - 1].dx[$-1] = chunk.x[0] - current_shed_filament.chunks[c_idx - 1].x[$-1];
					current_shed_filament.chunks[c_idx - 1].dy[$-1] = chunk.y[0] - current_shed_filament.chunks[c_idx - 1].y[$-1];
					current_shed_filament.chunks[c_idx - 1].dz[$-1] = chunk.z[0] - current_shed_filament.chunks[c_idx - 1].z[$-1];
				}
			}

			foreach(shed_idx, ref shed_filament; wake_history.history[0].rotor_wakes[rotor_idx].shed_vortices[blade_idx].shed_filaments[1..$]) {

				auto last_shed_filament = wake_history.history[1].rotor_wakes[rotor_idx].shed_vortices[blade_idx].shed_filaments.ptr_at_idx(shed_idx);
				/+static if(AC == ArrayContainer.none) {
					auto last_shed_filament = &wake_history.history[1].rotor_wakes[rotor_idx].shed_vortices[blade_idx].shed_filaments[shed_idx];
				} else {
					auto last_shed_filament = wake_history.history[1].rotor_wakes[rotor_idx].shed_vortices[blade_idx].shed_filaments[shed_idx];
				}+/

				foreach(c_idx, ref chunk; last_shed_filament.chunks) {

					Chunk lambda_i_z = 0;
					Chunk lambda_i_x = 0;
					foreach(i_rotor_idx, ref inflow; inflows) {
						immutable i_sin_aoa =  ac_input_state.rotor_inputs[i_rotor_idx].sin_aoa;
						immutable i_cos_aoa =  ac_input_state.rotor_inputs[i_rotor_idx].cos_aoa;

						immutable Chunk x_rel = chunk.x[];// - ac.rotors[rotor_idx].origin[0];
						immutable Chunk y_rel = -(chunk.y[] - ac.rotors[i_rotor_idx].origin[1]);
						immutable Chunk z_rel = chunk.z[];// - ac.rotors[rotor_idx].origin[2];

						immutable Chunk x_tppv = x_rel[]*i_cos_aoa - z_rel[]*i_sin_aoa - ac.rotors[i_rotor_idx].origin[0];
						immutable Chunk z_tppv = -x_rel[]*i_sin_aoa - z_rel[]*i_cos_aoa + ac.rotors[i_rotor_idx].origin[2];

						immutable i_omega = std.math.abs(ac_input_state.rotor_inputs[i_rotor_idx].angular_velocity);
						//immutable i_omega = abs(free_stream_x);
						immutable Chunk lambda_i = i_omega*inflows[i_rotor_idx].inflow_at(x_tppv, y_rel, z_tppv, chunk.x_e, ac_input_state.rotor_inputs[i_rotor_idx].angle_of_attack)[];
						lambda_i_z[] += -lambda_i[]*i_cos_aoa;
						lambda_i_x[] += -lambda_i[]*i_sin_aoa;
					}

					immutable Chunk v_x = free_stream_x[] + lambda_i_x[];
					immutable Chunk v_y = 0;
					immutable Chunk v_z = lambda_i_z[];
					
					shed_filament.chunks[c_idx].gamma[] = chunk.gamma[];
					shed_filament.chunks[c_idx].l_0[] = chunk.l_0[];
					shed_filament.chunks[c_idx].phi[] = chunk.phi[];
					shed_filament.chunks[c_idx].r_0[] = chunk.r_0[];
					shed_filament.chunks[c_idx].r_c[] = chunk.r_c[];
					shed_filament.chunks[c_idx].x[] = chunk.x[] + dt*v_x[];
					shed_filament.chunks[c_idx].y[] = chunk.y[] + dt*v_y[];
					shed_filament.chunks[c_idx].z[] = chunk.z[] + dt*v_z[];
					shed_filament.chunks[c_idx].v_z[] = v_z[];

					shed_filament.chunks[c_idx].dx[0..$-1] = shed_filament.chunks[c_idx].x[1..$] - shed_filament.chunks[c_idx].x[0..$-1];
					shed_filament.chunks[c_idx].dy[0..$-1] = shed_filament.chunks[c_idx].y[1..$] - shed_filament.chunks[c_idx].y[0..$-1];
					shed_filament.chunks[c_idx].dz[0..$-1] = shed_filament.chunks[c_idx].z[1..$] - shed_filament.chunks[c_idx].z[0..$-1];

					if(c_idx > 0) {
						shed_filament.chunks[c_idx - 1].dx[$-1] = shed_filament.chunks[c_idx].x[0] - shed_filament.chunks[c_idx - 1].x[$-1];
						shed_filament.chunks[c_idx - 1].dy[$-1] = shed_filament.chunks[c_idx].y[0] - shed_filament.chunks[c_idx - 1].y[$-1];
						shed_filament.chunks[c_idx - 1].dz[$-1] = shed_filament.chunks[c_idx].z[0] - shed_filament.chunks[c_idx - 1].z[$-1];
					}
				}

				shed_filament.chunks[$-1].r_c[$-1] = 0;
				shed_filament.chunks[$-1].phi[$-1] = 0;
			}

			/+immutable max_gamma = ac_state
				.rotor_states[rotor_idx]
				.blade_states[blade_idx]
				.chunks[]
				.map!(a => a.gamma[].map!((g) { /+writeln(g);+/ return tuple(g, std.math.abs(g));} ))
				.map!(c => c.maxElement!"a[1]")
				.maxElement!"a[1]";
				//.map!(c => maxElement(abs(c.gamma)[]))
				//.maxElement;
			+/

			double max_gamma = 0;
			foreach(ref c; ac_state
				.rotor_states[rotor_idx]
				.blade_states[blade_idx]
				.chunks) {
				foreach(ref g; c.gamma) {
					if(abs(g) > abs(max_gamma)) {
						max_gamma = g;
					}
				}
			}
			//writeln(max_gamma);
			current_tip_filament.chunks[0].y[0] = y;
			current_tip_filament.chunks[0].z[0] = z;
			current_tip_filament.chunks[0].x[0] = x;
			current_tip_filament.chunks[0].gamma[0] = max_gamma;
			current_tip_filament.chunks[0].l_0[0] = 0;
			current_tip_filament.chunks[0].r_0[0] = ac_input_state.rotor_inputs[rotor_idx].r_0[blade_idx];
			current_tip_filament.chunks[0].r_c[0] = current_tip_filament.chunks[0].r_0[0];
			current_tip_filament.chunks[0].x_e[0] = std.math.sqrt(1.0 - y_tpp*y_tpp);
			current_tip_filament.chunks[0].phi[0] = 0;

			foreach(c_idx, ref chunk; wake_history.history[1].rotor_wakes[rotor_idx].tip_vortices[blade_idx].chunks) {

				Chunk lambda_i_z = 0;
				Chunk lambda_i_x = 0;

				foreach(i_rotor_idx, ref inflow; inflows) {

					immutable i_aoa = ac_input_state.rotor_inputs[i_rotor_idx].angle_of_attack;
					immutable i_sin_aoa =  ac_input_state.rotor_inputs[i_rotor_idx].sin_aoa;
					immutable i_cos_aoa =  ac_input_state.rotor_inputs[i_rotor_idx].cos_aoa;

					immutable Chunk x_rel = chunk.x[];// - ac.rotors[i_rotor_idx].origin[0];
					immutable Chunk y_rel = -(chunk.y[] - ac.rotors[i_rotor_idx].origin[1]);
					immutable Chunk z_rel = chunk.z[];// - ac.rotors[i_rotor_idx].origin[2];

					immutable Chunk x_tppc = (x_rel[]*i_cos_aoa - z_rel[]*i_sin_aoa)[] - ac.rotors[i_rotor_idx].origin[0];
					immutable Chunk z_tppc = (-x_rel[]*i_sin_aoa - z_rel[]*i_cos_aoa)[] + ac.rotors[i_rotor_idx].origin[2];
					
					immutable i_omega = std.math.abs(ac_input_state.rotor_inputs[i_rotor_idx].angular_velocity);
					//immutable i_omega = abs(free_stream_x);
					immutable Chunk lambda_i = i_omega*inflow.inflow_at(x_tppc, y_rel, z_tppc, chunk.x_e, i_aoa)[];
					
					lambda_i_z[] += -lambda_i[]*i_cos_aoa;
					lambda_i_x[] += -lambda_i[]*i_sin_aoa;
				}

				immutable Chunk v_x = free_stream_x[] + lambda_i_x[];
				immutable Chunk v_y = 0;
				immutable Chunk v_z = lambda_i_z[];
				
				immutable Chunk tmp_x = chunk.x[] + dt*v_x[];
				immutable Chunk tmp_y = chunk.y[] + dt*v_y[];
				immutable Chunk tmp_z = chunk.z[] + dt*v_z[];

				current_tip_filament.chunks[c_idx].x[1..$] = tmp_x[0..$-1];
				current_tip_filament.chunks[c_idx].y[1..$] = tmp_y[0..$-1];
				current_tip_filament.chunks[c_idx].z[1..$] = tmp_z[0..$-1];
				current_tip_filament.chunks[c_idx].v_z[1..$] = v_z[0..$-1];
				current_tip_filament.chunks[c_idx].x_e[1..$] = chunk.x_e[0..$-1];
				current_tip_filament.chunks[c_idx].gamma[1..$] = chunk.gamma[0..$-1];
				current_tip_filament.chunks[c_idx].l_0[1..$] = chunk.l_0[0..$-1];

				current_tip_filament.chunks[c_idx].dx[0..$-1] = current_tip_filament.chunks[c_idx].x[0..$-1] - current_tip_filament.chunks[c_idx].x[1..$];
				current_tip_filament.chunks[c_idx].dy[0..$-1] = current_tip_filament.chunks[c_idx].y[0..$-1] - current_tip_filament.chunks[c_idx].y[1..$];
				current_tip_filament.chunks[c_idx].dz[0..$-1] = current_tip_filament.chunks[c_idx].z[0..$-1] - current_tip_filament.chunks[c_idx].z[1..$];

				current_tip_filament.chunks[c_idx].dx[$-1] = current_tip_filament.chunks[c_idx].x[$-1] - tmp_x[$-1];
				current_tip_filament.chunks[c_idx].dy[$-1] = current_tip_filament.chunks[c_idx].y[$-1] - tmp_y[$-1];
				current_tip_filament.chunks[c_idx].dz[$-1] = current_tip_filament.chunks[c_idx].z[$-1] - tmp_z[$-1];

				immutable Chunk l2 = current_tip_filament.chunks[c_idx].dx[]*current_tip_filament.chunks[c_idx].dx[] + current_tip_filament.chunks[c_idx].dy[]*current_tip_filament.chunks[c_idx].dy[] + current_tip_filament.chunks[c_idx].dz[]*current_tip_filament.chunks[c_idx].dz[];
				immutable Chunk l = sqrt(l2);
				immutable Chunk Rev = current_tip_filament.chunks[c_idx].gamma[]*omega*ac.rotors[rotor_idx].radius*ac.rotors[rotor_idx].blades[blade_idx].average_chord/atmo.kinematic_viscosity;


				immutable Chunk tmp_phi = chunk.phi[] + dt*(current_tip_filament.chunks[c_idx].l_0[]/l[]);

				immutable Chunk r_0 = current_tip_filament.chunks[c_idx].r_0;
				immutable Chunk delta = 1 + a1*abs(Rev)[];

				Chunk r_c_rad = r_0[]*r_0[];
				r_c_rad[] += 4.0*alpha_l*atmo.kinematic_viscosity*delta[]*tmp_phi[];
				current_tip_filament.chunks[c_idx].r_c = sqrt(r_c_rad)[]/ac.rotors[rotor_idx].radius;
				current_tip_filament.chunks[c_idx].volume = current_tip_filament.chunks[c_idx].r_c[]*current_tip_filament.chunks[c_idx].r_c[]*l[];
				current_tip_filament.chunks[c_idx].d_volume[1..$] = chunk.volume[0..$-1] - current_tip_filament.chunks[c_idx].volume[1..$];
				current_tip_filament.chunks[c_idx].phi[1..$] = tmp_phi[0..$-1];
				current_tip_filament.chunks[c_idx].r_0[1..$] = chunk.r_0[0..$-1];

				if(c_idx > 0) {
					current_tip_filament.chunks[c_idx].d_volume[0] = wake_history.history[1].rotor_wakes[rotor_idx].tip_vortices[blade_idx].chunks[c_idx - 1].volume[$ - 1] - current_tip_filament.chunks[c_idx].volume[0];
				}

				// Only update index 0 if this is not the last chunk. Index 0 of the first chunk is
				// the most recent and is updated later.
				if(c_idx != wake_history.history[1].rotor_wakes[rotor_idx].tip_vortices[blade_idx].chunks.length - 1) {
					current_tip_filament.chunks[c_idx + 1].x[0] = tmp_x[$-1];
					current_tip_filament.chunks[c_idx + 1].y[0] = tmp_y[$-1];
					current_tip_filament.chunks[c_idx + 1].z[0] = tmp_z[$-1];
					current_tip_filament.chunks[c_idx + 1].v_z[0] = v_z[$-1];
					current_tip_filament.chunks[c_idx + 1].x_e[0] = chunk.x_e[$-1];
					current_tip_filament.chunks[c_idx + 1].gamma[0] = chunk.gamma[$-1];
					current_tip_filament.chunks[c_idx + 1].l_0[0] = chunk.l_0[$-1];
					current_tip_filament.chunks[c_idx + 1].r_0[0] = chunk.r_0[$-1];
					current_tip_filament.chunks[c_idx + 1].phi[0] = tmp_phi[$-1];
				} else {
					current_tip_filament.chunks[c_idx].phi[$-1] = current_tip_filament.chunks[c_idx].phi[$-2];
					current_tip_filament.chunks[c_idx].r_c[$-1] = current_tip_filament.chunks[c_idx].r_c[$-2];
				}
			}

			double dx = current_tip_filament.chunks[0].dx[0];
			double dy = current_tip_filament.chunks[0].dy[0];
			double dz = current_tip_filament.chunks[0].dz[0];
			double l = sqrt(dx*dx + dy*dy + dz*dz);

			current_tip_filament.chunks[0].l_0[0] = l;
		}

		if(wake_history.history[0].rotor_wakes[rotor_idx].current_shed_idx < wake_history.history[0].rotor_wakes[rotor_idx].shed_vortices[0].shed_wake_length) {
			wake_history.history[0].rotor_wakes[rotor_idx].current_shed_idx++;
		}
	}
}
