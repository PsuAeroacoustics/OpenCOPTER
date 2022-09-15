module opencopter.io;

import opencopter.math;
import opencopter.wake;

import numd.linearalgebra.matrix : Vector;

import std.conv;
import std.stdio;

immutable ulong WAKE_MAGIC = 0xb4aca933bc40db78;
immutable ulong WAKE_VERSION = 1;

struct WakeFile {
	File file;
	double[] tip_buffer;
	double[] shed_buffer;

	this(string filename, size_t wake_history_length, size_t elements) {
		file = File(filename, "wb");
		tip_buffer = new double[wake_history_length];
		shed_buffer = new double[elements];
	}
}

auto start_wake_timehistory_file(AG, AIS)(string filename, double dt, auto ref AG aircraft_geometry, auto ref AIS aircraft_input_state, size_t wake_history_length, size_t shed_filaments, double aoa, size_t elements) {
	auto wake_file = WakeFile(filename, wake_history_length, elements);

	ulong[] meta_1_buff = [
		WAKE_MAGIC,
		WAKE_VERSION,
		aircraft_geometry.rotors.length.to!ulong,
		wake_history_length.to!ulong,
		shed_filaments.to!ulong
	];
	double[] meta_2_buff = [dt];
	wake_file.file.rawWrite(meta_1_buff);
	wake_file.file.rawWrite(meta_2_buff);

	foreach(r_idx; 0..aircraft_geometry.rotors.length) {
		meta_1_buff = [aircraft_geometry.rotors[r_idx].blades.length.to!ulong];
		wake_file.file.rawWrite(meta_1_buff);
		meta_2_buff = [
			aoa,
			aircraft_geometry.rotors[r_idx].origin[0],
			aircraft_geometry.rotors[r_idx].origin[1],
			aircraft_geometry.rotors[r_idx].origin[2],
			aircraft_input_state.rotor_inputs[r_idx].angular_velocity
		];
		wake_file.file.rawWrite(meta_2_buff);
	}

	return wake_file;
}

void write_wake_timestep(W)(ref WakeFile wake_file, auto ref W wake, size_t timestep) {
	
	ulong[] timestep_buff = [timestep.to!ulong];
	wake_file.file.rawWrite(timestep_buff);

	foreach(r_idx, ref rotor_wake; wake.rotor_wakes) {
		foreach(b_idx, ref shed_vortices; rotor_wake.shed_vortices) {
			foreach(ref shed_vortex; shed_vortices.shed_filaments) {
				shed_vortex.get_wake_component!"x"(wake_file.shed_buffer);
				wake_file.file.rawWrite(wake_file.shed_buffer);
				shed_vortex.get_wake_component!"y"(wake_file.shed_buffer);
				wake_file.file.rawWrite(wake_file.shed_buffer);
				shed_vortex.get_wake_component!"z"(wake_file.shed_buffer);
				wake_file.file.rawWrite(wake_file.shed_buffer);
				shed_vortex.get_wake_component!"v_z"(wake_file.shed_buffer);
				wake_file.file.rawWrite(wake_file.shed_buffer);
				shed_vortex.get_wake_component!"gamma"(wake_file.shed_buffer);
				wake_file.file.rawWrite(wake_file.shed_buffer);
			}
		}

		foreach(b_idx, ref tip_vortex; rotor_wake.tip_vortices) {
			tip_vortex.get_wake_component!"x"(wake_file.tip_buffer);
			wake_file.file.rawWrite(wake_file.tip_buffer);
			tip_vortex.get_wake_component!"y"(wake_file.tip_buffer);
			wake_file.file.rawWrite(wake_file.tip_buffer);
			tip_vortex.get_wake_component!"z"(wake_file.tip_buffer);
			wake_file.file.rawWrite(wake_file.tip_buffer);
			tip_vortex.get_wake_component!"v_z"(wake_file.tip_buffer);
			wake_file.file.rawWrite(wake_file.tip_buffer);
			tip_vortex.get_wake_component!"gamma"(wake_file.tip_buffer);
			wake_file.file.rawWrite(wake_file.tip_buffer);
			tip_vortex.get_wake_component!"r_c"(wake_file.tip_buffer);
			wake_file.file.rawWrite(wake_file.tip_buffer);
			tip_vortex.get_wake_component!"d_volume"(wake_file.tip_buffer);
			wake_file.file.rawWrite(wake_file.tip_buffer);
		}
	}
}
