module opencopter.inflow;

import opencopter.aircraft;
import opencopter.memory;
import opencopter.math;

public import opencopter.inflow.beddos;
public import opencopter.inflow.huangpeters;
public import opencopter.inflow.simplewing;
public import opencopter.inflow.winginflow;

import std.range;
import std.algorithm;
import std.math;
import std.conv;
enum Direction {
	clockwise,
	counter_clockwise
}

interface Inflow {
	void update(ArrayContainer AC = ArrayContainer.None)(ref AircraftInputStateT!AC ac_input , Inflow[] inflows, double freestream_velocity, double advance_ratio, double axial_advance_ratio, double dt);
	void update(ArrayContainer AC = ArrayContainer.None)(AircraftInputStateT!AC* ac_input , Inflow[] inflows, double freestream_velocity, double advance_ratio, double axial_advance_ratio, double dt);
	Chunk inflow_at(immutable Chunk x, immutable Chunk y, immutable Chunk z, immutable Chunk x_e, double angle_of_attack);
	void update_wing_circulation();
	//@nogc Chunk inflow_at(immutable Chunk r, immutable double cos_azimuth, immutable double sin_azimuth);
	//@nogc double wake_skew();
}

void get_ind_vel_on_rotor(RS,RG,RIS,WIS,WS)(auto ref RS rotor_states,auto ref RG rotors, auto ref RIS rotor_inputs,auto ref WIS wing_inputs, auto ref WS wing_states, Inflow[] inflows){
	
   	size_t num_rotors = rotors.length; 
   	double[] v_z = new double[num_rotors];
   	double[] v_x = new double[num_rotors];
   	Chunk x_rotor; Chunk y_rotor; Chunk z_rotor;
	immutable Chunk chunk_of_zeros = 0.0;
   	Chunk psi_i = iota(0.0,8.0).map!(x => (x*PI/8.0).to!double).array;

   	foreach(rotor_idx, rotor; rotors){
       	double cos_aoa = rotor_inputs[rotor_idx].cos_aoa;
   		double sin_aoa = rotor_inputs[rotor_idx].sin_aoa;

       	double half_rotor_radius = rotor.blades[0].chunks[$-1].r[$-1]/2;
       	// check the relative distance of the x, y, zcoordinates for the rotors
       	x_rotor[] = half_rotor_radius*cos(psi_i)[]*cos_aoa - rotor.origin[0];
       	y_rotor[] = -(half_rotor_radius*sin(psi_i)[] + rotor.origin[1]);
       	z_rotor[] = -half_rotor_radius*cos(psi_i)[]*sin_aoa + rotor.origin[2];
       	foreach(inflow_idx,inflow; inflows){
           	if(inflow_idx != rotor_idx){
               	x_rotor[] = (x_rotor[] - rotors[inflow_idx].origin[0])*rotor_inputs[inflow_idx].cos_aoa;
               	y_rotor[] = (-y_rotor[] - rotors[inflow_idx].origin[1]);
               	z_rotor[] = (-z_rotor[] - rotors[inflow_idx].origin[2])*rotor_inputs[inflow_idx].sin_aoa;
               	auto ind_vel = inflow.inflow_at(x_rotor,y_rotor,z_rotor,chunk_of_zeros,rotor_inputs[rotor_idx].angle_of_attack);
               	if(inflow_idx < num_rotors){
                   	v_z[rotor_idx] += -ind_vel[].mean*rotor_inputs[inflow_idx].cos_aoa*rotor_inputs[inflow_idx].angular_velocity*rotor_inputs[inflow_idx].r_0[$];
                   	v_x[rotor_idx] += ind_vel[].mean*rotor_inputs[inflow_idx].sin_aoa*rotor_inputs[inflow_idx].angular_velocity*rotor_inputs[inflow_idx].r_0[$];
               	} else {
                   	v_z[rotor_idx] += ind_vel[].mean*wing_inputs[inflow_idx-num_rotors].cos_aoa;
                   	v_x[rotor_idx] += -ind_vel[].mean*wing_inputs[inflow_idx-num_rotors].sin_aoa;
               	}
               
               	rotor_states[rotor_idx].advance_ratio = (rotor_inputs[0].freestream_velocity + v_x[rotor_idx])/rotor_inputs[inflow_idx].angular_velocity*rotor_inputs[inflow_idx].r_0[$];
               	rotor_states[rotor_idx].axial_advance_ratio = v_z[rotor_idx]/rotor_inputs[inflow_idx].angular_velocity*rotor_inputs[inflow_idx].r_0[$];  
           	}
   		}
   	}
}
