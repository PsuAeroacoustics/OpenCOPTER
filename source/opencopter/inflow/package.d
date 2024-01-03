module opencopter.inflow;

import opencopter.aircraft;
import opencopter.memory;
import opencopter.math;
import opencopter.wake;

public import opencopter.inflow.beddos;
public import opencopter.inflow.huangpeters;
public import opencopter.inflow.simplewing;
public import opencopter.inflow.winginflow;

import std.stdio;
import std.range;
import std.algorithm;
import std.math;
import std.conv;
enum Direction {
	clockwise,
	counter_clockwise
}

alias Inflow = InflowT!(ArrayContainer.none);


interface InflowT(ArrayContainer AC = ArrayContainer.none) {
	void update(ref AircraftInputStateT!(AC) ac_input , ref AircraftT!(AC) aircraft,  Inflow[] inflows, double freestream_velocity, double advance_ratio, double axial_advance_ratio, double dt);
	void update(AircraftInputStateT!(AC)* ac_input , AircraftT!(AC)* aircraft, Inflow[] inflows, double freestream_velocity, double advance_ratio, double axial_advance_ratio, double dt);
	Chunk inflow_at(immutable Chunk x, immutable Chunk y, immutable Chunk z, immutable Chunk x_e, double angle_of_attack);
	void update_wing_circulation();
	void update_wing_dC_L();
	InducedVelocities compute_wing_induced_vel_on_blade(immutable Chunk x, immutable Chunk y, immutable Chunk z);
	//@nogc Chunk inflow_at(immutable Chunk r, immutable double cos_azimuth, immutable double sin_azimuth);
	//@nogc double wake_skew();
}

void get_ind_vel_on_rotor(RS,RG,RIS,WG, WIS,WS)(auto ref RS[] rotor_states, auto ref RG[] rotors, auto ref RIS[] rotor_inputs, auto ref WG[] wings, auto ref WIS[] wing_inputs, auto ref WS[] wing_states, Inflow[] inflows){
	
   	size_t num_rotors = rotors.length; 
   	double[] v_z = new double[num_rotors];
   	double[] v_x = new double[num_rotors];
   	Chunk x_rotor; Chunk y_rotor; Chunk z_rotor;
	Chunk x_rotor_op; Chunk y_rotor_op; Chunk z_rotor_op;
	Chunk x_inflow; Chunk y_inflow; Chunk z_inflow;
	immutable Chunk chunk_of_zeros = 0.0;
   	Chunk psi_i = iota(0.0,8.0).map!(x => (x*PI/8.0).to!double).array;

   	foreach(rotor_idx, rotor; rotors){
		v_z[rotor_idx] = 0.0;
		v_x[rotor_idx] = 0.0;
       	double cos_aoa = rotor_inputs[rotor_idx].cos_aoa;
   		double sin_aoa = rotor_inputs[rotor_idx].sin_aoa;

       	double half_rotor_radius = rotor.blades[0].chunks[$-1].r[$-1]/2;
       	// check the relative distance of the x, y, zcoordinates for the rotors
       	x_rotor[] = -half_rotor_radius*cos(psi_i)[];
       	y_rotor[] = half_rotor_radius*sin(psi_i)[];
       	z_rotor[] = chunk_of_zeros;
		
       	foreach(inflow_idx,inflow; inflows){
           	if(inflow_idx != rotor_idx){
               	if(inflow_idx < num_rotors){
					x_rotor_op[] = x_rotor[]*cos_aoa + z_rotor[]*sin_aoa + rotors[rotor_idx].origin[0] - rotors[inflow_idx].origin[0];
					y_rotor_op[] = -y_rotor[] + rotors[rotor_idx].origin[1] - rotors[inflow_idx].origin[1];
					z_rotor_op[] = -z_rotor[]*cos_aoa + x_rotor[]*sin_aoa + rotors[rotor_idx].origin[2] - rotors[inflow_idx].origin[2];

					x_inflow[] = x_rotor_op[]*rotor_inputs[inflow_idx].cos_aoa + z_rotor_op[]*rotor_inputs[inflow_idx].sin_aoa;
					y_inflow[] = -y_rotor_op[];
					z_inflow[] = x_rotor_op[]*rotor_inputs[inflow_idx].sin_aoa - z_rotor_op[]*rotor_inputs[inflow_idx].cos_aoa;

					immutable Chunk ind_vel = inflow.inflow_at(x_inflow,y_inflow,z_inflow,chunk_of_zeros,rotor_inputs[inflow_idx].angle_of_attack);
					//writeln("ind_vel_rotor: ", ind_vel);

                   	v_z[rotor_idx] += -ind_vel[].mean*rotor_inputs[inflow_idx].cos_aoa*abs(rotor_inputs[rotor_idx].angular_velocity*rotors[rotor_idx].radius);
                   	v_x[rotor_idx] += -ind_vel[].mean*rotor_inputs[inflow_idx].sin_aoa*abs(rotor_inputs[rotor_idx].angular_velocity*rotors[rotor_idx].radius);
					//debug writeln("v_x:", v_x);
					//debug writeln("v_z:", v_z);
               	} else {

					x_rotor_op[] = x_rotor[]*cos_aoa + z_rotor[]*sin_aoa + rotors[rotor_idx].origin[0] - wings[inflow_idx-num_rotors].origin[0];
					y_rotor_op[] = -y_rotor[] + rotors[rotor_idx].origin[1] - wings[inflow_idx-num_rotors].origin[1];
					z_rotor_op[] = -z_rotor[]*cos_aoa + x_rotor[]*sin_aoa + rotors[rotor_idx].origin[2] - wings[inflow_idx- num_rotors].origin[2];

					x_inflow[] = -x_rotor_op[]*wing_inputs[inflow_idx-num_rotors].cos_aoa - z_rotor_op[]*wing_inputs[inflow_idx-num_rotors].sin_aoa;
					y_inflow[] = -y_rotor_op[];
					z_inflow[] = -x_rotor_op[]*wing_inputs[inflow_idx-num_rotors].sin_aoa + z_rotor_op[]*wing_inputs[inflow_idx-num_rotors].cos_aoa;
					
					//writeln("x: ", wing_inputs[inflow_idx-num_rotors].cos_aoa);
					//writeln("y: ", wing_inputs[inflow_idx-num_rotors].sin_aoa);
					//writeln("z: ", z_rotor_op);

					immutable Chunk ind_vel = inflow.inflow_at(x_inflow,y_inflow,z_inflow,chunk_of_zeros,wing_inputs[inflow_idx - num_rotors].angle_of_attack);
					//writeln("ind_vel_wing: ", ind_vel);

                   	v_z[rotor_idx] += ind_vel[].mean*wing_inputs[inflow_idx-num_rotors].cos_aoa;
                   	v_x[rotor_idx] += -ind_vel[].mean*wing_inputs[inflow_idx-num_rotors].sin_aoa;
					
               	}               
               	//rotor_states[rotor_idx].advance_ratio = (rotor_inputs[0].freestream_velocity + v_x[rotor_idx])/rotor_inputs[inflow_idx].angular_velocity*rotor_inputs[inflow_idx].r_0[$-1];
               	//rotor_states[rotor_idx].axial_advance_ratio = v_z[rotor_idx]/rotor_inputs[inflow_idx].angular_velocity*rotor_inputs[inflow_idx].r_0[$-1];  
           	}
   		}
		
		rotor_states[rotor_idx].advance_ratio = (rotor_inputs[0].freestream_velocity*rotor_inputs[rotor_idx].cos_aoa + v_x[rotor_idx])/abs(rotor_inputs[rotor_idx].angular_velocity * rotors[rotor_idx].radius);
        rotor_states[rotor_idx].axial_advance_ratio = (rotor_inputs[0].freestream_velocity*rotor_inputs[rotor_idx].sin_aoa - v_z[rotor_idx])/abs(rotor_inputs[rotor_idx].angular_velocity * rotors[rotor_idx].radius);

		immutable adv_ratio_t1 = rotor_inputs[0].freestream_velocity*rotor_inputs[rotor_idx].cos_aoa/abs(rotor_inputs[rotor_idx].angular_velocity*rotors[rotor_idx].radius);
		immutable adv_ratio_t2 = -v_x[rotor_idx]/abs(rotor_inputs[rotor_idx].angular_velocity*rotors[rotor_idx].radius);

		immutable axial_adv_ratio_t1 = rotor_inputs[0].freestream_velocity*rotor_inputs[rotor_idx].sin_aoa/abs(rotor_inputs[rotor_idx].angular_velocity*rotors[rotor_idx].radius);
		immutable axial_adv_ratio_t2 = -v_z[rotor_idx]/abs(rotor_inputs[rotor_idx].angular_velocity*rotors[rotor_idx].radius);
		
		/*writeln("adv_ratio_t1 = ", adv_ratio_t1, "adv_ratio_t2 = ", adv_ratio_t2);
		writeln("axial_adv_ratio_t1 = ", axial_adv_ratio_t1, "axial_adv_ratio_t2 = ", axial_adv_ratio_t2);
		writeln("adv_ratio: ", rotor_states[rotor_idx].advance_ratio);
		writeln("axial_adv_ratio: ", rotor_states[rotor_idx].axial_advance_ratio);*/
   	}
}
