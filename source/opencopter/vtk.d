module opencopter.vtk;

import opencopter.aircraft;
import opencopter.config;
import opencopter.math;
import opencopter.memory;
import opencopter.wake;

import numd.linearalgebra.matrix : Matrix, Vector;

import std.conv;
import std.math;
import std.stdio;

double[2][] naca0012 = [
	[1.000000,  0.001260],
	[0.998459,  0.001476],
	[0.993844,  0.002120],
	[0.986185,  0.003182],
	[0.975528,  0.004642],
	[0.961940,  0.006478],
	[0.945503,  0.008658],
	[0.926320,  0.011149],
	[0.904508,  0.013914],
	[0.880203,  0.016914],
	[0.853553,  0.020107],
	[0.824724,  0.023452],
	[0.793893,  0.026905],
	[0.761249,  0.030423],
	[0.726995,  0.033962],
	[0.691342,  0.037476],
	[0.654508,  0.040917],
	[0.616723,  0.044237],
	[0.578217,  0.047383],
	[0.539230,  0.050302],
	[0.500000,  0.052940],
	[0.460770,  0.055241],
	[0.421783,  0.057148],
	[0.383277,  0.058609],
	[0.345492,  0.059575],
	[0.308658,  0.060000],
	[0.273005,  0.059848],
	[0.238751,  0.059092],
	[0.206107,  0.057714],
	[0.175276,  0.055709],
	[0.146447,  0.053083],
	[0.119797,  0.049854],
	[0.095492,  0.046049],
	[0.073680,  0.041705],
	[0.054497,  0.036867],
	[0.038060,  0.031580],
	[0.024472,  0.025893],
	[0.013815,  0.019854],
	[0.006156,  0.013503],
	[0.001541,  0.006877],
	[0.000000,  0.000000],
	[0.001541, -0.006877],
	[0.006156, -0.013503],
	[0.013815, -0.019854],
	[0.024472, -0.025893],
	[0.038060, -0.031580],
	[0.054497, -0.036867],
	[0.073680, -0.041705],
	[0.095492, -0.046049],
	[0.119797, -0.049854],
	[0.146447, -0.053083],
	[0.175276, -0.055709],
	[0.206107, -0.057714],
	[0.238751, -0.059092],
	[0.273005, -0.059848],
	[0.308658, -0.060000],
	[0.345492, -0.059575],
	[0.383277, -0.058609],
	[0.421783, -0.057148],
	[0.460770, -0.055241],
	[0.500000, -0.052940],
	[0.539230, -0.050302],
	[0.578217, -0.047383],
	[0.616723, -0.044237],
	[0.654508, -0.040917],
	[0.691342, -0.037476],
	[0.726995, -0.033962],
	[0.761249, -0.030423],
	[0.793893, -0.026905],
	[0.824724, -0.023452],
	[0.853553, -0.020107],
	[0.880203, -0.016914],
	[0.904508, -0.013914],
	[0.926320, -0.011149],
	[0.945503, -0.008658],
	[0.961940, -0.006478],
	[0.975528, -0.004642],
	[0.986185, -0.003182],
	[0.993844, -0.002120],
	[0.998459, -0.001476],
	[1.000000, -0.001260]
];

alias PosVector = Vector!(3, size_t);

version(Have_vtkd) {
	import vtkd.vtkd;
}

class VtkRotor {

	this() {
	}

	version(Have_vtkd) {
		private vtkUnstructuredGrid* grid;
		private vtkPoints* points;
		private double[] azimuth_offsets;
		private Vec3 origin;
		private Vec3[vtkIdType][] base_points;
		private vtkIdType[][] r_to_point_map;
		private vtkDoubleArray* loads;
		private vtkDoubleArray* aoa;
		private vtkDoubleArray* aoa_eff;
		private vtkDoubleArray* u_p;
		private vtkDoubleArray* shed_u_p;
		private vtkDoubleArray* u_t;
		private vtkDoubleArray* inflow_angle;
		private vtkDoubleArray* gamma;
		private vtkDoubleArray* d_gamma;
		private vtkDoubleArray* theta;
		private vtkDoubleArray* dC_T_dot;
		private vtkDoubleArray* dC_L_dot;
		private vtkDoubleArray* dC_N;
		private vtkDoubleArray* dC_c;
		private vtkDoubleArray* dC_T;
		private vtkDoubleArray* dC_Q;
		private vtkDoubleArray* dC_D;
		private vtkDoubleArray* af_norm;
		private vtkDoubleArray* blade_local_vel;
		private vtkDoubleArray* projected_vel;

		private this(size_t num_blades) {
			azimuth_offsets = new double[num_blades];
			base_points = new Vec3[vtkIdType][num_blades];
			loads = vtkDoubleArray.New;
			dC_T = vtkDoubleArray.New;
			dC_Q = vtkDoubleArray.New;
			dC_D = vtkDoubleArray.New;
			dC_T_dot = vtkDoubleArray.New;
			dC_L_dot = vtkDoubleArray.New;
			dC_N = vtkDoubleArray.New;
			dC_c = vtkDoubleArray.New;
			u_p = vtkDoubleArray.New;
			shed_u_p = vtkDoubleArray.New;
			u_t = vtkDoubleArray.New;
			inflow_angle = vtkDoubleArray.New;
			aoa = vtkDoubleArray.New;
			aoa_eff = vtkDoubleArray.New;
			gamma = vtkDoubleArray.New;
			d_gamma = vtkDoubleArray.New;
			theta = vtkDoubleArray.New;
			af_norm = vtkDoubleArray.New;
			blade_local_vel = vtkDoubleArray.New;
			projected_vel = vtkDoubleArray.New;
			grid = vtkUnstructuredGrid.New;
			points = vtkPoints.New;
		}
	}
}

class VtkWing {

	this() {
	}

	version(Have_vtkd) {
		private vtkUnstructuredGrid* grid;
		private vtkPoints* points;
		private Vec3 origin;
		private Vec3[vtkIdType][] base_points;
		private vtkIdType[][][] base_point_ids;
		//private vtkDoubleArray* loads;
		private vtkDoubleArray* aoa;
		private vtkDoubleArray* u_p;
		private vtkDoubleArray* u_t;
		//private vtkDoubleArray* gamma;
		private vtkDoubleArray* dC_T;

		private this(size_t num_wing_parts) {
			base_points = new Vec3[vtkIdType][num_wing_parts];
			//loads = vtkDoubleArray.New;
			u_p = vtkDoubleArray.New;
			u_t = vtkDoubleArray.New;
			aoa = vtkDoubleArray.New;
			//gamma = vtkDoubleArray.New;
			dC_T = vtkDoubleArray.New;
			grid = vtkUnstructuredGrid.New;
			points = vtkPoints.New;
		}
	}
}

void write_wing_vtu(WS, WIS)(string base_filename, size_t iteration, size_t wing_idx, ref VtkWing wing, auto ref WS wing_state, auto ref WIS wing_input){

	version(Have_vtkd) {
		writeln("writing wing vtu");
		immutable spanwise_elements = wing_state.wing_part_states[0].chunks.length*chunk_size;
		immutable chordwise_elements = wing_state.wing_part_states[0].ctrl_chunks.length*chunk_size / spanwise_elements;

		auto wing_writer = vtkXMLUnstructuredGridWriter.New;

		auto aoa_rotation =
			Mat3(
				-std.math.cos(wing_input.angle_of_attack /++ flip_angle+/), 0, -std.math.sin(wing_input.angle_of_attack/+ + flip_angle+/),
				0, -1, 0,
				-std.math.sin(wing_input.angle_of_attack /++ flip_angle+/), 0, std.math.cos(wing_input.angle_of_attack /++ flip_angle+/)
			);
		
		wing.grid.SetPoints(wing.points);

		foreach(wp_idx, wing_part; wing_state.wing_part_states) {
			foreach(pi; wing.base_points[wp_idx].byKeyValue){
				vtkIdType id = pi.key;
				Vec3 point = pi.value;
				auto origin = wing.origin;

				auto aoa_rot = aoa_rotation*point ;
				auto final_p = aoa_rot + origin;
				writeln("final_p = ", final_p);

				wing.points.SetPoint(id, final_p[0], final_p[1], final_p[2]);
			}	
		}

		foreach(wp_idx, wp_state; wing_state.wing_part_states){
			foreach(ch_idx, chunk ; wp_state.ctrl_chunks){
				foreach(c_idx; 0..chunk_size){
					immutable span_chunk_length = wp_state.chunks.length;
					immutable span_chunk_idx = (ch_idx%span_chunk_length);
					immutable span_idx = span_chunk_idx*chunk_size + c_idx;
					immutable chord_idx = (ch_idx - ch_idx%span_chunk_length)/span_chunk_length;

					auto id = wing.base_point_ids[wp_idx][span_idx][chord_idx];
					wing.aoa.SetTuple1(id, wp_state.ctrl_chunks[ch_idx].ctrl_pt_aoa[c_idx]);
					wing.u_p.SetTuple1(id, wp_state.ctrl_chunks[ch_idx].ctrl_pt_up[c_idx]);
					wing.u_t.SetTuple1(id, wp_state.ctrl_chunks[ch_idx].ctrl_pt_ut[c_idx]);
					wing.dC_T.SetTuple1(id, wp_state.chunks[span_chunk_idx].dC_L[c_idx]);
				}
			}
		}

		import std.string : toStringz;

		auto filename = base_filename~"_"~wing_idx.to!string~"_"~iteration.to!string~".vtu";

		wing_writer.SetFileName(filename.toStringz);
		wing_writer.SetInputData(wing.grid);
		wing_writer.Write;
	}
}

void write_rotor_vtu(RS, RG)(string base_filename, size_t iteration, size_t rotor_idx, ref VtkRotor rotor, auto ref RS rotor_state, auto ref RG rotor_geom) {
	version(Have_vtkd) {
		immutable elements = rotor_state.blade_states[0].chunks.length*chunk_size;

		auto rotor_writer = vtkXMLUnstructuredGridWriter.New;

		rotor.grid.SetPoints(rotor.points);

		foreach(blade_idx, blade; rotor_state.blade_states) {
			foreach(pi; rotor.base_points[blade_idx].byKeyValue) {
				vtkIdType id = pi.key;
				auto point = Vec4(pi.value[0], pi.value[1], pi.value[2], 1.0/rotor_geom.radius)*rotor_geom.radius;

				auto final_p = rotor_geom.blades[blade_idx].frame.global_matrix * point;
				
				rotor.points.SetPoint(id, final_p[0], final_p[1], final_p[2]);
			}

			foreach(radial_idx, loop; rotor.r_to_point_map[blade_idx*elements..elements*(blade_idx + 1)]) {

				auto chunk_idx = radial_idx/chunk_size;
				auto inner_idx = radial_idx%chunk_size;

				auto af_norm = rotor_geom.blades[blade_idx].frame.global_matrix*rotor_geom.blades[blade_idx].chunks[chunk_idx].af_norm;

				auto blade_local_vel = blade.chunks[chunk_idx].blade_local_vel;
				auto projected_vel = blade.chunks[chunk_idx].projected_vel;

				foreach(l_idx, id; loop) {
					rotor.loads.SetTuple1(id, blade.chunks[chunk_idx].dC_L[inner_idx]);
					rotor.dC_T.SetTuple1(id, blade.chunks[chunk_idx].dC_T[inner_idx]);
					rotor.dC_Q.SetTuple1(id, blade.chunks[chunk_idx].dC_Q[inner_idx]);
					rotor.dC_D.SetTuple1(id, blade.chunks[chunk_idx].dC_D[inner_idx]);
					rotor.dC_L_dot.SetTuple1(id, blade.chunks[chunk_idx].dC_L_dot[inner_idx]);
					rotor.dC_T_dot.SetTuple1(id, blade.chunks[chunk_idx].dC_T_dot[inner_idx]);
					rotor.dC_N.SetTuple1(id, blade.chunks[chunk_idx].dC_N[inner_idx]);
					rotor.dC_c.SetTuple1(id, blade.chunks[chunk_idx].dC_c[inner_idx]);
					rotor.u_p.SetTuple1(id, blade.chunks[chunk_idx].u_p[inner_idx]);
					rotor.shed_u_p.SetTuple1(id, blade.chunks[chunk_idx].shed_u_p[inner_idx]);
					rotor.u_t.SetTuple1(id, blade.chunks[chunk_idx].u_t[inner_idx]);
					rotor.aoa.SetTuple1(id, blade.chunks[chunk_idx].aoa[inner_idx]*(180.0/PI));
					rotor.aoa_eff.SetTuple1(id, blade.chunks[chunk_idx].aoa_eff[inner_idx]*(180.0/PI));
					rotor.inflow_angle.SetTuple1(id, blade.chunks[chunk_idx].inflow_angle[inner_idx]*(180.0/PI));
					rotor.gamma.SetTuple1(id, blade.chunks[chunk_idx].gamma[inner_idx]);
					rotor.d_gamma.SetTuple1(id, blade.chunks[chunk_idx].d_gamma[inner_idx]);
					rotor.theta.SetTuple1(id, blade.chunks[chunk_idx].theta[inner_idx]*(180.0/PI));
					rotor.af_norm.SetTuple3(id, af_norm[0][inner_idx], af_norm[1][inner_idx], af_norm[2][inner_idx]);

					rotor.blade_local_vel.SetTuple3(id, blade_local_vel[0][inner_idx], blade_local_vel[1][inner_idx], blade_local_vel[2][inner_idx]);
					rotor.projected_vel.SetTuple3(id, projected_vel[0][inner_idx], projected_vel[1][inner_idx], projected_vel[2][inner_idx]);

				}
			}
		}

		import std.string : toStringz;

		auto filename = base_filename~"_"~rotor_idx.to!string~"_"~iteration.to!string~".vtu";

		rotor_writer.SetFileName(filename.toStringz);
		rotor_writer.SetInputData(rotor.grid);
		rotor_writer.Write;
	}
}

alias Mat3 = Matrix!(3, 3, double);

VtkWing build_base_vtu_wing(WG)(auto ref WG wing_geo) {
	version(Have_vtkd) {
		//immutable num_wing_parts = wing_geo.wing_parts.length;
		immutable span_elements = wing_geo.wing_parts[0].chunks.length*chunk_size;
		immutable chord_elements = wing_geo.wing_parts[0].ctrl_chunks.length*chunk_size/span_elements;

		//std.stdio.writeln("build_vtu_num_wing_parts = ",wing_geo.wing_parts.length );
		auto vtk_wing = new VtkWing(wing_geo.wing_parts.length);

		vtk_wing.origin = wing_geo.origin;

		/*vtk_wing.loads.SetNumberOfComponents(1);
		vtk_wing.loads.SetNumberOfTuples(span_elements*chord_elements*wing_geo.wing_parts.length);
		vtk_wing.loads.SetName("l");*/

		vtk_wing.aoa.SetNumberOfComponents(1);
		vtk_wing.aoa.SetNumberOfTuples(span_elements*chord_elements*wing_geo.wing_parts.length);
		vtk_wing.aoa.SetName("aoa");

		vtk_wing.u_p.SetNumberOfComponents(1);
		vtk_wing.u_p.SetNumberOfTuples(span_elements*chord_elements*wing_geo.wing_parts.length);
		vtk_wing.u_p.SetName("u_p");

		vtk_wing.u_t.SetNumberOfComponents(1);
		vtk_wing.u_t.SetNumberOfTuples(span_elements*chord_elements*wing_geo.wing_parts.length);
		vtk_wing.u_t.SetName("u_t");

		/*vtk_wing.gamma.SetNumberOfComponents(1);
		vtk_wing.gamma.SetNumberOfTuples(span_elements*chord_elements*wing_geo.wing_parts.length);
		vtk_wing.gamma.SetName("gamma");*/

		vtk_wing.dC_T.SetNumberOfComponents(1);
		vtk_wing.dC_T.SetNumberOfTuples(span_elements*chord_elements*wing_geo.wing_parts.length);
		vtk_wing.dC_T.SetName("dC_T");

		vtk_wing.base_point_ids = new vtkIdType[][][](wing_geo.wing_parts.length,span_elements,chord_elements);

		foreach(wp_idx, ref wp_geom; wing_geo.wing_parts){
			size_t num_span_elements = wp_geom.chunks.length;
			size_t num_chord_elements = wp_geom.ctrl_chunks.length/num_span_elements;

			foreach(idx, ref ctrl_chunk; wp_geom.ctrl_chunks) {
				foreach(c_idx; 0..chunk_size) {
					immutable span_chunk_idx = idx%wp_geom.chunks.length;
					immutable span_idx = span_chunk_idx*chunk_size + c_idx;
					immutable chord_idx = (idx-span_chunk_idx)/(wp_geom.chunks.length);
					writeln("span_idx =", span_idx, "\tchunk_idx = ", chord_idx);

					double twist = wp_geom.chunks[span_chunk_idx].twist[c_idx];
					double xp0 = ctrl_chunk.ctrl_pt_x[c_idx]; // chordwise location
					double y0 = ctrl_chunk.ctrl_pt_y[c_idx]; // spanwise location
					double zp0 = ctrl_chunk.ctrl_pt_z[c_idx]; // camber(thickness) location

					double x0 = xp0*std.math.cos(twist) + zp0*std.math.sin(twist);
					double z0 = -xp0*std.math.sin(twist) + zp0*std.math.cos(twist);

					auto p0 = Vec3(x0,y0,z0);
					auto id0 = vtk_wing.points.InsertNextPoint(p0[0], p0[1], p0[2]);

					//writeln("p0 = ", p0);

					vtk_wing.base_points[wp_idx][id0] = p0;
					vtk_wing.base_point_ids[wp_idx][span_idx][chord_idx] = id0;
				}
			}

			foreach(sp_idx; 0..num_span_elements-1){
				foreach(crd_idx; 0..num_chord_elements-1){
					vtkIdType[4] ids = [
						vtk_wing.base_point_ids[wp_idx][sp_idx][crd_idx],
						vtk_wing.base_point_ids[wp_idx][sp_idx + 1][crd_idx],
						vtk_wing.base_point_ids[wp_idx][sp_idx + 1][crd_idx + 1],
						vtk_wing.base_point_ids[wp_idx][sp_idx][crd_idx + 1]
					];
					vtk_wing.grid.InsertNextCell(VTK__POLYGON, ids.length, ids.ptr);
				}
			}
		}

		auto point_data = vtk_wing.grid.GetPointData;
		//point_data.AddArray(vtk_wing.loads);
		point_data.AddArray(vtk_wing.aoa);
		point_data.AddArray(vtk_wing.u_p);
		point_data.AddArray(vtk_wing.u_t);
		//point_data.AddArray(vtk_wing.gamma);
		point_data.AddArray(vtk_wing.dC_T);

		return vtk_wing;
	} else {
		auto vtk_wing = new VtkWing;
		return vtk_wing;
	}
}

VtkRotor build_base_vtu_rotor(RG)(auto ref RG rotor_geo) {
	version(Have_vtkd) {
		immutable elements = rotor_geo.blades[0].chunks.length*chunk_size;

		auto vtk_rotor = new VtkRotor(rotor_geo.blades.length);

		vtk_rotor.origin = rotor_geo.origin;

		vtk_rotor.loads.SetNumberOfComponents(1);
		vtk_rotor.loads.SetNumberOfTuples(elements*rotor_geo.blades.length*naca0012.length);
		vtk_rotor.loads.SetName("dC_L");
		
		vtk_rotor.dC_L_dot.SetNumberOfComponents(1);
		vtk_rotor.dC_L_dot.SetNumberOfTuples(elements*rotor_geo.blades.length*naca0012.length);
		vtk_rotor.dC_L_dot.SetName("dC_L_dot");

		vtk_rotor.dC_T_dot.SetNumberOfComponents(1);
		vtk_rotor.dC_T_dot.SetNumberOfTuples(elements*rotor_geo.blades.length*naca0012.length);
		vtk_rotor.dC_T_dot.SetName("dC_T_dot");

		vtk_rotor.dC_T.SetNumberOfComponents(1);
		vtk_rotor.dC_T.SetNumberOfTuples(elements*rotor_geo.blades.length*naca0012.length);
		vtk_rotor.dC_T.SetName("dC_T");

		vtk_rotor.dC_Q.SetNumberOfComponents(1);
		vtk_rotor.dC_Q.SetNumberOfTuples(elements*rotor_geo.blades.length*naca0012.length);
		vtk_rotor.dC_Q.SetName("dC_Q");

		vtk_rotor.dC_D.SetNumberOfComponents(1);
		vtk_rotor.dC_D.SetNumberOfTuples(elements*rotor_geo.blades.length*naca0012.length);
		vtk_rotor.dC_D.SetName("dC_D");

		vtk_rotor.dC_N.SetNumberOfComponents(1);
		vtk_rotor.dC_N.SetNumberOfTuples(elements*rotor_geo.blades.length*naca0012.length);
		vtk_rotor.dC_N.SetName("dC_N");

		vtk_rotor.dC_c.SetNumberOfComponents(1);
		vtk_rotor.dC_c.SetNumberOfTuples(elements*rotor_geo.blades.length*naca0012.length);
		vtk_rotor.dC_c.SetName("dC_c");

		vtk_rotor.aoa.SetNumberOfComponents(1);
		vtk_rotor.aoa.SetNumberOfTuples(elements*rotor_geo.blades.length*naca0012.length);
		vtk_rotor.aoa.SetName("aoa");

		vtk_rotor.aoa_eff.SetNumberOfComponents(1);
		vtk_rotor.aoa_eff.SetNumberOfTuples(elements*rotor_geo.blades.length*naca0012.length);
		vtk_rotor.aoa_eff.SetName("aoa_eff");

		vtk_rotor.u_p.SetNumberOfComponents(1);
		vtk_rotor.u_p.SetNumberOfTuples(elements*rotor_geo.blades.length*naca0012.length);
		vtk_rotor.u_p.SetName("u_p");

		vtk_rotor.shed_u_p.SetNumberOfComponents(1);
		vtk_rotor.shed_u_p.SetNumberOfTuples(elements*rotor_geo.blades.length*naca0012.length);
		vtk_rotor.shed_u_p.SetName("shed_u_p");

		vtk_rotor.u_t.SetNumberOfComponents(1);
		vtk_rotor.u_t.SetNumberOfTuples(elements*rotor_geo.blades.length*naca0012.length);
		vtk_rotor.u_t.SetName("u_t");

		vtk_rotor.inflow_angle.SetNumberOfComponents(1);
		vtk_rotor.inflow_angle.SetNumberOfTuples(elements*rotor_geo.blades.length*naca0012.length);
		vtk_rotor.inflow_angle.SetName("inflow angle");

		vtk_rotor.gamma.SetNumberOfComponents(1);
		vtk_rotor.gamma.SetNumberOfTuples(elements*rotor_geo.blades.length*naca0012.length);
		vtk_rotor.gamma.SetName("gamma");

		vtk_rotor.d_gamma.SetNumberOfComponents(1);
		vtk_rotor.d_gamma.SetNumberOfTuples(elements*rotor_geo.blades.length*naca0012.length);
		vtk_rotor.d_gamma.SetName("d_gamma");

		vtk_rotor.theta.SetNumberOfComponents(1);
		vtk_rotor.theta.SetNumberOfTuples(elements*rotor_geo.blades.length*naca0012.length);
		vtk_rotor.theta.SetName("theta");

		vtk_rotor.af_norm.SetNumberOfComponents(3);
		vtk_rotor.af_norm.SetNumberOfTuples(elements*rotor_geo.blades.length*naca0012.length);
		vtk_rotor.af_norm.SetName("af_norm");

		vtk_rotor.blade_local_vel.SetNumberOfComponents(3);
		vtk_rotor.blade_local_vel.SetNumberOfTuples(elements*rotor_geo.blades.length*naca0012.length);
		vtk_rotor.blade_local_vel.SetName("blade_local_vel");

		vtk_rotor.projected_vel.SetNumberOfComponents(3);
		vtk_rotor.projected_vel.SetNumberOfTuples(elements*rotor_geo.blades.length*naca0012.length);
		vtk_rotor.projected_vel.SetName("projected_vel");

		foreach(blade_idx, ref blade_geo; rotor_geo.blades) {

			vtk_rotor.azimuth_offsets[blade_idx] = blade_geo.azimuth_offset;

			vtkIdType last_id = 0;
			immutable double r_c = blade_geo.r_c;
			
			foreach(idx, ref chunk; blade_geo.chunks) {
				foreach(c_idx, ref r; chunk.r) {
					vtk_rotor.r_to_point_map ~= new vtkIdType[1];

					auto af_local_point = Vec3((r - r_c), chunk.xi[c_idx], 0);
					//auto af_local_point = Vec3(r, chunk.xi[c_idx], 0);

					//auto p0 = blade_geo.frame.global_matrix*af_local_point;

					auto id0 = vtk_rotor.points.InsertNextPoint(af_local_point[0], af_local_point[1], af_local_point[2]);
					
					vtk_rotor.base_points[blade_idx][id0] = Vec3(af_local_point[0], af_local_point[1], af_local_point[2]);
					vtk_rotor.r_to_point_map[$-1][0] = id0;

					if(c_idx != 0 || idx != 0) {
						vtkIdType[2] ids = [last_id, id0];
						vtk_rotor.grid.InsertNextCell(VTK__POLY_LINE, ids.length, ids.ptr);
					}

					last_id = id0;
				}
			}
		}

		auto point_data = vtk_rotor.grid.GetPointData;
		point_data.AddArray(vtk_rotor.loads);
		point_data.AddArray(vtk_rotor.dC_L_dot);
		point_data.AddArray(vtk_rotor.dC_T_dot);
		point_data.AddArray(vtk_rotor.dC_T);
		point_data.AddArray(vtk_rotor.dC_D);
		point_data.AddArray(vtk_rotor.dC_Q);
		point_data.AddArray(vtk_rotor.dC_N);
		point_data.AddArray(vtk_rotor.dC_c);
		point_data.AddArray(vtk_rotor.aoa);
		point_data.AddArray(vtk_rotor.aoa_eff);
		point_data.AddArray(vtk_rotor.u_p);
		point_data.AddArray(vtk_rotor.shed_u_p);
		point_data.AddArray(vtk_rotor.u_t);
		point_data.AddArray(vtk_rotor.inflow_angle);
		point_data.AddArray(vtk_rotor.gamma);
		point_data.AddArray(vtk_rotor.d_gamma);
		point_data.AddArray(vtk_rotor.theta);
		point_data.AddArray(vtk_rotor.af_norm);
		point_data.AddArray(vtk_rotor.blade_local_vel);
		point_data.AddArray(vtk_rotor.projected_vel);

		return vtk_rotor;
	} else {
		auto vtk_rotor = new VtkRotor();
		return vtk_rotor;
	}
	
}

class VtkWake {
	this() {

	}

	VtkRotorWake[] rotor_wakes;

	version(Have_vtkd) {
		private this(size_t num_rotors, size_t num_blades, size_t[] shed_length, size_t elements) {
			rotor_wakes = new VtkRotorWake[num_rotors];

			foreach(r_idx, ref rotor_wake; rotor_wakes) {
				rotor_wake = new VtkRotorWake(num_blades, shed_length[r_idx], elements);
			}
		}
	}
}

class VtkWingWake {
	this(){

	}

	version(Have_vtkd) {
		private vtkIdType[][][] basepoint_ids;
		private vtkPoints* filament_points;
		private vtkUnstructuredGrid* grid;
		private Vec3[vtkIdType][] base_points;
		private vtkDoubleArray* gamma;
		private vtkDoubleArray* A_kl;

		private this(size_t num_wing_parts) {
			base_points = new Vec3[vtkIdType][num_wing_parts];
			gamma = vtkDoubleArray.New;
			A_kl = vtkDoubleArray.New;
			grid = vtkUnstructuredGrid.New;
			filament_points = vtkPoints.New;
		}		
	}
}

VtkWingWake build_base_vtu_wing_wake(W, WLS)(auto ref W wing_geom, auto ref WLS wing_lift_surf) {

	version(Have_vtkd) {
		immutable num_wing_part = wing_geom.wing_parts.length;
		immutable num_span_elements = wing_geom.wing_parts[0].chunks.length * chunk_size;
		immutable num_chord_elements = wing_lift_surf.wing_part_lift_surf[0].spanwise_filaments.length;

		auto vtk_wing_wake = new VtkWingWake(num_wing_part);

		vtk_wing_wake.gamma.SetNumberOfComponents(1);
		vtk_wing_wake.gamma.SetNumberOfTuples(num_span_elements*num_chord_elements*num_wing_part);
		vtk_wing_wake.gamma.SetName("gamma");


		vtk_wing_wake.A_kl.SetNumberOfComponents(1);
		vtk_wing_wake.A_kl.SetNumberOfTuples(num_span_elements*num_chord_elements*num_wing_part);
		vtk_wing_wake.A_kl.SetName("A_kl");

		vtk_wing_wake.basepoint_ids = new vtkIdType[][][](num_wing_part,num_span_elements,num_chord_elements);

		foreach(wp_idx, ref wing_part_ls ; wing_lift_surf.wing_part_lift_surf) {
			foreach(fl_idx, fil; wing_part_ls.spanwise_filaments){
				foreach(ch_idx, chunk; fil.chunks){
					foreach(c_idx; 0..chunk_size){
						double twist = wing_geom.wing_parts[wp_idx].chunks[ch_idx].twist[c_idx];
						double xp0 = chunk.x[c_idx];
						double y0 = chunk.y[c_idx];
						double zp0 = chunk.z[c_idx];

						double x0 = xp0*std.math.cos(twist) + zp0*std.math.sin(twist);
						double z0 = -xp0*std.math.sin(twist) + zp0*std.math.cos(twist);

						auto p0 = Vec3(x0,y0,z0);
						auto id0 = vtk_wing_wake.filament_points.InsertNextPoint(p0[0], p0[1], p0[2]);

						vtk_wing_wake.base_points[wp_idx][id0] = p0;
						vtk_wing_wake.basepoint_ids[wp_idx][ch_idx][fl_idx] = id0;
					}
				}
			}

			foreach(crd_idx; 0..num_chord_elements){
				foreach(sp_idx; 0..num_span_elements-1){
					vtkIdType[2] ids = [
						vtk_wing_wake.basepoint_ids[wp_idx][sp_idx][crd_idx],
						vtk_wing_wake.basepoint_ids[wp_idx][sp_idx+1][crd_idx]
					];
					vtk_wing_wake.grid.InsertNextCell(VTK__POLY_LINE, ids.length, ids.ptr);
				}
			}
		}

		auto point_data = vtk_wing_wake.grid.GetPointData;
		point_data.AddArray(vtk_wing_wake.gamma);
		point_data.AddArray(vtk_wing_wake.A_kl);

		return vtk_wing_wake;
	} else {
		auto vtk_wing_wake = new VtkWingWake;
		return vtk_wing;
	}
}

void write_wing_wake_vtu(WG, WLS, WIS)(string base_filename, size_t iteration, size_t wing_idx, ref VtkWingWake wing_wake, auto ref WG wing_geom, auto ref WLS wing_lift_surf, auto ref WIS wing_input) {
	
	version(Have_vtkd) {
		writeln("writing wing wake vtu");
		immutable spanwise_elements = wing_lift_surf.wing_part_lift_surf[0].spanwise_filaments[0].length;
		immutable chordwise_elements = wing_lift_surf.wing_part_lift_surf[0].spanwise_filaments.length;

		auto wing_wake_writer = vtkXMLUnstructuredGridWriter.New;

		auto aoa_rotation = 
			Mat3(
				-std.math.cos(wing_input.angle_of_attack /++ flip_angle+/), 0, -std.math.sin(wing_input.angle_of_attack/+ + flip_angle+/),
				0, -1, 0,
				-std.math.sin(wing_input.angle_of_attack /++ flip_angle+/), 0, std.math.cos(wing_input.angle_of_attack /++ flip_angle+/)
			);

		wing_wake.grid.SetPoints(wing_wake.filament_points);

		foreach(wp_idx, wing_part_ls; wing_lift_surf.wing_part_lift_surf){
			foreach(pi; wing_wake.base_points[wp_idx].byKeyValue){
				vtkIdType id = pi.key;
				Vec3 point = pi.value;
				auto origin = wing_geom.origin;

				auto aoa_rot = aoa_rotation*point;
				auto final_p =  aoa_rot + origin;
				writeln("final_p = ", final_p);

				wing_wake.filament_points.SetPoint(id, final_p[0], final_p[1], final_p[2]);				
			}

			foreach(fl_idx, fil; wing_part_ls.spanwise_filaments){
				foreach(ch_idx, chunk; fil.chunks){
					foreach(c_idx; 0..chunk_size){

						auto id = wing_wake.basepoint_ids[wp_idx][ch_idx][fl_idx];

						wing_wake.gamma.SetTuple1(id, chunk.gamma[c_idx]);
						wing_wake.A_kl.SetTuple1(id, chunk.A_kl[c_idx]);
					}
				}
			}
		}

		import std.string : toStringz;

		auto filename = base_filename~"_"~wing_idx.to!string~"_"~iteration.to!string~".vtu";

		wing_wake_writer.SetFileName(filename.toStringz);
		wing_wake_writer.SetInputData(wing_wake.grid);
		wing_wake_writer.Write;
	}
}

class VtkRotorWake {
	this() {
	}

	version(Have_vtkd) {
		private vtkIdType[] shed_point_ids;
		private vtkIdType[] shed_cell_ids;
		private vtkXMLUnstructuredGridWriter* writer;
		private vtkAppendFilter*[] shed_appenders;
		private vtkUnstructuredGrid*[] shed_grids;
		private vtkPoints*[] shed_points;
		private vtkDoubleArray*[] shed_induced;
		private vtkDoubleArray*[] shed_circ;
		private vtkIdType[] tip_point_ids;
		private vtkIdType[] tip_cell_ids;
		private vtkUnstructuredGrid*[] tip_grids;
		private vtkPoints*[] tip_points;
		private vtkDoubleArray*[] tip_induced;
		private vtkDoubleArray*[] tip_circ;
		private vtkDoubleArray*[] tip_core_size;
		private vtkDoubleArray*[] tip_d_volume;

		this(size_t num_blades, size_t shed_length, size_t elements) {

			writer = vtkXMLUnstructuredGridWriter.New;

			shed_appenders = new vtkAppendFilter*[num_blades];
			foreach(ref shed_appender; shed_appenders) {
				shed_appender = vtkAppendFilter.New;
			}
			
			shed_grids = new vtkUnstructuredGrid*[num_blades*shed_length];
			foreach(ref shed_grid; shed_grids) {
				shed_grid = vtkUnstructuredGrid.New;
			}

			shed_points = new vtkPoints*[num_blades*shed_length];
			foreach(ref shed_point; shed_points) {
				shed_point = vtkPoints.New(VTK__DOUBLE);
			}

			shed_induced = new vtkDoubleArray*[num_blades*shed_length];
			foreach(ref induced; shed_induced) {
				induced = vtkDoubleArray.New;
			}

			shed_circ = new vtkDoubleArray*[num_blades*shed_length];
			foreach(ref circ; shed_circ) {
				circ = vtkDoubleArray.New;
			}
	
			tip_grids = new vtkUnstructuredGrid*[num_blades];
			foreach(ref tip_grid; tip_grids) {
				tip_grid = vtkUnstructuredGrid.New;
			}

			tip_points = new vtkPoints*[num_blades];
			foreach(ref tip_point; tip_points) {
				tip_point = vtkPoints.New(VTK__DOUBLE);
			}

			tip_induced = new vtkDoubleArray*[num_blades];
			foreach(ref induced; tip_induced) {
				induced = vtkDoubleArray.New;
			}

			tip_circ = new vtkDoubleArray*[num_blades];
			foreach(ref circ; tip_circ) {
				circ = vtkDoubleArray.New;
			}

			tip_core_size = new vtkDoubleArray*[num_blades];
			foreach(ref core_size; tip_core_size) {
				core_size = vtkDoubleArray.New;
			}

			tip_d_volume = new vtkDoubleArray*[num_blades];
			foreach(ref d_volume; tip_d_volume) {
				d_volume = vtkDoubleArray.New;
			}
		}

		~this() {
			writer.Delete;

			foreach(ref shed_appender; shed_appenders) {
				shed_appender.Delete;
			}
			
			foreach(ref shed_grid; shed_grids) {
				shed_grid.Delete;
			}

			foreach(ref shed_point; shed_points) {
				shed_point.Delete;
			}

			foreach(ref induced; shed_induced) {
				induced.Delete;
			}

			foreach(ref circ; shed_circ) {
				circ.Delete;
			}
	
			foreach(ref tip_grid; tip_grids) {
				tip_grid.Delete;
			}

			foreach(ref tip_point; tip_points) {
				tip_point.Delete;
			}

			foreach(ref induced; tip_induced) {
				induced.Delete;
			}

			foreach(ref circ; tip_circ) {
				circ.Delete;
			}

			foreach(ref core_size; tip_core_size) {
				core_size.Delete;
			}

			foreach(ref d_volume; tip_d_volume) {
				d_volume.Delete;
			}
		}
	}
}

VtkWake build_base_vtu_wake(W)(auto ref W wake) {
	
	version(Have_vtkd) {
		import std.algorithm : map;
		import std.array : array;

		size_t[] shed_length = wake.rotor_wakes.map!(r => r.shed_vortices[0].shed_filaments.length).array;

		immutable elements = wake.rotor_wakes[0].shed_vortices[0].shed_filaments[0].length*chunk_size;
		
		auto vtk_wake = new VtkWake(wake.rotor_wakes.length, wake.rotor_wakes[0].tip_vortices.length, shed_length, elements);

		foreach(rotor_idx, rotor_wake; wake.rotor_wakes) {
			vtkIdType last_point_id;

			size_t shed_idx = 0;
			foreach(blade_idx, shed_wake; rotor_wake.shed_vortices) {
				
				foreach(shed_filament; shed_wake.shed_filaments) {
					
					auto shed_wake_len = shed_filament.chunks.length*chunk_size;

					vtk_wake.rotor_wakes[rotor_idx].shed_grids[shed_idx].Allocate(shed_wake_len);
					
					vtk_wake.rotor_wakes[rotor_idx].shed_points[shed_idx].Allocate(shed_wake_len);

					vtk_wake.rotor_wakes[rotor_idx].shed_induced[shed_idx].SetNumberOfComponents(1);
					vtk_wake.rotor_wakes[rotor_idx].shed_induced[shed_idx].SetNumberOfTuples(shed_wake_len);
					vtk_wake.rotor_wakes[rotor_idx].shed_induced[shed_idx].SetName("induced");

					vtk_wake.rotor_wakes[rotor_idx].shed_circ[shed_idx].SetNumberOfComponents(1);
					vtk_wake.rotor_wakes[rotor_idx].shed_circ[shed_idx].SetNumberOfTuples(shed_wake_len - 1);
					vtk_wake.rotor_wakes[rotor_idx].shed_circ[shed_idx].SetName("d_circulation");

					foreach(idx, shed_chunk; shed_filament.chunks) {
						foreach(c_idx; 0..chunk_size) {
							auto point_id = vtk_wake.rotor_wakes[rotor_idx].shed_points[shed_idx].InsertNextPoint(0, 0, 0);

							vtk_wake.rotor_wakes[rotor_idx].shed_point_ids ~= point_id;
							vtk_wake.rotor_wakes[rotor_idx].shed_induced[shed_idx].SetTuple1(point_id, 0);
							if((idx > 0) || ((idx == 0) && (c_idx > 0))) {
								vtkIdType[2] ids = [last_point_id, point_id];

								auto cell_id = vtk_wake.rotor_wakes[rotor_idx].shed_grids[shed_idx].InsertNextCell(VTK__POLY_LINE, ids.length, ids.ptr);
								vtk_wake.rotor_wakes[rotor_idx].shed_cell_ids ~= cell_id;

								vtk_wake.rotor_wakes[rotor_idx].shed_circ[shed_idx].SetTuple1(cell_id, 0);
							}

							last_point_id = point_id;
						}
					}

					vtk_wake.rotor_wakes[rotor_idx].shed_grids[shed_idx].SetPoints(vtk_wake.rotor_wakes[rotor_idx].shed_points[shed_idx]);

					auto point_data = vtk_wake.rotor_wakes[rotor_idx].shed_grids[shed_idx].GetPointData;
					point_data.AddArray(vtk_wake.rotor_wakes[rotor_idx].shed_induced[shed_idx]);

					auto cell_data = vtk_wake.rotor_wakes[rotor_idx].shed_grids[shed_idx].GetCellData;
					cell_data.AddArray(vtk_wake.rotor_wakes[rotor_idx].shed_circ[shed_idx]);

					vtk_wake.rotor_wakes[rotor_idx].shed_appenders[blade_idx].AddInputData(vtk_wake.rotor_wakes[rotor_idx].shed_grids[shed_idx]);
					vtk_wake.rotor_wakes[rotor_idx].shed_appenders[blade_idx].Update;

					shed_idx++;
				}
			}

			foreach(blade_idx, tip_vortex; rotor_wake.tip_vortices) {

				immutable wake_length = tip_vortex.chunks.length*chunk_size;

				vtk_wake.rotor_wakes[rotor_idx].tip_grids[blade_idx].Allocate(wake_length);

				vtk_wake.rotor_wakes[rotor_idx].tip_points[blade_idx].Allocate(wake_length);

				vtk_wake.rotor_wakes[rotor_idx].tip_induced[blade_idx].SetNumberOfComponents(1);
				vtk_wake.rotor_wakes[rotor_idx].tip_induced[blade_idx].SetNumberOfTuples(wake_length);
				vtk_wake.rotor_wakes[rotor_idx].tip_induced[blade_idx].SetName("induced");

				vtk_wake.rotor_wakes[rotor_idx].tip_d_volume[blade_idx].SetNumberOfComponents(1);
				vtk_wake.rotor_wakes[rotor_idx].tip_d_volume[blade_idx].SetNumberOfTuples(wake_length - 1);
				vtk_wake.rotor_wakes[rotor_idx].tip_d_volume[blade_idx].SetName("d volume");

				vtk_wake.rotor_wakes[rotor_idx].tip_core_size[blade_idx].SetNumberOfComponents(1);
				vtk_wake.rotor_wakes[rotor_idx].tip_core_size[blade_idx].SetNumberOfTuples(wake_length - 1);
				vtk_wake.rotor_wakes[rotor_idx].tip_core_size[blade_idx].SetName("core size");

				vtk_wake.rotor_wakes[rotor_idx].tip_circ[blade_idx].SetNumberOfComponents(1);
				vtk_wake.rotor_wakes[rotor_idx].tip_circ[blade_idx].SetNumberOfTuples(wake_length - 1);
				vtk_wake.rotor_wakes[rotor_idx].tip_circ[blade_idx].SetName("circulation");

				foreach(idx, tip_chunk; tip_vortex.chunks) {
					foreach(c_idx; 0..chunk_size) {
						auto point_id =  vtk_wake.rotor_wakes[rotor_idx].tip_points[blade_idx].InsertNextPoint(0, 0, 0);

						vtk_wake.rotor_wakes[rotor_idx].tip_point_ids ~= point_id;

						vtk_wake.rotor_wakes[rotor_idx].tip_induced[blade_idx].SetTuple1(point_id, 0);
						if((idx > 0) || ((idx == 0) && (c_idx > 0))) {
							vtkIdType[2] ids = [last_point_id, point_id];

							auto cell_id = vtk_wake.rotor_wakes[rotor_idx].tip_grids[blade_idx].InsertNextCell(VTK__POLY_LINE, ids.length, ids.ptr);

							vtk_wake.rotor_wakes[rotor_idx].tip_cell_ids ~= cell_id;

							vtk_wake.rotor_wakes[rotor_idx].tip_core_size[blade_idx].SetTuple1(cell_id, 0);
							vtk_wake.rotor_wakes[rotor_idx].tip_d_volume[blade_idx].SetTuple1(cell_id, 0);
							vtk_wake.rotor_wakes[rotor_idx].tip_circ[blade_idx].SetTuple1(cell_id, 0);
						}

						last_point_id = point_id;
					}
				}

				vtk_wake.rotor_wakes[rotor_idx].tip_grids[blade_idx].SetPoints(vtk_wake.rotor_wakes[rotor_idx].tip_points[blade_idx]);

				auto point_data = vtk_wake.rotor_wakes[rotor_idx].tip_grids[blade_idx].GetPointData;
				point_data.AddArray(vtk_wake.rotor_wakes[rotor_idx].tip_induced[blade_idx]);

				auto cell_data = vtk_wake.rotor_wakes[rotor_idx].tip_grids[blade_idx].GetCellData;
				cell_data.AddArray(vtk_wake.rotor_wakes[rotor_idx].tip_circ[blade_idx]);
				cell_data.AddArray(vtk_wake.rotor_wakes[rotor_idx].tip_core_size[blade_idx]);
				cell_data.AddArray(vtk_wake.rotor_wakes[rotor_idx].tip_d_volume[blade_idx]);
			}
		}

		return vtk_wake;
	} else {
		auto vtk_wake = new VtkWake();
		return vtk_wake;
	}
}

void write_wake_vtu(W)(string base_filename, size_t iteration, VtkWake vtk_wake, auto ref W wake) {

	version(Have_vtkd) {
		foreach(rotor_idx, rotor_wake; wake.rotor_wakes) {
			vtkIdType last_point_id;

			size_t shed_point_idx = 0;
			size_t shed_cell_idx = 0;
			size_t shed_idx = 0;
			foreach(blade_idx, shed_wake; rotor_wake.shed_vortices) {
				
				vtk_wake.rotor_wakes[rotor_idx].shed_appenders[blade_idx].RemoveAllInputs;

				auto shed_wake_len = shed_wake.shed_filaments[0].chunks.length*chunk_size;
				
				foreach(shed_filament; shed_wake.shed_filaments) {

					double last_circ;
					double last_core_size;
					foreach(idx, shed_chunk; shed_filament.chunks) {
						foreach(c_idx; 0..chunk_size) {
							auto point_id = vtk_wake.rotor_wakes[rotor_idx].shed_point_ids[shed_point_idx];

							vtk_wake.rotor_wakes[rotor_idx].shed_points[shed_idx].SetPoint(point_id, shed_chunk.x[c_idx], shed_chunk.y[c_idx], shed_chunk.z[c_idx]);

							last_circ = shed_chunk.gamma[c_idx];

							vtk_wake.rotor_wakes[rotor_idx].shed_induced[shed_idx].SetTuple1(point_id, shed_chunk.v_z[c_idx]);

							if((idx > 0) || ((idx == 0) && (c_idx > 0))) {
								auto cell_id = vtk_wake.rotor_wakes[rotor_idx].shed_cell_ids[shed_cell_idx];
								vtk_wake.rotor_wakes[rotor_idx].shed_circ[shed_idx].SetTuple1(cell_id, last_circ);

								shed_cell_idx++;
							}

							last_point_id = point_id;

							shed_point_idx++;
						}
					}

					vtk_wake.rotor_wakes[rotor_idx].shed_appenders[blade_idx].AddInputData(vtk_wake.rotor_wakes[rotor_idx].shed_grids[shed_idx]);
					vtk_wake.rotor_wakes[rotor_idx].shed_appenders[blade_idx].Update;

					shed_idx++;
				}

				import std.string : toStringz;
				auto filename = base_filename~"_shed_"~rotor_idx.to!string~"_"~blade_idx.to!string~"_"~iteration.to!string~".vtu";
				vtk_wake.rotor_wakes[rotor_idx].writer.SetFileName(filename.toStringz);
				vtk_wake.rotor_wakes[rotor_idx].writer.SetInputData(vtk_wake.rotor_wakes[rotor_idx].shed_appenders[blade_idx].GetOutput);
				vtk_wake.rotor_wakes[rotor_idx].writer.Write;
			}

			size_t tip_point_idx = 0;
			size_t tip_cell_idx = 0;

			foreach(blade_idx, tip_vortex; rotor_wake.tip_vortices) {

				double last_circ;
				double last_core_size;
				double last_d_volume;
				foreach(idx, tip_chunk; tip_vortex.chunks) {
					foreach(c_idx; 0..chunk_size) {
						auto point_id = vtk_wake.rotor_wakes[rotor_idx].tip_point_ids[tip_point_idx];
						vtk_wake.rotor_wakes[rotor_idx].tip_points[blade_idx].SetPoint(point_id, tip_chunk.x[c_idx], tip_chunk.y[c_idx], tip_chunk.z[c_idx]);

						last_circ = tip_chunk.gamma[c_idx];
						last_core_size = tip_chunk.r_c[c_idx];
						last_d_volume = tip_chunk.d_volume[c_idx];
						vtk_wake.rotor_wakes[rotor_idx].tip_induced[blade_idx].SetTuple1(point_id, tip_chunk.v_z[c_idx]);

						if((idx != tip_vortex.chunks.length - 1) || ((idx == tip_vortex.chunks.length - 1) && (c_idx < (chunk_size - 1)) )) {

							auto cell_id = vtk_wake.rotor_wakes[rotor_idx].tip_cell_ids[tip_cell_idx];
							vtk_wake.rotor_wakes[rotor_idx].tip_core_size[blade_idx].SetTuple1(cell_id, last_core_size);
							vtk_wake.rotor_wakes[rotor_idx].tip_d_volume[blade_idx].SetTuple1(cell_id, last_d_volume);
							vtk_wake.rotor_wakes[rotor_idx].tip_circ[blade_idx].SetTuple1(cell_id, last_circ);

							tip_cell_idx++;
						}

						last_point_id = point_id;
						tip_point_idx++;
					}
				}

				import std.string : toStringz;
				auto filename = base_filename~"_"~rotor_idx.to!string~"_"~blade_idx.to!string~"_"~iteration.to!string~".vtu";
				vtk_wake.rotor_wakes[rotor_idx].writer.SetFileName(filename.toStringz);
				vtk_wake.rotor_wakes[rotor_idx].writer.SetInputData(vtk_wake.rotor_wakes[rotor_idx].tip_grids[blade_idx]);
				vtk_wake.rotor_wakes[rotor_idx].writer.Write;
			}
		}
	}
}

void write_inflow_vtu(I, RGA)(string filename, I[] inflows, Vec3 delta, Vec3 starts, size_t num_x, size_t num_y, size_t num_z, double aoa, double[] omegas, RGA rotors) {

	version(Have_vtkd) {
		vtkIdType[PosVector] node_id_map;
		auto points = vtkPoints.New(VTK__DOUBLE);
		auto induced_writer = vtkXMLUnstructuredGridWriter.New;

		auto induced_grid = vtkUnstructuredGrid.New;

		points.Allocate(num_x*num_y*num_z);
		induced_grid.Allocate(num_x*num_y*num_z);

		auto induced = vtkDoubleArray.New();
		induced.SetNumberOfComponents(3);
		induced.SetNumberOfTuples(num_x*num_y*num_z);
		induced.SetName("induced");

		import std.stdio : writeln;
		writeln("Generating vel field");
		foreach(z_idx; 0..num_z) {
			foreach(y_idx; 0..num_y) {
				foreach(x_chunk; 0..num_x/chunk_size) {
					vtkIdType[chunk_size] ids;

					Chunk x_e;
					auto g_pos = Vector!(4, Chunk)(1.0);
					g_pos[1][] = y_idx*delta[1] + starts[1];
					g_pos[2][] = z_idx*delta[2] + starts[2];

					foreach(x_c; 0..chunk_size) {
						size_t x_idx = x_chunk*chunk_size + x_c;
						g_pos[0][x_c] = x_idx*delta[0] + starts[0];
						auto int_pos = PosVector(x_idx, y_idx, z_idx);
						ids[x_c] = points.InsertNextPoint(g_pos[0][x_c], g_pos[1][x_c], g_pos[2][x_c]);
						
						node_id_map[int_pos] = ids[x_c];
					}

					auto global_inflow = Vector!(4, Chunk)(1);

					foreach(rotor_idx, ref inflow; inflows) {

						g_pos[3][] = 1;
						
						auto l_pos = inflow.frame.global_matrix.inverse().get() * (g_pos);

						l_pos /= rotors[rotor_idx].radius;

						immutable Chunk infl = rotors[rotor_idx].radius*omegas[rotor_idx]*inflow.inflow_at(l_pos[0], l_pos[1], l_pos[2], x_e, aoa)[];

						auto local_inflow = Vector!(4, Chunk)(zero, zero, infl, zero);

						local_inflow[2][] = infl[];
						global_inflow += inflow.frame.global_matrix * local_inflow;
					}

					foreach(c_idx; 0..chunk_size) {
						induced.SetTuple3(ids[c_idx], global_inflow[0][c_idx], global_inflow[1][c_idx], global_inflow[2][c_idx]);
					}
				}
			}
		}

		writeln("Done generating vel field");

		foreach(z; 0..num_z) {
			foreach(y; 0..num_y) {
				foreach(x; 0..num_x) {
					auto pos1 = PosVector(x, y, z);
					auto pos2 = PosVector(x + 1, y, z);
					auto pos3 = PosVector(x, y + 1, z);
					auto pos4 = PosVector(x + 1, y + 1, z);
					
					auto pos5 = PosVector(x, y, z + 1);
					auto pos6 = PosVector(x + 1, y, z + 1);
					auto pos7 = PosVector(x, y + 1, z + 1);
					auto pos8 = PosVector(x + 1, y + 1, z + 1);
					
					if((pos1 in node_id_map) && (pos2 in node_id_map) &&
						(pos3 in node_id_map) && (pos4 in node_id_map) &&
						(pos5 in node_id_map) && (pos6 in node_id_map) &&
						(pos7 in node_id_map) && (pos8 in node_id_map)

					) {
						vtkIdType[8] ids = [
							node_id_map[pos1],
							node_id_map[pos2],
							node_id_map[pos3],
							node_id_map[pos4],
							node_id_map[pos5],
							node_id_map[pos6],
							node_id_map[pos7],
							node_id_map[pos8]
						];

						induced_grid.InsertNextCell(VTK__VOXEL, ids.length, ids.ptr);
					}
				}
			}
		}

		writeln("Done generating grid");

		induced_grid.SetPoints(points);

		auto point_data = induced_grid.GetPointData;
		point_data.AddArray(induced);

		import std.string : toStringz;

		induced_writer.SetFileName(filename.toStringz);
		induced_writer.SetInputData(induced_grid);
		induced_writer.Write;

		writeln("Done writing induced flow field");
	}
}

void write_wake_field_vtu(ACS, W)(string filename, auto ref ACS ac_state, auto ref W wake, Vec3 delta, Vec3 starts, size_t num_x, size_t num_y, size_t num_z) {

	version(Have_vtkd) {
		vtkIdType[PosVector] node_id_map;
		auto points = vtkPoints.New(VTK__DOUBLE);
		auto induced_writer = vtkXMLUnstructuredGridWriter.New;

		auto induced_grid = vtkUnstructuredGrid.New;

		points.Allocate(num_x*num_y*num_z);
		induced_grid.Allocate(num_x*num_y*num_z);
		
		auto induced = vtkDoubleArray.New();
		induced.SetNumberOfComponents(3);
		induced.SetNumberOfTuples(num_x*num_y*num_z);
		induced.SetName("v");

		import std.stdio : writeln;
		writeln("Generating vel field");
		foreach(z_idx; 0..num_z) {
			foreach(y_idx; 0..num_y) {
				foreach(x_chunk; 0..num_x/chunk_size) {
					vtkIdType[chunk_size] ids;
					Chunk y = y_idx*delta[1] + starts[1];
					Chunk z = z_idx*delta[2] + starts[2];
					Chunk x;
					Chunk x_e;
					foreach(x_c; 0..chunk_size) {
						size_t x_idx = x_chunk*chunk_size + x_c;
						x[x_c] = x_idx*delta[0] + starts[0];
						auto int_pos = PosVector(x_idx, y_idx, z_idx);
						ids[x_c] = points.InsertNextPoint(x[x_c], y[x_c], z[x_c]);
						
						node_id_map[int_pos] = ids[x_c];
					}

					auto lambda = compute_wake_induced_velocities(wake, x, y, z, ac_state, size_t.max, size_t.max);

					foreach(c_idx; 0..chunk_size) {
						induced.SetTuple3(ids[c_idx], lambda.v_x[c_idx], lambda.v_y[c_idx], lambda.v_z[c_idx]);
					}
				}
			}
		}

		writeln("Done generating vel field");

		foreach(z; 0..num_z) {
			foreach(y; 0..num_y) {
				foreach(x; 0..num_x) {
					auto pos1 = PosVector(x, y, z);
					auto pos2 = PosVector(x + 1, y, z);
					auto pos3 = PosVector(x, y + 1, z);
					auto pos4 = PosVector(x + 1, y + 1, z);
					
					auto pos5 = PosVector(x, y, z + 1);
					auto pos6 = PosVector(x + 1, y, z + 1);
					auto pos7 = PosVector(x, y + 1, z + 1);
					auto pos8 = PosVector(x + 1, y + 1, z + 1);
					
					if((pos1 in node_id_map) && (pos2 in node_id_map) &&
						(pos3 in node_id_map) && (pos4 in node_id_map) &&
						(pos5 in node_id_map) && (pos6 in node_id_map) &&
						(pos7 in node_id_map) && (pos8 in node_id_map)

					) {
						vtkIdType[8] ids = [
							node_id_map[pos1],
							node_id_map[pos2],
							node_id_map[pos3],
							node_id_map[pos4],
							node_id_map[pos5],
							node_id_map[pos6],
							node_id_map[pos7],
							node_id_map[pos8]
						];

						induced_grid.InsertNextCell(VTK__VOXEL, ids.length, ids.ptr);
					}
				}
			}
		}

		writeln("Done generating grid");

		induced_grid.SetPoints(points);

		auto point_data = induced_grid.GetPointData;
		point_data.AddArray(induced);

		import std.string : toStringz;

		induced_writer.SetFileName(filename.toStringz);
		induced_writer.SetInputData(induced_grid);
		induced_writer.Write;

		writeln("Done writing induced flow field");
	}
}