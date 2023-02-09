module opencopter.vtk;

import opencopter.aircraft;
import opencopter.config;
import opencopter.math;
import opencopter.memory;
import opencopter.wake;

import numd.linearalgebra.matrix : Matrix, Vector;

import std.conv;
import std.math;

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
		private vtkDoubleArray* u_p;
		private vtkDoubleArray* u_t;
		private vtkDoubleArray* inflow_angle;
		private vtkDoubleArray* gamma;
		private vtkDoubleArray* dC_T_dot;
		private vtkDoubleArray* dC_L_dot;

		private this(size_t num_blades) {
			azimuth_offsets = new double[num_blades];
			base_points = new Vec3[vtkIdType][num_blades];
			loads = vtkDoubleArray.New;
			dC_T_dot = vtkDoubleArray.New;
			dC_L_dot = vtkDoubleArray.New;
			u_p = vtkDoubleArray.New;
			u_t = vtkDoubleArray.New;
			inflow_angle = vtkDoubleArray.New;
			aoa = vtkDoubleArray.New;
			gamma = vtkDoubleArray.New;
			grid = vtkUnstructuredGrid.New;
			points = vtkPoints.New;
		}
	}
}

void write_rotor_vtu(RS, RIS)(string base_filename, size_t iteration, size_t rotor_idx, ref VtkRotor rotor, auto ref RS rotor_state, auto ref RIS rotor_input) {

	version(Have_vtkd) {
		immutable elements = rotor_state.blade_states[0].chunks.length*chunk_size;

		auto rotor_writer = vtkXMLUnstructuredGridWriter.New;

		auto rotor_sgn = sgn(rotor_input.angular_velocity);
		double flip_angle = 0;
		if(rotor_input.angular_velocity < 0) {
			flip_angle = PI;
		}
		auto aoa_rotation =
			Mat3(
				std.math.cos(rotor_input.angle_of_attack + flip_angle), 0, std.math.sin(rotor_input.angle_of_attack + flip_angle),
				0, 1, 0,
				-std.math.sin(rotor_input.angle_of_attack + flip_angle), 0, std.math.cos(rotor_input.angle_of_attack + flip_angle)
			);

		rotor.grid.SetPoints(rotor.points);

		foreach(b_idx, blade; rotor_state.blade_states) {

			auto azimuth_rotation =
				Mat3(
					cos(rotor_sgn*(rotor_input.azimuth + rotor.azimuth_offsets[b_idx])), -sin(rotor_sgn*(rotor_input.azimuth + rotor.azimuth_offsets[b_idx])), 0,
					sin(rotor_sgn*(rotor_input.azimuth + rotor.azimuth_offsets[b_idx])), cos(rotor_sgn*(rotor_input.azimuth + rotor.azimuth_offsets[b_idx])), 0,
					0, 0, 1
				);

			auto pitch_rotation =
				Mat3(
					1, 0, 0,
					0, cos(rotor_input.blade_pitches[b_idx]), -sin(rotor_input.blade_pitches[b_idx]),
					0, sin(rotor_input.blade_pitches[b_idx]), cos(rotor_input.blade_pitches[b_idx])
				);
			
			auto flap_rotation =
				Mat3(
					cos(rotor_input.blade_flapping[b_idx]), 0, sin(rotor_input.blade_flapping[b_idx]),
					0, 1, 0,
					-sin(rotor_input.blade_flapping[b_idx]), 0, cos(rotor_input.blade_flapping[b_idx])
				);

			foreach(pi; rotor.base_points[b_idx].byKeyValue) {
				vtkIdType id = pi.key;
				Vec3 point = pi.value;

				auto origin = rotor.origin;
				if(rotor_input.angular_velocity < 0) {
					origin[0] = -origin[0];
					origin[2] = -origin[2];
				}
				//auto pitch_rotated = pitch_rotation*point;
				//auto az_rotated = azimuth_rotation*pitch_rotated;
				//auto pitch_az_rot = azimuth_rotation*pitch_rotation;
				auto post_az_rot = azimuth_rotation*flap_rotation*pitch_rotation*point + origin;
				//auto post_az_rot = az_rotated + origin;
				auto final_p = aoa_rotation*post_az_rot;

				rotor.points.SetPoint(id, final_p[0], final_p[1], final_p[2]);
			}

			foreach(r_idx, loop; rotor.r_to_point_map[b_idx*elements..elements*(b_idx + 1)]) {
				auto actual_idx = r_idx + b_idx*elements;

				auto chunk_idx = r_idx/chunk_size;
				auto inner_idx = r_idx%chunk_size;

				foreach(l_idx, id; loop) {
					rotor.loads.SetTuple1(id, blade.chunks[chunk_idx].dC_T[inner_idx]);
					rotor.dC_L_dot.SetTuple1(id, blade.chunks[chunk_idx].dC_L_dot[inner_idx]);
					rotor.dC_T_dot.SetTuple1(id, blade.chunks[chunk_idx].dC_T_dot[inner_idx]);
					rotor.u_p.SetTuple1(id, blade.chunks[chunk_idx].u_p[inner_idx]);
					rotor.u_t.SetTuple1(id, blade.chunks[chunk_idx].u_t[inner_idx]);
					rotor.aoa.SetTuple1(id, blade.chunks[chunk_idx].aoa[inner_idx]);
					rotor.inflow_angle.SetTuple1(id, blade.chunks[chunk_idx].inflow_angle[inner_idx]);
					rotor.gamma.SetTuple1(id, blade.chunks[chunk_idx].gamma[inner_idx]);
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

VtkRotor build_base_vtu_rotor(RG)(auto ref RG rotor_geo) {
	version(Have_vtkd) {
		immutable elements = rotor_geo.blades[0].chunks.length*chunk_size;

		auto vtk_rotor = new VtkRotor(rotor_geo.blades.length);

		vtk_rotor.origin = rotor_geo.origin;

		vtk_rotor.loads.SetNumberOfComponents(1);
		vtk_rotor.loads.SetNumberOfTuples(elements*rotor_geo.blades.length*naca0012.length);
		vtk_rotor.loads.SetName("l");
		
		vtk_rotor.dC_L_dot.SetNumberOfComponents(1);
		vtk_rotor.dC_L_dot.SetNumberOfTuples(elements*rotor_geo.blades.length*naca0012.length);
		vtk_rotor.dC_L_dot.SetName("dC_L_dot");

		vtk_rotor.dC_T_dot.SetNumberOfComponents(1);
		vtk_rotor.dC_T_dot.SetNumberOfTuples(elements*rotor_geo.blades.length*naca0012.length);
		vtk_rotor.dC_T_dot.SetName("dC_T_dot");

		vtk_rotor.aoa.SetNumberOfComponents(1);
		vtk_rotor.aoa.SetNumberOfTuples(elements*rotor_geo.blades.length*naca0012.length);
		vtk_rotor.aoa.SetName("aoa");

		vtk_rotor.u_p.SetNumberOfComponents(1);
		vtk_rotor.u_p.SetNumberOfTuples(elements*rotor_geo.blades.length*naca0012.length);
		vtk_rotor.u_p.SetName("u_p");

		vtk_rotor.u_t.SetNumberOfComponents(1);
		vtk_rotor.u_t.SetNumberOfTuples(elements*rotor_geo.blades.length*naca0012.length);
		vtk_rotor.u_t.SetName("u_t");

		vtk_rotor.inflow_angle.SetNumberOfComponents(1);
		vtk_rotor.inflow_angle.SetNumberOfTuples(elements*rotor_geo.blades.length*naca0012.length);
		vtk_rotor.inflow_angle.SetName("inflow angle");

		vtk_rotor.gamma.SetNumberOfComponents(1);
		vtk_rotor.gamma.SetNumberOfTuples(elements*rotor_geo.blades.length*naca0012.length);
		vtk_rotor.gamma.SetName("gamma");


		foreach(b_idx, ref blade_geo; rotor_geo.blades) {

			vtk_rotor.azimuth_offsets[b_idx] = blade_geo.azimuth_offset;
			foreach(idx, ref chunk; blade_geo.chunks) {
				foreach(c_idx, ref r; chunk.r) {
					double x0 = r;

					double c0 = chunk.chord[c_idx];
					double t0 = chunk.twist[c_idx] + PI;

					vtk_rotor.r_to_point_map ~= new vtkIdType[naca0012.length];

					foreach(yz_idx, ref yz_point; naca0012) {
						auto yp0 = c0*(yz_point[0] - 0.25);// + chunk.xi[c_idx];
						auto zp0 = c0*yz_point[1];

						auto y0 = yp0*std.math.cos(t0) - zp0*std.math.sin(t0);
						auto z0 = yp0*std.math.sin(t0) + zp0*std.math.cos(t0);

						y0 -= chunk.xi[c_idx];
						auto p0 = Vec3(x0, y0, z0);
						
						auto id0 = vtk_rotor.points.InsertNextPoint(p0[0], p0[1], p0[2]);
						
						vtk_rotor.base_points[b_idx][id0] = p0;
						vtk_rotor.r_to_point_map[$-1][yz_idx] = id0;
						
					}
				}
			}

			foreach(r_idx, loop; vtk_rotor.r_to_point_map[b_idx*elements..elements*(b_idx + 1)]) {
				auto actual_idx = r_idx + b_idx*elements;
				if(r_idx < (elements - 1)) {
					auto next_loop = vtk_rotor.r_to_point_map[actual_idx + 1];

					foreach(l_idx, id; loop[0..$-1]) {

						vtkIdType[4] ids = [
							id,
							next_loop[l_idx],
							next_loop[l_idx + 1],
							loop[l_idx + 1]
						];
						vtk_rotor.grid.InsertNextCell(VTK__POLYGON, ids.length, ids.ptr);
					}
				}
			}
		}

		auto point_data = vtk_rotor.grid.GetPointData;
		point_data.AddArray(vtk_rotor.loads);
		point_data.AddArray(vtk_rotor.dC_L_dot);
		point_data.AddArray(vtk_rotor.dC_T_dot);
		point_data.AddArray(vtk_rotor.aoa);
		point_data.AddArray(vtk_rotor.u_p);
		point_data.AddArray(vtk_rotor.u_t);
		point_data.AddArray(vtk_rotor.inflow_angle);
		point_data.AddArray(vtk_rotor.gamma);

		return vtk_rotor;
	} else {
		auto vtk_rotor = new VtkRotor();
		return vtk_rotor;
	}
	
}

class VtkWake {
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
	}
}

VtkWake build_base_vtu_wake(W)(auto ref W wake) {
	
	version(Have_vtkd) {
		immutable wake_length = wake.rotor_wakes[0].tip_vortices[0].chunks.length*chunk_size;
		immutable shed_length = wake.rotor_wakes[0].shed_vortices[0].shed_filaments.length;
		immutable elements = wake.rotor_wakes[0].shed_vortices[0].shed_filaments[0].length*chunk_size;

		auto vtk_wake = new VtkWake(wake.rotor_wakes[0].tip_vortices.length, shed_length, elements);

		foreach(r_idx, rotor_wake; wake.rotor_wakes) {
			vtkIdType last_point_id;

			size_t shed_idx = 0;
			foreach(b_idx, shed_wake; rotor_wake.shed_vortices) {
				
				auto shed_wake_len = shed_wake.shed_filaments[0].chunks.length*chunk_size;
				
				foreach(shed_filament; shed_wake.shed_filaments) {

					vtk_wake.shed_grids[shed_idx].Allocate(shed_wake_len);
					
					vtk_wake.shed_points[shed_idx].Allocate(shed_wake_len);

					vtk_wake.shed_induced[shed_idx].SetNumberOfComponents(1);
					vtk_wake.shed_induced[shed_idx].SetNumberOfTuples(shed_wake_len);
					vtk_wake.shed_induced[shed_idx].SetName("induced");

					vtk_wake.shed_circ[shed_idx].SetNumberOfComponents(1);
					vtk_wake.shed_circ[shed_idx].SetNumberOfTuples(shed_wake_len - 1);
					vtk_wake.shed_circ[shed_idx].SetName("d_circulation");

					foreach(idx, shed_chunk; shed_filament.chunks) {
						foreach(c_idx; 0..chunk_size) {
							auto point_id = vtk_wake.shed_points[shed_idx].InsertNextPoint(0, 0, 0);

							vtk_wake.shed_point_ids ~= point_id;
							vtk_wake.shed_induced[shed_idx].SetTuple1(point_id, 0);
							if((idx > 0) || ((idx == 0) && (c_idx > 0))) {
								vtkIdType[2] ids = [last_point_id, point_id];

								auto cell_id = vtk_wake.shed_grids[shed_idx].InsertNextCell(VTK__POLY_LINE, ids.length, ids.ptr);
								vtk_wake.shed_cell_ids ~= cell_id;

								vtk_wake.shed_circ[shed_idx].SetTuple1(cell_id, 0);
							}

							last_point_id = point_id;
						}
					}

					vtk_wake.shed_grids[shed_idx].SetPoints(vtk_wake.shed_points[shed_idx]);

					auto point_data = vtk_wake.shed_grids[shed_idx].GetPointData;
					point_data.AddArray(vtk_wake.shed_induced[shed_idx]);

					auto cell_data = vtk_wake.shed_grids[shed_idx].GetCellData;
					cell_data.AddArray(vtk_wake.shed_circ[shed_idx]);

					vtk_wake.shed_appenders[b_idx].AddInputData(vtk_wake.shed_grids[shed_idx]);
					vtk_wake.shed_appenders[b_idx].Update;

					shed_idx++;
				}
			}

			foreach(b_idx, tip_vortex; rotor_wake.tip_vortices) {

				vtk_wake.tip_grids[b_idx].Allocate(wake_length);

				//auto tip_points = vtkPoints.New(VTK__DOUBLE);
				vtk_wake.tip_points[b_idx].Allocate(wake_length);

				vtk_wake.tip_induced[b_idx].SetNumberOfComponents(1);
				vtk_wake.tip_induced[b_idx].SetNumberOfTuples(wake_length);
				vtk_wake.tip_induced[b_idx].SetName("induced");

				vtk_wake.tip_d_volume[b_idx].SetNumberOfComponents(1);
				vtk_wake.tip_d_volume[b_idx].SetNumberOfTuples(wake_length - 1);
				vtk_wake.tip_d_volume[b_idx].SetName("d volume");

				vtk_wake.tip_core_size[b_idx].SetNumberOfComponents(1);
				vtk_wake.tip_core_size[b_idx].SetNumberOfTuples(wake_length - 1);
				vtk_wake.tip_core_size[b_idx].SetName("core size");

				vtk_wake.tip_circ[b_idx].SetNumberOfComponents(1);
				vtk_wake.tip_circ[b_idx].SetNumberOfTuples(wake_length - 1);
				vtk_wake.tip_circ[b_idx].SetName("circulation");

				foreach(idx, tip_chunk; tip_vortex.chunks) {
					foreach(c_idx; 0..chunk_size) {
						auto point_id =  vtk_wake.tip_points[b_idx].InsertNextPoint(0, 0, 0);

						vtk_wake.tip_point_ids ~= point_id;

						vtk_wake.tip_induced[b_idx].SetTuple1(point_id, 0);
						if((idx > 0) || ((idx == 0) && (c_idx > 0))) {
							vtkIdType[2] ids = [last_point_id, point_id];

							auto cell_id = vtk_wake.tip_grids[b_idx].InsertNextCell(VTK__POLY_LINE, ids.length, ids.ptr);

							vtk_wake.tip_cell_ids ~= cell_id;

							vtk_wake.tip_core_size[b_idx].SetTuple1(cell_id, 0);
							vtk_wake.tip_d_volume[b_idx].SetTuple1(cell_id, 0);
							vtk_wake.tip_circ[b_idx].SetTuple1(cell_id, 0);
						}

						last_point_id = point_id;
					}
				}

				vtk_wake.tip_grids[b_idx].SetPoints(vtk_wake.tip_points[b_idx]);

				auto point_data = vtk_wake.tip_grids[b_idx].GetPointData;
				point_data.AddArray(vtk_wake.tip_induced[b_idx]);

				auto cell_data = vtk_wake.tip_grids[b_idx].GetCellData;
				cell_data.AddArray(vtk_wake.tip_circ[b_idx]);
				cell_data.AddArray(vtk_wake.tip_core_size[b_idx]);
				cell_data.AddArray(vtk_wake.tip_d_volume[b_idx]);

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
		foreach(r_idx, rotor_wake; wake.rotor_wakes) {
			vtkIdType last_point_id;

			size_t shed_point_idx = 0;
			size_t shed_cell_idx = 0;
			size_t shed_idx = 0;
			foreach(b_idx, shed_wake; rotor_wake.shed_vortices) {
				
				vtk_wake.shed_appenders[b_idx].RemoveAllInputs;

				auto shed_wake_len = shed_wake.shed_filaments[0].chunks.length*chunk_size;
				
				foreach(shed_filament; shed_wake.shed_filaments) {

					double last_circ;
					double last_core_size;
					foreach(idx, shed_chunk; shed_filament.chunks) {
						foreach(c_idx; 0..chunk_size) {
							auto point_id = vtk_wake.shed_point_ids[shed_point_idx];

							vtk_wake.shed_points[shed_idx].SetPoint(point_id, shed_chunk.x[c_idx], shed_chunk.y[c_idx], shed_chunk.z[c_idx]);

							last_circ = shed_chunk.gamma[c_idx];

							vtk_wake.shed_induced[shed_idx].SetTuple1(point_id, shed_chunk.v_z[c_idx]);

							if((idx > 0) || ((idx == 0) && (c_idx > 0))) {
								vtkIdType[2] ids = [last_point_id, point_id];

								//auto cell_id = wake_grid.InsertNextCell(VTK__POLY_LINE, ids.length, ids.ptr);
								auto cell_id = vtk_wake.shed_cell_ids[shed_cell_idx];
								vtk_wake.shed_circ[shed_idx].SetTuple1(cell_id, last_circ);

								shed_cell_idx++;
							}

							last_point_id = point_id;

							shed_point_idx++;
						}
					}

					vtk_wake.shed_appenders[b_idx].AddInputData(vtk_wake.shed_grids[shed_idx]);
					vtk_wake.shed_appenders[b_idx].Update;

					shed_idx++;
				}

				import std.string : toStringz;
				auto filename = base_filename~"_shed_"~r_idx.to!string~"_"~b_idx.to!string~"_"~iteration.to!string~".vtu";
				vtk_wake.writer.SetFileName(filename.toStringz);
				vtk_wake.writer.SetInputData(vtk_wake.shed_appenders[b_idx].GetOutput);
				vtk_wake.writer.Write;
			}

			size_t tip_point_idx = 0;
			size_t tip_cell_idx = 0;

			foreach(b_idx, tip_vortex; rotor_wake.tip_vortices) {

				double last_circ;
				double last_core_size;
				double last_d_volume;
				foreach(idx, tip_chunk; tip_vortex.chunks) {
					foreach(c_idx; 0..chunk_size) {
						auto point_id = vtk_wake.tip_point_ids[tip_point_idx];
						vtk_wake.tip_points[b_idx].SetPoint(point_id, tip_chunk.x[c_idx], tip_chunk.y[c_idx], tip_chunk.z[c_idx]);

						last_circ = tip_chunk.gamma[c_idx];
						last_core_size = tip_chunk.r_c[c_idx];
						last_d_volume = tip_chunk.d_volume[c_idx];
						vtk_wake.tip_induced[b_idx].SetTuple1(point_id, tip_chunk.v_z[c_idx]);

						if((idx != tip_vortex.chunks.length - 1) || ((idx == tip_vortex.chunks.length - 1) && (c_idx < (chunk_size - 1)) )) {
						//if((idx > 0) || ((idx == 0) && (c_idx > 0))) {

							auto cell_id = vtk_wake.tip_cell_ids[tip_cell_idx];
							vtk_wake.tip_core_size[b_idx].SetTuple1(cell_id, last_core_size);
							vtk_wake.tip_d_volume[b_idx].SetTuple1(cell_id, last_d_volume);
							vtk_wake.tip_circ[b_idx].SetTuple1(cell_id, last_circ);

							tip_cell_idx++;
						}

						last_point_id = point_id;
						tip_point_idx++;
					}
				}

				import std.string : toStringz;
				auto filename = base_filename~"_"~r_idx.to!string~"_"~b_idx.to!string~"_"~iteration.to!string~".vtu";
				vtk_wake.writer.SetFileName(filename.toStringz);
				vtk_wake.writer.SetInputData(vtk_wake.tip_grids[b_idx]);
				vtk_wake.writer.Write;
			}
		}
	}
}

void write_inflow_vtu(I)(string filename, I[] inflows, Vec3 delta, Vec3 starts, size_t num_x, size_t num_y, size_t num_z, Vec3[] origins, double aoa, double omega) {

	version(Have_vtkd) {
		vtkIdType[PosVector] node_id_map;
		auto points = vtkPoints.New(VTK__DOUBLE);
		auto induced_writer = vtkXMLUnstructuredGridWriter.New;
		auto wake_writer = vtkXMLUnstructuredGridWriter.New;

		auto induced_grid = vtkUnstructuredGrid.New;

		points.Allocate(num_x*num_y*num_z);
		induced_grid.Allocate(num_x*num_y*num_z);

		auto induced_z = vtkDoubleArray.New();
		induced_z.SetNumberOfComponents(1);
		induced_z.SetNumberOfTuples(num_x*num_y*num_z);
		induced_z.SetName("induced");

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

					Chunk lambda = 0;

					foreach(r_idx, ref inflow; inflows) {
						immutable Chunk x_rel = x[];// - origins[r_idx][0];
						immutable Chunk y_rel = -(y[] - origins[r_idx][1]);
						immutable Chunk z_rel = z[];// - origins[r_idx][2];

						//immutable Chunk x_rel = x[] - origin[0];
						//immutable Chunk y_rel = y[] - origin[1];
						//immutable Chunk z_rel = z[] - origin[2];

						immutable Chunk x_m = x_rel[]*cos(aoa) - z_rel[]*sin(aoa) - origins[r_idx][0];
						immutable Chunk z_m = -x_rel[]*sin(aoa) - z_rel[]*cos(aoa) + origins[r_idx][2];

						auto infl = inflow.inflow_at(x_m, y_rel, z_m, x_e, aoa);

						lambda[] += infl[]*cos(aoa);
					}

					foreach(c_idx; 0..chunk_size) {
						induced_z.SetTuple1(ids[c_idx], -omega*lambda[c_idx]);
					}
				}
			}
			//writeln("z_idx: ", z_idx);
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
			//writeln("z_idx: ", z);
		}

		writeln("Done generating grid");

		induced_grid.SetPoints(points);

		auto point_data = induced_grid.GetPointData;
		point_data.SetScalars(induced_z);

		import std.string : toStringz;

		induced_writer.SetFileName(filename.toStringz);
		induced_writer.SetInputData(induced_grid);
		induced_writer.Write;

		writeln("Done writing induced flow field");
	}
}