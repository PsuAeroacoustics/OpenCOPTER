module opencopter.io;

import opencopter.aircraft;
import opencopter.airfoilmodels;
import opencopter.math;
import opencopter.memory;
import opencopter.wake;

import numd.linearalgebra.matrix : Vector;

import std.algorithm;
import std.array;
import std.conv;
import std.exception : enforce;
import std.file : readText;
import std.math;
import std.range;
import std.stdio;
import std.typecons	: Tuple, tuple;

import kxml.xml;

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

enum PlanarSymetry {
	none = 0,
	xy = 1,
	xz = 2,
	xy_xz = 3,
	yz = 4,
	xy_yz = 5,
	xz_yz = 6,
	xy_xz_yz = 7
}

immutable Mat4[PlanarSymetry] symmetry_mats;
immutable InterpMap interp_map;
shared static this() {
	symmetry_mats = [
		PlanarSymetry.xy: Mat4(
			1, 0, 0, 0,
			0, 1, 0, 0,
			0, 0, -1, 0,
			0, 0, 0, 1),
		PlanarSymetry.xz: Mat4(
			1, 0, 0, 0,
			0, -1, 0, 0,
			0, 0, 1, 0,
			0, 0, 0, 1),
		PlanarSymetry.yz: Mat4(
			-1, 0, 0, 0,
			0, 1, 0, 0,
			0, 0, 1, 0,
			0, 0, 0, 1)
	];

	interp_map[CurveType.linear] = &linear;
	interp_map[CurveType.spline_pchip] = &spline_pchip;
	interp_map[CurveType.cubic_bezier] = &cubic_bezier;
	interp_map[CurveType.cubic_bezier_approx] = &cubic_bezier_approx;
}

alias VspFrameData = Tuple!(Frame*, "frame", string, "parent_id", string[], "child_ids", PlanarSymetry, "symmetry", int, "symmetry_parent", XmlNode, "xml_node");
alias InterpFunc = double[] function(double[], double[], double[]);
alias InterpMap = InterpFunc[CurveType];


enum CurveType {
	linear = 0,
	spline_pchip = 1,
	cubic_bezier = 2,
	cubic_bezier_approx = 3
}

double[] cubic_bezier(double[] control_points, double[] control_values, double[] radial_points) {
	double[] values = new double[radial_points.length];

	foreach(idx, ref r; radial_points) {

		size_t max_idx = size_t.max;
		for(size_t c_idx = 3; c_idx < control_points.length; c_idx += 3) {
			if(r < control_points[c_idx]) {
				max_idx = c_idx;
				break;
			}
		}

		immutable double min_control_point = control_points[max_idx - 3];
		immutable double max_control_point = control_points[max_idx];
		immutable double interval = max_control_point - min_control_point;

		immutable double t = (r - min_control_point)/interval;

		immutable double p0 = control_values[max_idx - 3];
		immutable double p1 = control_values[max_idx - 2];
		immutable double p2 = control_values[max_idx - 1];
		immutable double p3 = control_values[max_idx];

		values[idx] = (1.0 - t)^^3.0*p0 + 3.0*(1.0 - t)^^2.0*t*p1 + 3.0*(1.0 - t)*t^^2.0*p2 + t^^3.0*p3;
	}

	return values;
}

double[] spline_pchip(double[] control_points, double[] control_values, double[] radial_points) {
	double[] values = new double[radial_points.length];

	// double[] m = new double[control_points.length];

	// foreach(m_idx, ref _m; m) {
	// 	// if(m_idx == 0) {
	// 	// 	_m = (control_values[m_idx + 1] - control_values[m_idx])/(control_points[m_idx + 1] - control_points[m_idx]);
	// 	// } else/+ if(m_idx == m.length - 1)+/ {
	// 	// 	_m = (control_values[m_idx] - control_values[m_idx - 1])/(control_points[m_idx] - control_points[m_idx - 1]);
	// 	// }/+ else {
	// 	// 	_m = 0.5*(
	// 	// 		(control_values[m_idx] - control_values[m_idx - 1])/(control_points[m_idx] - control_points[m_idx - 1]) + 
	// 	// 		(control_values[m_idx] - control_values[m_idx + 1])/(control_points[m_idx] - control_points[m_idx + 1])
	// 	// 	);
	// 	// }+/

	// 	if(m_idx == m.length - 1) {
	// 		_m = (control_values[m_idx] - control_values[m_idx - 1])/(control_points[m_idx] - control_points[m_idx - 1]);
	// 	} else {
	// 		_m = (control_values[m_idx + 1] - control_values[m_idx])/(control_points[m_idx + 1] - control_points[m_idx]);
	// 	}
	// }

	// writeln("m = ", m, ";");
	// writeln("control_points = ", control_points, ";");
	// writeln("control_values = ", control_values, ";");

	// foreach(idx, ref r; radial_points) {

	// 	size_t max_idx = size_t.max;
	// 	for(size_t c_idx = 0; c_idx < control_points.length; c_idx++) {
	// 		if(r < control_points[c_idx]) {
	// 			max_idx = c_idx;
	// 			break;
	// 		}
	// 	}

	// 	immutable double min_control_point = control_points[max_idx - 1];
	// 	immutable double max_control_point = control_points[max_idx];
	// 	immutable double interval = max_control_point - min_control_point;

	// 	immutable double p0 = control_values[max_idx - 1];
	// 	immutable double p1 = control_values[max_idx];

	// 	immutable double m0 = m[max_idx - 1];
	// 	immutable double m1 = m[max_idx];

	// 	immutable double t = (r - min_control_point)/interval;

	// 	values[idx] = 
	// 		(2.0*t^^3 - 3.0*t^^2.0 + 1.0)*p0 + 
	// 		(t^^3 - 2.0*t^^2.0 + t)*m0 + 
	// 		(-2.0*t^^3 + 3.0*t^^2.0)*p1 + 
	// 		(t^^3 - t^^2.0)*m1;
	// }
	enforce(false, "Spline (PCHIP) not supported by OpenCOPTER");
	return values;
}

double[] linear(double[] control_points, double[] control_values, double[] radial_points) {
	double[] values = new double[radial_points.length];

	foreach(idx, ref r; radial_points) {

		size_t max_idx = size_t.max;
		for(size_t c_idx = 1; c_idx < control_points.length; c_idx++) {
			if(r < control_points[c_idx]) {
				max_idx = c_idx;
				break;
			}
		}

		immutable double min_control_point = control_points[max_idx - 1];
		immutable double max_control_point = control_points[max_idx];
		immutable double interval = max_control_point - min_control_point;

		immutable double p1 = control_values[max_idx];
		immutable double p0 = control_values[max_idx - 1];

		immutable double m = (p1 - p0)/interval;
		immutable double b = p1 - m*max_control_point;

		values[idx] = m*r + b;
	}

	return values;
}

double[] cubic_bezier_approx(double[], double[], double[]) {
	enforce(false, "Approximate cubic bezier not supported by OpenCOPTER");
	return new double[0];
}

AircraftT!AC create_aircraft_from_vsp(ArrayContainer AC)(string filename, size_t elements = 48) {

	AircraftT!AC ac;

	auto vsp_lines = readText(filename);

	auto root_node = vsp_lines.readDocument;

	auto file_version_nodes = root_node.parseXPath("/Vsp_Geometry/Version");
	enforce(file_version_nodes.length == 1, "Incorrect number of Vsp_Geometry/Version nodes. Expected 1 got "~file_version_nodes.length.to!string);
	auto file_version = file_version_nodes[0].getCData.to!int;

	enforce(file_version == 5, "Unsupported vsp3 file version. Version 5 supported, file is version "~file_version.to!string);

	auto search_list = root_node.parseXPath("/Vsp_Geometry/Vehicle/Geom");

	VspFrameData[string] geom_dict;

	foreach(geom; search_list) {

		auto parm_nodes = geom.parseXPath("ParmContainer");
		auto geom_base_nodes = geom.parseXPath("GeomBase");
		enforce(parm_nodes.length == 1, "Incorrect number of ParmContainer nodes. Expected 1 got "~parm_nodes.length.to!string);
		auto parm_node = parm_nodes[0];

		enforce(geom_base_nodes.length == 1, "Incorrect number of GeomBase nodes. Expected 1 got "~geom_base_nodes.length.to!string);
		auto geom_base_node = geom_base_nodes[0];

		auto geom_id_nodes = parm_node.parseXPath("ID");
		enforce(geom_id_nodes.length == 1, "Incorrect number of Geom ID nodes. Expected 1 got "~geom_id_nodes.length.to!string);
		string geom_id = geom_id_nodes[0].getCData;

		auto name_nodes = parm_node.parseXPath("Name");
		enforce(name_nodes.length == 1, "Incorrect number of Geom ID nodes. Expected 1 got "~name_nodes.length.to!string);
		string name = name_nodes[0].getCData;

		auto sym_nodes = parm_node.parseXPath("Sym");
		enforce(sym_nodes.length == 1, "Incorrect number of Sym nodes. Expected 1 got "~sym_nodes.length.to!string);
		auto sym_node = sym_nodes[0];

		auto axial_sym_nodes = sym_node.parseXPath("Sym_Axial_Flag");
		enforce(axial_sym_nodes.length == 1, "Incorrect number of Sym_Axial_Flag nodes. Expected 1 got "~axial_sym_nodes.length.to!string);
		auto axial_sym = axial_sym_nodes[0].getAttribute("Value").to!double.to!int;
		enforce(axial_sym == 0, "OpenCOPTER currently does not support axial symmetry");

		auto sym_ancestor_nodes = sym_node.parseXPath("Sym_Ancestor");
		enforce(sym_ancestor_nodes.length == 1, "Incorrect number of Sym_Ancestor nodes. Expected 1 got "~sym_ancestor_nodes.length.to!string);
		auto sym_ancestor = sym_ancestor_nodes[0].getAttribute("Value").to!double.to!int;

		auto planar_sym_nodes = sym_node.parseXPath("Sym_Planar_Flag");
		enforce(planar_sym_nodes.length == 1, "Incorrect number of Sym_Planar_Flag nodes. Expected 1 got "~planar_sym_nodes.length.to!string);
		auto planar_sym = planar_sym_nodes[0].getAttribute("Value").to!double.to!int.to!PlanarSymetry;
		enforce(
			(planar_sym == PlanarSymetry.xy) || (planar_sym == PlanarSymetry.xz) || (planar_sym == PlanarSymetry.yz) || (planar_sym == PlanarSymetry.none),
			"OpenCOPTER currently supports only a single plane of symmetry. Component: "~name
		);

		auto origin_flag_nodes = sym_node.parseXPath("Sym_Ancestor_Origin_Flag");
		enforce(origin_flag_nodes.length == 1, "Incorrect number of Sym_Ancestor_Origin_Flag nodes. Expected 1 got "~origin_flag_nodes.length.to!string);
		auto origin_flag = origin_flag_nodes[0].getAttribute("Value").to!double.to!bool;

		enforce(!origin_flag || (planar_sym == PlanarSymetry.none), "XForm:Symmetry must be set to \"Object\" not \"Attach\". Component: "~name);

		auto translation_attach_nodes = parm_node.parseXPath("Attach/Trans_Attach_Flag");
		enforce(translation_attach_nodes.length == 1, "Incorrect number of Attach/Trans_Attach_Flag nodes. Expected 1 got "~translation_attach_nodes.length.to!string);
		auto translation_attach_flag = translation_attach_nodes[0].getAttribute("Value").to!double.to!int;

		auto rotation_attach_nodes = parm_node.parseXPath("Attach/Rots_Attach_Flag");
		enforce(rotation_attach_nodes.length == 1, "Incorrect number of Attach/Rots_Attach_Flag nodes. Expected 1 got "~rotation_attach_nodes.length.to!string);
		auto rotation_attach_flag = rotation_attach_nodes[0].getAttribute("Value").to!double.to!int;

		auto xform_nodes = parm_node.parseXPath("XForm");
		enforce(xform_nodes.length == 1, "Incorrect number of XForm nodes. Expected 1 got "~xform_nodes.length.to!string);
		auto xform_node = xform_nodes[0];

		auto x_rel_loc_nodes = xform_node.parseXPath("X_Rel_Location");
		enforce(x_rel_loc_nodes.length == 1, "Incorrect number of X_Rel_Location nodes. Expected 1 got "~x_rel_loc_nodes.length.to!string);

		auto y_rel_loc_nodes = xform_node.parseXPath("Y_Rel_Location");
		enforce(y_rel_loc_nodes.length == 1, "Incorrect number of Y_Rel_Location nodes. Expected 1 got "~y_rel_loc_nodes.length.to!string);

		auto z_rel_loc_nodes = xform_node.parseXPath("Z_Rel_Location");
		enforce(z_rel_loc_nodes.length == 1, "Incorrect number of Z_Rel_Location nodes. Expected 1 got "~z_rel_loc_nodes.length.to!string);

		auto x_rel_rot_nodes = xform_node.parseXPath("X_Rel_Rotation");
		enforce(x_rel_rot_nodes.length == 1, "Incorrect number of X_Rel_Rotation nodes. Expected 1 got "~x_rel_rot_nodes.length.to!string);
		
		auto y_rel_rot_nodes = xform_node.parseXPath("Y_Rel_Rotation");
		enforce(y_rel_rot_nodes.length == 1, "Incorrect number of Y_Rel_Rotation nodes. Expected 1 got "~y_rel_rot_nodes.length.to!string);
		
		auto z_rel_rot_nodes = xform_node.parseXPath("Z_Rel_Rotation");
		enforce(z_rel_rot_nodes.length == 1, "Incorrect number of Z_Rel_Rotation nodes. Expected 1 got "~z_rel_rot_nodes.length.to!string);

		auto position = Vec3(
			x_rel_loc_nodes[0].getAttribute("Value").to!double,
			y_rel_loc_nodes[0].getAttribute("Value").to!double,
			z_rel_loc_nodes[0].getAttribute("Value").to!double
		);

		auto children_id_nodes = geom_base_node.parseXPath("Child_List/Child/ID");

		string[] children_ids = children_id_nodes.map!(n => n.getCData).array;

		auto type_name_nodes = geom_base_node.parseXPath("TypeName");
		enforce(type_name_nodes.length == 1, "Incorrect number of TypeName nodes. Expected 1 got "~type_name_nodes.length.to!string);
		string type_name = type_name_nodes[0].getCData;

		auto parent_id_nodes = geom_base_node.parseXPath("ParentID");
		enforce(parent_id_nodes.length == 1, "Incorrect number of ParentID nodes. Expected 1 got "~parent_id_nodes.length.to!string);
		string parent_id = parent_id_nodes[0].getCData;

		if(parent_id != "NONE") {
			enforce(
				(translation_attach_flag == 1),
				"Settings XForm:Attach to Parent:Translate must be set to \"Comp\" when component is a child of another component. Component: "~name
			);
			enforce(
				(rotation_attach_flag == 1),
				"Settings XForm:Attach to Parent:Rotate must be set to \"Comp\" when component is a child of another component. Component: "~name
			);
		}
		
		FrameType frame_type;

		if(type_name != "Propeller") {
			frame_type = FrameType.connection;
		} else {
			frame_type = FrameType.rotor;
		}

		auto frame = new Frame(Vec3(1, 0, 0), 0, position, null, name, frame_type);

		if(frame_type == FrameType.rotor) {
			frame.rotate(Vec3(1, 0, 0), x_rel_rot_nodes[0].getAttribute("Value").to!double*(PI/180.0));
			frame.rotate(Vec3(0, 1, 0), y_rel_rot_nodes[0].getAttribute("Value").to!double*(PI/180.0) - PI/2.0);
			frame.rotate(Vec3(0, 0, 1), z_rel_rot_nodes[0].getAttribute("Value").to!double*(PI/180.0));
		} else {
			frame.rotate(Vec3(1, 0, 0), x_rel_rot_nodes[0].getAttribute("Value").to!double*(PI/180.0));
			frame.rotate(Vec3(0, 1, 0), y_rel_rot_nodes[0].getAttribute("Value").to!double*(PI/180.0));
			frame.rotate(Vec3(0, 0, 1), z_rel_rot_nodes[0].getAttribute("Value").to!double*(PI/180.0));
		}

		geom_dict[geom_id] = VspFrameData(frame, parent_id, children_ids, planar_sym, sym_ancestor, geom);
	}

	auto root_components = geom_dict.byPair.filter!(a => a.value.parent_id == "NONE").assocArray;

	RotorGeometryT!AC*[] build_oc_frame(VspFrameData* frame_data, VspFrameData* parent_frame_data, bool symmetry_applied, bool symmetry_parent) {
		
		RotorGeometryT!AC*[] rotors;

		if(parent_frame_data != null) {
			writeln("Setting ", parent_frame_data.frame.name, " parent to ", frame_data.frame.name);
			frame_data.frame.parent = parent_frame_data.frame;
			writeln("parent_frame_data.frame.children: ", parent_frame_data.frame.children);
		}

		foreach(child_id; frame_data.child_ids) {
			if((frame_data.symmetry != PlanarSymetry.none) && !symmetry_applied) {
				auto frame = geom_dict[child_id].frame;

				auto sym_frame = new Frame(frame.axis, frame.angle, frame.local_position, frame.parent, frame.name~"_symmetry", frame.frame_type);
				sym_frame.local_matrix = symmetry_mats[frame_data.symmetry]*sym_frame.local_matrix;

				frame_data.frame.children ~= [frame, sym_frame];

				auto child1_frame_data = new VspFrameData(
					geom_dict[child_id].frame,
					geom_dict[child_id].parent_id,
					geom_dict[child_id].child_ids,
					geom_dict[child_id].symmetry,
					geom_dict[child_id].symmetry_parent,
					geom_dict[child_id].xml_node
				);

				auto child2_frame_data = new VspFrameData(
					sym_frame,
					geom_dict[child_id].parent_id,
					geom_dict[child_id].child_ids,
					geom_dict[child_id].symmetry,
					geom_dict[child_id].symmetry_parent,
					geom_dict[child_id].xml_node
				);

				rotors ~= build_oc_frame(child1_frame_data, frame_data, true, false);
				rotors ~= build_oc_frame(child2_frame_data, frame_data, true, true);

			} else {
				writeln("geom_dict[child_id].frame.name: ", geom_dict[child_id].frame.name);
				
				string postfix = symmetry_parent ? "_symmetry" : "";

				auto new_frame = new VspFrameData(
					new Frame(geom_dict[child_id].frame, geom_dict[child_id].frame.name~postfix),
					geom_dict[child_id].parent_id,
					geom_dict[child_id].child_ids,
					geom_dict[child_id].symmetry,
					geom_dict[child_id].symmetry_parent,
					geom_dict[child_id].xml_node
				);
				
				frame_data.frame.children ~= new_frame.frame;

				rotors ~= build_oc_frame(new_frame, frame_data, symmetry_applied, symmetry_parent);
			}
		}

		if(frame_data.frame.frame_type == FrameType.rotor) {
			// add blades
			auto num_blades = frame_data.xml_node.parseXPath("ParmContainer/Design/NumBlade")[0].getAttribute("Value").to!double.to!size_t;
			auto R = 0.5*frame_data.xml_node.parseXPath("ParmContainer/Design/Diameter")[0].getAttribute("Value").to!double;
			auto c_ave = frame_data.xml_node.parseXPath("ParmContainer/Design/Chord")[0].getAttribute("Value").to!double;
			auto solidity = frame_data.xml_node.parseXPath("ParmContainer/Design/Solidity")[0].getAttribute("Value").to!double;

			auto r_c = frame_data.xml_node.parseXPath("PropellerGeom/Chord/ParmContainer/Chord/r_0")[0].getAttribute("Value").to!double;
			auto r = generate_radius_points(elements, r_c);

			auto chord_curve_points = frame_data.xml_node.parseXPath("PropellerGeom/Chord/PCurve/NumPts")[0].getCData().to!double.to!size_t;
			auto chord_curve_type = frame_data.xml_node.parseXPath("PropellerGeom/Chord/ParmContainer/Chord/CrvType")[0].getAttribute("Value").to!double.to!int.to!CurveType;

			double[] chord_values = 
				iota(0, chord_curve_points).
				map!(a => frame_data.xml_node.parseXPath("PropellerGeom/Chord/ParmContainer/Chord/crd_"~a.to!string)[0].getAttribute("Value").to!double).
				array;

			double[] chord_points = 
				iota(0, chord_curve_points).
				map!(a => frame_data.xml_node.parseXPath("PropellerGeom/Chord/ParmContainer/Chord/r_"~a.to!string)[0].getAttribute("Value").to!double).
				array;

			double[] chord = interp_map[chord_curve_type](chord_points, chord_values, r);

			auto twist_curve_points = frame_data.xml_node.parseXPath("PropellerGeom/Twist/PCurve/NumPts")[0].getCData().to!double.to!size_t;
			auto twist_curve_type = frame_data.xml_node.parseXPath("PropellerGeom/Twist/ParmContainer/Twist/CrvType")[0].getAttribute("Value").to!double.to!int.to!CurveType;

			double[] twist_values = 
				iota(0, twist_curve_points).
				map!(a => frame_data.xml_node.parseXPath("PropellerGeom/Twist/ParmContainer/Twist/tw_"~a.to!string)[0].getAttribute("Value").to!double*(PI/180.0).to!double).
				array;

			double[] twist_points = 
				iota(0, twist_curve_points).
				map!(a => frame_data.xml_node.parseXPath("PropellerGeom/Twist/ParmContainer/Twist/r_"~a.to!string)[0].getAttribute("Value").to!double).
				array;

			double[] twist = interp_map[twist_curve_type](twist_points, twist_values, r);

			auto sweep_curve_points = frame_data.xml_node.parseXPath("PropellerGeom/Sweep/PCurve/NumPts")[0].getCData().to!double.to!size_t;
			auto sweep_curve_type = frame_data.xml_node.parseXPath("PropellerGeom/Sweep/ParmContainer/Sweep/CrvType")[0].getAttribute("Value").to!double.to!int.to!CurveType;

			double[] sweep_values = 
				iota(0, sweep_curve_points).
				map!(a => -frame_data.xml_node.parseXPath("PropellerGeom/Sweep/ParmContainer/Sweep/sw_"~a.to!string)[0].getAttribute("Value").to!double*(PI/180.0).to!double).
				array;

			double[] sweep_points = 
				iota(0, sweep_curve_points).
				map!(a => frame_data.xml_node.parseXPath("PropellerGeom/Sweep/ParmContainer/Sweep/r_"~a.to!string)[0].getAttribute("Value").to!double).
				array;

			double[] sweep = interp_map[sweep_curve_type](sweep_points, sweep_values, r);
			double[] xi_p = sweep.map!(s => tan(s)).array;
			writeln("r = ", r, ";");
			writeln("xi_p = ", xi_p, ";");
			writeln("sweep = ", sweep, ";");
			writeln("twist = ", twist, ";");
			writeln("chord = ", chord, ";");

			double[] r_edges = new double[elements + 1];

			r_edges[0] = r_c;

			foreach(r_idx; 0..elements) {
				immutable d = r[r_idx] - r_edges[r_idx];
				r_edges[r_idx + 1] = r[r_idx] + d;
			}

			double[] r_delta = iota(0, elements).map!(i => r_edges[i + 1] - r_edges[i]).array;

			double[] xi = new double[elements];

			xi[] = r_delta[]*xi_p[];
			foreach(x_idx, ref x; xi) {
				if(x_idx == 0) {
					x = 0;//xi_p[x_idx]*r_delta[x_idx];
				} else {
					x = xi_p[x_idx]*(r[x_idx] - r[x_idx - 1]) + xi[x_idx - 1];
				}
				// auto tmp = new double[x_idx];
				// tmp[] = r_delta[0..x_idx]*xi_p[0..x_idx];
				// x = sum(tmp);
			}

			writeln("xi = ", xi, ";");

			auto rotor = new RotorGeometryT!AC(
				num_blades,
				Vec3(0),
				R,
				solidity
			);

			auto fixed_frame = frame_data.frame;
			auto rotor_name = fixed_frame.name;
			fixed_frame.name = fixed_frame.name~"_fixed";
			fixed_frame.frame_type = FrameType.connection;
			rotor.frame = new Frame(Vec3(1, 0, 0), 0, Vec3(0), fixed_frame, rotor_name, FrameType.rotor);

			fixed_frame.children ~= rotor.frame;

			foreach(b_idx; 0..num_blades) {
				double azimuth = (b_idx.to!double/num_blades.to!double)*2.0*PI;

				auto blade_frame = new Frame(Vec3(0, 0, 1), azimuth, Vec3(0), rotor.frame, rotor.frame.name~"_blade_"~b_idx.to!string, FrameType.blade);

				auto af_model = new ThinAirfoil(0.0);
				size_t[2] extent = [0, elements - 1];

				auto blade_airfoil = new BladeAirfoil([af_model], [extent]);

				rotor.blades[b_idx] = BladeGeometryT!AC(
					elements,
					azimuth,
					c_ave,
					blade_airfoil,
					0
				);

				rotor.blades[b_idx].frame = blade_frame;

				rotor.blades[b_idx].set_geometry_array!"r"(r);
				rotor.blades[b_idx].set_geometry_array!"twist"(twist);
				rotor.blades[b_idx].set_geometry_array!"chord"(chord);
				rotor.blades[b_idx].set_geometry_array!"sweep"(sweep);
				rotor.blades[b_idx].set_geometry_array!"xi"(xi);
				rotor.blades[b_idx].set_geometry_array!"xi_p"(xi_p);
				rotor.blades[b_idx].blade_length = R*(1.0 - r_c);
				compute_blade_vectors(rotor.blades[b_idx]);
				
				writeln("blade_geom.chunks.length: ", rotor.blades[b_idx].chunks.length);

				writeln("rotor.blades[b_idx].chunks.length: ", rotor.blades[b_idx].chunks.length);

				rotor.frame.children.length++;
				rotor.frame.children[$-1] = blade_frame;
			}

			rotors ~= rotor;
		}

		writeln("rotors[$-1].blades[$-1].chunks.length: ", rotors[$-1].blades[$-1].chunks.length);
		return rotors;
	}

	ac.root_frame = new Frame(Vec3(0, 0, 1), 0, Vec3(0, 0, 0), null, "aircraft", "connection");

	RotorGeometryT!AC*[] rotors;
	foreach(comp_id, ref root_component; root_components) {
		ac.root_frame.children ~= root_component.frame;
		root_component.frame.parent = ac.root_frame;
		rotors ~= build_oc_frame(&root_component, null, false, false);
	}

	writeln("rotors.length: ", rotors.length);
	writeln("rotors[$-1].blades[$-1].chunks.length: ", rotors[$-1].blades[$-1].chunks.length);

	ac.rotors = new RotorGeometryT!AC[rotors.length];
	foreach(r_idx; 0..rotors.length) {
		ac.rotors[r_idx] = *rotors[r_idx];
	}

	ac.root_frame.update(Mat4.identity);

	return ac;
}

unittest {
	auto ac = create_aircraft_from_vsp!(ArrayContainer.none)("./oc_fly/example/test_evtol.vsp3");



	print_frame(ac.root_frame);


	import opencopter.vtk;
	import opencopter.aircraft.state;

// alias AircraftState = AircraftStateT!(ArrayContainer.none);

//  struct AircraftStateT(ArrayContainer _AC) {
// 	alias AC = _AC;
// 	mixin ArrayDeclMixin!(AC, RotorStateT!(AC), "rotor_states");

// 	Vec4 freestream;

// 	this(size_t num_rotors, size_t num_blades, size_t num_elements, ref AircraftT!AC ac) {

	writeln("ac.rotors[0].blades[0].chunks.length: ", ac.rotors[0].blades[0].chunks.length);

	auto ac_state = AircraftState(ac.rotors.length, ac.rotors.map!(r => r.blades.length).array, 48, ac);

	foreach(rotor_idx; 0..ac.rotors.length) {
	
	
		foreach(blade_idx, ref blade; ac.rotors[rotor_idx].blades) {

			auto local_blade_pos = Vector!(4, Chunk)(0);
			
			foreach(chunk_idx, ref state_chunk; ac_state.rotor_states[rotor_idx].blade_states[blade_idx].chunks) {

				immutable Chunk adjusted_r = blade.chunks[chunk_idx].r[] - blade.r_c;
				local_blade_pos[0][] = adjusted_r[]*ac.rotors[rotor_idx].radius;
				local_blade_pos[1][] = blade.chunks[chunk_idx].xi[]*ac.rotors[rotor_idx].radius;
				local_blade_pos[1][] = 0;
				local_blade_pos[3][] = 1;

				auto global_blade_pos = blade.frame.global_matrix * local_blade_pos;

				state_chunk.x[] = global_blade_pos[0][];
				state_chunk.y[] = global_blade_pos[1][];
				state_chunk.z[] = global_blade_pos[2][];
			}
		}
	}

	foreach(r_idx, ref rotor; ac.rotors) {
		auto vtk_rotor = build_base_vtu_rotor(rotor);
		write_rotor_vtu("vsp_rotor_"~r_idx.to!string, 0, r_idx, vtk_rotor, ac_state.rotor_states[r_idx], rotor);
	}
}