import opencopter.aircraft;
import opencopter.config;
import opencopter.atmosphere;
import opencopter.bladeelement;
import opencopter.inflow;
import opencopter.liftmodels;
import opencopter.math;
import opencopter.memory;
import opencopter.trim;
import opencopter.vtk;
import opencopter.wake;

import plt = matplotlibd.pyplot;

import std.algorithm;
import std.array;
import std.conv;
import std.exception;
import std.getopt;
import std.math;
import std.range;
import std.stdio;

//import mpid;
import wopwopd;

import numd.linearalgebra.matrix : Matrix, Vector;
import numd.utility : linspace, meshgrid;

enum double RAD_TO_HZ = 0.1591549;

double[] generate_radius_points(size_t n_sections) {
    /+auto i = linspace(0., (n_sections - 1).to!double, n_sections);
	auto r_array = new double[n_sections];
    r_array[] = (1./n_sections)*(i[] + 0.5);+/
    //return r_array;

    //return iota(1.0, n_sections + 1.0).map!(n => 0.5*(cos(n*PI/(n_sections.to!double + 1.0)) + 1.0).to!double).retro.array;
    return iota(1.0, n_sections + 1.0).map!((n) {
    	immutable psi = n*PI/(n_sections.to!double + 1.0);
    	auto r = 0.5*(cos(psi) + 1.0).to!double;
    	import std.stdio : writeln;
    	//writeln("n: ", n);
    	//writeln("psi: ", psi);
    	//writeln("r: ", r);
    	//writeln;

    	return r;
    }).retro.array;
}

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

auto generate_blade_geom(double[] radial_stations, double[] twist, double radius, double[] chord) {
	// radial_stations are in the y direction and the airfoil points are in the x-z plane.

	size_t num_nodes = radial_stations.length*naca0012.length;
	auto geom = ConstantGeometryData(num_nodes);

	size_t node_idx = 0;
	foreach(p_idx, p; naca0012) {
		foreach(r_idx, rs; radial_stations) {
			static import std.math;
			immutable xp = p[0]*chord[r_idx];
			immutable zp = p[1]*chord[r_idx];
			geom.y_nodes[node_idx] = xp*std.math.cos(twist[r_idx]) - zp*std.math.sin(twist[r_idx]);
			geom.x_nodes[node_idx] = rs*radius;
			geom.z_nodes[node_idx] = xp*std.math.sin(twist[r_idx]) + zp*std.math.cos(twist[r_idx]);

			node_idx++;
		}
	}

	node_idx = 0;
	foreach(p_idx, p; naca0012) {
		foreach(r_idx, rs; radial_stations) {
			immutable n1 = p_idx*radial_stations.length + r_idx - 1;
			immutable n2 = p_idx*radial_stations.length + r_idx + 1;
			immutable n3 = (p_idx - 1)*radial_stations.length + r_idx;
			immutable n4 = (p_idx + 1)*radial_stations.length + r_idx;

			if(r_idx == 0) {
				if(p_idx == 0) {
					// only use n2, n4, and node_idx

					immutable n2n = Vec3(
						geom.x_nodes[n2] - geom.x_nodes[node_idx],
						geom.y_nodes[n2] - geom.y_nodes[node_idx],
						geom.z_nodes[n2] - geom.z_nodes[node_idx]
					);

					immutable n4n = Vec3(
						geom.x_nodes[n4] - geom.x_nodes[node_idx],
						geom.y_nodes[n4] - geom.y_nodes[node_idx],
						geom.z_nodes[n4] - geom.z_nodes[node_idx]
					);

					immutable normal = n2n.cross(n4n).normalize;

					geom.x_normals[node_idx] = normal[0];
					geom.y_normals[node_idx] = normal[1];
					geom.z_normals[node_idx] = normal[2];

				} else if(p_idx == naca0012.length - 1) {
					// only use n1, n4, and node_idx

					immutable n1n = Vec3(
						geom.x_nodes[n3] - geom.x_nodes[node_idx],
						geom.y_nodes[n3] - geom.y_nodes[node_idx],
						geom.z_nodes[n3] - geom.z_nodes[node_idx]
					);

					immutable n4n = Vec3(
						geom.x_nodes[n2] - geom.x_nodes[node_idx],
						geom.y_nodes[n2] - geom.y_nodes[node_idx],
						geom.z_nodes[n2] - geom.z_nodes[node_idx]
					);

					immutable normal = n4n.cross(n1n).normalize;

					geom.x_normals[node_idx] = normal[0];
					geom.y_normals[node_idx] = normal[1];
					geom.z_normals[node_idx] = normal[2];

				} else {
					// only use n1, n2, n4, and node_idx

					immutable n1n = Vec3(
						geom.x_nodes[n1] - geom.x_nodes[node_idx],
						geom.y_nodes[n1] - geom.y_nodes[node_idx],
						geom.z_nodes[n1] - geom.z_nodes[node_idx]
					);

					immutable n4n = Vec3(
						geom.x_nodes[n4] - geom.x_nodes[node_idx],
						geom.y_nodes[n4] - geom.y_nodes[node_idx],
						geom.z_nodes[n4] - geom.z_nodes[node_idx]
					);

					immutable n2n = Vec3(
						geom.x_nodes[n2] - geom.x_nodes[node_idx],
						geom.y_nodes[n2] - geom.y_nodes[node_idx],
						geom.z_nodes[n2] - geom.z_nodes[node_idx]
					);

					immutable normal1 = n4n.cross(n1n).normalize;
					immutable normal2 = n2n.cross(n4n).normalize;

					// Average the 2 local face normals for the vert normal
					immutable ave_norm = 0.5*(normal1 + normal2).normalize;
					geom.x_normals[node_idx] = ave_norm[0];
					geom.y_normals[node_idx] = ave_norm[1];
					geom.z_normals[node_idx] = ave_norm[2];
				}
			} else if(r_idx == radial_stations.length - 1) {
				if(p_idx == 0) {
					// only use n2, n3, and node_idx

					immutable n2n = Vec3(
						geom.x_nodes[n2] - geom.x_nodes[node_idx],
						geom.y_nodes[n2] - geom.y_nodes[node_idx],
						geom.z_nodes[n2] - geom.z_nodes[node_idx]
					);

					immutable n3n = Vec3(
						geom.x_nodes[n1] - geom.x_nodes[node_idx],
						geom.y_nodes[n1] - geom.y_nodes[node_idx],
						geom.z_nodes[n1] - geom.z_nodes[node_idx]
					);

					immutable normal = n3n.cross(n2n).normalize;

					geom.x_normals[node_idx] = normal[0];
					geom.y_normals[node_idx] = normal[1];
					geom.z_normals[node_idx] = normal[2];

				} else if(p_idx == naca0012.length - 1) {
					// only use n1, n3, and node_idx

					immutable n1n = Vec3(
						geom.x_nodes[n1] - geom.x_nodes[node_idx],
						geom.y_nodes[n1] - geom.y_nodes[node_idx],
						geom.z_nodes[n1] - geom.z_nodes[node_idx]
					);

					immutable n3n = Vec3(
						geom.x_nodes[n3] - geom.x_nodes[node_idx],
						geom.y_nodes[n3] - geom.y_nodes[node_idx],
						geom.z_nodes[n3] - geom.z_nodes[node_idx]
					);

					immutable normal = n1n.cross(n3n).normalize;

					geom.x_normals[node_idx] = normal[0];
					geom.y_normals[node_idx] = normal[1];
					geom.z_normals[node_idx] = normal[2];

				} else {
					// only use n1, n2, n3, and node_idx

					immutable n1n = Vec3(
						geom.x_nodes[n1] - geom.x_nodes[node_idx],
						geom.y_nodes[n1] - geom.y_nodes[node_idx],
						geom.z_nodes[n1] - geom.z_nodes[node_idx]
					);

					immutable n3n = Vec3(
						geom.x_nodes[n3] - geom.x_nodes[node_idx],
						geom.y_nodes[n3] - geom.y_nodes[node_idx],
						geom.z_nodes[n3] - geom.z_nodes[node_idx]
					);

					immutable n2n = Vec3(
						geom.x_nodes[n2] - geom.x_nodes[node_idx],
						geom.y_nodes[n2] - geom.y_nodes[node_idx],
						geom.z_nodes[n2] - geom.z_nodes[node_idx]
					);

					immutable normal1 = n1n.cross(n3n).normalize;
					immutable normal2 = n3n.cross(n2n).normalize;

					// Average the 2 local face normals for the vert normal
					immutable ave_norm = 0.5*(normal1 + normal2).normalize;
					geom.x_normals[node_idx] = ave_norm[0];
					geom.y_normals[node_idx] = ave_norm[1];
					geom.z_normals[node_idx] = ave_norm[2];
				}
			} else {
				if(p_idx == 0) {
					// only use n2, n3, n4, and node_idx

					immutable n4n = Vec3(
						geom.x_nodes[n4] - geom.x_nodes[node_idx],
						geom.y_nodes[n4] - geom.y_nodes[node_idx],
						geom.z_nodes[n4] - geom.z_nodes[node_idx]
					);

					immutable n3n = Vec3(
						geom.x_nodes[n1] - geom.x_nodes[node_idx],
						geom.y_nodes[n1] - geom.y_nodes[node_idx],
						geom.z_nodes[n1] - geom.z_nodes[node_idx]
					);

					immutable n2n = Vec3(
						geom.x_nodes[n2] - geom.x_nodes[node_idx],
						geom.y_nodes[n2] - geom.y_nodes[node_idx],
						geom.z_nodes[n2] - geom.z_nodes[node_idx]
					);

					immutable normal1 = n2n.cross(n4n).normalize;
					immutable normal2 = n3n.cross(n2n).normalize;

					// Average the 2 local face normals for the vert normal
					immutable ave_norm = 0.5*(normal1 + normal2).normalize;
					geom.x_normals[node_idx] = ave_norm[0];
					geom.y_normals[node_idx] = ave_norm[1];
					geom.z_normals[node_idx] = ave_norm[2];
					
				} else if(p_idx == naca0012.length - 1) {
					// only use n1, n3, n4, and node_idx
					immutable n4n = Vec3(
						geom.x_nodes[n2] - geom.x_nodes[node_idx],
						geom.y_nodes[n2] - geom.y_nodes[node_idx],
						geom.z_nodes[n2] - geom.z_nodes[node_idx]
					);

					immutable n3n = Vec3(
						geom.x_nodes[n1] - geom.x_nodes[node_idx],
						geom.y_nodes[n1] - geom.y_nodes[node_idx],
						geom.z_nodes[n1] - geom.z_nodes[node_idx]
					);

					immutable n1n = Vec3(
						geom.x_nodes[n3] - geom.x_nodes[node_idx],
						geom.y_nodes[n3] - geom.y_nodes[node_idx],
						geom.z_nodes[n3] - geom.z_nodes[node_idx]
					);

					immutable normal1 = n4n.cross(n1n).normalize;
					immutable normal2 = n1n.cross(n3n).normalize;

					// Average the 2 local face normals for the vert normal
					immutable ave_norm = 0.5*(normal1 + normal2).normalize;
					geom.x_normals[node_idx] = ave_norm[0];
					geom.y_normals[node_idx] = ave_norm[1];
					geom.z_normals[node_idx] = ave_norm[2];
				} else {
					// only use all
					immutable n4n = Vec3(
						geom.x_nodes[n4] - geom.x_nodes[node_idx],
						geom.y_nodes[n4] - geom.y_nodes[node_idx],
						geom.z_nodes[n4] - geom.z_nodes[node_idx]
					);

					immutable n3n = Vec3(
						geom.x_nodes[n3] - geom.x_nodes[node_idx],
						geom.y_nodes[n3] - geom.y_nodes[node_idx],
						geom.z_nodes[n3] - geom.z_nodes[node_idx]
					);

					immutable n2n = Vec3(
						geom.x_nodes[n2] - geom.x_nodes[node_idx],
						geom.y_nodes[n2] - geom.y_nodes[node_idx],
						geom.z_nodes[n2] - geom.z_nodes[node_idx]
					);

					immutable n1n = Vec3(
						geom.x_nodes[n1] - geom.x_nodes[node_idx],
						geom.y_nodes[n1] - geom.y_nodes[node_idx],
						geom.z_nodes[n1] - geom.z_nodes[node_idx]
					);

					immutable normal1 = n1n.cross(n3n).normalize;
					immutable normal2 = n3n.cross(n2n).normalize;
					immutable normal3 = n2n.cross(n4n).normalize;
					immutable normal4 = n4n.cross(n1n).normalize;
					
					// Average the 4 local face normals for the vert normal
					immutable ave_norm = 0.25*(normal1 + normal2 + normal3 + normal4).normalize;
					geom.x_normals[node_idx] = ave_norm[0];
					geom.y_normals[node_idx] = ave_norm[1];
					geom.z_normals[node_idx] = ave_norm[2];
				}
			}

			//geom.x_normals[node_idx] = -geom.x_normals[node_idx];
			//geom.y_normals[node_idx] = -geom.y_normals[node_idx];
			//geom.z_normals[node_idx] = -geom.z_normals[node_idx];
			
			node_idx++;
		}
	}

	return geom;
}

struct BladeFlappingState {
	double h;
	double h_star;
}

immutable double w = 0.9912;
immutable double a0 = -2.143;
immutable double a1 = 0.4677;
immutable double b1 = -1.442;
immutable double a2 = -0.2405;
immutable double b2 = 0.1289;
immutable double a3 = -0.1355;
immutable double b3 = 0.2638;
immutable double a4 = 0.02706;
immutable double b4 = 0.07691;
immutable double a5 = 0.009568;
immutable double b5 = 0.01391;
immutable double a6 = -0.0005588;
immutable double b6 = 0.01027;
immutable double a7 = -0.002984;
immutable double b7 = 0.004347;
immutable double a8 = -0.0008832;
immutable double b8 = 0.01659;

BladeFlappingState get_blade_flapping_state(double azimuth) {
	BladeFlappingState bfs;

	immutable double cos_1 = cos(w*azimuth);
	immutable double sin_1 = sin(w*azimuth);
	immutable double cos_2 = cos(2.0*w*azimuth);
	immutable double sin_2 = sin(2.0*w*azimuth);
	immutable double cos_3 = cos(3.0*w*azimuth);
	immutable double sin_3 = sin(3.0*w*azimuth);
	immutable double cos_4 = cos(4.0*w*azimuth);
	immutable double sin_4 = sin(4.0*w*azimuth);
	immutable double cos_5 = cos(5.0*w*azimuth);
	immutable double sin_5 = sin(5.0*w*azimuth);
	immutable double cos_6 = cos(6.0*w*azimuth);
	immutable double sin_6 = sin(6.0*w*azimuth);
	immutable double cos_7 = cos(7.0*w*azimuth);
	immutable double sin_7 = sin(7.0*w*azimuth);
	immutable double cos_8 = cos(8.0*w*azimuth);
	immutable double sin_8 = sin(8.0*w*azimuth);

	bfs.h = a0
		+ a1*cos_1 + b1*sin_1
		+ a2*cos_2 + b2*sin_2
		+ a3*cos_3 + b3*sin_3
		+ a4*cos_4 + b4*sin_4
		+ a5*cos_5 + b5*sin_5
		+ a6*cos_6 + b6*sin_6
		+ a7*cos_7 + b7*sin_7
		+ a8*cos_8 + b8*sin_8;

	bfs.h_star = 
		- w*a1*sin_1 + w*b1*cos_1
		- 2.0*w*a2*sin_2 + 2.0*w*b2*cos_2
		- 3.0*w*a3*sin_3 + 3.0*w*b3*cos_3
		- 4.0*w*a4*sin_4 + 4.0*w*b4*cos_4
		- 5.0*w*a5*sin_5 + 5.0*w*b5*cos_5
		- 6.0*w*a6*sin_6 + 6.0*w*b6*cos_6
		- 7.0*w*a7*sin_7 + 7.0*w*b7*cos_7
		- 8.0*w*a8*sin_8 + 8.0*w*b8*cos_8;

	bfs.h /= 100.0;
	bfs.h_star /= 100.0;
	return bfs;
}

auto charm_wake_hist = [0.000, 5.000, 10.00, 15.00, 20.00, 25.00, 30.00, 35.00, 40.00, 45.00, 50.00, 55.00, 60.00, 65.00, 70.00, 75.00, 80.00, 85.00, 90.00, 95.00, 100.0, 105.0, 110.0, 115.0, 120.0, 125.0, 130.0, 135.0, 140.0, 145.0, 150.0, 155.0, 160.0, 165.0, 170.0, 175.0, 180.0, 185.0, 190.0, 195.0, 200.0, 205.0, 210.0, 215.0, 220.0, 225.0, 230.0, 235.0, 240.0, 245.0, 250.0, 255.0, 260.0, 265.0, 270.0, 275.0, 280.0, 285.0, 290.0, 295.0, 300.0, 305.0, 310.0, 315.0, 320.0, 325.0, 330.0, 335.0, 340.0, 345.0, 350.0, 355.0, 360.0, 365.0, 370.0, 375.0, 380.0, 385.0, 390.0, 395.0, 400.0, 405.0, 410.0, 415.0, 420.0, 425.0, 430.0, 435.0, 440.0, 445.0, 450.0, 455.0, 460.0, 465.0, 470.0, 475.0, 480.0, 485.0, 490.0, 495.0, 500.0, 505.0, 510.0, 515.0, 520.0, 525.0, 530.0, 535.0, 540.0, 545.0, 550.0, 555.0, 560.0, 565.0, 570.0, 575.0, 580.0, 585.0, 590.0, 595.0, 600.0, 605.0, 610.0, 615.0, 620.0, 625.0, 630.0, 635.0, 640.0, 645.0, 650.0, 655.0, 660.0, 665.0, 670.0, 675.0, 680.0, 685.0, 690.0, 695.0, 700.0, 705.0, 710.0, 715.0, 720.0, 725.0, 730.0, 735.0, 740.0, 745.0, 750.0, 755.0, 760.0, 765.0, 770.0, 775.0, 780.0, 785.0, 790.0, 795.0, 800.0, 805.0, 810.0, 815.0, 820.0, 825.0, 830.0, 835.0, 840.0, 845.0, 850.0, 855.0, 860.0, 865.0, 870.0, 875.0, 880.0, 885.0, 890.0, 895.0, 900.0, 905.0, 910.0, 915.0, 920.0, 925.0, 930.0, 935.0, 940.0, 945.0, 950.0, 955.0, 960.0, 965.0, 970.0, 975.0, 980.0, 985.0, 990.0, 995.0, 1000., 1005., 1010., 1015., 1020., 1025., 1030., 1035., 1040., 1045., 1050., 1055., 1060., 1065., 1070., 1075.];
auto charm_wake_z = [0.1415E-01,  0.5659E-02,  0.4794E-02,  0.1652E-02, -0.1370E-02,  0.1667E-02,  0.1129E-02,  0.8581E-03, -0.2589E-04,  0.4925E-02,  0.3237E-02,  0.1251E-02, -0.5850E-03, -0.2450E-02, -0.5382E-02, -0.1752E-02, -0.7786E-02, -0.2334E-01,  0.6998E-02, -0.1770E-03, -0.1488E-02,  0.2561E-02,  0.4327E-01,  0.3913E-01,  0.3534E-01,  0.2946E-01,  0.2354E-01,  0.2341E-01,  0.2880E-01,  0.2052E-01,  0.2353E-01,  0.1710E-01, -0.1070E-01,  0.7732E-02,  0.1565E-01, -0.1559E-02, -0.8739E-01, -0.9637E-01,  0.2002    ,  0.6928    ,  0.6689    ,  0.6733    ,  0.7187    ,  0.7764    ,  0.7878    ,  0.7504    ,  0.7236    ,  0.7069    ,  0.6999    ,  0.7019    ,  0.6979    ,  0.6863    ,  0.6672    ,  0.6431    ,  0.6292    ,  0.6194    ,  0.6213    ,  0.6354    ,  0.6510    ,  0.6686    ,  0.6829    ,  0.6995    ,  0.7162    ,  0.7093    ,  0.5629    ,  0.6536    ,  0.5814    ,  0.6826    ,  0.1803    , -0.1344E-01,  0.1350    ,  0.1907    ,  0.2542    ,  0.2006    ,  0.2340    ,  0.2027    ,  0.2644    ,  0.3104    ,  0.3305    ,  0.3325    ,  0.3397    ,  0.3620    ,  0.3681    ,  0.3819    ,  0.4452    ,  0.4119    ,  0.4244    ,  0.4659    ,  0.5299    ,  0.5015    ,  0.5512    ,  0.6258    ,  0.7240    ,  0.8323    ,  0.8067    ,  0.8007    ,  0.7754    ,  0.7831    ,  0.7525    ,  0.6768    ,  0.4041    ,  0.2550    ,  0.2524    ,  0.2256    ,  0.1804    ,  0.6244E-01,  0.3343    ,  0.6143    ,  0.4831    ,  0.6057    ,   1.090    ,  0.7927    ,   1.005    ,   1.150    ,   1.563    ,   1.689    ,   1.834    ,   1.886    ,   1.902    ,   1.893    ,   1.885    ,   1.861    ,   1.820    ,   1.752    ,   1.660    ,   1.552    ,   1.454    ,   1.392    ,   1.366    ,   1.379    ,   1.428    ,   1.458    ,   1.488    ,   1.460    ,   1.287    ,   1.113    ,  0.9391    ,  0.8104    ,  0.9867    ,  0.9818    ,  0.1565    ,  0.7422    ,  0.7045    ,  0.3969    ,  0.3591E-01,  0.3729E-01,  0.2666    ,  0.1002    ,  0.8929E-01,  0.7928E-02,  0.3439    ,  0.6203    ,  0.8306    ,  0.9004    ,  0.8291    ,  0.8372    ,  0.8246    ,  0.8180    ,  0.7586    ,  0.7847    ,  0.8753    ,  0.9501    ,  0.9981    ,   1.146    ,   1.397    ,   1.582    ,   1.509    ,   1.508    ,   1.485    ,   1.420    ,   1.125    ,  0.5892    ,  0.7355    ,  0.5649    ,  0.7976    ,  0.8740    ,   1.035    ,  0.9153    ,   1.137    ,  0.6561    ,  0.8456    ,  0.7585    ,   1.079    ,   1.322    ,   1.211    ,   1.477    ,   1.648    ,   2.191    ,   2.618    ,   3.192    ,   2.927    ,   2.732    ,   2.633    ,   2.605    ,   2.584    ,   2.503    ,   2.348    ,   2.191    ,   2.072    ,   1.920    ,   1.889    ,   1.862    ,   1.788    ,   1.886    ,   2.125    ,   2.126    ,   1.898    ,   1.536    ,  0.9970    ,  0.9774    ,   1.181    ,  0.9170    ,   1.066    ,   1.142    ,  0.9814    ,  0.8361];

alias LoadingFileDef = AperiodicStructuredLoadingFile!(LoadingType.surface_loading);
alias LoadingFile = LoadingFileT!(LoadingType.surface_loading);

int main(string[] args) {
	import std.conv : to;
	import std.math : PI;
	import std.stdio : File, writeln, write;

	size_t elements = 1*48;

	string output_dir = "wopwop_cases";

	double AR = 16.5;
	double sos = 343;
	double density = 1.125;
	double omega = 109.12;

	size_t iterations = 8*1024;
	size_t wake_history_length = 1*1024;

	// 1 degree per revolution
	double d_psi = 1;
	double dt = d_psi*(PI/180.0)/abs(omega);

	double R = 2;
	bool plot = false;
	bool save = false;
	
	double theta_75 = 3;

	float tMin = 1.1;

	immutable num_rotors = 1;
	immutable num_blades = 4;
	
	immutable theta_tw_1 = -8.0*(PI/180.0);
	//immutable theta_tw_1 = -12.0*(PI/180.0);
	
	double theta_1c = 0;
	double theta_1s = 0;

	double shed_history_angle = 45;
	size_t shed_history = round(shed_history_angle/d_psi).to!size_t;

	double aoa = -5.3;
	double V_inf = 32.9;
	bool plot_wake = false;

	arraySep = ",";
	auto help_information = getopt(
		args,
		std.getopt.config.bundling,
		"iterations|i", "Number of iterations to simulate. Default = "~iterations.to!string, &iterations,
		"output_dir|o", "Base output directory. Default = "~output_dir.to!string, &output_dir,
		"plot|p", "Plot results. Default = "~plot.to!string, &plot,
		"save|s", "Save plot to files. Default = "~save.to!string, &save,
		"collective|c", "Collective pitch in degrees. Default = "~theta_75.to!string, &theta_75,
		"tMin|t", "Minimum observer time", &tMin,
		"aoa|a", "Angle of attack (degrees). Default = "~aoa.to!string, &aoa,
		"vinf|v", "Free streem velocity (m/s). Default = "~V_inf.to!string, &V_inf,
		"theta_1s|b", "Theta 1s", &theta_1s,
		"theta_1c|n", "Theta 1s", &theta_1c,
		"plot_wake|w", "Save wake to vtu", &plot_wake
	);

	if(help_information.helpWanted) {
		defaultGetoptPrinter("Program options:", help_information.options);
		import core.stdc.stdlib : exit;
		exit(-1);
	}

	aoa *= (PI/180.0);
	theta_75 *= (PI/180.0);

	writeln("dt = ", dt);
	writeln("dt*omega = ", dt*omega);
	
	size_t num_cases = 1;
	CaseList case_list = {
		cases: new Casename[num_cases]
	};

	auto num_chunks = elements/chunk_size;
	if(elements % chunk_size != 0) {
		num_chunks++;
	}

	double z_mic_loc = -1.1075*R;
	double[] x_mic_loc = linspace(-2.0*R, 2.0*R, 23);
	double[] y_mic_loc = linspace(-1.25*R, 1.25*R, 17);

	immutable total_mics = x_mic_loc.length*y_mic_loc.length;

	double[] x_obs = new double[total_mics];
	double[] y_obs = new double[total_mics];
	double[] z_obs = new double[total_mics];

	size_t tot_idx = 0;
	foreach(y_idx, ref y_mic; y_mic_loc) {
		foreach(x_idx, ref x_mic; x_mic_loc) {
			x_obs[tot_idx] = x_mic;
			y_obs[tot_idx] = y_mic;
			z_obs[tot_idx] = z_mic_loc;
			tot_idx++;
		}
	}


	string case_directory = "./";

	import std.file : mkdirRecurse, getcwd, exists, write;
	import std.path : buildPath, dirSeparator;

	string full_output_dir = getcwd.buildPath(output_dir, case_directory);

	if(!full_output_dir.exists) {
		full_output_dir.mkdirRecurse;
	}

	string observer_filename = "observer_grid.inp";
	auto obs_file = File(full_output_dir.buildPath(observer_filename), "w");

	obs_file.writeln(x_mic_loc.length, " ", y_mic_loc.length, " ", 1);
	obs_file.writeln(x_obs[0], x_obs[1..$].map!(a => " "~a.to!string).join);
	obs_file.writeln(y_obs[0], y_obs[1..$].map!(a => " "~a.to!string).join);
	obs_file.writeln(z_obs[0], z_obs[1..$].map!(a => " "~a.to!string).join);
	
	obs_file.close;

	string geom_filename = "blade_geom.dat";

	immutable atmo = Atmosphere(density, 18.03e-6);

	writeln("aoa: ", aoa*(180.0/PI));
	writeln("V_inf: ", V_inf);
	writeln("mu = ", V_inf*cos(aoa)/(R*omega));
	writeln("mu_charm = ", V_inf/(R*omega));
	auto r = generate_radius_points(elements);
	double[] C_l_alpha = new double[elements];
	C_l_alpha[] = 2.0*PI;
	immutable double c_r = 1.0/AR;

	double[] c = new double[elements];

	c[] = c_r;
	double[] twist = new double[elements];
	double[] sweep = new double[elements];
	sweep[] = 0;

	twist[] = (r[] - 0.75)*theta_tw_1;

	writeln("Twist: ", twist.map!(a => a*(180.0/PI)));
	writeln("c: ", c.map!(a => a*R));
	writeln("r: ", r);

	double[] alpha_0 = new double[elements];
	alpha_0[] = -2.5*(PI/180.0);

	auto aircraft = Aircraft(num_rotors);

	immutable double d_azimuth = (2.0*PI)/num_blades;

	CB aircraft_cob = {
		Title: "Forward Velocity",
		TranslationType: TranslationType.known_function,
		AH: FVec3(0, 0, 0),
		VH: FVec3(V_inf/+*cos(flight_path_angle)+/, 0, 0/+V_inf*sin(flight_path_angle)+/),
		Y0: FVec3(0, 0, 0)
	};

	ContainerIn wopwop_aircraft = {
		Title: "Tilt Rotor Aircraft",
		//tauMin: 0,
		dTau: dt,
		tauMax: dt*iterations.to!float,
		//nTau: iterations,
		children: new ContainerIn[num_rotors],
		cobs: [aircraft_cob]
	};

	Vec3[] origins = [Vec3(0, 0, 0), Vec3(-2.3, 0, 0), Vec3(4, 0, 0)];

	foreach(r_idx, ref rotor; aircraft.rotors) {
		rotor = RotorGeometry(
			num_blades,
			origins[r_idx],
			R,
			0
		);

		wopwop_aircraft.children[r_idx].Title = "Nacelle "~r_idx.to!string;
		CB nacelle_tilt = {
			Title: "Nacelle "~r_idx.to!string~" tilt angle",
			AxisType: AxisType.time_independent,
			AngleType: AngleType.time_independent,
			AxisValue: FVec3(0, 1, 0),
			AngleValue: -aoa//0*(PI/180.0)
		};

		CB nacelle_offset = {
			Title: "Necelle "~r_idx.to!string~" wing offset",
			TranslationType: TranslationType.time_independent,
			TranslationValue: FVec3(0, R*origins[r_idx][1], 0)
		};

		wopwop_aircraft.children[r_idx].cobs ~= nacelle_tilt;
		wopwop_aircraft.children[r_idx].cobs ~= nacelle_offset;

		ContainerIn wopwop_rotor = {
			Title: "Rotor "~r_idx.to!string,
			cobs: [
				{
					Title: "Rotor "~r_idx.to!string~" verticle offset",
					TranslationType: TranslationType.time_independent,
					TranslationValue: FVec3(0, 0, R*origins[r_idx][2])
				},
				{
					Title: "Rotor "~r_idx.to!string~" rotation",
					AxisType: AxisType.time_independent,
					AxisValue: FVec3(0, 0, 1),
					AngleType: AngleType.known_function,
				}
			]
		};

		foreach(b_idx, ref blade; rotor.blades) {
			blade.chunks = new BladeGeometryChunk[num_chunks];
			blade.set_geometry_array!"r"(r);
			blade.set_geometry_array!"twist"(twist);
			blade.set_geometry_array!"sweep"(sweep);
			blade.set_geometry_array!"alpha_0"(alpha_0);
			blade.set_geometry_array!"chord"(c);
			blade.set_geometry_array!"C_l_alpha"(C_l_alpha);
			blade.azimuth_offset = b_idx*d_azimuth;
			import std.algorithm : sum;
			blade.average_chord = R*c.sum/c.length.to!double;

			ContainerIn wopwop_blade = {
				Title: "blade "~b_idx.to!string~" container",
				patchGeometryFile: geom_filename,
				patchLoadingFile: "blade_"~r_idx.to!string~"_"~b_idx.to!string~"_loading.dat",
				cobs: [
					{
						Title: "blade "~b_idx.to!string~" azimuth offset",
						AxisType: AxisType.time_independent,
						AngleType: AngleType.time_independent,
						AxisValue: FVec3(0, 0, 1),
						AngleValue: blade.azimuth_offset
					}
				]
			};

			wopwop_rotor.children ~= wopwop_blade;
		}

		wopwop_aircraft.children[r_idx].children ~= wopwop_rotor;

		rotor.solidity = num_blades.to!double*rotor.blades[0].average_chord/(PI*rotor.radius);
	}

	auto ac_timehistory = AircraftTimehistory(aircraft, 2);


	auto input_state = AircraftInputState(num_rotors, num_blades);
	input_state.rotor_inputs[0].angle_of_attack = aoa;
	input_state.rotor_inputs[0].angular_velocity = omega;
	input_state.rotor_inputs[0].angular_accel = 0;
	input_state.rotor_inputs[0].freestream_velocity = V_inf;
	input_state.rotor_inputs[0].azimuth = 0;
	input_state.rotor_inputs[0].r_0[] = 0.0000001*aircraft.rotors[0].blades[0].average_chord;
	//ac_timehistory.aircraft_history[0].rotor_states[0].C_T = 0.005;

	if(num_rotors > 1) {
		input_state.rotor_inputs[1].angle_of_attack = aoa;
		input_state.rotor_inputs[1].angular_velocity = -omega;
		input_state.rotor_inputs[1].angular_accel = 0;
		input_state.rotor_inputs[1].freestream_velocity = V_inf;
		input_state.rotor_inputs[1].azimuth = 0;
		input_state.rotor_inputs[1].r_0[] = 0.00001*aircraft.rotors[0].blades[0].average_chord;
		//ac_timehistory.aircraft_history[0].rotor_states[1].C_T = 0.005;
	}

	foreach(r_idx, ref inp; input_state.rotor_inputs) {
		wopwop_aircraft.children[r_idx].children[0].cobs[1].Omega = inp.angular_velocity;
	}

	string namelist_filename = full_output_dir.buildPath("case.nam");
	
	Casename casename = {
		globalFolderName: ".".buildPath(case_directory)~dirSeparator,
		caseNameFile: "case.nam"
	};

	case_list.cases[0] = casename;

	double[] real_chord = new double[r.length];
	real_chord[] = (R/AR);
	auto blade_geom = generate_blade_geom(r, twist, R, real_chord);

	auto lifting_line_geometry_data = ConstantGeometryData(r.length);

	lifting_line_geometry_data.x_nodes[] = r.map!(a => (R*a).to!float).array;
	lifting_line_geometry_data.y_nodes[] = 0;
	lifting_line_geometry_data.z_nodes[] = 0;

	lifting_line_geometry_data.x_normals[] = 0;
	lifting_line_geometry_data.y_normals[] = 0;
	lifting_line_geometry_data.z_normals[] = 1;

	auto geometry = ConstantStructuredGeometryFile(
		"Tilt rotor blade geometry",
		"Pa",
		DataAlignment.node_centered,
		[
			ConstantStructuredGeometryFile.HeaderType(
				"blade",
				r.length.to!int,
				naca0012.length.to!int,
			),
			ConstantStructuredGeometryFile.HeaderType(
				"lifting line",
				r.length.to!int,
				1
			)
		]
	);


	immutable blade_passing_freq = num_blades.to!double*abs(omega)*RAD_TO_HZ;
	immutable lower_mid_freq = 6*blade_passing_freq;
	immutable upper_mid_freq = 40*blade_passing_freq;

	size_t observer_nt = 4*(1.0*(2.0*PI/abs(omega))/dt).to!size_t;

	import std.algorithm : minElement, maxElement;
	Namelist namelist = {
		environment_in: {
			pressureFolderName: "pressure",
			SPLFolderName: "spl",
			sigmaFolderName: "sigma",
			debugLevel: 12,
			ASCIIOutputFlag: false,
			OASPLdBFlag: true,
			OASPLdBAFlag: false,
			spectrumFlag: true,
			SPLdBFlag: true,
			SPLdBAFlag: false,
			pressureGradient1AFlag: false,
			acousticPressureFlag: true,
			thicknessNoiseFlag: true,
			loadingNoiseFlag: true,
			totalNoiseFlag: true,
			sigmaFlag: true
		},
		environment_constants: {
			rho: density,
			c: sos,
			nu: 1.8e-5
		},
		observers: [
			{
				Title: "Mic",
				nt: observer_nt,
				tMin: tMin,
				tMax: tMin + 2.0*(2.0*PI/abs(omega)),
				nbx: 23,
				nby: 17,
				nbz: 1,
				xMin: -2.0*R,
				xMax: 2.0*R,
				yMin: -1.345*R,
				yMax: 1.345*R,
				zMin: -1.1075*R,
				zMax: -1.1075*R,
				highPassFrequency: lower_mid_freq,
				lowPassFrequency: upper_mid_freq,
				cobs: [aircraft_cob]
			}
		],
		containers: [wopwop_aircraft]
	};

	namelist.write_namelist(namelist_filename);

	auto geometry_file = create_geometry_file(geometry, full_output_dir.buildPath(geom_filename));
	geometry_file.append_geometry_data(blade_geom, 0);
	geometry_file.append_geometry_data(lifting_line_geometry_data, 1);
	geometry_file.close_geometry_file();

	LoadingFile[][] loading_files = new LoadingFile[][](num_rotors, num_blades);

	foreach(r_idx; 0..num_rotors) {
		foreach(b_idx; 0..num_blades) {
			auto loading = LoadingFileDef(
				"Tilt rotor blade loading",
				ReferenceFrame.patch_fixed,
				DataAlignment.node_centered,
				[
					LoadingFileDef.HeaderType(
						"Dummy blade loading",
						iterations.to!int,
						r.length.to!int,
						naca0012.length.to!int,
						1, // zone # (can't wait for named arguments)
						true, // Do compute thickness noise for this patch,
						false // No loading data
					),
					LoadingFileDef.HeaderType(
						"lifting line loading",
						iterations.to!int,
						r.length.to!int,
						1,
						2, // zone # (can't wait for named arguments)
						false, // Do not compute thickness noise for this patch,
						true // Has loading data
					)
				]
			);

			immutable loading_filename = full_output_dir.buildPath("blade_"~r_idx.to!string~"_"~b_idx.to!string~"_loading.dat");
			loading_files[r_idx][b_idx] = create_loading_file(loading, loading_filename);
		}
	}

	import std.datetime.stopwatch : benchmark;

	HuangPetersInflow[] inflows = new HuangPetersInflow[num_rotors];

	foreach(i, ref inflow; inflows) {
		//inflow = new HuangPetersInflow(4, 2, &aircraft.rotors[i], dt);
		inflow = new HuangPetersInflow(6, 4, &aircraft.rotors[i], dt);
		//inflow = new HuangPetersInflow(6, 4, &aircraft.rotors[i], dt);
	}
	
	auto wake_history = WakeHistory(num_rotors, num_blades, wake_history_length, 2, elements, shed_history);

	double[] C_T_history = new double[iterations];
	double[] station_loading_1 = new double[iterations];
	double[] station_loading_2 = new double[iterations];

	double[] u_t_1 = new double[iterations];
	double[] u_t_2 = new double[iterations];

	double[] inflow_1 = new double[iterations];
	double[] inflow_2 = new double[iterations];

	double[] inflow_angle_1 = new double[iterations];
	double[] inflow_angle_2 = new double[iterations];

	double[] v_z_1 = new double[iterations];
	double[] v_z_2 = new double[iterations];

	double[] v_y_1 = new double[iterations];
	double[] v_y_2 = new double[iterations];

	double[] v_x_1 = new double[iterations];
	double[] v_x_2 = new double[iterations];


	double[] azimuth_history = new double[iterations];
	double[] time = new double[iterations];

	import std.algorithm : map, sum;

	ZoneLoadingData loading_data;
	loading_data.x_loading = new float[r.length];
	loading_data.y_loading = new float[r.length];
	loading_data.z_loading = new float[r.length];

	double[] z_loading = new double[r.length];

	loading_data.x_loading[] = 0;
	loading_data.y_loading[] = 0;

	import std.datetime.stopwatch : StopWatch, AutoStart;

	auto s = StopWatch(AutoStart.yes);
	auto start_time = s.peek;
	auto total_start_time = start_time;

	double[] U_p = new double[elements];
	double[] U_t = new double[elements];
	double[] c_l = new double[elements];
	double[] alpha = new double[elements];
	double[] gamma = new double[elements];
	double[] d_gamma = new double[elements];

	double label_font_size = 13;
	double tick_font_size = 12;
	double marker_size = 0.5;
	double line_width_wake = 0.5;
	double line_width_tpp = 0.7;

	auto colors = ["r*-", "g*-", "b*-", "c*-", "m*-", "y*-"];

	double max_dim_1 = -double.infinity;
	double min_dim_1 = double.infinity;

	double max_dim_2 = -double.infinity;
	double min_dim_2 = double.infinity;

	double max_dim_3 = -double.infinity;
	double min_dim_3 = double.infinity;

	auto vtk_rotor = build_base_vtu_rotor(aircraft.rotors[0]);
	auto vtk_wake = build_base_vtu_wake(wake_history[0]);

	size_t C_T_len = round(2.0*PI/(dt*abs(omega))).to!size_t;
	auto average_C_T_array = new double[C_T_len];
	average_C_T_array[] = 0;
	double average_C_T = 0;

	auto average_C_My_array = new double[C_T_len];
	average_C_My_array[] = 0;
	double average_C_My = 0;

	auto average_C_Mx_array = new double[C_T_len];
	average_C_Mx_array[] = 0;
	double average_C_Mx = 0;

	size_t plot_blade = 2;


	double x_min = -7.0;
	double x_max = 7.0;
	double z_min = -4.2;
	double z_max = 1.2;
	double y_min = -1.2;
	double y_max = 1.2;

	size_t x_res = 512;
	size_t y_res = 512;
	size_t z_res = 512;
	double[] x_mesh = linspace(x_min, x_max, x_res);
	//auto x_mesh = linspace(-1.0, 1.0, x_res);
	double[] y_mesh = linspace(y_min, y_max, y_res);
	double[] z_mesh = linspace(z_min, z_max, z_res);
	//auto z_mesh = linspace(0.0, 2.0, z_res);
	double[][] inflow_xz = new double[][](z_res, x_res);
	double[][] inflow_xy = new double[][](y_res, x_res);
	double[][] inflow_yz = new double[][](z_res, y_res);

	Chunk y_slice = 0.96;
	Chunk z_slice = 0.0017;
	Chunk x_slice = -0.25;
	Chunk x_e = 0;

	foreach(i; 0..iterations) {
		average_C_T = average_C_T_array.sum/C_T_len.to!double;
		average_C_My = average_C_My_array.sum/C_T_len.to!double;
		average_C_Mx = average_C_Mx_array.sum/C_T_len.to!double;

		if(i % 100 == 0) {

			//get_state_array!"aoa"(ac_timehistory.aircraft_history[0].rotor_states[0].blade_states[0], alpha);
			//get_state_array!"gamma"(ac_timehistory.aircraft_history[0].rotor_states[0].blade_states[0], gamma);
			//get_state_array!"d_gamma"(ac_timehistory.aircraft_history[0].rotor_states[0].blade_states[0], d_gamma);
			auto now = s.peek;

			write(i, ", ");
			write("𝜓: ", input_state.rotor_inputs[0].azimuth, ", ");
			write("χ_1: ", inflows[0].wake_skew*(180.0/PI));
			write(", v_z_1: ", inflows[0].get_average_inflow, "; ");
			write(", C_T: ", average_C_T, "; ");
			write(", C_Mx: ", average_C_Mx, "; ");
			write(", C_My: ", average_C_My, "; ");
			if(num_rotors > 1) {
				write(" χ_2: ", inflows[1].wake_skew*(180.0/PI));
				write(", v_z_2: ", inflows[1].get_average_inflow, "; ");
			}
			writeln(now - start_time);
			//writeln(inflows[0].lambda);
			//writeln(inflows[0].lambda_s);

			start_time = now;
		}

		// Super simple accelerating rotor dynamics. No consideration for the actual inertia of the rotor. A dynamics code
		// can/should be coupled in here.
		foreach(r_idx, ref rotor; input_state.rotor_inputs) {
			rotor.azimuth += rotor.angular_velocity*dt + rotor.angular_accel*dt*dt;

			// Keep the azimuth between 0 and 2*PI so we don't
			// lose fp precicion as the sim marches in time and
			// the azimuth grows unbounded.
			if(rotor.azimuth > 2.0*PI) {
				rotor.azimuth = fmod(abs(rotor.azimuth), 2.0*PI);
			}

			foreach(b_idx, ref blade; aircraft.rotors[r_idx].blades) {
				auto blade_azimuth = rotor.azimuth + blade.azimuth_offset;
				auto sin_azimuth = sin(blade_azimuth);
				auto cos_azimuth = cos(blade_azimuth);

				rotor.blade_pitches[b_idx] = theta_75 + theta_1c*cos_azimuth + theta_1s*sin_azimuth;

				auto flapping_state = get_blade_flapping_state(abs(blade_azimuth));
				rotor.blade_flapping[b_idx] = -flapping_state.h/R - 2.0*PI/180.0; // Adding on a 2 degree static coning angle
				rotor.blade_flapping_rate[b_idx] = abs(omega)*flapping_state.h_star/R;
			}
		}

		step(ac_timehistory.aircraft_history[0], aircraft, input_state, inflows, wake_history, atmo, i, dt);

		loading_data.time = i.to!double*dt;
		foreach(r_idx, ref rotor; ac_timehistory.aircraft_history[0].rotor_states[0..1]) {
			foreach(b_idx, ref blade; rotor.blade_states) {
				blade.get_state_array!"dC_T"(z_loading);

				foreach(z_idx, ref z_load; z_loading) {
					loading_data.z_loading[z_idx] = 0.5*z_load.to!float*density*PI*R*R*R*R*abs(omega)*abs(omega);
				}

				loading_files[r_idx][b_idx].append_loading_data(loading_data);
			}
		}

		//u_t_1[i] = ac_timehistory.aircraft_history[0].rotor_states[0].blade_states[0].chunks[$-1].u_t[0];
		//u_t_2[i] = ac_timehistory.aircraft_history[0].rotor_states[1].blade_states[0].chunks[$-1].u_t[0];

		//inflow_1[i] = ac_timehistory.aircraft_history[0].rotor_states[0].blade_states[0].chunks[$-1].inflow[0];
		//inflow_2[i] = ac_timehistory.aircraft_history[0].rotor_states[1].blade_states[0].chunks[$-1].inflow[0];

		//inflow_angle_1[i] = ac_timehistory.aircraft_history[0].rotor_states[0].blade_states[0].chunks[$-1].inflow_angle[0]*(180.0/PI);
		//inflow_angle_2[i] = ac_timehistory.aircraft_history[0].rotor_states[1].blade_states[0].chunks[$-1].inflow_angle[0]*(180.0/PI);

		//v_x_1[i] = ac_timehistory.aircraft_history[0].rotor_states[0].blade_states[0].chunks[$-1].v_x[0];
		//v_x_2[i] = ac_timehistory.aircraft_history[0].rotor_states[1].blade_states[0].chunks[$-1].v_x[0];

		//v_y_1[i] = ac_timehistory.aircraft_history[0].rotor_states[0].blade_states[0].chunks[$-1].v_y[0];
		//v_y_2[i] = ac_timehistory.aircraft_history[0].rotor_states[1].blade_states[0].chunks[$-1].v_y[0];

		//v_z_1[i] = ac_timehistory.aircraft_history[0].rotor_states[0].blade_states[0].chunks[$-1].v_z[0];
		//v_z_2[i] = ac_timehistory.aircraft_history[0].rotor_states[1].blade_states[0].chunks[$-1].v_z[0];

		//station_loading_1[i] = ac_timehistory.aircraft_history[0].rotor_states[0].blade_states[0].chunks[$-1].dC_T[0];
		//station_loading_2[i] = ac_timehistory.aircraft_history[0].rotor_states[1].blade_states[0].chunks[$-1].dC_T[0];

		station_loading_1[i] = ac_timehistory.aircraft_history[0].rotor_states[0].blade_states[plot_blade].chunks[$-2].dC_T[5];
		azimuth_history[i] = ac_timehistory.aircraft_history[0].rotor_states[0].blade_states[plot_blade].azimuth*180.0/PI;

		average_C_T_array[i % C_T_len] = ac_timehistory.aircraft_history[0].rotor_states[0].C_T;
		average_C_Mx_array[i % C_T_len] = ac_timehistory.aircraft_history[0].rotor_states[0].C_My;
		average_C_My_array[i % C_T_len] = ac_timehistory.aircraft_history[0].rotor_states[0].C_Mx;

		C_T_history[i] = ac_timehistory.aircraft_history[0].rotor_states[0].C_T;
		
		time[i] = i.to!double*dt;

		
		if(plot_wake) {
			if(i > (iterations - (360/d_psi).round.to!size_t)) {
				write_wake_vtu("boxwell_wake_", i, vtk_wake, wake_history[0]);
				write_rotor_vtu("rotor", i, 0, vtk_rotor, ac_timehistory.aircraft_history[0].rotor_states[0], input_state.rotor_inputs[0]);
			}
		}

		
		if(i >= (iterations - 360)) {
			if((i % 5) == 0) {
				//Chunk x_e = 0;
				double max_z = 0.20;//-double.infinity;
				foreach(r_idx; 0..1) {
					//inflows[r_idx].contraction_mapping = true;
					inflows[r_idx].debug_coords = true;
					foreach(z_idx, _z_m; z_mesh) {
						foreach(y_idx, _y_m; y_mesh.chunks(chunk_size).enumerate) {

							
							immutable Chunk x_rel = x_slice[] - origins[r_idx][0];
							immutable Chunk y_rel = _y_m[] - origins[r_idx][1];
							immutable Chunk z_rel = _z_m - origins[r_idx][2];

							immutable Chunk x_m = z_rel[];//x_rel[]*cos(aoa) - z_rel[]*sin(aoa);
							immutable Chunk z_m = 0.17;//-x_rel[]*sin(aoa) - z_rel[]*cos(aoa);

							immutable Chunk infl = inflows[r_idx].inflow_at(x_m, y_rel, z_m, x_e, 0)[];

							//max_z = max(max_z, infl[].maxElement);

							inflow_yz[z_idx][y_idx*chunk_size..y_idx*chunk_size + chunk_size] = infl[];
						}
						//writeln("--------------------------------------------------------------------------------------------------");
					}
					//inflows[r_idx].contraction_mapping = false;
					inflows[r_idx].debug_coords = false;
					//writeln("========================================================================================================================================================");
				}

				//writeln("max_z:", max_z);
				auto levels = linspace(-max_z, max_z, 18);
				static import matplotlibd.pyplot;
				matplotlibd.pyplot.clear();
				plt.figure;
				plt.title("YZ Plane Vertical Inflow");
				//plt.contourf(y_mesh, z_mesh, inflow_yz, ["levels": 100], ["cmap": "bwr"], ["vmin": -max_z], ["vmax": max_z]);
				//plt.contourf(y_mesh, z_mesh, inflow_yz, ["levels": 100], ["vmin": -max_z], ["vmax": max_z]);
				plt.contourf(y_mesh, z_mesh, inflow_yz, ["levels": levels], ["cmap": "bwr"]);
				plt.colorbar;
				plt.axis("equal");
				plt.savefig("vert_yz_"~i.to!string~".png", ["dpi": 500]);
			}
		}
	}

	auto total_time = s.peek - total_start_time;

	writeln("Sim took ", total_time);

	//size_t num_x = 512;
	//size_t num_y = 512;
	//size_t num_z = 512;
	//auto starts = Vec3(-2, -2, -2);
	//auto ends = Vec3(2, 2, 2);
	//auto deltas = Vec3((ends[0] - starts[0])/num_x, (ends[1] - starts[1])/num_y, (ends[2] - starts[2])/num_z);

	x_res = 2048;
	y_res = 2048;
	z_res = 2048;
	x_mesh = linspace(x_min, x_max, x_res);
	y_mesh = linspace(y_min, y_max, y_res);
	z_mesh = linspace(z_min, z_max, z_res);
	inflow_xz = new double[][](z_res, x_res);
	inflow_xy = new double[][](y_res, x_res);
	inflow_yz = new double[][](z_res, y_res);

	foreach(r_idx, ref rotor; ac_timehistory.aircraft_history[0].rotor_states) {
		foreach(b_idx, ref _blade; rotor.blade_states) {
			loading_files[r_idx][b_idx].close_loading_file;
		}
	}

	writeln("rotor 0 C_T: ", average_C_T);

	if(plot) {
		/+
		plt.figure;

		//plt.figure;
		plt.gca(["projection": "3d"]);
		max_dim_1 = -double.infinity;
		min_dim_1 = double.infinity;

		max_dim_2 = -double.infinity;
		min_dim_2 = double.infinity;

		max_dim_3 = -double.infinity;
		min_dim_3 = double.infinity;
		foreach(r_idx, ref rotor; ac_timehistory.aircraft_history[0].rotor_states) {
			foreach(b_idx, ref _blade; rotor.blade_states) {

				foreach(ref shed_filament; wake_history[0].rotor_wakes[r_idx].shed_vortices[b_idx].shed_filaments) {
					auto x = shed_filament.get_wake_component!"x";
					auto y = shed_filament.get_wake_component!"y";
					auto z = shed_filament.get_wake_component!"z";
					plt.plot(x, y, z, colors[b_idx], ["linewidth": line_width_wake]);
				}

				auto x = wake_history[0].rotor_wakes[r_idx][b_idx].get_wake_component!"x";
				auto y = wake_history[0].rotor_wakes[r_idx][b_idx].get_wake_component!"y";
				auto z = wake_history[0].rotor_wakes[r_idx][b_idx].get_wake_component!"z";

				max_dim_3 = max(max_dim_3, z.maxElement);
				min_dim_3 = min(min_dim_3, z.minElement);

				max_dim_2 = max(max_dim_2, y.maxElement);
				min_dim_2 = min(min_dim_2, y.minElement);

				max_dim_1 = max(max_dim_1, x.maxElement);
				min_dim_1 = min(min_dim_1, x.minElement);

				plt.plot(x, y, z, colors[b_idx], ["linewidth": 0.5]);
			}
		}


		plt.title("Top View: V = "~V_inf.to!string~", $\\gamma$ = "~_flight_path_angle.to!string~", Tilt angle = "~_tilt_angle.to!string);
		//plt.axis("square");
		double span_1 = max_dim_1 - min_dim_1;
		double span_2 = max_dim_2 - min_dim_2;
		double span_3 = max_dim_3 - min_dim_3;

		double mid_x = 0.5*(max_dim_1 + min_dim_1);
		double mid_y = 0.5*(max_dim_2 + min_dim_2);
		double mid_z = 0.5*(max_dim_3 + min_dim_3);
		double max_span = max(span_1, span_2, span_3);

		plt.xlim(mid_x - max_span, mid_x + max_span);
		plt.ylim(mid_y - max_span, mid_y + max_span);
		plt.zlim(mid_z - max_span, mid_z + max_span);
		if(save) plt.savefig("tilt_rotor_v"~V_inf.to!string~"_fpa"~_flight_path_angle.to!string~"_ta"~_tilt_angle.to!string~"_full_wake_3d.png", ["dpi": 500]);

		plt.figure;
		plt.gca(["projection": "3d"]);
		max_dim_1 = -double.infinity;
		min_dim_1 = double.infinity;

		max_dim_2 = -double.infinity;
		min_dim_2 = double.infinity;

		max_dim_3 = -double.infinity;
		min_dim_3 = double.infinity;
		foreach(r_idx, ref rotor; ac_timehistory.aircraft_history[0].rotor_states) {

			foreach(ref shed_filament; wake_history[0].rotor_wakes[r_idx].shed_vortices[0].shed_filaments) {
				auto x = shed_filament.get_wake_component!"x";
				auto y = shed_filament.get_wake_component!"y";
				auto z = shed_filament.get_wake_component!"z";
				plt.plot(x, y, z, colors[0], ["linewidth": line_width_wake]);
			}

			auto x = wake_history[0].rotor_wakes[r_idx][0].get_wake_component!"x";
			auto y = wake_history[0].rotor_wakes[r_idx][0].get_wake_component!"y";
			auto z = wake_history[0].rotor_wakes[r_idx][0].get_wake_component!"z";

			max_dim_3 = max(max_dim_3, z.maxElement);
			min_dim_3 = min(min_dim_3, z.minElement);

			max_dim_2 = max(max_dim_2, y.maxElement);
			min_dim_2 = min(min_dim_2, y.minElement);

			max_dim_1 = max(max_dim_1, x.maxElement);
			min_dim_1 = min(min_dim_1, x.minElement);

			plt.plot(x, y, z, colors[0], ["linewidth": 0.5]);
		}

		plt.title("Top View: V = "~V_inf.to!string~", $\\gamma$ = "~_flight_path_angle.to!string~", Tilt angle = "~_tilt_angle.to!string);
		//plt.axis("square");
		span_1 = max_dim_1 - min_dim_1;
		span_2 = max_dim_2 - min_dim_2;
		span_3 = max_dim_3 - min_dim_3;

		mid_x = 0.5*(max_dim_1 + min_dim_1);
		mid_y = 0.5*(max_dim_2 + min_dim_2);
		mid_z = 0.5*(max_dim_3 + min_dim_3);
		max_span = max(span_1, span_2, span_3);

		plt.xlim(mid_x - max_span, mid_x + max_span);
		plt.ylim(mid_y - max_span, mid_y + max_span);
		plt.zlim(mid_z - max_span, mid_z + max_span);
		if(save) plt.savefig("tilt_rotor_v"~V_inf.to!string~"_fpa"~_flight_path_angle.to!string~"_ta"~_tilt_angle.to!string~"_single_wake_3d.png", ["dpi": 500]);

		//plt.gca(["projection": "3d"]);

		//azimuth_history[] = -(azimuth_history[] - azimuth_history.maxElement);
		+/


		foreach(y_idx, _y_m; y_mesh) {
			inflow_xy[y_idx][] = 0;
		}

		foreach(z_idx, _z_m; z_mesh) {
			inflow_xz[z_idx][] = 0;
			inflow_yz[z_idx][] = 0;
		}

		double max_infl = 0;

		foreach(r_idx; 0..num_rotors) {
			//writeln(inflows[r_idx].contraction_array);
			//inflows[r_idx].contraction_mapping = true;
			foreach(z_idx, _z_m; z_mesh) {
				foreach(x_idx, _x_m; x_mesh.chunks(chunk_size).enumerate) {

					immutable Chunk x_rel = _x_m[] - origins[r_idx][0];
					immutable Chunk y_rel = y_slice[] - origins[r_idx][1];
					//immutable Chunk z_rel = _z_m - origins[r_idx][2];

					//immutable Chunk x_m = x_rel[]*cos(aoa) - z_rel[]*sin(aoa);
					//immutable Chunk z_m = -x_rel[]*sin(aoa) - z_rel[]*cos(aoa);

					immutable Chunk x_m = x_rel[];
					immutable Chunk z_m = -_z_m;
					auto infl = inflows[r_idx].inflow_at(x_m, y_rel, z_m, x_e, 0);

					max_infl = (max_infl~infl[]).map!(a => abs(a)).maxElement;

					inflow_xz[z_idx][x_idx*chunk_size..x_idx*chunk_size + chunk_size] += infl[];
				}
			}
		}

		
		foreach(r_idx; 0..num_rotors) {
			foreach(y_idx, _y_m; y_mesh) {
				foreach(x_idx, _x_m; x_mesh.chunks(chunk_size).enumerate) {

					immutable Chunk x_rel = _x_m[] - origins[r_idx][0];
					immutable Chunk y_rel = _y_m - origins[r_idx][1];
					immutable Chunk z_rel = z_slice[] - origins[r_idx][2];

					immutable Chunk x_m = x_rel[]*cos(aoa) - z_rel[]*sin(aoa);
					//immutable Chunk z_m = -x_rel[]*sin(aoa) - z_rel[]*cos(aoa);

					immutable Chunk infl = inflows[r_idx].inflow_at(x_m, y_rel, z_slice, x_e, 0)[]/+*cos(aoa).to!double+/;

					max_infl = (max_infl~infl[]).map!(a => abs(a)).maxElement;

					inflow_xy[y_idx][x_idx*chunk_size..x_idx*chunk_size + chunk_size] += infl[];
				}
			}
		}

		foreach(r_idx; 0..num_rotors) {
			foreach(z_idx, _z_m; z_mesh) {
				foreach(y_idx, _y_m; y_mesh.chunks(chunk_size).enumerate) {

					immutable Chunk x_rel = x_slice[] - origins[r_idx][0];
					immutable Chunk y_rel = _y_m[] - origins[r_idx][1];
					immutable Chunk z_rel = _z_m - origins[r_idx][2];

					//immutable Chunk x_m = x_rel[]*cos(aoa) - z_rel[]*sin(aoa);
					immutable Chunk z_m = -z_rel[];//-x_rel[]*sin(aoa) - z_rel[]*cos(aoa);

					immutable Chunk infl = inflows[r_idx].inflow_at(x_slice, y_rel, z_m, x_e, 0)[]/+*cos(aoa).to!double+/;

					max_infl = (max_infl~infl[]).map!(a => abs(a)).maxElement;

					inflow_yz[z_idx][y_idx*chunk_size..y_idx*chunk_size + chunk_size] += infl[];
				}
			}
		}
		
		auto levels = linspace(-max_infl, max_infl, 50);
		//writeln(inflow_xz);

		plt.gca(["projection": "rectilinear"]);
		plt.figure;
		plt.title("XZ Plane Vertical Inflow");
		plt.contourf(x_mesh, z_mesh, inflow_xz, ["levels": levels], ["cmap": "bwr"]);
		plt.colorbar;
		//if(save) plt.savefig("vert_xz.png", ["dpi": 500]);

		/+
		max_dim_1 = -double.infinity;
		min_dim_1 = double.infinity;

		max_dim_2 = -double.infinity;
		min_dim_2 = double.infinity;
		foreach(r_idx, ref rotor; ac_timehistory.aircraft_history[0].rotor_states) {
			foreach(b_idx, ref _blade; rotor.blade_states[0..1]) {

				foreach(ref shed_filament; wake_history[0].rotor_wakes[r_idx].shed_vortices[b_idx].shed_filaments) {
					auto x = shed_filament.get_wake_component!"x";
					auto z = shed_filament.get_wake_component!"z";
					plt.plot(x, z, colors[b_idx], ["linewidth": line_width_wake], ["markersize": marker_size]/+colors[b_idx]+/);
				}

				auto x = wake_history[0].rotor_wakes[r_idx][b_idx].get_wake_component!"x";
				//auto y = wake_history[0].rotor_wakes[r_idx][b_idx].get_wake_component!"y";
				auto z = wake_history[0].rotor_wakes[r_idx][b_idx].get_wake_component!"z";

				//writeln("x: ", x);
				//writeln("z: ", z);
				max_dim_2 = max(max_dim_2, z.maxElement);
				min_dim_2 = min(min_dim_2, z.minElement);

				max_dim_1 = max(max_dim_1, x.maxElement);
				min_dim_1 = min(min_dim_1, x.minElement);
				plt.plot(x, z, colors[b_idx], ["linewidth": line_width_wake], ["markersize": marker_size]/+colors[b_idx]+/);
				//plt.plot(x, y, z, /+colors[b_idx]+/);
				//plt.plot(azimuth_history, z.retro, /+colors[b_idx]+/);
				//plt.plot(x, z, y, /+colors[b_idx]+/);
				//plt.plot(y, z, x/+colors[b_idx]+/);

				static import std.math;
			
				auto x_r = -std.math.cos(_blade.azimuth);
				auto x_b = aircraft.rotors[r_idx].origin[0] + 1.0*std.math.cos(input_state.rotor_inputs[r_idx].angle_of_attack);
				//auto y_b = aircraft.rotors[r_idx].origin[1] + std.math.sin(blade.azimuth);
				auto z_b = aircraft.rotors[r_idx].origin[2] - 1.0*std.math.sin(input_state.rotor_inputs[r_idx].angle_of_attack);

				auto x_b_1 = aircraft.rotors[r_idx].origin[0] - cos(aoa);
				auto z_b_1 = aircraft.rotors[r_idx].origin[2] + sin(aoa);

				auto x_b_2 = aircraft.rotors[r_idx].origin[0] + cos(aoa);
				auto z_b_2 = aircraft.rotors[r_idx].origin[2] - sin(aoa);

				plt.plot([x_b_1, x_b_2], [z_b_1, z_b_2], "k", ["linewidth": line_width_tpp], ["label": "Tip path plane"]);
				//plt.plot([aircraft.rotors[r_idx].origin[0], x_b], [aircraft.rotors[r_idx].origin[1], y_b], [aircraft.rotors[r_idx].origin[2], z_b], "k-", ["linewidth": 1.5]);
				//plt.plot([aircraft.rotors[r_idx].origin[0], x_b], [aircraft.rotors[r_idx].origin[2], z_b], "k-", ["linewidth": line_width_tpp]);
				//plt.plot([aircraft.rotors[r_idx].origin[0], x_b], [aircraft.rotors[r_idx].origin[2], z_b], "k-");
				//plt.plot([aircraft.rotors[r_idx].origin[0], x_b], [aircraft.rotors[r_idx].origin[2], z_b], [aircraft.rotors[r_idx].origin[1], y_b], "k-");
				//plt.plot([aircraft.rotors[r_idx].origin[1], y_b], [aircraft.rotors[r_idx].origin[2], z_b], [aircraft.rotors[r_idx].origin[0], x_b], "k-");
			}
		}
		
		plt.title("Side View: V = "~V_inf.to!string~", $\\gamma$ = "~_flight_path_angle.to!string~", Tilt angle = "~_tilt_angle.to!string);
		
		
		plt.axis("square");

		span_1 = max_dim_1 - min_dim_1;
		span_2 = max_dim_2 - min_dim_2;
		plt.xlim(["left": min_dim_1 - .1*span_1], ["right": max_dim_1 + .1*span_1]);
		plt.ylim(["bottom": min_dim_2 - .9*span_2], ["top": max_dim_2 + .9*span_2]);
		plt.xlabel("x/R", ["fontsize": label_font_size]);
		plt.ylabel("z/R", ["fontsize": label_font_size]);
		plt.xticks(["fontsize": tick_font_size]);
		plt.yticks(["fontsize": tick_font_size]);
		if(save) plt.savefig("tilt_rotor_v"~V_inf.to!string~"_fpa"~_flight_path_angle.to!string~"_ta"~_tilt_angle.to!string~"_side.png", ["dpi": 500]);
		+/

		max_dim_1 = -double.infinity;
		min_dim_1 = double.infinity;

		max_dim_2 = -double.infinity;
		min_dim_2 = double.infinity;
		plt.figure;

		//plt.figure;
		plt.title("XY Plane Vertical Inflow");
		plt.contourf(y_mesh, x_mesh, inflow_xy, ["levels": levels], ["cmap": "bwr"]);
		plt.colorbar;
		//plt.axis("square");
		if(save) plt.savefig("vert_xy.png", ["dpi": 500]);

		/+
		max_dim_1 = -double.infinity;
		min_dim_1 = double.infinity;

		max_dim_2 = -double.infinity;
		min_dim_2 = double.infinity;
		foreach(r_idx, ref rotor; ac_timehistory.aircraft_history[0].rotor_states) {
			foreach(b_idx, ref _blade; rotor.blade_states) {

				foreach(ref shed_filament; wake_history[0].rotor_wakes[r_idx].shed_vortices[b_idx].shed_filaments) {
					auto x = shed_filament.get_wake_component!"x";
					auto y = shed_filament.get_wake_component!"y";
					plt.plot(x, y, colors[b_idx], ["linewidth": line_width_wake], ["markersize": marker_size]/+colors[b_idx]+/);
				}

				auto x = wake_history[0].rotor_wakes[r_idx][b_idx].get_wake_component!"x";
				auto y = wake_history[0].rotor_wakes[r_idx][b_idx].get_wake_component!"y";
				//auto z = wake_history[0].rotor_wakes[r_idx][b_idx].get_wake_component!"z";

				max_dim_2 = max(max_dim_2, y.maxElement);
				min_dim_2 = min(min_dim_2, y.minElement);

				max_dim_1 = max(max_dim_1, x.maxElement);
				min_dim_1 = min(min_dim_1, x.minElement);

				plt.plot(x, y, colors[b_idx], ["linewidth": line_width_wake], ["markersize": marker_size]);
				//plt.plot(x, y, z, /+colors[b_idx]+/);
				//plt.plot(azimuth_history, z.retro, /+colors[b_idx]+/);
				//plt.plot(x, z, y, /+colors[b_idx]+/);
				//plt.plot(y, z, x/+colors[b_idx]+/);

				static import std.math;
			
				auto x_r = -std.math.cos(_blade.azimuth);
				auto x_b = aircraft.rotors[r_idx].origin[0] + x_r*std.math.cos(input_state.rotor_inputs[r_idx].angle_of_attack);
				auto y_b = aircraft.rotors[r_idx].origin[1] + std.math.sin(_blade.azimuth);
				//auto z_b = aircraft.rotors[r_idx].origin[2] + x_r*std.math.sin(input_state.rotor_inputs[r_idx].angle_of_attack);

				//auto y_b_1 = /+rotor.origin[1] ++/ -1;
				//auto z_b_1 = 0;

				//auto y_b_2 = /+rotor.origin[1] ++/ 1;
				//auto z_b_2 = 0;

				//plt.plot([y_b_1, y_b_2], [z_b_1, z_b_2], "k", ["linewidth": line_width_tpp], ["label": "Tip path plane"]);
				//plt.plot([aircraft.rotors[r_idx].origin[0], x_b], [aircraft.rotors[r_idx].origin[1], y_b], [aircraft.rotors[r_idx].origin[2], z_b], "k-", ["linewidth": 1.5]);
				plt.plot([aircraft.rotors[r_idx].origin[0], x_b], [aircraft.rotors[r_idx].origin[1], y_b], "k-", ["linewidth": line_width_tpp]);
				//plt.plot([aircraft.rotors[r_idx].origin[0], x_b], [aircraft.rotors[r_idx].origin[2], z_b], "k-");
				//plt.plot([aircraft.rotors[r_idx].origin[0], x_b], [aircraft.rotors[r_idx].origin[2], z_b], [aircraft.rotors[r_idx].origin[1], y_b], "k-");
				//plt.plot([aircraft.rotors[r_idx].origin[1], y_b], [aircraft.rotors[r_idx].origin[2], z_b], [aircraft.rotors[r_idx].origin[0], x_b], "k-");
			}
		}
		

		plt.title("Top View: V = "~V_inf.to!string~", $\\gamma$ = "~_flight_path_angle.to!string~", Tilt angle = "~_tilt_angle.to!string);
		plt.axis("square");
		span_1 = max_dim_1 - min_dim_1;
		span_2 = max_dim_2 - min_dim_2;
		plt.xlim(["left": min_dim_1 - .1*span_1], ["right": max_dim_1 + .1*span_1]);
		plt.ylim(["top": min_dim_2 - .1*span_2], ["bottom": max_dim_2 + .1*span_2]);
		plt.xlabel("x/R", ["fontsize": label_font_size]);
		plt.ylabel("y/R", ["fontsize": label_font_size]);
		plt.xticks(["fontsize": tick_font_size]);
		plt.yticks(["fontsize": tick_font_size]);
		if(save) plt.savefig("tilt_rotor_v"~V_inf.to!string~"_fpa"~_flight_path_angle.to!string~"_ta"~_tilt_angle.to!string~"_top.png", ["dpi": 500]);
		
		+/
		plt.figure;

		plt.contourf(y_mesh, z_mesh, inflow_yz, ["levels": levels], ["cmap": "bwr"]);
		plt.colorbar;
		//if(save) plt.savefig("tilt_rotor_v"~V_inf.to!string~"_fpa"~_flight_path_angle.to!string~"_ta"~_tilt_angle.to!string~"_front.png", ["dpi": 500]);
		/+
		max_dim_1 = -double.infinity;
		min_dim_1 = double.infinity;

		max_dim_2 = -double.infinity;
		min_dim_2 = double.infinity;
		foreach(r_idx, ref rotor; ac_timehistory.aircraft_history[0].rotor_states) {
			foreach(b_idx, ref _blade; rotor.blade_states[0..1]) {

				foreach(ref shed_filament; wake_history[0].rotor_wakes[r_idx].shed_vortices[b_idx].shed_filaments) {
					auto y = shed_filament.get_wake_component!"y";
					auto z = shed_filament.get_wake_component!"z";
					plt.plot(y, z, colors[b_idx], ["linewidth": line_width_wake], ["markersize": marker_size]/+colors[b_idx]+/);
				}

				auto x = wake_history[0].rotor_wakes[r_idx][b_idx].get_wake_component!"x";
				auto y = wake_history[0].rotor_wakes[r_idx][b_idx].get_wake_component!"y";
				auto z = wake_history[0].rotor_wakes[r_idx][b_idx].get_wake_component!"z";

				max_dim_2 = max(max_dim_2, z.maxElement);
				min_dim_2 = min(min_dim_2, z.minElement);

				max_dim_1 = max(max_dim_1, y.maxElement);
				min_dim_1 = min(min_dim_1, y.minElement);
				plt.plot(y, z, colors[b_idx], ["linewidth": line_width_wake], ["markersize": marker_size]);
				//plt.plot(x, y, z, /+colors[b_idx]+/);
				//plt.plot(azimuth_history, z.retro, /+colors[b_idx]+/);
				//plt.plot(x, z, y, /+colors[b_idx]+/);
				//plt.plot(y, z, x/+colors[b_idx]+/);

				static import std.math;
			
				auto x_r = -std.math.cos(_blade.azimuth);
				//auto x_b = aircraft.rotors[r_idx].origin[0] + x_r*std.math.cos(input_state.rotor_inputs[r_idx].angle_of_attack);
				auto y_b = aircraft.rotors[r_idx].origin[1] + std.math.sin(_blade.azimuth);
				auto z_b = aircraft.rotors[r_idx].origin[2] - 1.0*std.math.sin(input_state.rotor_inputs[r_idx].angle_of_attack);

				auto y_b_1 = /+rotor.origin[1] ++/ -1;
				auto z_b_1 = 0;

				auto y_b_2 = /+rotor.origin[1] ++/ 1;
				auto z_b_2 = 0;

				plt.plot([y_b_1, y_b_2], [z_b_1, z_b_2], "k", ["linewidth": line_width_tpp], ["label": "Tip path plane"]);

				//plt.plot([aircraft.rotors[r_idx].origin[0], x_b], [aircraft.rotors[r_idx].origin[1], y_b], [aircraft.rotors[r_idx].origin[2], z_b], "k-", ["linewidth": 1.5]);
				//plt.plot([aircraft.rotors[r_idx].origin[1], y_b], [aircraft.rotors[r_idx].origin[2], z_b], "k-", ["linewidth": line_width_tpp], ["markersize": marker_size]);
				//plt.plot([aircraft.rotors[r_idx].origin[0], x_b], [aircraft.rotors[r_idx].origin[2], z_b], "k-");
				//plt.plot([aircraft.rotors[r_idx].origin[0], x_b], [aircraft.rotors[r_idx].origin[2], z_b], [aircraft.rotors[r_idx].origin[1], y_b], "k-");
				//plt.plot([aircraft.rotors[r_idx].origin[1], y_b], [aircraft.rotors[r_idx].origin[2], z_b], [aircraft.rotors[r_idx].origin[0], x_b], "k-");
			}
		}

		span_1 = max_dim_1 - min_dim_1;
		span_2 = max_dim_2 - min_dim_2;
		plt.title("Front View: V = "~V_inf.to!string~", $\\gamma$ = "~_flight_path_angle.to!string~", Tilt angle = "~_tilt_angle.to!string);
		plt.axis("square");
		plt.xlim(["left": min_dim_1 - .1*span_1], ["right": max_dim_1 + .1*span_1]);
		plt.ylim(["bottom": min_dim_2 - 0.2*span_2], ["top": max_dim_2 + 0.2*span_2]);
		plt.xlabel("y/R", ["fontsize": label_font_size]);
		plt.ylabel("z/R", ["fontsize": label_font_size]);
		plt.xticks(["fontsize": tick_font_size]);
		plt.yticks(["fontsize": tick_font_size]);
		if(save) plt.savefig("tilt_rotor_v"~V_inf.to!string~"_fpa"~_flight_path_angle.to!string~"_ta"~_tilt_angle.to!string~"_front.png", ["dpi": 500]);
		+/

		//auto psis = linspace(0, 360.0, 360);
		auto psis = linspace(0.0, 2.0*PI, (360/d_psi).to!size_t);
		//auto psis = linspace(2.0*PI, 0, 360);
		double[] radial_points = linspace(0., 1., elements);
		//auto psis = linspace(0, 2.0*PI, 8*360);
		//double[] radial_points = linspace(0., 1., elements*6);

		auto u_t = new double[][](radial_points.length, psis.length);
		auto wake_u_t = new double[][](radial_points.length, psis.length);
		//auto wake_u_t = new double[][](r.length, phis.length);
		auto b_load = new double[][](radial_points.length, psis.length);
		auto b_load_dot = new double[][](radial_points.length, psis.length);
		auto u_p = new double[][](radial_points.length, psis.length);
		auto u_p_dot = new double[][](radial_points.length, psis.length);

		auto cos_alpha = cos(aoa);
		auto sin_alpha = sin(aoa);

		foreach(p_idx, psi; psis) {
			auto cos_azimuth = cos(psi/+*(PI/180)+/);
			auto sin_azimuth = sin(psi/+*(PI/180)+/);

			size_t old_idx;
			if(p_idx == 0) {
				old_idx = psis.length - 1;
			} else {
				old_idx = p_idx - 1;
			}
			auto old_cos_azimuth = cos(psis[old_idx]);
			auto old_sin_azimuth = sin(psis[old_idx]);

			//foreach(chunk_idx; 0..elements/chunk_size) {
			foreach(chunk_idx, radial_chunk; radial_points.chunks(chunk_size).enumerate) {
				immutable Chunk x = origins[0][0] - radial_chunk[]*cos_azimuth*cos_alpha;
				immutable Chunk y = (origins[0][1] + radial_chunk[]*sin_azimuth);
				immutable Chunk z = origins[0][2] - radial_chunk[]*cos_azimuth*sin_alpha;

				immutable Chunk old_x = origins[0][0] - radial_chunk[]*old_cos_azimuth*cos_alpha;
				immutable Chunk old_y = (origins[0][1] + radial_chunk[]*old_sin_azimuth);
				immutable Chunk old_z = origins[0][2] - radial_chunk[]*old_cos_azimuth*sin_alpha;
				//immutable Chunk x = origins[0][0] - aircraft.rotors[0].blades[0].chunks[chunk_idx].r[]*sin_azimuth*cos_alpha;
				//immutable Chunk y = origins[0][1] + aircraft.rotors[0].blades[0].chunks[chunk_idx].r[]*cos_azimuth;
				//immutable Chunk z = origins[0][2] - aircraft.rotors[0].blades[0].chunks[chunk_idx].r[]*sin_azimuth*sin_alpha;

				auto old_wake_velocities = wake_history[1].compute_wake_induced_velocities(old_x, old_y, old_z, ac_timehistory.aircraft_history[0], input_state.rotor_inputs[0].angular_velocity, 0, 0);
				auto wake_velocities = wake_history[0].compute_wake_induced_velocities(x, y, z, ac_timehistory.aircraft_history[0], input_state.rotor_inputs[0].angular_velocity, 0, 0);

				immutable Chunk wake_z = wake_velocities.v_x[]*sin_alpha - wake_velocities.v_z[]*cos_alpha;
				immutable Chunk old_wake_z = old_wake_velocities.v_x[]*sin_alpha - old_wake_velocities.v_z[]*cos_alpha;

				immutable Chunk wake_x = wake_velocities.v_x[]*cos_alpha - wake_velocities.v_z[]*sin_alpha;
				immutable Chunk old_wake_x = old_wake_velocities.v_x[]*cos_alpha - old_wake_velocities.v_z[]*sin_alpha;
				immutable Chunk wake_y = wake_velocities.v_y[];
				immutable Chunk old_wake_y = old_wake_velocities.v_y[];

				immutable Chunk _wake_u_t = sin_azimuth*wake_x[] - cos_azimuth*wake_y[];
				immutable Chunk old_wake_u_t = old_sin_azimuth*old_wake_x[] - old_cos_azimuth*old_wake_y[];

				immutable mu_cos_azimuth = ac_timehistory.aircraft_history[0].rotor_states[0].advance_ratio*sin_azimuth;
				immutable old_mu_cos_azimuth = ac_timehistory.aircraft_history[0].rotor_states[0].advance_ratio*old_sin_azimuth;

				
				immutable Chunk _u_t = _wake_u_t[] + radial_chunk[] + mu_cos_azimuth;
				immutable Chunk old_u_t = old_wake_u_t[] + radial_chunk[] + old_mu_cos_azimuth;
				immutable Chunk _u_p = -wake_z[];
				immutable Chunk old_u_p = -old_wake_z[];
				import opencopter.math.trigonometry : atan2;
				immutable Chunk inflow_angle = atan2(_u_p, _u_t);
				immutable Chunk old_inflow_angle = atan2(old_u_p, old_u_t);

				immutable Chunk dC_L = steady_lift_model(_u_p, R, aircraft.rotors[0].blades[0].chunks[chunk_idx], ac_timehistory.aircraft_history[0].rotor_states[0].blade_states[0].chunks[chunk_idx], _u_t, inflow_angle, 0);
				immutable Chunk old_dC_L = steady_lift_model(old_u_p, R, aircraft.rotors[0].blades[0].chunks[chunk_idx], ac_timehistory.aircraft_history[0].rotor_states[0].blade_states[0].chunks[chunk_idx], old_u_t, old_inflow_angle, 0);

				//writeln("p_idx: ", p_idx);
				//writeln("wake_u_t.length: ", wake_u_t.length);
				//writeln("chunk_idx*chunk_size: ", chunk_idx*chunk_size);
				//writeln("chunk_idx*chunk_size + chunk_size: ", chunk_idx*chunk_size + chunk_size);
				foreach(c_idx; 0..chunk_size) {
					wake_u_t[chunk_idx*chunk_size + c_idx][p_idx] = _wake_u_t[c_idx];
					u_t[chunk_idx*chunk_size + c_idx][p_idx] = _wake_u_t[c_idx] + radial_chunk[c_idx] + mu_cos_azimuth;
					u_p[chunk_idx*chunk_size + c_idx][p_idx] = -wake_z[c_idx];//wake_velocities.v_x[c_idx];
					b_load[chunk_idx*chunk_size + c_idx][p_idx] = dC_L[c_idx];
					b_load_dot[chunk_idx*chunk_size + c_idx][p_idx] = (dC_L[c_idx] - old_dC_L[c_idx])/dt;
					u_p_dot[chunk_idx*chunk_size + c_idx][p_idx] = -(wake_z[c_idx] - old_wake_z[c_idx])/dt;
					//u_p_dot[chunk_idx*chunk_size + c_idx][p_idx] = (wake_velocities.v_z[c_idx] - old_wake_velocities.v_z[c_idx])/dt;
				}
			}
		}

		double[] wake_age = linspace(0., dt*abs(omega)*wake_history_length.to!double*(180.0/PI.to!double), wake_history_length);
		/+
		foreach(r_idx, ref rotor; ac_timehistory.aircraft_history[0].rotor_states) {
			plt.figure;
			auto z = wake_history[0].rotor_wakes[r_idx][2].get_wake_component!"z";
			//auto x = wake_history[0].rotor_wakes[r_idx][0].get_wake_component!"x";

			writeln("Rotor ", r_idx, " azimuth: ", input_state.rotor_inputs[0].azimuth*(180/PI));

			plt.plot(wake_age, z, charm_wake_hist, charm_wake_z.map!(a => -a/(R*3.281)));
			plt.xlabel("Wake age [rad]");
			plt.ylabel("z/R");
		}+/

		//auto dC_ = ac_timehistory.aircraft_history[0].rotor_states[0].blade_states[0].get_state_array!"dC_T";
		auto blade_loads = ac_timehistory.aircraft_history[0].rotor_states[0].blade_states[0].get_state_array!"dC_T";
		writeln(blade_loads);
		plt.figure;
		plt.plot(r, blade_loads, ["linewidth": 0.5]);
		plt.title("Blade spanwise loading");
		plt.xlabel("r");
		plt.ylabel("$dC_T$");
		if(save) plt.savefig("blade_loading.png", ["dpi": 500]);

		writeln("aircraft.rotors[0].blades[plot_blade].azimuth_offset: ", aircraft.rotors[0].blades[plot_blade].azimuth_offset);
		plt.figure;

		auto az = azimuth_history[$ - (360/d_psi).round.to!size_t..$].map!(a => a - 2*aircraft.rotors[0].blades[plot_blade].azimuth_offset*(PI/180.0));
		auto st = station_loading_1[$ - (360/d_psi).round.to!size_t..$];
		plt.plot(az, st/+, time, station_loading_2+/, ["linewidth": 0.5]);
		plt.legend(["Rotor 1", "Rotor 2"]);

		plt.title("Outer radial station loading");
		plt.xlabel("Time [s]");
		plt.ylabel("$dC_T$");
		if(save) plt.savefig("station_loading_1.png", ["dpi": 500]);

		/+plt.figure;
		plt.plot(time, inflow_1/+, time, inflow_2+/, ["linewidth": 0.5]);
		plt.legend(["Rotor 1", "Rotor 2"]);
		plt.title("Outer radial station inflow");
		plt.xlabel("Time [s]");
		plt.ylabel("$dC_T$");
		if(save) plt.savefig("station_loading_1.png", ["dpi": 500]);

		plt.figure;
		plt.plot(time, inflow_angle_1/+, time, inflow_angle_2+/, ["linewidth": 0.5]);
		plt.legend(["Rotor 1", "Rotor 2"]);
		plt.title("Outer radial station inflow angle");
		plt.xlabel("Time [s]");
		plt.ylabel("$dC_T$");
		if(save) plt.savefig("station_loading_1.png", ["dpi": 500]);

		plt.figure;
		plt.plot(time, u_t_1/+, time, u_t_2+/, ["linewidth": 0.5]);
		plt.legend(["Rotor 1", "Rotor 2"]);
		plt.title("Outer radial station tangential inflow");
		plt.xlabel("Time [s]");
		plt.ylabel("$dC_T$");
		if(save) plt.savefig("station_loading_1.png", ["dpi": 500]);

		plt.figure;
		plt.plot(time, v_x_1/+, time, v_x_2+/, ["linewidth": 0.5]);
		plt.legend(["Rotor 1", "Rotor 2"]);
		plt.title("Outer radial station wake x");
		plt.xlabel("Time [s]");
		plt.ylabel("$dC_T$");
		if(save) plt.savefig("station_loading_1.png", ["dpi": 500]);

		plt.figure;
		plt.plot(time, v_y_1/+, time, v_y_2+/, ["linewidth": 0.5]);
		plt.legend(["Rotor 1", "Rotor 2"]);
		plt.title("Outer radial station wake y");
		plt.xlabel("Time [s]");
		plt.ylabel("$dC_T$");
		if(save) plt.savefig("station_loading_1.png", ["dpi": 500]);

		plt.figure;
		plt.plot(time, v_z_1/+, time, v_z_2+/, ["linewidth": 0.5]);
		plt.legend(["Rotor 1", "Rotor 2"]);
		plt.title("Outer radial station wake z");
		plt.xlabel("Time [s]");
		plt.ylabel("$dC_T$");
		if(save) plt.savefig("station_loading_1.png", ["dpi": 500]);+/

		/+plt.figure;
		plt.gca(["projection": "polar"]);
		plt.title("Tangential Inflow");
		plt.contourf(psis, radial_points, u_t, ["levels": 500]);
		plt.colorbar;
		if(save) plt.savefig("tangential.png", ["dpi": 500]);

		plt.figure;
		plt.gca(["projection": "polar"]);
		plt.title("Tangential Inflow wake only");
		plt.contourf(psis, radial_points, wake_u_t, ["levels": 500]);
		plt.colorbar;
		if(save)  plt.savefig("tangential_wake_only.png", ["dpi": 500]);+/

		plt.figure;
		plt.gca(["projection": "polar"]);
		plt.title("Perpendicular Inflow");
		plt.contourf(psis, radial_points, u_p, ["levels": 100]);
		plt.colorbar;
		if(save)  plt.savefig("perp.png", ["dpi": 500]);

		plt.figure;
		plt.gca(["projection": "polar"]);
		plt.title("Perpendicular Inflow time derivative");
		plt.contourf(psis, radial_points, u_p_dot, ["levels": 100]);
		plt.colorbar;
		if(save)  plt.savefig("perp_dot.png", ["dpi": 500]);

		plt.figure;
		plt.gca(["projection": "polar"]);
		plt.title("Blade loading");
		plt.contourf(psis, radial_points, b_load, ["levels": 100]);
		plt.colorbar;
		if(save)  plt.savefig("b_load.png", ["dpi": 500]);

		plt.figure;
		plt.gca(["projection": "polar"]);
		plt.title("Blade loading time derivative");
		plt.contourf(psis, radial_points, b_load_dot, ["levels": 100]);
		plt.colorbar;
		if(save)  plt.savefig("b_load_dot.png", ["dpi": 500]);

		auto core_size = wake_history[0].rotor_wakes[0][0].get_wake_component!"r_c";
		auto strain_integral = wake_history[0].rotor_wakes[0][0].get_wake_component!"phi";
		//double[] wake_age = linspace(0., dt*omega*core_size.length.to!double, core_size.length);
		plt.figure;
		plt.plot(wake_age, core_size);
		//plt.ylabel("$r_c$");
		//plt.xlabel(`$\xi$`);
		plt.title("Core size");
		if(save)  plt.savefig("core_size.png", ["dpi": 500]);

		plt.figure;
		plt.plot(wake_age, strain_integral);
		//plt.ylabel(`$\phi$`);
		//plt.xlabel(`$\xi$`);
		plt.title("Strain integral");
		if(save)  plt.savefig("strain_integral.png", ["dpi": 500]);

		//plt.xlim(["left": -1], ["right": 3]);
		//plt.xlim(["left": 0], ["right": 1080]);
		//plt.ylim(["bottom": -1], ["top": .2]);

		//plt.xlim(["left": min_dim - 0.05*span], ["right": max_dim + 0.05*span]);
		//plt.ylim(["bottom": min_dim - 0.05*span], ["top": max_dim + 0.05*span]);
		//plt.set_box_aspect([span, span, span]);
		//plt.zlim(["bottom": min_dim - 0.05*span], ["top": max_dim + 0.05*span]);
		//plt.xlim()
		//plt.axis("equal");


		/+plt.figure;

		import std.range : iota;
		foreach(b_idx, ref blade; ac_timehistory.aircraft_history[0].rotor_states[0].blade_states) {
			auto gamma = wake_history[0].rotor_wakes[0][b_idx].get_wake_component!"gamma";
			
			auto age = iota(0, gamma.length);

			plt.plot(age, gamma);
		}+/

		/+plt.figure;
		string[] legend;
		//auto age = iota(0, C_T_history.length);
		//foreach(idx, ref c_t; C_T_history[20..21]) {
			plt.plot(time, C_T_history);
			//legend ~= idx.to!string;
		//}
		+/

		/+foreach(idx, ref c_t; wake_x_history[20..21]) {
			plt.plot(azimuth_history, c_t);
			legend ~= idx.to!string;
		}

		foreach(idx, ref c_t; wake_y_history[20..21]) {
			plt.plot(azimuth_history, c_t);
			legend ~= idx.to!string;
		}

		foreach(idx, ref c_t; wake_u_t_history[20..21]) {
			plt.plot(azimuth_history, c_t);
			legend ~= idx.to!string;
		}+/

		//plt.legend(["C_T", "wake x", "wake y", "wake u_p"]);
		//plt.legend(legend);

		//plt.figure;
		//plt.plot(azimuth_history, inflow_history);

		if(!save) plt.show;
	}

	//write_inflow_vtu("boxwell.vtu", inflows[0], deltas, starts, num_x, num_y, num_z, origins[0], aoa, omega);

	case_list.write_caselist(output_dir);

	return 0;
}
