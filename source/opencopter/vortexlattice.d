module opencopter.vortexlattice;

import opencopter.aircraft;
import opencopter.config;
import opencopter.math;
import opencopter.math.blas;
import opencopter.math.lapacke;
import opencopter.memory;
import opencopter.wake;

import numd.utility;

import std.algorithm;
import std.conv;
import std.math;
import std.math;
import std.range;
import std.stdio : writeln;


double[] generate_spanwise_vortex_nodes(size_t span_n_sections) {
	import std.algorithm : map;
	import std.array : array;
	import std.math : cos, PI;
	import std.range : iota, retro;
	// Spanwise Votex nodes
	immutable num_points = span_n_sections%chunk_size == 0 ? span_n_sections : span_n_sections + (chunk_size - span_n_sections%chunk_size);
	return iota(1.0,num_points + 1.0).map!((l){
		immutable phi = (2*l-1)*PI/(2.0*(num_points.to!double));
		auto span_vortex_pt = 0.25*(1-cos(phi)).to!double;
		return span_vortex_pt;
	}).array;
}


double[] generate_chordwise_votex_nodes(size_t chord_n_sections) {
	import std.algorithm : map;
	import std.array : array;
	import std.math : cos, PI;
	import std.range : iota, retro;
	// Spanwise Votex nodes
	immutable num_points = chord_n_sections%chunk_size == 0 ? chord_n_sections : chord_n_sections + (chunk_size - chord_n_sections%chunk_size);
	return iota(1.0,num_points + 1.0).map!((k){
		immutable theta = (2*k-1)*PI/(2.0*(num_points.to!double));
		auto chord_vortex_pt = 0.5*(1 - cos(theta)).to!double;
		return chord_vortex_pt;
	}).array;
}

extern (C++) struct WingFilamentChunk {
	Chunk x;
	Chunk y;
	Chunk z;
	Chunk x_e;
	Chunk gamma; // Vortex circulation
    Chunk A_kl; // d_gamma/dy|_x
	//Chunk phi; // Time marched wake strain integral
	//Chunk r_c; // vortex core radius.
	//Chunk r_0; // vortex core radius.
	Chunk dx;
	Chunk dy;
	Chunk dz;
    Chunk trail_end;
	//Chunk l_0;
	//Chunk v_z;
	//Chunk volume;
	//Chunk d_volume;
}

extern (C++) struct WingVortexFilamentT(ArrayContainer AC) {

	mixin ArrayDeclMixin!(AC, WingFilamentChunk, "chunks");

	size_t length;

	/*this(size_t wake_history) {

		assert(wake_history % chunk_size == 0);

		immutable num_chunks = wake_history/chunk_size;

		mixin(array_ctor_mixin!(AC, "FilamentChunk", "chunks", "num_chunks"));

		foreach(ref chunk; chunks) {
			chunk.gamma[] = 0;
			chunk.phi[] = 0;
			chunk.r_c = 0.00001;
		}

		length = 0;
	}*/

	this(WingFilamentChunk* slice, size_t slice_length) {
		mixin(array_ctor_mixin_slice!(AC, "WingFilamentChunk", "chunks", "cast(WingFilamentChunk[])slice[0..slice_length]"));

		foreach(ref chunk; chunks) {
			chunk.gamma[] = 0;
			//chunk.phi[] = 0;
			//chunk.r_c = 0.000001;
		}

		length = 0;
	}
}

template is_wing_lift_surface(A) {
	enum bool is_wing_lifting_surface = {
		static if(isPointer!(A)) {
			return isInstanceOf!(WingLiftSurfT, PointerTarget!A);
		} else {
			return isInstanceOf!(WingLiftSurfT, A);
		}
	}();
}

alias WingPartLiftingSurf = WingPartLiftingSurfT!(ArrayContainer.none);

extern(C++) struct WingPartLiftingSurfT(ArrayContainer AC) {
	mixin ArrayDeclMixin!(AC, WingVortexFilamentT!(AC), "spanwise_filaments");
	//mixin ArrayDeclMixin!(AC, WingVortexFilamentT!(AC), "chordwise_filaments");
	
	private size_t spanwise_nodes;
    private size_t chordwise_nodes;

	this(size_t _spanwise_nodes, size_t _chordwise_nodes) {
		spanwise_nodes = _spanwise_nodes;
        chordwise_nodes = _chordwise_nodes;
        immutable actual_num_span_nodes = spanwise_nodes%chunk_size == 0 ? spanwise_nodes : spanwise_nodes + (chunk_size - spanwise_nodes%chunk_size);
        immutable spanwise_chunks = actual_num_span_nodes/chunk_size;
        //immutable chordwise_chunks = chordwise_nodes/chunk_size;
        
		//immutable radial_chunks = radial_elements/chunk_size;

		auto spanwise_filament_mem = new WingFilamentChunk[spanwise_chunks*(chordwise_nodes)];
        //auto chordwise_filament_mem = new FilamentChunk[(spanwise_nodes)*chordwise_chunks];

		mixin(array_ctor_mixin!(AC, "WingVortexFilamentT!(AC)", "spanwise_filaments", "chordwise_nodes"));
        //mixin(array_ctor_mixin!(AC, "VortexFilamentT!(AC)", "chordwise_filaments", "spanwise_nodes"));

		foreach(l, ref spanwise_filament; spanwise_filaments) {
			spanwise_filament = WingVortexFilamentT!AC(&spanwise_filament_mem[l*spanwise_chunks], spanwise_chunks);
		}

        /*foreach(k, ref chordwise_filament; chordwise_filaments) {
			chordwise_filament = VortexFilamentT!AC(&chordwise_filament_mem[k*chordwise_chunks], chordwise_chunks);
		}*/
	}
}

alias WingLiftSurf = WingLiftSurfT!(ArrayContainer.none);

extern(C++) struct WingLiftSurfT(ArrayContainer AC){
    mixin ArrayDeclMixin!(AC, WingPartLiftingSurfT!(AC), "wing_part_lift_surf");

    this(size_t num_wing_parts){
        mixin(array_ctor_mixin!(AC, "WingPartLiftingSurfT!(AC)", "wing_part_lift_surf", "num_wing_parts"));
    }

    ref typeof(this) opAssign(typeof(this) wing_lift_surf) {
		this.wing_part_lift_surf = wing_lift_surf.wing_part_lift_surf;
		return this;
	}

	ref typeof(this) opAssign(ref typeof(this) wing_lift_surf) {
		this.wing_part_lift_surf = wing_lift_surf.wing_part_lift_surf;
		return this;
	}

	ref typeof(this) opAssign(typeof(this)* wing_lift_surf) {
		this.wing_part_lift_surf = wing_lift_surf.wing_part_lift_surf;
		return this;
	}
}

void set_wing_vortex_geometry(WLS,WG)(auto ref WLS wing_lifting_surf, auto ref WG wing_geometry, size_t _spanwise_chunks, size_t _chordwise_nodes){
    
    import std.math : cos,tan, PI;
    import std.math : abs;
    Chunk random_mult_by_minus_one = -1.0;
    double wing_span = wing_geometry.wing_parts[0].wing_span;
    Chunk span = wing_span;
	double root_chord = wing_geometry.wing_parts[0].wing_root_chord;
    auto span_vr_nodes = generate_spanwise_vortex_nodes(_spanwise_chunks*chunk_size);
    //auto span_vr_nodes_left = generate_spanwise_vortex_nodes(_spanwise_chunks*chunk_size);
    //auto chord_vr_nodes = generate_chordwise_votex_nodes(_chordwise_nodes);

    foreach(wp_idx, ref wp_lift_surf; wing_lifting_surf.wing_part_lift_surf){
        foreach(sf_idx, ref spanwise_filament; wp_lift_surf.spanwise_filaments){
            foreach(c_idx;0.._spanwise_chunks){
                if(wp_idx%2 == 0){
                    spanwise_filament.chunks[c_idx].y[]= random_mult_by_minus_one[]*span_vr_nodes[c_idx*chunk_size..c_idx*chunk_size + chunk_size]*span[];
                    writeln("sf_idx = ", sf_idx, "y = ", spanwise_filament.chunks[c_idx].y[]);
                    spanwise_filament.chunks[c_idx].x[]= 0.5*wing_geometry.wing_parts[0].chunks[c_idx].chord[]*(1-cos((2*(sf_idx+1) - 1)*PI/(2.0*_chordwise_nodes))) - spanwise_filament.chunks[c_idx].y[]*tan(wing_geometry.wing_parts[0].le_sweep_angle*PI/180.0);
                    writeln("sf_idx = ", sf_idx, "x = ", spanwise_filament.chunks[c_idx].x[]);
                    spanwise_filament.chunks[c_idx].trail_end[]= spanwise_filament.chunks[c_idx].x[] + span[] + span[];
                } else {
                    spanwise_filament.chunks[c_idx].y[]= span_vr_nodes[c_idx*chunk_size..c_idx*chunk_size + chunk_size]*span[];
                    spanwise_filament.chunks[c_idx].x[]= 0.5*wing_geometry.wing_parts[1].chunks[c_idx].chord[]*(1-cos((2*(sf_idx+1) - 1)*PI/(2.0*_chordwise_nodes))) + spanwise_filament.chunks[c_idx].y[]*tan(wing_geometry.wing_parts[1].le_sweep_angle*PI/180.0);
                    spanwise_filament.chunks[c_idx].trail_end[]= spanwise_filament.chunks[c_idx].x[] + span[] + span[];
                }            
                spanwise_filament.chunks[c_idx].z[]=0.0;
            }        
        }
    }
}

//import std.math : PI;
immutable Chunk one_over_four_pi = 1.0/(4.0*PI);

InducedVelocities compute_horseshoe_vortex_induce_vel(WFC)(auto ref WFC wing_filament_chunks, immutable Chunk x, immutable Chunk y, immutable Chunk z){
    InducedVelocities ret;
    Chunk[chunk_size] v_x, v_y, v_z;

    foreach(n_idx; 0..Chunk.length){
        v_x[n_idx][] = 0;
        v_y[n_idx][] = 0;
        v_z[n_idx][] = 0;
    }

    foreach(i_c_idx, ref chunk_i; wing_filament_chunks){
        
        Chunk x_a; Chunk y_a; Chunk z_a;
        Chunk x_b; Chunk y_b; Chunk z_b;

        Chunk x_c; Chunk x_d;

        Chunk x_amb; Chunk y_amb; Chunk z_amb;

        Chunk x_1ma; Chunk y_1ma; Chunk z_1ma;
        Chunk x_1mb; Chunk y_1mb; Chunk z_1mb;

        Chunk x_dma; Chunk x_cmb;
        Chunk x_1md; Chunk x_1mc;
        
        x_b[] = chunk_i.x[];
        y_b[] = chunk_i.y[];
        z_b[] = chunk_i.z[];

        x_c[] = chunk_i.trail_end[];
        x_d[1..$] = chunk_i.trail_end[0..$-1]; 

        x_a[1..$] = chunk_i.x[0..$-1];
        y_a[1..$] = chunk_i.y[0..$-1];
        z_a[1..$] = chunk_i.z[0..$-1];

        if(i_c_idx == 0){
            x_a[0] = x_a[1] - ((x_a[2]-x_a[1])/(y_a[2]-y_a[1]))*y_a[1];
            y_a[0] = 0;
            z_a[0] = z_a[1] - ((z_a[2]-z_a[1])/(y_a[2]-y_a[1]))*y_a[1];

            x_d[0] = x_d[1];
        } else {
            x_a[0] = wing_filament_chunks[i_c_idx-1].x[$-1];
            y_a[0] = wing_filament_chunks[i_c_idx-1].y[$-1];
            z_a[0] = wing_filament_chunks[i_c_idx-1].z[$-1];

            x_d[0] = wing_filament_chunks[i_c_idx-1].trail_end[$-1];
        }
        /*
        writeln("x_a = ", x_a);
        writeln("x_b = ", x_b);
        writeln("y_a = ", y_a);
        writeln("y_b = ", y_b);
        writeln("z_a = ", z_a);
        writeln("z_b = ", z_b);
        writeln("x_d = ", x_d);
        writeln("x_c = ", x_c);*/

        x_amb = x_a[] - x_b[];
        y_amb = y_a[] - y_b[];
        z_amb = z_a[] - z_b[];

        x_dma = x_d[] - x_a[];
        x_cmb = x_c[] - x_b[];

        Chunk gamma = chunk_i.gamma[];

        foreach(n_idx; 0..chunk_size){
            x_1ma[] = x[n_idx] - x_a[];
            //writeln("x_1ma = ", x_1ma);
            y_1ma[] = y[n_idx] - y_a[];
            z_1ma[] = z[n_idx] - z_a[];

            x_1mb[] = x[n_idx] - x_b[];
            y_1mb[] = y[n_idx] - y_b[];
            z_1mb[] = z[n_idx] - z_b[];

            x_1md[] = x[n_idx] - x_d[];
            //writeln("x_1md = ", x_1md);
            x_1mc[] = x[n_idx] - x_c[];

            //writeln("y_1ma = ", y_1ma);
            //writeln("y_1mb = ", y_1mb);
            //induced velocity due to bound vortex AB
            immutable Chunk fac2_t1_den_square = x_1ma[]*x_1ma[] + y_1ma[]*y_1ma[] + z_1ma[]*z_1ma[];
            immutable Chunk fac2_t2_den_square = x_1mb[]*x_1mb[] + y_1mb[]*y_1mb[] + z_1mb[]*z_1mb[];
            immutable Chunk fac2_t1_den = sqrt(fac2_t1_den_square);
            immutable Chunk fac2_t2_den = sqrt(fac2_t2_den_square);
            //writeln("fac2_t1_den : ", fac2_t1_den[]);
            //writeln("fac2_t2_den : ", fac2_t2_den[]);
            immutable Chunk fac2_t1 = -(x_amb[]*x_1ma[] + y_amb[]*y_1ma[] + z_amb[]*z_1ma[])/fac2_t1_den[];
            immutable Chunk fac2_t2 = -(x_amb[]*x_1mb[] + y_amb[]*y_1mb[] + z_amb[]*z_1mb[])/fac2_t2_den[];
            immutable Chunk fac2_ab = fac2_t1[] - fac2_t2[];
  
            immutable Chunk AB_cross_mag_1 = (y_1ma[]*z_1mb[] - y_1mb[]*z_1ma[])*(y_1ma[]*z_1mb[] - y_1mb[]*z_1ma[]);
            immutable Chunk AB_cross_mag_2 = (x_1ma[]*z_1mb[] - x_1mb[]*z_1ma[])*(x_1ma[]*z_1mb[] - x_1mb[]*z_1ma[]);
            immutable Chunk AB_cross_mag_3 = (x_1ma[]*y_1mb[] - x_1mb[]*y_1ma[])*(x_1ma[]*y_1mb[] - x_1mb[]*y_1ma[]);
            immutable Chunk AB_cross_mag_suare = AB_cross_mag_1[] + AB_cross_mag_2[] + AB_cross_mag_3[];
            //writeln("AB_cross_mag_1 : ", AB_cross_mag_1);
            //writeln("AB_cross_mag_2 : ", AB_cross_mag_2);
            //writeln("AB_cross_mag_3 : ", AB_cross_mag_3);

            immutable Chunk AB_cross_x = (y_1ma[]*z_1mb[] - y_1mb[]*z_1ma[])/AB_cross_mag_suare[];
            immutable Chunk AB_cross_y = -(x_1ma[]*z_1mb[] - x_1mb[]*z_1ma[])/AB_cross_mag_suare[];
            immutable Chunk AB_cross_z = (x_1ma[]*y_1mb[] - x_1mb[]*y_1ma[])/AB_cross_mag_suare[];
            //writeln("AB_cross_mag_x : ", AB_cross_x);
            //writeln("AB_cross_mag_y : ", AB_cross_y);
            //writeln("AB_cross_mag_z : ", AB_cross_z);

            immutable Chunk AB_v_x = AB_cross_x[]*fac2_ab[];
            immutable Chunk AB_v_y = AB_cross_y[]*fac2_ab[];
            immutable Chunk AB_v_z = AB_cross_z[]*fac2_ab[];

            //induced velocity due to vortex AD
            immutable Chunk AD_fac2_t1_den_square = x_1md[]*x_1md[] + y_1ma[]*y_1ma[] + z_1ma[]*z_1ma[];
            immutable Chunk AD_fac2_t2_den_square = x_1ma[]*x_1ma[] + y_1ma[]*y_1ma[] + z_1ma[]*z_1ma[];
            immutable Chunk AD_fac2_t1_den = sqrt(AD_fac2_t1_den_square);
            immutable Chunk AD_fac2_t2_den = sqrt(AD_fac2_t2_den_square);
            //writeln("point_idx : ", n_idx,"\tAD_fac2_t1_den : ", AD_fac2_t1_den[]);
            //writeln("point_idx : ", n_idx,"\tAD_fac2_t2_den : ", AD_fac2_t2_den[]);
            immutable Chunk AD_fac2_t1 = -(x_1md[])/AD_fac2_t1_den[];
            immutable Chunk AD_fac2_t2 = x_1ma[]/AD_fac2_t2_den[];
            immutable Chunk AD_fac2 = x_dma[]*(AD_fac2_t1[] + AD_fac2_t2[]);
            //writeln("point_idx : ", n_idx, "\tAD_fac2_t1 : ", AD_fac2_t1, "\tAD_fac_2_t2 : ", AD_fac2_t2);

            immutable Chunk AD_fac1_den = (z_1ma[]*z_1ma[] + y_1ma[]*y_1ma[])*x_dma[];
            immutable Chunk AD_cross_y = z_1ma[]/AD_fac1_den[];
            immutable Chunk AD_cross_z = -y_1ma[]/AD_fac1_den[];

                // AD_v_x = 0
            immutable Chunk AD_v_y = AD_cross_y[]*AD_fac2[];
            immutable Chunk AD_v_z = AD_cross_z[]*AD_fac2[];

            //induced velocity due to vortex BC
            immutable Chunk BC_fac2_t1_den_square = x_1mc[]*x_1mc[] + y_1mb[]*y_1mb[] + z_1mb[]*z_1mb[];
            immutable Chunk BC_fac2_t2_den_square = x_1mb[]*x_1mb[] + y_1mb[]*y_1mb[] + z_1mb[]*z_1mb[];
            immutable Chunk BC_fac2_t1_den = sqrt(BC_fac2_t1_den_square);
            immutable Chunk BC_fac2_t2_den = sqrt(BC_fac2_t2_den_square);
            immutable Chunk BC_fac2_t1 = -x_1mc[]/BC_fac2_t1_den[];
            immutable Chunk BC_fac2_t2 = x_1mb[]/BC_fac2_t2_den[];
            immutable Chunk BC_fac2 = x_cmb[]*(BC_fac2_t1[] + BC_fac2_t2[]);

            immutable Chunk BC_fac1_den = (z_1mb[]*z_1mb[] + y_1mb[]*y_1mb[])*x_cmb[];
            immutable Chunk BC_cross_y = (z_1mb[])/BC_fac1_den[];
            immutable Chunk BC_cross_z = -(y_1mb[])/BC_fac1_den[];

                //BC_v_x = 0
            immutable Chunk BC_v_y = BC_cross_y[]*BC_fac2[];
            immutable Chunk BC_v_z = BC_cross_z[]*BC_fac2[];

            immutable Chunk circulation = gamma[]*one_over_four_pi[];

            immutable Chunk temp_v_x = circulation[]*AB_v_x[];
            immutable Chunk temp_v_y = circulation[]*(AB_v_y[] + AD_v_y[] - BC_v_y[]);
            immutable Chunk temp_v_z = circulation[]*(AB_v_z[] + AD_v_z[] - BC_v_z[]);

            v_x[n_idx][] += temp_v_x[];
            v_y[n_idx][] += temp_v_y[];
            v_z[n_idx][] += temp_v_z[];

            /*writeln("v_x = ", v_x);
            writeln("v_y = ", v_y);
            writeln("v_z = ", v_z);*/
        }
    }

    foreach(n_idx; 0..Chunk.length){
        ret.v_x[n_idx] = v_x[n_idx].sum;
        ret.v_y[n_idx] = v_y[n_idx].sum;
        ret.v_z[n_idx] = v_z[n_idx].sum; 
    }

    return ret;

}

InducedVelocities compute_wing_induced_vel(WLS)(auto ref WLS wing_lift_surface, Chunk x, Chunk y, Chunk z){
    static import std.math;

    InducedVelocities ret;
    ret.v_x[] = 0.0;
    ret.v_y[] = 0.0;
    ret.v_z[] = 0.0;

    foreach(wp_idx, ref wp_lift_surf; wing_lift_surface.wing_part_lift_surf){
        foreach(fil_idx, ref horseshoe_vortex; wp_lift_surf.spanwise_filaments) {

            auto ind_vel = compute_horseshoe_vortex_induce_vel(horseshoe_vortex.chunks, x, y, z);
            if(wp_idx%2 == 0){
                ret.v_x[] -= ind_vel.v_x[];
                ret.v_y[] -= ind_vel.v_y[];
                ret.v_z[] -= ind_vel.v_z[];
            }else{
                ret.v_x[] += ind_vel.v_x[];
                ret.v_y[] += ind_vel.v_y[];
                ret.v_z[] += ind_vel.v_z[];
            }
        }
    }

    return ret;
}

alias VortexLattice = VortexLatticeT!(ArrayContainer.none);

struct VortexLatticeT(ArrayContainer AC) {

    private Chunk[][] influence_inv;
    private size_t span_elements;
    private size_t chord_elements;
    this(size_t _span_elements, size_t _chord_elements, ref WingPartGeometryT!AC wing) {
        

        span_elements = _span_elements;
        chord_elements = _chord_elements;

        immutable wing_AR = 2*wing.wing_span/(wing.wing_root_chord + wing.wing_tip_chord);
        immutable lamda = wing.wing_tip_chord/wing.wing_root_chord;
        writeln("lamda : ", lamda);
        immutable delta_m = 4*(lamda - 1)/(wing_AR * (lamda +1));
        writeln("delta_m : ", delta_m);
        immutable m_LE = tan(wing.le_sweep_angle*PI/180.0);
        writeln("m_LE : ", m_LE);

        immutable chunks = span_elements*chord_elements/chunk_size;

        auto influence = allocate_dense(span_elements*chord_elements , span_elements*chord_elements);
        double[][] _influence_inv = allocate_dense(span_elements*chord_elements , span_elements*chord_elements);
        //writeln(_influence_inv.length);
        double theta_i; double theta_k;
        double phi_j; double phi_l;
        double W_R; double W_L;
        double W_RC; double W_LC; double W_RC_plus_LC;
        double del_ij;

        foreach(ch1; 0..chunks) {
            foreach(c1; 0..chunk_size) {
                immutable v = ch1*chunk_size + c1;
                immutable j = (v%span_elements); // index of spanwise control point
                immutable i = (v - v%span_elements)/span_elements +1; // index of chordwise control point
                assert(i<= chord_elements, "i is out of bound");
                assert(j<= span_elements, "j is out of bound");

                theta_i = i*PI/chord_elements;
                if(j<span_elements){
                    phi_j = j*PI/span_elements;    
                }
                
                foreach(ch2 ; 0..chunks) {
                    foreach(c2 ; 0..chunk_size) {
                        immutable n = ch2*chunk_size + c2;
                        //writeln(n);
                        immutable l = (n%span_elements)+1; // index of span vortex node
                        //writeln(l);
                        immutable k = (n - n%span_elements)/span_elements +1; // index of chord vortex node

                        assert(k <= chord_elements, "k is out of bound");
                        assert(l <= span_elements, "l is out of bound");

                        theta_k = (2*k - 1)*PI/(2*chord_elements);
                        phi_l = (2*l - 1)*PI/(2*span_elements);

                        double A1 = m_LE + (delta_m/2)*(1-cos(theta_k));
                        double B1 = (cos(theta_k) - cos(theta_i))*(2.0/(wing_AR*(1+lamda)) + delta_m*(1.0-cos(phi_j))/4);
                        double C1 = PI/(wing_AR*(1.0+lamda)*span_elements*chord_elements);

                        double K_RC = 1.0 + sqrt((0.25*(A1*A1 + 1.0)*(1.0-cos(phi_j))*(1.0-cos(phi_j))) + A1*B1*(1.0-cos(phi_j)) + B1*B1) / B1;

                        double K_R =  1.0 + sqrt((0.25*(A1*A1 + 1.0)*(cos(phi_l)-cos(phi_j))*(cos(phi_l)-cos(phi_j))) + A1*B1*(cos(phi_l)-cos(phi_j)) + B1*B1) / B1;

                        double A2 = -A1;
                        double B2 = (m_LE + (delta_m/4)*(2.0-cos(theta_k)-cos(theta_i)))*(1.0-cos(phi_j)) + (2.0/(wing_AR*(1.0+lamda)))*(cos(theta_k) - cos(theta_i));

                        double K_LC = 1.0 + sqrt((0.25*(A2*A2 + 1.0)*(1.0-cos(phi_j))*(1.0-cos(phi_j))) + A2*B2*(1.0-cos(phi_j)) + B2*B2) / B2;
                        double K_L = 1.0 + sqrt((0.25*(A2*A2 + 1.0)*(2.0-cos(phi_l)-cos(phi_j))*(2.0-cos(phi_l)-cos(phi_j))) + A2*B2*(2.0-cos(phi_l)-cos(phi_j)) + B2*B2) / B2;

                        W_R = -C1 * K_R * sin(theta_k)/(cos(phi_l)-cos(phi_j));
                        W_L = C1 * K_L * sin(theta_k)/(2.0-cos(phi_l)-cos(phi_j));


                        if(j == 0){
                            if(cos(theta_k) - cos(theta_i) > 0)
                                del_ij = 1.0;
                                else{
                                    del_ij = -1.0;
                                }
                            double test_term1 = (PI/(4.0*span_elements*chord_elements))* del_ij * sin(theta_k);
                            double test_term2 = (m_LE + (delta_m/2)*(1.0-cos(theta_i)));
                            W_RC_plus_LC = (PI/(4.0*span_elements*chord_elements))* del_ij * sin(theta_k) * (2.0*(m_LE + (delta_m/2)*(1.0-cos(theta_i)))/(cos(theta_k)-cos(theta_i)) - delta_m);
                            //writeln("test_term_1 : ", test_term1);
                            //writeln("test_term_2 : ", test_term2);
                            //writeln("k = ", k, "\tl = ", l, "\ti = ", i, "\tW_RC+LC : ", W_RC_plus_LC);
                            influence[v][n] = W_R + W_L + W_RC_plus_LC;
                        }
                        else{
                            W_RC = C1 * K_RC * sin(theta_k)/(1.0-cos(phi_j));
                            W_LC = -C1 * K_LC * sin(theta_k)/(1.0-cos(phi_j));
                            influence[v][n] = W_R + W_RC + W_L + W_LC;
                        }

                        /*if(i==2 && j==7 && k==2 && l==8){
                            writeln("A1 =", A1);
                            writeln("B1 =", B1);
                            writeln("B2 =", B2);
                            writeln("C1 =", C1);
                            writeln("K_L =", K_L);
                            writeln("K_R =", K_R);
                            writeln("W_R =", W_R);
                            writeln("W_L =", W_L);
                            writeln("W_LC+RC = ", W_RC_plus_LC);
                            writeln("Phi_l = ", phi_l);
                            writeln("inf = ", influence[v][n]);
                        }*/
                    }
                }
                //writeln("inf = ", influence[v][]);
            }
        }
        
        // check this carefully on wednesday (today)!!!!
        immutable total_elements = span_elements*chord_elements;
        foreach(r_idx; 0..total_elements) {
			_influence_inv[r_idx][] = influence[r_idx][];
		}

		openblas_set_num_threads(1);

		int info = 0;
		auto ipiv = new int[total_elements];
		info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, total_elements.to!int, total_elements.to!int, _influence_inv[0].ptr, total_elements.to!int, ipiv.ptr);
		assert(info == 0, "Failed to invert influence matrix");
		info = LAPACKE_dgetri(LAPACK_ROW_MAJOR, total_elements.to!int, _influence_inv[0].ptr, total_elements.to!int, ipiv.ptr);

		assert(info == 0, "Failed to invert influence matrix");


		influence_inv = allocate_dense_chunk_aliased(total_elements, total_elements);

		foreach(r; 0..total_elements) {
			foreach(ch; 0..chunks) {
				influence_inv[r][ch][] = _influence_inv[r][ch*chunk_size..ch*chunk_size+chunk_size];
			}
		}
        //writeln("inf = ", influence);
    }

    void compute_d_gamma_coefficients(WLS,WS)(auto ref WLS wing_lift_surface, auto ref WS wing_part_state, size_t wp_idx, size_t span_chunk_idx, size_t chord_node_idx, immutable Chunk u){ 
        foreach(c1; 0..chunk_size){
            
            double A_kl = 0;
            size_t span_chunks_length = wing_part_state.chunks.length;
            //size_t num_half_filaments = wing_part_state.ctrl_chunks.length/span_chunks_length;
            immutable r = span_chunk_idx*chunk_size + c1 + span_chunks_length*chunk_size*chord_node_idx;
            //writeln("span_chunk_length = ",span_chunks_length);
            //immutable Chunk tmp;
            foreach(ch, ref inf; influence_inv[r]) {
                //writeln("inf = ", inf[]);
                Chunk tmp = -inf[]*sin(wing_part_state.ctrl_chunks[ch].ctrl_pt_aoa)[];
                //writeln("Chunk_tmp = ", tmp);
                A_kl += tmp.sum;
            }
            //writeln("span_idx= ", span_chunk_idx*chunk_size + c1, "\tChord_node_idx = ", chord_node_idx, "\tA_kl = ", A_kl);
            A_kl *= u[c1];
                wing_lift_surface.wing_part_lift_surf[wp_idx].spanwise_filaments[chord_node_idx].chunks[span_chunk_idx].A_kl[c1] = A_kl;
        }
    }

    void compute_bound_circulation(WLS,WPG)(auto ref WLS wing_lift_surface, auto ref WPG wing_part, size_t wp_idx, size_t span_chunk_idx, size_t chord_node_idx){
        double gamma = 0.0;
        size_t num_span_chunks = wing_part.chunks.length;
        //size_t num_half_filaments = wing_part.ctrl_chunks.length/num_span_chunks;
        //writeln(num_span_chunks);
        foreach(c1; 0..chunk_size){
            double tmp_gamma = wing_lift_surface.wing_part_lift_surf[wp_idx].spanwise_filaments[chord_node_idx].chunks[span_chunk_idx+1..$]
                .map!(c => c.A_kl.sum).sum;
            writeln("wing_part = ", wp_idx, "\tchord_idx = ", chord_node_idx, "\tspan_chunk =", span_chunk_idx ,"temp_gamma = ", tmp_gamma);    
            gamma = tmp_gamma + wing_lift_surface.wing_part_lift_surf[wp_idx].spanwise_filaments[chord_node_idx].chunks[span_chunk_idx].A_kl[c1..$].sum;
            //writeln("span_idx= ", span_chunk_idx*chunk_size + c1, "\tChord_node_idx = ", chord_node_idx, "\tgamma = ", gamma);
            wing_lift_surface.wing_part_lift_surf[wp_idx].spanwise_filaments[chord_node_idx].chunks[span_chunk_idx].gamma[c1] = -PI*gamma/(num_span_chunks*chunk_size*wing_part.chunks[span_chunk_idx].chord[c1]);  
            gamma = 0.0;
        }         
    }

    void compute_dCl(WLS,WPS)(auto ref WLS wing_lift_surface, auto ref WPS wing_part_state, size_t wp_idx, size_t span_chunk_idx, Chunk u1){
        double dCl = 0;
        size_t num_span_chunks = wing_part_state.chunks.length;
        //size_t num_half_filaments = wing_part_state.ctrl_chunks.length/num_span_chunks;
        size_t num_chord_pt = wing_part_state.ctrl_chunks.length/num_span_chunks;
        
        foreach(c1; 0..chunk_size){
            foreach(n_idx; 0..num_chord_pt){
                double theta_n = (2*n_idx +1)*PI/(2*num_chord_pt);
                immutable dCl_inter = PI*(wing_lift_surface.wing_part_lift_surf[wp_idx].spanwise_filaments[n_idx].chunks[span_chunk_idx].gamma[c1])*sin(theta_n)/num_chord_pt;
                dCl += dCl_inter;
            }
            wing_part_state.chunks[span_chunk_idx].dC_L[c1]= dCl/u1[c1];
            dCl = 0.0;            
        }
    }
}


