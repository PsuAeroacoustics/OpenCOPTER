module opencopter.vortexlattice;

import opencopter.aircraft;
import opencopter.config;
import opencopter.math;
import opencopter.math.blas;
import opencopter.math.lapacke;
import opencopter.memory;

import numd.utility;

import std.algorithm;
import std.conv;
import std.math;
import std.math : sin, cos, PI;
import std.range;


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

	this(FilamentChunk* slice, size_t slice_length) {
		mixin(array_ctor_mixin_slice!(AC, "WingFilamentChunk", "chunks", "cast(FilamentChunk[])slice[0..slice_length]"));

		foreach(ref chunk; chunks) {
			chunk.gamma[] = 0;
			//chunk.phi[] = 0;
			//chunk.r_c = 0.000001;
		}

		length = 0;
	}
}

extern(C++) struct WingLiftingSurfT(ArrayContainer AC) {
	mixin ArrayDeclMixin!(AC, WingVortexFilamentT!(AC), "spanwise_filaments");
	//mixin ArrayDeclMixin!(AC, WingVortexFilamentT!(AC), "chordwise_filaments");
	
	private size_t spanwise_nodes;
    private size_t chordwise_nodes;

	this(size_t _spanwise_nodes, size_t _chordwise_nodes) {
		spanwise_nodes = _spanwise_nodes;
        chordwise_nodes = _chordwise_nodes;
        immutable spanwise_chunks = spanwise_nodes/chunk_size;
        //immutable chordwise_chunks = chordwise_nodes/chunk_size;
        
		//immutable radial_chunks = radial_elements/chunk_size;

		auto spanwise_filament_mem = new WingFilamentChunk[2*spanwise_chunks*(chordwise_nodes)];
        //auto chordwise_filament_mem = new FilamentChunk[(spanwise_nodes)*chordwise_chunks];

		mixin(array_ctor_mixin!(AC, "WingVortexFilamentT!(AC)", "spanwise_filaments", "2*chordwise_nodes"));
        //mixin(array_ctor_mixin!(AC, "VortexFilamentT!(AC)", "chordwise_filaments", "spanwise_nodes"));

		foreach(l, ref spanwise_filament; spanwise_filaments) {
			spanwise_filament = VortexFilamentT!AC(&spanwise_filament_mem[l*spanwise_chunks], spanwise_chunks);
		}

        /*foreach(k, ref chordwise_filament; chordwise_filaments) {
			chordwise_filament = VortexFilamentT!AC(&chordwise_filament_mem[k*chordwise_chunks], chordwise_chunks);
		}*/
	}
}

void set_wing_vortex_geometry(WLS,WG)(auto ref WLS wing_lifting_surf, auto ref WG wing_geometry, size_t _spanwise_chunks, size_t _chordwise_nodes){
    
    import std.math : cos,tan, PI;
    import std.math : abs;
    auto span_vr_nodes_right = generate_spanwise_vortex_nodes(_spanwise_chunks*chunk_size);
    auto span_vr_nodes_left = -generate_spanwise_vortex_nodes(_spanwise_chunks*chunk_size);
    span_vr_nodes_left = span_vr_nodes_left.retro;
    //auto chord_vr_nodes = generate_chordwise_votex_nodes(_chordwise_nodes);

    foreach(sf_idx, ref spanwise_filament; wing_lifting_surf.spanwise_filaments){
        foreach(c_idx;_spanwise_chunks){
            if(sf_idx < spnawise_chunks*_chordwise_nodes){
                spanwise_filament.y[]=wing_geometry.wing_span*span_vr_nodes_left[c_idx*chunk_size..c_idx*chunk_size + chunk_size];
                spanwise_filament.x[]=0.5*wing_geometry.chunks[c_idx].chord[]*(1-cos((2*(sf_idx+1) - 1)*PI/_chordwise_nodes)) + abs(spanwise_filament.y)[]*tan(wing_geometry.le_sweep_angle*PI/180.0);
            } else {
                spanwise_filament.y[]=wing_geometry.wing_span*span_vr_nodes_right[c_idx*chunk_size..c_idx*chunk_size + chunk_size];
                spanwise_filament.x[]=0.5*wing_geometry.chunks[c_idx].chord[]*(1-cos((2*(sf_idx+1 -_chordwise_nodes) - 1)*PI/_chordwise_nodes)) + abs(spanwise_filament.y)[]*tan(wing_geometry.le_sweep_angle*PI/180.0);
            }            
            spanwise_filament.z[]=0.0;
        }        
    }
}

import std.math : PI;
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
        
        Chunk x_a, Chunk y_a, Chunk z_a;
        Chunk x_b, Chunk y_b, Chunk z_b;

        Chunk x_c, Chunk x_d;

        Chunk x_amb, Chunk y_amb, Chunk z_xmb;

        Chunk x_1ma, Chunk y_1ma, Chunk z_1ma;
        Chunk x_1mb, Chunk y_1mb, Chunk z_1mb;

        Chunk x_dma, Chunk x_cmb;
        
        x_a[] = chunk_i.x[];
        y_a[] = chunk_i.y[];
        z_a[] = chunk_i.z[];

        x_d[] = chunk_i.trail_end[];
        x_c[0..$-1] = chunk_i.trail_end[1..$]; 

        x_b[0..$-1] = chunk_i.x[1..$];
        y_b[0..$-1] = chunk_i.y[1..$];
        z_b[0..$-1] = chunk_i.z[1..$];

        if(i_c_idx == chunks.length - 1){
            x_b[$-1] = x_b[$-2];
            y_b[$-1] = y_b[$-2];
            z_b[$-1] = z_b[$-2];

            x_c[$-1] = x_c[$-2];
        } else {
            x_b[$-1] = chunks[i_c_idx+1].x[0];
            y_b[$-1] = chunks[i_c_idx+1].y[0];
            z_b[$-1] = chunks[i_c_idx+1].z[0];

            x_c[$-1] = chunks.[i_c_idx+1].trail_end[0];
        }

        x_amb = x_a[] - x_b[];
        y_amb = y_a[] - y_b[];
        z_amb = z_a[] - z_b[];

        x_dma = x_d[] - x_a[];
        x_cmb = x_c[] - x_b[];

        Chunk gamma = chunk_i.gamma[];

        foreach(n_idx; 0..chunk_size){
            x_1ma[] = x[n_idx] - x_a[];
            y_1ma[] = y[n_idx] - y_a[];
            z_1ma[] = z[n_idx] - z_a[];

            x_1mb[] = x[n_idx] - x_b[];
            y_1mb[] = y[n_idx] - y_b[];
            z_1mb[] = z[n_idx] - z_b[];

            x_1md[] = x[n_idx] - x_d[];
            x_1mc[] = x[n_idx] - x_c[];

            //induced velocity due to bound vortex AB
            immutable Chunk fac2_t1_den = sqrt(x_1ma[]*x_1ma[] + y_1ma[]*y_1ma[] + z_1ma[]*z_1ma[]);
            immutable Chunk fac2_t2_den = sqrt(x_1mb[]*x_1mb[] + y_1mb[]*y_1mb[] + z_1mb[]*z_1mb[]);
            immutable Chunk fac2_t1 = (x_amb[]*x_1ma[] + y_amb[]*y_1ma[] + z_amb[]*z_1ma[])/fac2_t1_den[];
            immutable Chunk fac2_t2 = (x_amb[]*x_1mb[] + y_amb[]*y_1mb[] + z_amb[]*z_1mb[])/fac2_t2_den[];
            immutable Chunk fac2_ab = fac2_t1[] - fac2_t2[];

            immutable Chunk AB_cross_mag_1 = (y_1ma[]*z_1mb[] - y_1mb[]*z_1ma[])*(y_1ma[]*z_1mb[] - y_1mb[]*z_1ma[]);
            immutable Chunk AB_cross_mag_2 = (x_1ma[]*z_1mb[] - x_1mb[]*z_1ma[])*(x_1ma[]*z_1mb[] - x_1mb[]*z_1ma[]);
            immutable Chunk AB_cross_mag_3 = (x_1ma[]*y_1mb[] - x_1mb[]*y_1ma[])*(x_1ma[]*y_1mb[] - x_1mb[]*y_1ma[]);
            immutable Chunk AB_cross_mag_suare = AB_cross_mag_1[] + AB_cross_mag_2[] + AB_cross_mag_3[];

            immutable Chunk AB_cross_x = (y_1ma[]*z_1mb[] - y_1mb[]*z_1ma[])/AB_cross_mag_suare[];
            immutable Chunk AB_cross_y = -(x_1ma[]*z_1mb[] - x_1mb[]*z_1ma[])/AB_cross_mag_suare[];
            immutable Chunk AB_cross_z = (x_1ma[]*y_1mb[] - x_1mb[]*y_1ma[])/AB_cross_mag_suare[];

            immutable Chunk AB_v_x = AB_cross_x[]*fac2_ab[];
            immutable Chunk AB_v_y = AB_cross_y[]*fac2_ab[];
            immutable Chunk AB_v_z = AB_cross_z[]*fac2_ab[];

            //induced velocity due to vortex AD
            immutable Chunk AD_fac2_t1_den = sqrt(x_1md[]*x_1md[] + y_1ma[]*x_1ma[] + z_1ma[]*z_1ma[]);
            immutable Chunk AD_fac2_t2_den = sqrt(x_1ma[]*x_1ma[] + y_1ma[]*x_1ma[] + z_1ma[]*z_1ma[]);
            immutable Chunk AD_fac2_t1 = -(x_1md[])/AD_fac2_t1_den[];
            immutable Chunk AD_fac2_t2 = x_1ma[]/AD_fac2_t2_den[];
            immutable Chunk AD_fac2 = x_dma[]*(AD_fac2_t1[] + AD_fac2_t2[]);

            immutable Chunk AD_fac1_den = (z_1ma[]*z_1ma[] + y_1ma[]*y_1ma[])*x_dma[];
            immutable Chunk AD_cross_y = z_1ma[]/AD_fac1_den[];
            immutable Chunk AD_cross_z = -y_1ma[]/AD_fac1_den[];

                // AD_v_x = 0
            immutable Chunk AD_v_y = AD_cross_y[]*AD_fac2[];
            immutable Chunk AD_v_z = AD_cross_z[]*AD_fac2[];

            //induced velocity due to vortex BC
            immutable Chunk BC_fac2_t1_den = sqrt(x_1mc[]*x_1mc[] + y_1mb[]*y_1mb[] + z_1mb[]*z_1mb[]);
            immutable Chunk BC_fac2_t2_den = sqrt(x_1mb[]*x_1mb[] + y_1mb[]*y_1mb[] + z_1mb[]*z_1mb[]);
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
        }
    }

    foreach(n_idx; 0..Chunk.lenght){
        ret.v_x[n_idx] = v_x[n_idx].sum;
        ret.v_y[n_idx] = v_y[n_idx].sum;
        ret.v_z[n_idx] = v_z[n_idx].sum; 
    }

    return ret;

}


struct VortexLattice(ArrayContainer AC){
    private Chunk[][] Inf_ceff_inv;
    private size_t span_elements;
    private size_t chord_elements;
    this(size_t _span_elements, size_t _chord_elements, ref WingGeometryT!AC wing) {
        import std.stdio : writeln;

        span_elements = _span_elements;
        chord_elements = _chord_elements;

        immutable wing_AR = 2* wing.wing_span/(wing.wing_root_chord + wing.wing_tip_chord);
        immutable lamda = wing.wing_tip_chord/wing.wing_root_chord;
        immutable delta_m = 4*(lamda - 1)/(wing_AR * (lamda +1));
        immutable m_LE = math.tan(wing.le_sweep_angle*PI/180.0);

        immutable chunks = span_elements*chord_elements/chunk_size;

        auto influence = allocate_dense(span_elements*chord_elements , span_elements*chord_elements);

        foreach(ch1; 0..chunks) {
            foreach(c1; 0..chunk_size) {
                immutable v = ch1*chunk_size + c1;
                immutable i = (v - v%span_elements)+1; // index of chordwise control point
                immutable j = v%span_elements; // index of spanwise control point
                assert(i<= chord_elements, "i is out of bound");
                assert(j<= span_elements, "j is out of bound");

                double theta_i = i*PI/chord_elements;
                if(j<spam_elements-1){
                    double phi_j = j*PI/span_elements;    
                }
                
                foreach(ch2 ; 0..chunks) {
                    foreach(c2 ; 0..chunk_size) {
                        immutable n = ch2*chunk_size + c2;
                        immutable k = (n - n%span_elements)+1; // index of chordwise vortex node
                        immutable l = n%span_elements+1; // index of spanwise vortex node

                        assert(k<= chord_elements, "k is out of bound");
                        assert(l<= span_elements, "l is out of bound");

                        double theta_k = (2*k - 1)*PI/(2*chord_elements);
                        double phi_l = (2*l - 1)*PI/(2*span_elements);

                        double A1 = m_LE + (delta_m/2)*(1-cos(theta_k));
                        double B1 = (cos(theta_k) - cos(theta_i))*( 2/(wing_AR*(1+lamda)) + delta_m*(1-cos(phi_j))/4 );
                        double C1 = PI/(wing_AR*(1+lamda)*span_elements*chord_elements);

                        double K_RC = 1 + (0.25*(A1^2 + 1)*(1-cos(phi_j))^2 + A1*B1*(1-cos(phi_j)) + B1^2)^0.5 / B1;

                        double K_R =  1 + (0.25*(A1^2 + 1)*(cos(phi_l)-cos(phi_j))^2 + A1*B1*(cos(phi_l)-cos(phi_j)) + B1^2)^0.5 / B1;

                        double A2 = -A1;
                        double B2 = (m_LE + (delta_m/4)*(2-cos(theta_k)-cos(theta_i)))*1-cos(phi_j) + (2/(wing_AR*(1+lamda)))*(cos(theta_k) - cos(theta_i));

                        double K_LC = 1 + (0.25*(A2^2 + 1)*(1-cos(phi_j))^2 + A2*B2*(1-cos(phi_j)) + B2^2)^0.5 / B2;
                        double K_L = 1 + (0.25*(A2^2 + 1)*(2-cos(phi_l)-cos(phi_j))^2 + A2*B2*(2-cos(phi_l)-cos(phi_j)) + B2^2)^0.5 / B2;

                        W_R = -C1 * K_R * sin(theta_k)/(cos(theta_l)-cos(theta_j));
                        W_L = C1 * K_L * sin(theta_k)/(1-cos(theta_l)-cos(theta_j));

                        if(j == 0){
                            if(cos(theta_k) - cos(theta_i) > 0)
                                double del_ij = 1;
                                else{
                                    double del_ij = -1;
                                }

                            double W_RC_plus_LC = PI/(4*span_elements*chord_elements)* del_ij * sin(theta_k) * ((2*(m_LE + (delta_m/2)*(1-cos(theta_i)))/(cos(theta_k)-cos(theta_i))) - delta_m);
                            influence[v][n] = W_R + W_L + W_RC_plus_LC;
                        }
                        else{
                            double W_RC = C1 * K_RC * sin(theta_k)/(1-cos(phi_j));
                            double W_LC = -C1 * K_LC * sin(theta_k)/(1-cos(phi_j));
                            influence[v][n] = W_R + W_RC + W_L + W_LC;
                        }
                    }
                }
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
    }

    void compute_d_gamma_coefficients(WLS,WS)(auto ref WLS wing_lift_surface, auto ref WS wing_state, size_t span_chunk_idx, size_t chord_node_idx){
        foreach(c1; 0..chunk_size){
            immutable r = span_chunk_idx*chunk_size + c1 + wing_lift_surface.spanwise_filaments[0].lenght*chunk_size*chord_node_idx;
            double A_kl = 0;
            double num_filaments = wing_lift_surface.spanwise_filaments.length;
            foreach(w_ch, ref wing_part_state; wing_state.wing_part_states){
                foreach(ch, ref inf; influence_inv[r]) {
                    Chunk tmp = inf[]*sin(wing_part_state.ctrl_chunk[ch].ctrl_pt_aoa)[];
                    A_kl += tmp.sum;
                }
                A_kl*= u[c1];
                if(mod(w_ch,2)==0){
                    wing_lift_surface.spanwise_filaments[chord_node_idx].chunks[wing_lift_surface.spanwise_filaments[0].length-span_chunk_idx].A_kl[chunk_size-c1] = A_kl;
                }else{
                    wing_lift_surface.spanwise_filaments[0.5*num_filaments+chord_node_idx].chunks[span_chunk_idx].A_kl[c1] = A_kl;
                }
                
            }
        }
    }

    void compute_bound_circulation(WLS)(auto ref WLS wing_lift_surface, size_t span_chunk_idx, size_t chord_node_idx){
        double gamma = 0;
        double num_filaments = wing_lift_surface.spanwise_filaments.length;
        double num_span_chunks = wing_lift_surface.spanwise_filaments[0].length;
        foreach(c1; 0..chunk_size){
            Chunk tmp_gamma_right = wing_lift_surface.spanwise_filaments[0.5*num_filaments+chord_node_idx].chunks[span_chunk_idx+1..$].sum;
            Chunk tmp_gamma_left = wing_lift_surface.spanwise_filaments[chord_node_idx].chunks[0..num_span_chunks-span_chunk_idx-1].sum;

            gamma_right = tmp_gamma_right.sum + wing_lift_surface.spanwise_filaments[0.5*num_filaments+chord_node_idx].chunks[span_chunk_idx+1..$].A_kl[c1..$].sum;
            gamma_left = tmp_gamma_left.sum + wing_lift_surface.spanwise_filaments[chord_node_idx].chunks[0..num_span_chunks-span_chunk_idx].A_kl[0..chunk_size-c1].sum;
            wing_lift_surface.spanwise_filaments[0.5*num_filaments+chord_node_idx].chunks[span_chunk_idx].gamma[c1] = gamma_right;
            wing_lfit_surface.spanwise_filaments[chord_node_idx].chunks[num_span_chunks-span_chunk_idx].gamma[chunk_size-c1] = gamma_left;
            gamma_left = 0;
            gamma_right = 0;
        }         
    }
}

