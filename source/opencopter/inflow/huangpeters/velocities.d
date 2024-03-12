module opencopter.inflow.huangpeters.velocities;

import opencopter.inflow.huangpeters;
import opencopter.inflow.huangpeters.iteration;
import opencopter.inflow.huangpeters.math;
import opencopter.inflow.huangpeters.system;

import opencopter.config;
import opencopter.math;
import opencopter.memory;

import std.algorithm;
import std.conv;
import std.math;
import std.range;
import std.stdio;

package auto compute_velocities_nh(ArrayContainer AC, T)(HuangPetersInflowT!AC infl, T[] coefficients, T[] coeffiecients_sin, immutable ElipticalCoords coords) {

	Chunk V = 0;

	size_t sin_idx = 0;
	auto _idx = iterate_odds!(
		(m, n, idx) {
			immutable Chunk cos_mpsi = infl.mpsi_buff[m];
			immutable Chunk Pmn = infl.P_nh_mn_bar[m][n];

			immutable Chunk alpha = coefficients[idx];
			if(m == 0) {
				immutable Chunk Qmn_cos_mpsi = infl.Qmn_bar[m][m+1][]*cos_mpsi[];
				V[] += alpha[]*Pmn[]*Qmn_cos_mpsi[];
			} else {
				immutable Chunk sin_mpsi = infl.sin_mpsi_buff[m];
				immutable Chunk beta = coeffiecients_sin[sin_idx];
				immutable Chunk Qmn_bar = infl.Qmn_bar[m][m+1][];
				sin_idx++;
				immutable Chunk a_cos = alpha[]*cos_mpsi[];
				immutable Chunk b_sin = beta[]*sin_mpsi[];

				V[] += Pmn[]*Qmn_bar[]*(a_cos[] + b_sin[]);
			}
		}
	)(infl.Mo, 0);

	iterate_evens!(
		(m, n, idx) {
			immutable Chunk cos_mpsi = infl.mpsi_buff[m];
			immutable Chunk Pmn = infl.P_md_mn_bar[m][n];

			immutable Chunk alpha = coefficients[idx];
			if(m == 0) {
				immutable Chunk Qmn_cos_mpsi = infl.Qmn_bar[m][m+1][]*cos_mpsi[];
				V[] += alpha[]*Pmn[]*Qmn_cos_mpsi[];
			} else {
				immutable Chunk sin_mpsi = infl.sin_mpsi_buff[m];
				immutable Chunk beta = coeffiecients_sin[sin_idx];
				immutable Chunk Qmn_bar = infl.Qmn_bar[m][m+1][];
				sin_idx++;
				immutable Chunk a_cos = alpha[]*cos_mpsi[];
				immutable Chunk b_sin = beta[]*sin_mpsi[];

				V[] += Pmn[]*Qmn_bar[]*(a_cos[] + b_sin[]);
			}
		}
	)(infl.Me, _idx);

	return V;
}

package auto compute_velocities_md(ArrayContainer AC, T)(HuangPetersInflowT!AC infl, T[] coefficients, T[] coeffiecients_sin, immutable ElipticalCoords coords) {

	Chunk V = 0;

	size_t sin_idx = 0;
	auto _idx = iterate_odds!(
		(m, n, idx) {
			immutable Chunk cos_mpsi = infl.mpsi_buff[m];
			immutable Chunk Pmn = infl.P_md_mn_bar[m][n];
			immutable Chunk alpha = coefficients[idx];

			if(m == 0) {
				immutable Chunk Qmn_cos_mpsi = infl.Qmn_bar[m][n][]*cos_mpsi[];
				V[] += alpha[]*Pmn[]*Qmn_cos_mpsi[];
			} else {
				immutable Chunk sin_mpsi = infl.sin_mpsi_buff[m];
				immutable Chunk beta = coeffiecients_sin[sin_idx];
				immutable Chunk Qmn_bar = infl.Qmn_bar[m][n][];
				sin_idx++;
				V[] += Pmn[]*Qmn_bar[]*(alpha[]*cos_mpsi[] + beta[]*sin_mpsi[]);
			}

		}
	)(infl.Mo, 0);

	iterate_evens!(
		(m, n, idx) {
			immutable Chunk cos_mpsi = infl.mpsi_buff[m];
			immutable Chunk Pmn = infl.P_md_mn_bar[m][n];
			immutable Chunk alpha = coefficients[idx];

			if(m == 0) {
				immutable Chunk Qmn_cos_mpsi = infl.Qmn_bar[m][n][]*cos_mpsi[];
				V[] += alpha[]*Pmn[]*Qmn_cos_mpsi[];
			} else {
				immutable Chunk sin_mpsi = infl.sin_mpsi_buff[m];
				immutable Chunk beta = coeffiecients_sin[sin_idx];
				immutable Chunk Qmn_bar = infl.Qmn_bar[m][n][];
				sin_idx++;
				V[] += Pmn[]*Qmn_bar[]*(alpha[]*cos_mpsi[] + beta[]*sin_mpsi[]);
			}
		}
	)(infl.Me, _idx);

	return V;
}

package auto compute_velocities_bl(ArrayContainer AC, T)(HuangPetersInflowT!AC infl, T[] coefficients_md, T[] coefficients_nh, T[] coefficients_md_sin, T[] coefficients_nh_sin, immutable CartisianCoords ccoords, immutable ElipticalCoords coords) {

	// Recomended buffer zone in Huang's dissertation
	immutable double eps = 0.01;

	immutable Chunk h =  zip(coords.eta[], eps.repeat).map!(a => a[0] < a[1] ? 0 : a[0] - a[1]).staticArray!Chunk;
	immutable abs_y = abs(ccoords.y);
	
	immutable Chunk b = zip(ccoords.x[], abs_y[], coords.eta[], infl.sin_chi.repeat).map!((xye) {
		if(xye[0].isNaN || xye[1].isNaN || xye[2].isNaN || xye[3].isNaN) {
			return 0;
		}
		if(xye[0] <= 0.0 && xye[1] <= 1.0) {
			return 20.0*(1.0 - xye[1]*xye[1]*xye[3]/(1.0 + xye[2]*xye[2]));
		} else if(xye[0] <= 0.0 && xye[1] > 1.0) {
			return 20.0*(1.0 - xye[1]*xye[1]*xye[3]/((1.0 + xye[2]*xye[2]) + 0.615*(xye[1]*xye[1] - 1.0)));
		} else if(xye[0] > 0.0 && xye[1] <= 1.0) {
			return 20.0*(1.0 - (xye[0]*xye[0] + xye[1]*xye[1])*xye[3]/(1.0 + xye[2]*xye[2]));
		} else if(xye[0] > 0.0 && xye[1] > 1.0) {
			return 20.0*(1.0 - (xye[0]*xye[0] + xye[1]*xye[1])*xye[3]/((1.0 + xye[2]*xye[2]) + 0.615*(xye[1]*xye[1] - 1.0)));
		}
		assert(false);
	}).staticArray!Chunk;

	immutable Chunk blend_denom = 1.0 + b[]*h[];

	immutable Chunk blend_nh = 1.0/blend_denom[];
	immutable Chunk blend_md = b[]*h[]/blend_denom[];
	
	immutable bool compute_nh = blend_nh[].map!(a => !a.isClose(0.0)).fold!((res, a) => res |= a)(false);
	immutable bool compute_md = blend_md[].map!(a => !a.isClose(0.0)).fold!((res, a) => res |= a)(false);

	auto _idx = iterate_odds!(
			(m, n, idx) {
				if(compute_md) {
					immutable Chunk Pmn_md = associated_legendre_polynomial(m, n, coords.nu, infl.P_coefficients[idx]);
					infl.P_md_mn_bar[m][n][] = Pmn_md[];
				}

				if(compute_nh) {
					immutable Chunk Pmn_nh = associated_legendre_polynomial_nh(m, n, coords.nu, infl.P_coefficients_nh[idx]);
					infl.P_nh_mn_bar[m][n][] = Pmn_nh[];
				}
			}
	)(infl.Mo, 0);

	iterate_evens!(
		(m, n, idx) {
			immutable Chunk Pmn_md = associated_legendre_polynomial(m, n, coords.nu, infl.P_coefficients[idx]);
			infl.P_md_mn_bar[m][n][] = Pmn_md[];
		}
	)(infl.Me, _idx);

	associated_legendre_function(coords.eta, infl.Qmn_bar, infl.K_table);

	foreach(m; 0..max(infl.Mo, infl.Me) + 1) {
		if(m == 0) {
			infl.mpsi_buff[m] = 1.0;
			infl.sin_mpsi_buff[m] = 0.0;
		} else {
			immutable Chunk mpsi = m.to!double*coords.psi[];
			immutable Chunk[2] sin_cos = sincos(mpsi);

			infl.mpsi_buff[m] = sin_cos[1];
			infl.sin_mpsi_buff[m] = sin_cos[0];
		}
	}

	immutable Chunk V_md = compute_md ? compute_velocities_md(infl, coefficients_md, coefficients_md_sin, coords) : Z;
	immutable Chunk V_nh = compute_nh ? compute_velocities_nh(infl, coefficients_nh, coefficients_nh_sin, coords) : Z;
	
	immutable Chunk V = blend_nh[]*V_nh[] + blend_md[]*V_md[];

	return V;
}

@nogc package auto final_blend(double cos_chi, double sin_chi, immutable Chunk y, immutable Chunk sigma, immutable Chunk z, immutable Chunk x, immutable Chunk s_0) {

	immutable g = 1.84*sqrt(cos_chi) - 4.06*cos_chi + 11.84*cos_chi^^(1.5);
	immutable double sin_xi_2 = sin_chi*sin_chi;

	immutable Chunk f =
		zip(sigma[], sin_xi_2.repeat, g.repeat, y[], z[])
		.map!(
			(sxgy) {
				if((sxgy[0] < 0.0) || (sxgy[1] <= 1.0e-14) || ((sxgy[0] < 0.0) && (abs(sxgy[4]) <= 1.0e-14))) {
					return 0.0;
				} else {
					if(abs(sxgy[3]) <= 1.0) {
						return sxgy[1]/(sxgy[1] + sxgy[0]*sxgy[2]);
					} else {
						return sxgy[1]/(sxgy[1] + (sxgy[0] + 1.5*sqrt(sxgy[3]*sxgy[3] - 1.0))*sxgy[2]);
					}
				}
				assert(false);
			}
		)
		.staticArray!Chunk;
	return f;
}

package auto compute_velocities_final(ArrayContainer AC, T)(HuangPetersInflowT!AC infl, T[] a, T[] alpha, T[] delta, T[] lambda, T[] b, T[] beta, T[] delta_s, T[] lambda_s, immutable CartisianCoords ccoords, Chunk t, immutable Chunk x_0, immutable Chunk z_0, bool above_disk) {

	immutable auto coords = to_eliptical(ccoords);

	immutable Chunk rho_axial = 1;

	immutable Chunk y2z2 = ccoords.y[]*ccoords.y[] + ccoords.z[]*ccoords.z[];
	immutable Chunk rho2 = rho_axial[]*rho_axial[];
	immutable Chunk s_0 = zip(y2z2[], rho2[]).map!(a => a[0] < a[1] ? sqrt(a[1] - a[0]) : 0).staticArray!Chunk;
	
	immutable Chunk neg_s_0 = -s_0[];
	immutable Chunk neg_y = -ccoords.y[];
	immutable Chunk sigma = -ccoords.x[] - s_0[];
	immutable Chunk sigma_p_s = sigma[] + s_0[];

	Chunk f = final_blend(infl.cos_chi, infl.sin_chi, ccoords.y, sigma, ccoords.z, ccoords.x, s_0);
	immutable bool compute_ds = f[].map!(a => !a.isClose(0.0)).fold!((res, a) => res |= a)(false);
	
	immutable Chunk V_bl = compute_velocities_bl(infl, a, alpha, b, beta, ccoords, coords);
	Chunk V_ds = 0.0;
	if(compute_ds) {

		double devisor;
		devisor = infl.advance_ratio;

		Chunk t_delay = t[] - sigma[]*infl.sin_chi;

		foreach(c_idx, ref t_d; t_delay) {
			if(t_d > infl.times[infl.get_circular_index(infl.curr_state)]) {
				while(t_d > infl.times[infl.get_circular_index(infl.curr_state)]) {
					if(abs(f[c_idx]) > 1.0e-14) {
						writeln("Dealing with the future: ", ccoords.x[c_idx], ", ", ccoords.y[c_idx], ", ", ccoords.z[c_idx], ", ", ccoords.x[c_idx]^^2.0+ccoords.y[c_idx]^^2.0 + ccoords.z[c_idx]^^2.0, ", f: ", f[c_idx]);
					}
					
					t_d -= 2.0*PI;
				}
			}
		}

		//t_delay[] = t_delay[]/devisor;

		infl.interpolate_to_time(t_delay, infl.time_delay_alpha_2, infl.time_delay_a_2, infl.time_delay_beta_2, infl.time_delay_b_2);

		immutable CartisianCoords[3] ccoords_ds = [
			CartisianCoords(neg_s_0, ccoords.y, ccoords.z),
			CartisianCoords(s_0, neg_y, ccoords.z),
			CartisianCoords(sigma_p_s, neg_y, ccoords.z)
		];

		immutable ElipticalCoords[3] coords_ds = [
			to_eliptical(ccoords_ds[0]),
			to_eliptical(ccoords_ds[1]),
			to_eliptical(ccoords_ds[2])
		];

		immutable Chunk v_ds_bl = compute_velocities_bl(infl,
			infl.time_delay_a_2[0..infl.total_states],
			infl.time_delay_alpha_2[0..infl.total_states],

			infl.time_delay_b_2[0..infl.total_sin_states],
			infl.time_delay_beta_2[0..infl.total_sin_states],
			
			ccoords_ds[0],
			coords_ds[0]
		);
		immutable Chunk v_ds_bl_a_1 = compute_velocities_bl(infl,
			infl.time_delay_a_2[infl.total_states..$],
			infl.time_delay_alpha_2[infl.total_states..$],

			infl.time_delay_b_2[infl.total_sin_states..$],
			infl.time_delay_beta_2[infl.total_sin_states..$],
			
			ccoords_ds[1],
			coords_ds[1]
		);
		immutable Chunk v_ds_bl_a_2 = compute_velocities_bl(infl,
			delta,
			lambda,
			delta_s,
			lambda_s,
			ccoords_ds[2],
			coords_ds[2]
		);

		V_ds[] = v_ds_bl[] + v_ds_bl_a_1[] - v_ds_bl_a_2[];
	}
	Chunk one_m_f = 1.0 - f[];

	immutable Chunk V_ds_f = V_ds[]*f[];
	immutable Chunk V = V_bl[]*one_m_f[] + V_ds_f[];

	return V;
}

package auto compute_velocities_final_adjoint(ArrayContainer AC, T)(HuangPetersInflowT!AC infl, T[] a, T[] alpha, T[] delta, T[] lambda, T[] b, T[] beta, T[] delta_s, T[] lambda_s, immutable CartisianCoords ccoords, Chunk t, immutable Chunk x_0, immutable Chunk z_0) {

	immutable auto coords = to_eliptical(ccoords);

	immutable Chunk rho_axial = 1;

	immutable Chunk y2z2 = ccoords.y[]*ccoords.y[] + ccoords.z[]*ccoords.z[];
	immutable Chunk rho2 = rho_axial[]*rho_axial[];
	immutable Chunk s_0 = zip(y2z2[], rho2[]).map!(a => a[0] < a[1] ? sqrt(a[1] - a[0]) : 0).staticArray!Chunk;

	immutable Chunk neg_s_0 = -s_0[];
	immutable Chunk neg_y = -ccoords.y[];
	immutable Chunk sigma = -ccoords.x[] - s_0[];
	immutable Chunk sigma_p_s = sigma[] + s_0[];

	Chunk f = final_blend(infl.cos_chi, infl.sin_chi, ccoords.y, sigma, ccoords.z, ccoords.x, s_0);

	immutable bool compute_ds = f[].map!(a => !a.isClose(0.0)).fold!((res, a) => res |= a)(false);

	immutable Chunk V_bl = compute_velocities_bl(infl, delta, lambda, delta_s, lambda_s, ccoords, coords);
	Chunk V_ds = 0.0;

	if(compute_ds) {
		double devisor;
		devisor = infl.advance_ratio;

		Chunk t_delay = t[] + sigma[]*infl.sin_chi;

		foreach(c_idx, ref t_d; t_delay) {
			if(t_d > infl.times[infl.get_circular_index(infl.curr_state)]) {
				while(t_d > infl.times[infl.get_circular_index(infl.curr_state)]) {
					t_d -= 2.0*PI;
				}
			}
		}

		infl.interpolate_to_time(t_delay, infl.time_delay_alpha_2, infl.time_delay_a_2, infl.time_delay_beta_2, infl.time_delay_b_2);

		immutable CartisianCoords[3] ccoords_ds = [
			CartisianCoords(neg_s_0, ccoords.y, ccoords.z),
			CartisianCoords(s_0, neg_y, ccoords.z),
			CartisianCoords(sigma_p_s, neg_y, ccoords.z)
		];

		immutable ElipticalCoords[3] coords_ds = [
			to_eliptical(ccoords_ds[0]),
			to_eliptical(ccoords_ds[1]),
			to_eliptical(ccoords_ds[2])
		];

		immutable Chunk v_ds_bl_a = compute_velocities_bl(infl,
			infl.time_delay_a_2[infl.total_states..$],
			infl.time_delay_alpha_2[infl.total_states..$],
			
			infl.time_delay_b_2[infl.total_sin_states..$],
			infl.time_delay_beta_2[infl.total_sin_states..$],

			ccoords_ds[0],
			coords_ds[0]
		);

		immutable Chunk v_ds_bl_1 = compute_velocities_bl(infl,
			infl.time_delay_a_2[0..infl.total_states],
			infl.time_delay_alpha_2[0..infl.total_states],
			
			infl.time_delay_b_2[0..infl.total_sin_states],
			infl.time_delay_beta_2[0..infl.total_sin_states],

			ccoords_ds[1],
			coords_ds[1]
		);
		immutable Chunk v_ds_bl_2 = compute_velocities_bl(infl,
			a,
			alpha,
			b,
			beta,
			ccoords_ds[2],
			coords_ds[2]
		);

		V_ds[] = v_ds_bl_a[] + v_ds_bl_1[] - v_ds_bl_2[];
	}
	
	Chunk one_m_f = 1.0 - f[];

	immutable Chunk V_ds_f = V_ds[]*f[];
	immutable Chunk V = V_bl[]*one_m_f[] + V_ds_f[];

	return V;
}

package auto inflow_at_impl(ArrayContainer AC, C)(HuangPetersInflowT!AC infl, auto ref C x, auto ref C y, auto ref C z) {
	
	static import opencopter.math;

	bool all_nan = x[].map!(a => a.isNaN).fold!((res, a) => res && a)(true);
	if(all_nan) {
		Chunk V = 0;
		return V;
	}

	immutable all_below_disk = z[].map!(a => a > 0).fold!((res, a) => a && res)(true);
	immutable all_above_or_on_disk = z[].map!(a => (a <= 0)).fold!((res, a) => a && res)(true);

	immutable all_somewhere = !all_above_or_on_disk && !all_below_disk;

	import core.stdc.stdio : printf;

	CartisianCoords coords = CartisianCoords(x, y, z);

	Chunk V_above = 0;
	Chunk V_below = 0;

	if(all_above_or_on_disk || all_somewhere) {
		Chunk t = infl.times[infl.get_circular_index(infl.curr_state)];

		V_above = compute_velocities_final(infl, infl.a, infl.alpha, infl.delta, infl.lambda, infl.b, infl.beta, infl.delta_s, infl.lambda_s, coords, t, x, z, true);
		if(all_above_or_on_disk && !all_somewhere) {
			return V_above;
		}
	}
	
	if(all_below_disk || all_somewhere) {	
		immutable tan_chi = tan(infl.chi);
		immutable Chunk x_offset = z[]*tan_chi;
		immutable Chunk x0_1 = x[] + x_offset[];
		immutable Chunk y0_1 = y[];
		immutable Chunk x0_2 = -x0_1[];
		immutable Chunk y0_2 = -y[];
		immutable Chunk zero = 0;

		Chunk s = z[]*z[] + x_offset[]*x_offset[];
		s = sqrt(s);
		immutable Chunk t = infl.times[infl.get_circular_index(infl.curr_state)];

		Chunk t_minus = s[];
		Chunk t_delay = t[] - t_minus[];

		foreach(c_idx, ref t_d; t_delay) {
			if(t_d > infl.times[infl.get_circular_index(infl.curr_state)]) {

				while(t_d > infl.times[infl.get_circular_index(infl.curr_state)]) {
					t_d -= 2.0*PI;
				}
			}
		}

		infl.interpolate_to_time(t_delay, infl.time_delay_alpha, infl.time_delay_a, infl.time_delay_beta, infl.time_delay_b);

		CartisianCoords coords1 = CartisianCoords(x0_1, y0_1, zero);
		CartisianCoords coords2 = CartisianCoords(x0_2, y0_2, zero);

		immutable Chunk v_f = compute_velocities_final(infl,
			infl.time_delay_a[0..infl.total_states],
			infl.time_delay_alpha[0..infl.total_states],
			infl.time_delay_a[infl.total_states..$],
			infl.time_delay_alpha[infl.total_states..$],

			infl.time_delay_b[0..infl.total_sin_states],
			infl.time_delay_beta[0..infl.total_sin_states],
			infl.time_delay_b[infl.total_sin_states..$],
			infl.time_delay_beta[infl.total_sin_states..$],

			coords1,
			t_delay,
			x,
			z,
			false
		);

		immutable Chunk v_f_a_1 = compute_velocities_final_adjoint(infl,
			infl.time_delay_a[0..infl.total_states],
			infl.time_delay_alpha[0..infl.total_states],
			infl.time_delay_a[infl.total_states..$],
			infl.time_delay_alpha[infl.total_states..$],
			
			infl.time_delay_b[0..infl.total_sin_states],
			infl.time_delay_beta[0..infl.total_sin_states],
			infl.time_delay_b[infl.total_sin_states..$],
			infl.time_delay_beta[infl.total_sin_states..$],

			coords2,
			t_delay,
			x,
			z
		);

		immutable Chunk neg_z = -z[];
		
		immutable Chunk x3 = -x[];

		CartisianCoords coords3 = CartisianCoords(x3, y0_2, neg_z);

		immutable Chunk v_f_a_2 = compute_velocities_final_adjoint(infl,
			infl.a,
			infl.alpha,
			infl.delta,
			infl.lambda,

			infl.b,
			infl.beta,
			infl.delta_s,
			infl.lambda_s,

			coords3,
			t,
			x,
			z
		);

		V_below = v_f[] + v_f_a_1[] - v_f_a_2[];
		
		if(all_below_disk && !all_somewhere) {
			return V_below;
		}
	}

	// Merge results
	if(all_somewhere) {
		Chunk V = 0;
		foreach(c_idx; 0..chunk_size) {
			if(z[c_idx] > 0) {
				V[c_idx] = V_below[c_idx];
			} else {
				V[c_idx] = V_above[c_idx];
			}
		}
		return V;
	}

	assert(false);
}
