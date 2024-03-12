module opencopter.inflow.huangpeters.iteration;

import opencopter.inflow.huangpeters;

import std.conv;

package size_t iterate_odds_sin(alias f, Args...)(long Mo, size_t idx, auto ref Args args) {
	foreach(long _m; 1..Mo + 1) {
		foreach(long _n; odd_states[Mo][_m]) {
			f(_m.to!long, _n.to!long, idx, args);
			idx++;
		}
	}

	return idx;
}

package size_t iterate_evens_sin(alias f, Args...)(long Me, size_t idx, auto ref Args args) {
	foreach(long _m; 1..Me + 1) {
		foreach(long _n; even_states[Me][_m]) {
			f(_m.to!long, _n.to!long, idx, args);
			idx++;
		}
	}

	return idx;
}

size_t iterate_even_odd_sin(alias f, Args...)(long Mo, long Me, size_t start_row_idx, size_t start_col_idx, auto ref Args args) {
	return iterate_evens_sin!((r, j, row_idx) {
		iterate_odds_sin!(
			(m, n, col_idx) => f(r, j, m, n, row_idx, col_idx, args)
		)(Mo, start_col_idx);
	})(Me, start_row_idx);
}

size_t iterate_odd_even_sin(alias f, Args...)(long Mo, long Me, size_t start_row_idx, size_t start_col_idx, auto ref Args args) {
	return iterate_odds_sin!((r, j, row_idx) {
		iterate_evens_sin!(
			(m, n, col_idx) => f(r, j, m, n, row_idx, col_idx, args)
		)(Mo, start_col_idx);
	})(Me, start_row_idx);
}

size_t iterate_odd_odd_sin(alias f, Args...)(long Mo, long start_row_idx, size_t start_col_idx, auto ref Args args) {
	return iterate_odds_sin!((r, j, row_idx) {
		iterate_odds_sin!(
			(m, n, col_idx) => f(r, j, m, n, row_idx, col_idx, args)
		)(Mo, start_col_idx);
	})(Mo, start_row_idx);
}

size_t iterate_even_even_sin(alias f, Args...)(long Me, long start_row_idx, size_t start_col_idx, auto ref Args args) {
	return iterate_evens_sin!((r, j, row_idx) {
		iterate_evens_sin!(
			(m, n, col_idx) => f(r, j, m, n, row_idx, col_idx, args)
		)(Me, start_col_idx);
	})(Me, start_row_idx);
}

void iterate_whole_matrix_sin(alias f_oo, alias f_oe, alias f_eo, alias f_ee, Args...)(long Mo, long Me, auto ref Args args) {

	auto _row_idx = iterate_odds_sin!((r, j, row_idx) {
		auto _col_idx = iterate_odds_sin!(
			(m, n, col_idx) => f_oo(r, j, m, n, row_idx, col_idx, args)
		)(Mo, 0);

		iterate_evens_sin!(
			(m, n, col_idx) => f_oe(r, j, m, n, row_idx, col_idx, args)
		)(Me, _col_idx);
	})(Mo, 0);

	iterate_evens_sin!((r, j, row_idx) {
		auto _col_idx = iterate_odds_sin!(
			(m, n, col_idx) => f_eo(r, j, m, n, row_idx, col_idx, args)
		)(Mo, 0);

		iterate_evens_sin!(
			(m, n, col_idx) => f_ee(r, j, m, n, row_idx, col_idx, args)
		)(Me, _col_idx);
	})(Me, _row_idx);
}

package size_t iterate_odds(alias f, Args...)(long Mo, size_t idx, auto ref Args args) {
	foreach(long _m; 0..Mo + 1) {
		foreach(long _n; odd_states[Mo][_m]) {
			f(_m.to!long, _n.to!long, idx, args);
			idx++;
		}
	}

	return idx;
}

package size_t iterate_evens(alias f, Args...)(long Me, size_t idx, auto ref Args args) {
	foreach(long _m; 0..Me + 1) {
		foreach(long _n; even_states[Me + 0][_m]) {
			f(_m.to!long, _n.to!long, idx, args);
			idx++;
		}
	}

	return idx;
}

size_t iterate_even_odd(alias f, Args...)(long Mo, long Me, size_t start_row_idx, size_t start_col_idx, auto ref Args args) {
	return iterate_evens!((r, j, row_idx) {
		iterate_odds!(
			(m, n, col_idx) => f(r, j, m, n, row_idx, col_idx, args)
		)(Mo, start_col_idx);
	})(Me, start_row_idx);
}

size_t iterate_odd_even(alias f, Args...)(long Mo, long Me, size_t start_row_idx, size_t start_col_idx, auto ref Args args) {
	return iterate_odds!((r, j, row_idx) {
		iterate_evens!(
			(m, n, col_idx) => f(r, j, m, n, row_idx, col_idx, args)
		)(Mo, start_col_idx);
	})(Me, start_row_idx);
}

size_t iterate_odd_odd(alias f, Args...)(long Mo, long start_row_idx, size_t start_col_idx, auto ref Args args) {
	return iterate_odds!((r, j, row_idx) {
		iterate_odds!(
			(m, n, col_idx) => f(r, j, m, n, row_idx, col_idx, args)
		)(Mo, start_col_idx);
	})(Mo, start_row_idx);
}

size_t iterate_even_even(alias f, Args...)(long Me, long start_row_idx, size_t start_col_idx, auto ref Args args) {
	return iterate_evens!((r, j, row_idx) {
		iterate_evens!(
			(m, n, col_idx) => f(r, j, m, n, row_idx, col_idx, args)
		)(Me, start_col_idx);
	})(Me, start_row_idx);
}

void iterate_whole_matrix(alias f_oo, alias f_oe, alias f_eo, alias f_ee, Args...)(long Mo, long Me, auto ref Args args) {

	auto _row_idx = iterate_odds!((r, j, row_idx) {
		auto _col_idx = iterate_odds!(
			(m, n, col_idx) => f_oo(r, j, m, n, row_idx, col_idx, args)
		)(Mo, 0);

		iterate_evens!(
			(m, n, col_idx) => f_oe(r, j, m, n, row_idx, col_idx, args)
		)(Me, _col_idx);
	})(Mo, 0);

	iterate_evens!((r, j, row_idx) {
		auto _col_idx = iterate_odds!(
			(m, n, col_idx) => f_eo(r, j, m, n, row_idx, col_idx, args)
		)(Mo, 0);

		iterate_evens!(
			(m, n, col_idx) => f_ee(r, j, m, n, row_idx, col_idx, args)
		)(Me, _col_idx);
	})(Me, _row_idx);
}
