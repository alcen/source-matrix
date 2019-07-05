/*
List of functions:
first_row(mx)               : matrix -> matrix
row_down(mx)                : matrix -> matrix
col_to_row(mx)              : matrix -> row
col_left(mx)                : matrix -> matrix
width(mx)                   : matrix -> integer
is_empty_matrix(mx)         : matrix -> boolean
is_square(mx)               : matrix -> boolean
empty_mat()                 : -> matrix
draw_matrix(mx)             : matrix -> matrix, displays matrix
draw_row(row)               : row -> matrix, displays matrix
scalar_pdt(s, mx)           : number, matrix -> matrix
matrix_pdt(mx1, mx2)        : matrix, matrix -> matrix
augment_matrix(mx1, mx2)    : matrix, matrix -> matrix
transpose(mx)               : matrix -> matrix
del_row(mx, i)              : matrix, integer -> matrix
del_col(mx, j)              : matrix, integer -> matrix
make_zero(r, c)             : integer, integer -> matrix
make_identity(n)            : integer -> matrix
build_row(n, fx)            : integer, function(j) -> matrix
build_col(n, fx)            : integer, function(i) -> matrix
build_matrix(r, c, fx)      :
	integer, integer, function(i, j) -> matrix

matrix_from_list(r, c, xs)  : integer, integer, list -> matrix
row_sum(row1, row2)         : row, row -> row
matrix_sum(mx1, mx2)        : matrix, matrix -> matrix
add_row(mx, row, n)         : matrix, row, integer -> matrix
extract_row(mx, i)          : matrix, integer -> row
extract_col(mx, j)          : matrix, integer -> matrix
replace_row(mx, i, row)     : matrix, integer, row -> matrix
matrix_r_swap(mx, i1, i2)   : matrix, integer, integer -> matrix
matrix_r_times(mx, k, i)    : matrix, number, integer -> matrix
matrix_r_add(mx, i1, k, i2) :
	matrix, integer, number, integer -> matrix

ref(mx)                     : matrix -> matrix
rref(mx)                    : matrix -> matrix
del_n_rows(mx, n)           : matrix, integer -> matrix
del_n_cols(mx, n)           : matrix, integer -> matrix
inverse(mx)                 : matrix -> matrix
submatrix(mx, i, j)         : matrix, integer, integer -> matrix
cofactor(mx, i, j)          : matrix, integer, integer -> number
determinant(mx)             : matrix -> number
cof_matrix(mx)              : matrix -> matrix
matrix_adj(mx)              : matrix -> matrix

arr_to_mx(r, c, arr)        : integer, integer, array -> matrix
row_vector(arr)             : array -> vector (row matrix)
col_vector(arr)             : array -> vector (column matrix)
least_squares(mx)           : matrix -> vector (column matrix)
eigen_3x3(mx)               : 3x3 matrix -> vector (row matrix)
eig3_arr(arr)               : array -> vector (row matrix)
*/

// Function expressions

const first_row = head;

const row_down = tail;

const col_to_row = mx => map(row => head(row), mx);
// col_to_row returns first column of a matrix mx as a row

const col_left = mx => map(row => tail(row), mx);

const width = mx => length(head(mx));


// Function statements

function is_empty_matrix(mx) {
	return is_empty_list(mx) || is_empty_list(head(mx));
}

function is_square(mx) {
	return length(mx) === width(mx);
}

function empty_mat() {
	return [];
}

function draw_matrix(mx) {
	function accu_last(f, g, ini, xs) {
		if (is_empty_list(xs)) {
			return ini;
		} else if (is_empty_list(tail(xs))) {
			return g(head(xs), ini);
		} else {
			return f(head(xs), accu_last(f, g, ini, tail(xs)));
		}
	}
	const row_to_string = (x, y) => x + "   " + y;
	const last = (x, y) => x + y;
	const draw_row = row => accu_last(row_to_string, last, "", row);
	if (is_empty_matrix(mx)) {
		display("\n");
	} else {
		display("| " + draw_row(first_row(mx)) + " |");
		draw_matrix(row_down(mx));
	}
	return mx;
}

// Test
// draw_matrix(list(list(1, 2, 3), list(4, 5, 6), list(7, 8, 9)));

function draw_row(row) {
	return draw_matrix(pair(row, []));
}

function scalar_pdt(s, mx) {
	const pdt_func = x => s * x;
	const row_pdt = row => map(pdt_func, row);
	return map(row_pdt, mx);
}

// Test
// const matrix1 = list(list(1, 4, 7), list(2, 5, 8), list(3, 6, 9));
// draw_matrix(scalar_pdt(2, matrix1));

function matrix_pdt(mx1, mx2) {
	const cols_m1 = width(mx1);
	const rows_m2 = length(mx2);
	const cols_m2 = width(mx2);

	function row_col_pdt(row, col_r) {
		if (is_empty_matrix(row) && is_empty_matrix(col_r)) {
			return 0;
		} else {
			const others = row_col_pdt(tail(row), tail(col_r));
			return head(row) * head(col_r) + others;
		}
	}

	function row_builder(row, mx, new_row, len_c) {
		if (len_c <= 0) {
			return reverse(new_row);
		} else {
			const entry = row_col_pdt(row, col_to_row(mx));
			return row_builder(row, col_left(mx),
				pair(entry, new_row), len_c - 1);
		}
	}

	if (cols_m1 !== rows_m2) {
		error("wrong matrix dimensions");
	} else if (cols_m1 === 0) {
		return [];
	} else {
		const new_rows = row => row_builder(row, mx2, [], cols_m2);
		return map(new_rows, mx1);
	}
}

// Test
// const matrix3 = list(list(1, 2, 3), list(4, 5, 6));
// const matrix4 = list(list(9), list(8), list(7));
// draw_matrix(matrix_pdt(matrix3, matrix4));

function augment_matrix(mx1, mx2) {
	if (is_empty_matrix(mx1)) {
		return mx2;
	} else if (is_empty_matrix(mx2)) {
		return mx1;
	} else {
		const row1 = first_row(mx1);
		const row2 = first_row(mx2);
		const my1 = row_down(mx1);
		const my2 = row_down(mx2);
		return pair(append(row1, row2), augment_matrix(my1, my2));
	}
}

function transpose(mx) {
	function trans_iter(orig, my) {
		if (is_empty_matrix(orig)) {
			return reverse(my);
		} else {
			const add_row = pair(col_to_row(orig), my);
			return trans_iter(col_left(orig), add_row);
		}
	}
	return trans_iter(mx, []);
}

function del_row(mx, i) {
	if (is_empty_matrix(mx)) {
		error("deleting row of empty matrix");
	} else if (is_empty_matrix(row_down(mx)) || i <= 0) {
		return i === 1 ? [] : error("row outside matrix");
	} else if (i === 1) {
		return row_down(mx);
	} else {
		return pair(first_row(mx), del_row(row_down(mx), i - 1));
	}
}

function del_col(mx, j) {
	function col_iter(my, z, len) {
		if (len <= 0) {
			return z <= 0 ? [] : error("column outside matrix");
		} else if (z === 1) {
			return col_iter(col_left(my), z - 1, len - 1);
		} else {
			const others = col_left(my);
			return pair(col_to_row(my),
				col_iter(others, z - 1, len - 1));
		}
	}
	if (is_empty_matrix(mx)) {
		error("deleting column of empty matrix");
	} else if (is_empty_matrix(col_left(mx)) || j <= 0) {
		return j === 1 ? [] : error("column outside matrix");
	} else if (j === 1) {
		return col_left(mx);
	} else {
		return transpose(col_iter(mx, j, width(mx)));
	}
}

function make_zero(r, c) {
	function zero_iter(mx, n_rows) {
		if (n_rows <= 0) {
			return mx;
		} else {
			const ith_row = build_list(c, x => 0);
			return zero_iter(pair(ith_row, mx), n_rows - 1);
		}
	}
	return zero_iter([], r);
}

function make_identity(n) {
	function iden_iter(mx, n_rows) {
		if (n_rows <= 0) {
			return mx;
		} else {
			const diagonal = x => x === (n_rows - 1) ? 1 : 0;
			const ith_row = build_list(n, diagonal);
			return iden_iter(pair(ith_row, mx), n_rows - 1);
		}
	}
	return iden_iter([], n);
}

function build_row(n, fx) {
	return list(build_list(n, x => fx(x + 1)));
}

function build_col(n, fx) {
	return transpose(list(build_list(n, x => fx(x + 1))));
}

function build_matrix(r, c, fx) {
	function pair_map(fy, mx) {
		const pair_f = p => fy(head(p), tail(p));
		const map_row = row => map(pair_f, row);
		return map(map_row, mx);
	}
	if (r <= 0 || c <= 0) {
		return empty_mat();
	} else {
		const cli = build_col(c, x => x);
		const row_f = s => first_row(
			build_row(r, n => pair(head(s), n)));
		const mx1 = map(row_f, cli);
		display("| building matrix...");
		draw_matrix(mx1);
		return pair_map(fx, mx1);
	}
}

function matrix_from_list(r, c, xs) {
	function row_maker(n, ys, row, mx) {
		if (is_empty_list(ys)) {
			const rev = reverse(row);
			return pair(rev, mx);
		} else if (n <= 0) {
			const rev = reverse(row);
			return row_maker(c, ys, [], pair(rev, mx));
		} else {
			return row_maker(n - 1, tail(ys),
				pair(head(ys), row), mx);
		}
	}
	if (is_empty_list(xs)) {
		return empty_mat();
	} else if (length(xs) !== r * c) {
		error("matrix from list failed: wrong dimensions");
	} else {
		return reverse(row_maker(c, xs, [], empty_mat()));
	}
}

function row_sum(row1, row2) {
	function double_add(xs, ys) {
		if (is_empty_list(ys)) {
			return xs;
		} else {
			const new_item = head(xs) + head(ys);
			const others = double_add(tail(xs), tail(ys));
			return pair(new_item, others);
		}
	}
	if (length(row1) !== length(row2)) {
		error("wrong row length");
	} else {
		return double_add(row1, row2);
	}
}

function matrix_sum(mx1, mx2) {
	function double_sum(my1, my2) {
		if (is_empty_matrix(my2)) {
			return [];
		} else {
			const new_row = row_sum(first_row(my1), first_row(my2));
			const others = double_sum(row_down(my1), row_down(my2));
			return pair(new_row, others);
		}
	}
	if (length(mx1) !== length(mx2)) {
		error("wrong number of rows");
	} else {
		return double_sum(mx1, mx2);
	}
}

function add_row(mx, row, n) {
	if (is_empty_matrix(mx) || n <= 0) {
		return pair(row, mx);
	} else {
		return pair(first_row(mx),
			add_row(row_down(mx), row, n - 1));
	}
}

function extract_row(mx, i) {
	if (is_empty_matrix(mx) || i <= 0) {
		return empty_mat();
	} else if (i === 1) {
		return first_row(mx);
	} else {
		return extract_row(row_down(mx), i - 1);
	}
}

function extract_col(mx, j) {
	if (is_empty_matrix(mx) || j <= 0) {
		return empty_mat();
	} else if (j === 1) {
		return transpose(pair(col_to_row(mx), []));
	} else {
		return extract_col(col_left(mx), j - 1);
	}
}

function replace_row(mx, i, row) {
	if (is_empty_matrix(mx) || i < 1) {
		error("cannot replace row outside of matrix");
	} else if (i === 1) {
		return add_row(row_down(mx), row, 0);
	} else {
		return pair(first_row(mx),
			replace_row(row_down(mx), i - 1, row));
	}
}

function matrix_r_swap(mx, i1, i2) {
	if (is_empty_matrix(mx) || i1 <= 0 || i2 <= 0) {
		error("row swap performed on empty matrix");
	} else {
		const row_i1 = extract_row(mx, i1);
		const row_i2 = extract_row(mx, i2);
		const rep_i1 = replace_row(mx, i1, row_i2);
		return replace_row(rep_i1, i2, row_i1);
	}
}

function matrix_r_times(mx, k, i) {
	const k_x = x => k * x;
	const k_row = row => map(k_x, row);
	if (is_empty_matrix(mx) || i <= 0) {
		return mx;
	} else if (i === 1) {
		return pair(k_row(first_row(mx)), row_down(mx));
	} else {
		return pair(first_row(mx),
			matrix_r_times(row_down(mx), k, i - 1));
	}
}

function matrix_r_add(mx, i1, k, i2) {
	const k_x = x => k * x;
	const k_row = row => map(k_x, row);
	const num_rows = length(mx);
	function m_add_r(mx, i1, row) {
		if (is_empty_matrix(mx)) {
			return mx;
		} else if (i1 === 1) {
			return pair(row_sum(first_row(mx), row), row_down(mx));
		} else {
			return pair(first_row(mx),
				m_add_r(row_down(mx), i1 - 1, row));
		}
	}
	if (is_empty_matrix(mx) || i1 <= 0 || i2 <= 0) {
		error("row add performed on empty matrix");
	} else if (i1 > num_rows || i2 > num_rows) {
		error("row add performed on empty matrix");
	} else if (i1 === i2) {
		error("row add performed on same row");
	} else {
		const row2 = k_row(extract_row(mx, i2));
		return m_add_r(mx, i1, row2);
	}
}

function ref(mx) {
	// to get rid of rounding errors
	function round_zero(x) {
		// t is threshhold for accuracy
		const t = 0.0000001;
		if (x < 0) {
			const rem = -x % 1;
			const int = x + rem;
			return rem <= t ? int : (rem >= 1 - t ? int - 1 : x);
		} else if (x > 0) {
			const rem = x % 1;
			const int = x - rem;
			return rem <= t ? int : (rem >= 1 - t ? int + 1 : x);
		} else {
			return 0;
		}
	}

	function find_top_row(my, mz, i) {
		if (is_empty_matrix(my)) {
			// zero matrix
			return 1;
		} else if (is_empty_matrix(mz)) {
			const next_cols = col_left(my);
			return find_top_row(next_cols, next_cols, 1);
		} else if (head(first_row(mz)) !== 0) {
			return i;
		} else {
			return find_top_row(my, row_down(mz), i + 1);
		}
	}

	function swap_to_top(ma) {
		const t = find_top_row(ma, ma, 1);
		return t === 1 ? ma : matrix_r_swap(ma, 1, t);
	}

	function row_value(mb) {
		if (is_empty_matrix(mb)) {
			// zero matrix
			return 0;
		} else if (head(first_row(mb)) === 0) {
			return row_value(col_left(mb));
		} else {
			return tail(col_to_row(mb));
		}
	}

	function first_non_zero(xs) {
		if (is_empty_list(xs)) {
			// zero matrix
			return 1;
		} else if (head(xs) !== 0) {
			return head(xs);
		} else {
			return first_non_zero(tail(xs));
		}
	}

	function zero_down(mc, ys, row1) {
		if (is_empty_matrix(mc) || ys === 0) {
			// zero matrix
			return mc;
		} else if (head(ys) === 0) {
			const others = row_down(mc);
			return pair(first_row(mc),
				zero_down(others, tail(ys), row1));
		} else {
			const othr = row_down(mc);
			const a1 = first_non_zero(row1);
			const ai = head(ys);
			const k_x = x => -(ai / a1) * x;
			const k_row = row => map(k_x, row);
			const new_row = row_sum(first_row(mc), k_row(row1));
			// may return values affected by round-off error
			// return pair(new_row, zero_down(othr, tail(ys), row1));
			const r_map = xs => map(round_zero, xs);
			return pair(r_map(new_row),
				zero_down(othr, tail(ys), row1));
		}
	}

	if (is_empty_matrix(mx)) {
		return mx;
	} else {
		const m_sorted = swap_to_top(mx);
		// we fix round off error from zero_down:
		const top = map(round_zero, first_row(mx));
		const add1 = pair(top, zero_down(row_down(m_sorted),
			row_value(m_sorted), top));
		return pair(first_row(add1), ref(row_down(add1)));
	}
}

function rref(mx) {
	// to get rid of rounding errors
	function round_zero(x) {
		// t is threshhold for accuracy
		const t = 0.0000001;
		if (x < 0) {
			const rem = -x % 1;
			const int = x + rem;
			return rem <= t ? int : (rem >= 1 - t ? int - 1 : x);
		} else if (x > 0) {
			const rem = x % 1;
			const int = x - rem;
			return rem <= t ? int : (rem >= 1 - t ? int + 1 : x);
		} else {
			return 0;
		}
	}

	function first_non_zero(row) {
		if (is_empty_list(row)) {
			// zero row
			return 1;
		} else if (head(row) !== 0) {
			return head(row);
		} else {
			return first_non_zero(tail(row));
		}
	}

	function one_down(mr) {
		if (is_empty_matrix(mr)) {
			// zero matrix
			return mr;
		} else {
			const row1 = first_row(mr);
			const entry = first_non_zero(row1);
			const others = row_down(mr);
			return entry === 1
				? pair(row1, one_down(others))
				: pair(map(x => x / entry, row1), one_down(others));
		}
	}

	function zero_down(mc, ys, row1) {
		if (is_empty_matrix(mc) || ys === 0) {
			return mc;
		} else if (head(ys) === 0) {
			const others = row_down(mc);
			return pair(first_row(mc),
				zero_down(others, tail(ys), row1));
		} else {
			const othr = row_down(mc);
			const int = head(ys);
			const k_x = x => -int * x;
			const k_row = row => map(k_x, row);
			const new_row = row_sum(first_row(mc), k_row(row1));
			// may return values affected by round-off error
			// return pair(new_row, zero_down(othr, tail(ys), row1));
			const r_map = xs => map(round_zero, xs);
			return pair(r_map(new_row),
				zero_down(othr, tail(ys), row1));
		}
	}

	function row_value(my, xs) {
		if (is_empty_matrix(my)) {
			return xs;
		} else if (head(first_row(my)) === 0) {
			return row_value(col_left(my), xs);
		} else {
			return row_value(col_left(my), col_to_row(my));
		}
	}

	function down_adder(ms) {
		if (is_empty_matrix(ms)) {
			return ms;
		} else {
			const ys = row_value(ms, 0);
			// we fix round off error from zero_down:
			const top = map(round_zero, first_row(ms));
			const add1 = ys === 0 ? ms // zero matrix
				: pair(top, zero_down(row_down(ms), tail(ys), top));
			return pair(first_row(add1), down_adder(row_down(add1)));
		}
	}

	if (is_empty_matrix(mx)) {
		return mx;
	} else {
		const all_ones = one_down(ref(mx));
		display("ref found, calculating rref...");
		const mirror = my => map(row => reverse(row), my);
		const rot = reverse(mirror(all_ones));
		return reverse(mirror(down_adder(rot)));
	}
}

function del_n_rows(mx, n) {
	if (is_empty_matrix(mx) || n <= 0) {
		return mx;
	} else {
		return del_n_rows(row_down(mx), n - 1);
	}
}

function del_n_cols(mx, n) {
	if (is_empty_matrix(mx)) {
		return empty_mat();
	} else if (n <= 0) {
		return mx;
	} else {
		return del_n_cols(col_left(mx), n - 1);
	}
}

function inverse(mx) {
	function is_iden(my) {
		if (is_empty_matrix(my)) {
			return true;
		} else {
			const a11 = head(first_row(my));
			const next = col_left(row_down(my));
			return a11 === 1 && is_iden(next);
		}
	}

	function print_error(mr) {
		display("rref of (I | mx):");
		draw_matrix(mr);
		error("inverse failed: matrix is singular");
	}

	if (is_empty_matrix(mx)) {
		error("finding inverse of empty matrix");
	} else if (length(mx) !== width(mx)) {
		error("finding inverse of non-square matrix");
	} else {
		const dim = length(mx);
		const iden = make_identity(dim);
		const i_a = augment_matrix(mx, iden);
		display("calculating inverse using rref...");
		const gje = rref(i_a);
		return is_iden(gje) ? del_n_cols(gje, dim)
			: print_error(gje);
	}
}

function submatrix(mx, i, j) {
	if (is_empty_matrix(mx)) {
		error("finding submatrix of empty matrix");
	} else {
		return del_col(del_row(mx, i), j);
	}
}

function cofactor(mx, i, j) {
	const mij = determinant(submatrix(mx, i, j));
	return math_pow(-1, i + j) * mij;
}

function determinant(mx) {
	const len = length(mx);
	const wid = width(mx);

	function print_size() {
		display("n = " + len + " - calculating determinant...");
	}
	const dis_limit = 5;
	const dis = (len >= dis_limit ? print_size() : 0);

	function accu_j(f, row, j) {
		if (j > wid) {
			return 0;
		} else {
			const a = head(row);
			return f(a, j) + accu_j(f, tail(row), j + 1);
		}
	}

	if (is_empty_matrix(mx)) {
		error("finding determinant of empty matrix");
	} else if (len !== wid) {
		error("finding determinant of non-square matrix");
	} else if (len === 1) {
		return head(first_row(mx));
	} else if (len === 2) {
		const row1 = first_row(mx);
		const row2 = first_row(row_down(mx));
		const a = head(row1);
		const b = head(tail(row1));
		const c = head(row2);
		const d = head(tail(row2));
		return a * d - b * c;
	} else {
		const exp_row1 = (x, j) => x * cofactor(mx, 1, j);
		return accu_j(exp_row1, first_row(mx), 1);
	}
}

function cof_matrix(mx) {
	const len = length(mx);
	const wid = width(mx);

	function print_size() {
		display("n = " + len +
			" - calculating cofactor matrix...\n");
	}
	const dis_limit = 4;

	if (is_empty_matrix(mx)) {
		error("finding cofactors of empty matrix");
	} else if (len !== wid) {
		error("finding cofactors of non-square matrix");
	} else {
		const func = (i, j) => cofactor(mx, i, j);
		const dis = (len >= dis_limit ? print_size() : 0);
		return build_matrix(len, wid, func);
	}
}

function matrix_adj(mx) {
	return transpose(cof_matrix(mx));
}

// Test
// const m1 = list(list(0, 0, 2, 4, 2), list(1, 2, 4, 5, 3),
//	list(-2, -4, -5, -4, 3));
// const r2 = 5;
// const c2 = 5;
// const xs2 = list(1, 0, 3, 1, 8, 2, 1, 1, 2, -5, 4, 1, 7, 5, 10,
//	0, -1, 5, 0, 22, 0, 0, 1, 0, 0);
// const m2 = matrix_from_list(r2, c2, xs2);
// const xs3 = list(1, -9, 4, 6, -4, -3, -14, 2, -1, -7, -3, 5,
//	-1, 7, -4, -5);
// const m3 = matrix_from_list(4, 4, xs3);

// draw_matrix(m2);
// draw_matrix(m3);

// draw_matrix(inverse(m2));
// draw_matrix(matrix_adj(m3));

function arr_to_mx(r, c, arr) {
	const len = array_length(arr);
	let mx = empty_mat();
	if (len !== r * c) {
		error("matrix from array failed: wrong dimensions");
	} else {
		let n = len - 1;
		for (let k = r - 1; k >= 0; k = k - 1) {
			let row = [];
			for (let l = c - 1; l >= 0; l = l - 1) {
				row = pair(arr[n], row);
				n = n - 1;
			}
			mx = pair(row, mx);	
		}
	}
	draw_matrix(mx);
	return mx;
}

function row_vector(arr) {
	const len = array_length(arr);
	let vec = [];
	for (let i = len - 1; i >= 0; i = i - 1) {
		vec = pair(arr[i], vec);
	}
	const mx = pair(vec, []);
	draw_matrix(mx);
	return mx;
}

function col_vector(arr) {
	const len = array_length(arr);
	let mx = [];
	for (let i = len - 1; i >= 0; i = i - 1) {
		const row = pair(arr[i], []);
		mx = pair(row, mx);
	}
	draw_matrix(mx);
	return mx;
}

function least_squares(mx) {
	const rows = length(mx);
	const cols = length(first_row(mx));
	let rep = [];
	let b = [];
	let j = cols;

	function extract_col(my) {
		if (j <= 1) {
			b = transpose(pair(col_to_row(my), []));
			return pair(head(build_row(rows, x => 1)), []);
		} else {
			j = j - 1;
			return pair(col_to_row(my), extract_col(col_left(my)));
		}
	}

	const AT = extract_col(mx);
	const A = transpose(AT);
	const ATA = matrix_pdt(AT, A);
	const ATb = matrix_pdt(AT, b);
	const result = del_n_cols(rref(augment_matrix(ATA, ATb)), cols);
	draw_matrix(result);
	return result;
}

function eigen_3x3(mx) {
	function tr(mz) {
		const a11 = head(first_row(mz));
		const down = row_down(mz);
		const a22 = head(tail(first_row(down)));
		const a33 = head(tail(tail(first_row(row_down(down)))));
		return a11 + a22 + a33;
	}

	const rows = length(mx);
	const cols = length(first_row(mx));
	let result = [];

	if (rows === 2 && cols === 2) {
		const entry22 = head(tail(first_row(row_down(mx))));
		const tr2 = head(first_row(mx)) + entry22;
		const det = determinant(mx);
		const disc = tr2 * tr2 - 4 * det;
		if (disc < 0) {
			error("no real eigenvalues found");
		} else { }
		const rt1 = math_sqrt(disc);
		const rt2 = -rt1;
		const eig1 = (tr2 + rt1) / 2;
		const eig2 = (tr2 + rt2) / 2;
		result = list(list(eig1, eig2));
	} else if (rows === 3 && cols === 3) {
		const q = tr(mx) / 3;
		const pmy = matrix_sum(mx, scalar_pdt(-q, make_identity(3)));
		const p1 = tr(matrix_pdt(pmy, pmy));
		const p = math_sqrt(p1 / 6);
		const my = scalar_pdt(1 / p, pmy);
		const det = determinant(my);
		const phi = k =>
			(math_acos(det / 2) / 3) + (2 * k * math_PI)/ 3;
		const eig_f = k => p * (2 * math_cos(phi(k))) + q;
		const eig1 = eig_f(0);
		const eig2 = eig_f(1);
		const eig3 = eig_f(2);
		result = list(list(eig1, eig2, eig3));
	} else if (rows !== cols) {
		error("eigen_3x3 failed: not square matrix");
	} else {
		error("eigen_3x3 failed: not order 2 or 3");
	}
	draw_matrix(result);
	return result;
}

function eig3_arr(arr) {
	return eigen_3x3(arr_to_mx(3, 3, arr));
}