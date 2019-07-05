// 1. "Hidden" functions
function display(x) {
	console.log(x);
}

function runtime() {
//	Source Academy version:
//	return Date.now();
	return performance.now();
}

function error(x) {
	throw new Error(x);
}

function is_number(x) {
	return typeof x === "number";
}

function is_boolean(x) {
	return typeof x === "boolean";
}

function is_string(x) {
	return typeof x === "string";
}

function is_function(x) {
	return typeof x === "function";
}

function is_object(x) {
	return typeof x === "object";
}

// 2. New functions
function disp_data(xs) {
	function arr_to_str(ys) {
		let str = "";
		if (is_array(ys)) {
			str = array_length(ys) === 0 ? "" : arr_to_str(ys[0]);
			for (let i = 1; i < ys.length; i = i + 1) {
				const next = arr_to_str(ys[i]);
				str = str + ", " + next;
			}
			str = "[" + str + "]";
		} else {
			str = str + ys;
		}
		return str;
	}
	const string = arr_to_str(xs);
	display(string);
	return xs;
}

// 3. Implementation available on Source Academy
function is_array(x) {
	if (Array.isArray === undefined) {
		return x instanceof Array;
	} else {
		return Array.isArray(x);
	}
}

function pair(a, b) {
	return [a, b];
}

function is_pair(x) {
	return is_array(x) && x.length === 2;
}

function head(xs) {
	if (is_pair(xs)) {
		return xs[0];
	} else {
		error("head encountered " + xs);
	}
}

function tail(xs) {
	if (is_pair(xs)) {
		return xs[1];
	} else {
		error("tail encountered " + xs);
	}
}

function is_empty_list(xs) {
	if (is_array(xs) && xs.length === 0) {
		return true;
	} else {
		return false;
	}
}

function is_list(xs) {
	return is_empty_list(xs) || (is_pair(xs) && is_list(tail(xs)));
}

function list() {
	let xs = [];
	for (let i = arguments.length - 1; i >= 0; i = i - 1) {
		xs = pair(arguments[i], xs);
	}
	return xs;
}

function draw_list(xs) {
	if (is_list(xs)) {
		disp_data(xs);
	} else {
		error("draw_list encountered " + xs);
	}
	return xs;
}

function equal(item1, item2) {
	return (is_pair(item1) && is_pair(item2))
		? (equal(head(item1), head(item2)) && equal(tail(item1), tail(item2)))
		: (is_empty_list(item1) && is_empty_list(item2)) || item1 === item2;
}

function length(xs) {
	return is_empty_list(xs) ? 0 : 1 + length(tail(xs));
}

function map(f, xs) {
	return is_empty_list(xs) ? []
		: pair(f(head(xs)), map(f, tail(xs)));
}

function build_list(n, fun) {
	function build(i, fun, already_built) {
		return i < 0 ? already_built
			: build(i - 1, fun, pair(fun(i), already_built));
	}
	return build(n - 1, fun, []);
}

function for_each(fun, xs) {
	if (is_empty_list(xs)) {
		return true;
	} else {
		fun(head(xs)); return for_each(fun, tail(xs));
	}
}

function to_string(x) {
	return x + "";
}

function list_to_string(xs) {
	return is_empty_list(xs) ? "[]" : (is_pair(xs)
			? "[" + list_to_string(head(xs)) + "," + list_to_string(tail(xs)) + "]"
			: to_string(xs));
}

function reverse(xs) {
	function rev(original, reversed) {
		return is_empty_list(original) ? reversed
			: rev(tail(original), pair(head(original), reversed));
	}
	return rev(xs, []);
}

function append(xs, ys) {
	return is_empty_list(xs) ? ys
		: pair(head(xs), append(tail(xs), ys));
}


function member(v, xs) {
	return is_empty_list(xs) ? []
		: ((v === head(xs)) ? xs : member(v, tail(xs)));
}

function remove(v, xs) {
	return is_empty_list(xs) ? [] :
		(v === head(xs) ? tail(xs) : pair(head(xs), remove(v, tail(xs))));
}

function remove_all(v, xs) {
	return is_empty_list(xs) ? []
		: (v === head(xs) ? remove_all(v, tail(xs))
			: pair(head(xs), remove_all(v, tail(xs))));
}

function filter(pred, xs) {
	return is_empty_list(xs) ? xs
		: (pred(head(xs)) ? pair(head(xs), filter(pred, tail(xs)))
			: filter(pred, tail(xs)));
}

function enum_list(start, end) {
	return start > end ? [] : pair(start, enum_list(start + 1, end));
}

function list_ref(xs, n) {
	return n === 0 ? head(xs) : list_ref(tail(xs), n - 1);
}

function accumulate(f, initial, xs) {
	return is_empty_list(xs) ? initial
		: f(head(xs), accumulate(f, initial, tail(xs)));
}

function set_head(xs,x) {
	if (is_pair(xs)) {
		xs[0] = x;
	} else {
		error("set_head encountered " + xs);
	}
}

function set_tail(xs,x) {
	if (is_pair(xs)) {
		xs[1] = x;
	} else {
		error("set_tail encountered " + xs);
	}
}

function array_length(x) {
	return x.length;
}