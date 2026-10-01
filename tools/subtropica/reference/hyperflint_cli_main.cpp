// HyperFLINT CLI.
//
// Single JSON-in / JSON-out entry point:
//   hyperflint eval-json  < request.json  > response.json
//
// Request schema (flat, one request per line):
//   {"op": <opname>, <op-specific fields...>, "vars": [...]}
//
// Supported ops (Phase 0 and 1a):
//   factor   : {"op":"factor", "expr":<str>}              -> {constant, factors[]}
//   add/sub/mul : {"op":"<op>", "a":<str>, "b":<str>}     -> {result}
//   neg      : {"op":"neg", "a":<str>}                    -> {result}
//   pow      : {"op":"pow", "a":<str>, "n":<int>}         -> {result}
//
// If "vars" is omitted, the expression(s) are scanned for every
// identifier and those become the polynomial ring variables.

#include "hyperflint/core/poly.hpp"
#include "hyperflint/core/rat.hpp"
#include "hyperflint/algebra/linear_factors.hpp"
#include "hyperflint/core/zw_table.hpp"  // iter-52 C0c.1: ZWTable for linear_factors transient
#include "hyperflint/algebra/algebraic_letters.hpp"
#include "hyperflint/convert/expr.hpp"
#include "hyperflint/convert/parse.hpp"
#include "hyperflint/convert/convert_hlog.hpp"
#include "hyperflint/algebra/partial_fractions.hpp"
#include "hyperflint/algebra/shuffle.hpp"
#include "hyperflint/algebra/convert.hpp"
#include "hyperflint/algebra/diff.hpp"
#include "hyperflint/series/expansions.hpp"
#include "hyperflint/series/mpl_sum.hpp"
#include "hyperflint/series/hlog_series.hpp"
#include "hyperflint/series/mpl_series.hpp"
#include "hyperflint/integrator/regularize.hpp"
#include "hyperflint/integrator/differentiate.hpp"
#include "hyperflint/integrator/transform.hpp"
#include "hyperflint/integrator/primitive.hpp"
#include "hyperflint/integrator/integration_step.hpp"
#include "hyperflint/integrator/hyper_int.hpp"
#include "hyperflint/integrator/lr_search.hpp"
#include "hyperflint/bridge/env_flags.hpp"  // iter-96 §T7 26th chunk Track-diagnostic-dump partial bridge portion
#include "hyperflint/bridge/handlers.hpp"
#include <cmath>
#include "hyperflint/series/laurent.hpp"
#include "hyperflint/reduce/mzv_reduce.hpp"
#include "hyperflint/reduce/periods.hpp"
#include "hyperflint/reduce/break_up_contour.hpp"
#include "hyperflint/runtime/scalar_rep.hpp"     // C0b.4 iter-42: scalar_rep_enabled() at break_up_contour_sym callsite
#include "hyperflint/instrumentation/dag_hashcons_probe.hpp"  // HF FF Phase 5 §A.1 iter-49: probe A1 init (env-gated; default-OFF no-op)
#include "hyperflint/instrumentation/mpz_pool_probe.hpp"      // HF FF Phase 5 REC-1 iter-83: mpz-pool tracker init (env-gated; default-OFF no-op)
#include "hyperflint/symbols/word.hpp"
#include "hyperflint/core/symcoef.hpp"

#include <cctype>
#include <chrono>
#include <climits>
#include <cstdlib>      // std::abort (iter-44 require_persistent assertion)
#include <iostream>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace hf = hyperflint;

namespace {

// ---------- minimal JSON helpers (flat schema) ----------

std::string json_escape(const std::string& s) {
    std::string out;
    out.reserve(s.size());
    for (char c : s) {
        switch (c) {
            case '"':  out += "\\\""; break;
            case '\\': out += "\\\\"; break;
            case '\n': out += "\\n";  break;
            case '\r': out += "\\r";  break;
            case '\t': out += "\\t";  break;
            default:   out += c;
        }
    }
    return out;
}

std::string json_str_field(const std::string& body, const std::string& key) {
    std::regex re("\"" + key + "\"\\s*:\\s*\"((?:[^\"\\\\]|\\\\.)*)\"");
    std::smatch m;
    if (!std::regex_search(body, m, re)) return {};
    std::string v = m[1];
    std::string out;
    out.reserve(v.size());
    for (size_t i = 0; i < v.size(); ++i) {
        if (v[i] == '\\' && i + 1 < v.size()) {
            char nx = v[++i];
            switch (nx) {
                case 'n': out += '\n'; break;
                case 't': out += '\t'; break;
                case 'r': out += '\r'; break;
                default:  out += nx;
            }
        } else {
            out += v[i];
        }
    }
    return out;
}

std::vector<std::string> json_str_array(const std::string& body,
                                        const std::string& key) {
    std::regex re("\"" + key + "\"\\s*:\\s*\\[([^\\]]*)\\]");
    std::smatch m;
    std::vector<std::string> out;
    if (!std::regex_search(body, m, re)) return out;
    std::string inner = m[1];
    std::regex re_item("\"((?:[^\"\\\\]|\\\\.)*)\"");
    for (auto it = std::sregex_iterator(inner.begin(), inner.end(), re_item);
         it != std::sregex_iterator(); ++it) {
        out.push_back((*it)[1]);
    }
    return out;
}

long json_int_field(const std::string& body, const std::string& key) {
    std::regex re("\"" + key + "\"\\s*:\\s*(-?\\d+)");
    std::smatch m;
    if (!std::regex_search(body, m, re)) return 0;
    return std::stol(m[1]);
}

// ---------- Variable auto-scan ----------

std::vector<std::string> autoscan_vars(
    std::initializer_list<const std::string*> exprs) {
    std::set<std::string> seen;
    std::regex re("[A-Za-z][A-Za-z0-9_]*");
    for (const std::string* ep : exprs) {
        for (auto it = std::sregex_iterator(ep->begin(), ep->end(), re);
             it != std::sregex_iterator(); ++it) {
            seen.insert((*it)[0]);
        }
    }
    return std::vector<std::string>(seen.begin(), seen.end());
}

// ---------- Op handlers ----------

std::string emit_factor_json(const std::vector<std::string>& vars,
                             const hf::Factored& fac) {
    std::ostringstream o;
    o << "{\"op\":\"factor\""
      << ",\"constant\":\"" << json_escape(fac.constant) << "\""
      << ",\"factors\":[";
    for (size_t i = 0; i < fac.factors.size(); ++i) {
        if (i) o << ",";
        o << "[\"" << json_escape(fac.factors[i].first) << "\","
          << fac.factors[i].second << "]";
    }
    o << "]";
    o << ",\"vars\":[";
    for (size_t i = 0; i < vars.size(); ++i) {
        if (i) o << ",";
        o << "\"" << json_escape(vars[i]) << "\"";
    }
    o << "]}";
    return o.str();
}

// Forward declaration — defined later in the file, used by
// handle_find_lr_orders (which sits above its definition).
static std::string extract_top_array(const std::string& body,
                                     const std::string& key);

std::string emit_result_json(const std::string& op,
                             const std::vector<std::string>& vars,
                             const std::string& result) {
    std::ostringstream o;
    o << "{\"op\":\"" << op << "\""
      << ",\"result\":\"" << json_escape(result) << "\""
      << ",\"vars\":[";
    for (size_t i = 0; i < vars.size(); ++i) {
        if (i) o << ",";
        o << "\"" << json_escape(vars[i]) << "\"";
    }
    o << "]}";
    return o.str();
}

int handle_factor(const std::string& body) {
    std::string expr = json_str_field(body, "expr");
    if (expr.empty()) { std::cerr << "factor: missing \"expr\"\n"; return 1; }
    auto vars = json_str_array(body, "vars");
    if (vars.empty()) vars = autoscan_vars({&expr});
    if (vars.empty()) vars.push_back("x");

    hf::PolyCtx ctx(vars);
    hf::Poly p(ctx, expr);
    auto fac = hf::factor(p);
    std::cout << emit_factor_json(vars, fac) << "\n";
    return 0;
}

// Phase β.2: LR-order search op.  Port of STFasterFubini2
// (SubTropica.wl:11050-11149).  Heuristic = LeafCountLinear,
// FindRoots = False.  Multi-group supported.
//
// Request (two equivalent input shapes; `groups` takes priority if
// both are given):
//   {"op":"find_lr_orders",
//    "polys":["<p1_str>", ..., "<pn_str>"],     (single-group: flat list)
//    OR
//    "groups":[["<p1>", ...], ["<q1>", ...]],   (multi-group: list of lists)
//    "xvars":["x1", ..., "xk"],                 (integration vars, in
//                                                user order; ctx uses
//                                                a stable mapping)
//    "coeff_vars":["MM", "s12", ...]            (optional symbolic
//                                                kinematic parameters)}
// Response:
//   {"op":"find_lr_orders",
//    "best_order":["x5", "x7", ...],            (empty if NOLR)
//    "score": 2996.9535...,                     (null if NOLR)
//    "nolr": bool,
//    "timing_compute_s": N.NNN,
//    "nXVars": K, "nGroups": G, "nPolys": [N_1, ..., N_G]}
int handle_find_lr_orders(const std::string& body) {
    // Track 6.5 wire-in landed at iter-40 (satisfies iter-35 rec 3 BINDING
    // before the iter-45 deadline): the real @ARCH:step-strategy marker
    // pair now wraps the pick_step_strategy() call in
    // src/bridge/handlers.cpp::find_lr_orders, since this CLI shim just
    // forwards the body to the handler and prints its JSON response.  The
    // response JSON now carries a "strategy" field whose value is one of
    // the four StepStrategy enum names; Mma callers dispatch on it.
    // See HyperFLINT/docs/ARCHITECTURE.md §vi for the single source of
    // truth and HyperFLINT development notes (iter40 track 6.5 close audit)
    // for the wire-in audit.
    std::cout << hyperflint::handlers::find_lr_orders(body) << "\n";
    return 0;
}

// Factor-prediction table (spec 2026-06-11): transport-neutral
// handlers.cpp::factor_table does all the work; this shim just wraps
// stdout + newline emission, mirroring handle_find_lr_orders.
int handle_factor_table(const std::string& body) {
    std::cout << hyperflint::handlers::factor_table(body) << "\n";
    return 0;
}

// Doppio-port phase 3 bridge op (2026-06-06): CLI shim for the
// projective Cheng-Wu gauge scan.  Pure forward to the transport-
// neutral handler (src/bridge/handlers.cpp::find_lr_orders_scan);
// request/response schema documented there and in handlers.hpp.
int handle_find_lr_orders_scan(const std::string& body) {
    std::cout << hyperflint::handlers::find_lr_orders_scan(body) << "\n";
    return 0;
}

// Phase β.1: discriminant of a polynomial with respect to a named
// variable.  Matches Mathematica's `Discriminant[p, x]` up to an
// integer sign `(-1)^{n(n-1)/2}` that doesn't matter for LR-search
// (dedup canonicalizes away the sign + numeric factor).
//
// Request:
//   {"op":"discriminant", "expr":"<poly>", "var":"<name>",
//    "vars":[...all symbols the poly contains...]}
// Response:
//   {"op":"discriminant", "var":"<name>", "result":"<poly>",
//    "vars":[...]}
int handle_discriminant(const std::string& body) {
    std::string expr = json_str_field(body, "expr");
    std::string var  = json_str_field(body, "var");
    if (expr.empty()) { std::cerr << "discriminant: missing \"expr\"\n"; return 1; }
    if (var.empty())  { std::cerr << "discriminant: missing \"var\"\n";  return 1; }
    auto vars = json_str_array(body, "vars");
    if (vars.empty()) vars = autoscan_vars({&expr});
    // Ensure `var` is in the ctx.
    {
        bool present = false;
        for (const auto& v : vars) if (v == var) { present = true; break; }
        if (!present) vars.push_back(var);
    }

    hf::PolyCtx ctx(vars);
    size_t var_idx = ctx.index_of(var);
    if (var_idx == SIZE_MAX) {
        std::cerr << "discriminant: var not in ctx: " << var << "\n";
        return 1;
    }

    hf::Poly p(ctx, expr);
    hf::Poly d = p.discriminant_in_var(var_idx);

    std::cout << emit_result_json("discriminant", vars, d.to_string()) << "\n";
    return 0;
}

int handle_binary(const std::string& op, const std::string& body) {
    std::string a = json_str_field(body, "a");
    std::string b = json_str_field(body, "b");
    if (a.empty() || b.empty()) {
        std::cerr << op << ": need \"a\" and \"b\" fields\n";
        return 1;
    }
    auto vars = json_str_array(body, "vars");
    if (vars.empty()) vars = autoscan_vars({&a, &b});
    if (vars.empty()) vars.push_back("x");

    hf::PolyCtx ctx(vars);
    hf::Poly pa(ctx, a);
    hf::Poly pb(ctx, b);
    hf::Poly r(ctx);
    if      (op == "add") r = pa + pb;
    else if (op == "sub") r = pa - pb;
    else if (op == "mul") r = pa * pb;
    else { std::cerr << "unknown binary op: " << op << "\n"; return 1; }

    std::cout << emit_result_json(op, vars, r.to_string()) << "\n";
    return 0;
}

int handle_neg(const std::string& body) {
    std::string a = json_str_field(body, "a");
    if (a.empty()) { std::cerr << "neg: missing \"a\"\n"; return 1; }
    auto vars = json_str_array(body, "vars");
    if (vars.empty()) vars = autoscan_vars({&a});
    if (vars.empty()) vars.push_back("x");

    hf::PolyCtx ctx(vars);
    hf::Poly pa(ctx, a);
    hf::Poly r = -pa;
    std::cout << emit_result_json("neg", vars, r.to_string()) << "\n";
    return 0;
}

int handle_pow(const std::string& body) {
    std::string a = json_str_field(body, "a");
    long n = json_int_field(body, "n");
    if (a.empty() || n < 0) {
        std::cerr << "pow: need \"a\" and non-negative \"n\"\n";
        return 1;
    }
    auto vars = json_str_array(body, "vars");
    if (vars.empty()) vars = autoscan_vars({&a});
    if (vars.empty()) vars.push_back("x");

    hf::PolyCtx ctx(vars);
    hf::Poly pa(ctx, a);
    hf::Poly r = pa.pow(static_cast<unsigned long>(n));
    std::cout << emit_result_json("pow", vars, r.to_string()) << "\n";
    return 0;
}

// ---------- Phase 1b: calculus ----------

int handle_derivative(const std::string& body) {
    std::string a = json_str_field(body, "a");
    std::string var = json_str_field(body, "var");
    if (a.empty() || var.empty()) {
        std::cerr << "derivative: need \"a\" and \"var\"\n"; return 1;
    }
    auto vars = json_str_array(body, "vars");
    if (vars.empty()) vars = autoscan_vars({&a});
    if (vars.empty()) vars.push_back(var);
    // Ensure `var` is in `vars`.
    bool present = false;
    for (const auto& v : vars) if (v == var) { present = true; break; }
    if (!present) vars.push_back(var);

    hf::PolyCtx ctx(vars);
    hf::Poly pa(ctx, a);
    size_t idx = 0;
    for (; idx < vars.size(); ++idx) if (vars[idx] == var) break;
    hf::Poly r = pa.derivative(idx);
    std::cout << emit_result_json("derivative", vars, r.to_string()) << "\n";
    return 0;
}

int handle_eval(const std::string& body) {
    // Request: {"op":"eval","a":<expr>,"vars":[v1,...,vn],"values":[q1,...,qn]}
    // Returns: {"op":"eval","result":<rational string>,"vars":[...]}
    std::string a = json_str_field(body, "a");
    if (a.empty()) { std::cerr << "eval: missing \"a\"\n"; return 1; }
    auto vars = json_str_array(body, "vars");
    auto values = json_str_array(body, "values");
    if (vars.size() != values.size()) {
        std::cerr << "eval: vars and values must have the same length\n";
        return 1;
    }
    hf::PolyCtx ctx(vars);
    hf::Poly pa(ctx, a);
    std::string r = pa.evaluate_all(values);
    std::cout << emit_result_json("eval", vars, r) << "\n";
    return 0;
}

// ---------- Phase 1c: division + GCD ----------

int handle_divexact(const std::string& body) {
    std::string a = json_str_field(body, "a");
    std::string b = json_str_field(body, "b");
    if (a.empty() || b.empty()) {
        std::cerr << "divexact: need \"a\" and \"b\"\n"; return 1;
    }
    auto vars = json_str_array(body, "vars");
    if (vars.empty()) vars = autoscan_vars({&a, &b});
    if (vars.empty()) vars.push_back("x");

    hf::PolyCtx ctx(vars);
    hf::Poly pa(ctx, a), pb(ctx, b);
    hf::Poly r = pa.divexact(pb);
    std::cout << emit_result_json("divexact", vars, r.to_string()) << "\n";
    return 0;
}

int handle_gcd(const std::string& body) {
    std::string a = json_str_field(body, "a");
    std::string b = json_str_field(body, "b");
    if (a.empty() || b.empty()) {
        std::cerr << "gcd: need \"a\" and \"b\"\n"; return 1;
    }
    auto vars = json_str_array(body, "vars");
    if (vars.empty()) vars = autoscan_vars({&a, &b});
    if (vars.empty()) vars.push_back("x");

    hf::PolyCtx ctx(vars);
    hf::Poly pa(ctx, a), pb(ctx, b);
    hf::Poly r = pa.gcd(pb);
    std::cout << emit_result_json("gcd", vars, r.to_string()) << "\n";
    return 0;
}

// ---------- Phase 1d: Rat ops ----------

int handle_rat_binary(const std::string& op, const std::string& body) {
    // op in {rat_add, rat_sub, rat_mul, rat_div}
    std::string a = json_str_field(body, "a");
    std::string b = json_str_field(body, "b");
    if (a.empty() || b.empty()) {
        std::cerr << op << ": need \"a\" and \"b\"\n"; return 1;
    }
    auto vars = json_str_array(body, "vars");
    if (vars.empty()) vars = autoscan_vars({&a, &b});
    if (vars.empty()) vars.push_back("x");

    hf::PolyCtx ctx(vars);
    hf::Rat ra = hf::Rat::parse(ctx, a);
    hf::Rat rb = hf::Rat::parse(ctx, b);
    hf::Rat rr(hf::Poly(ctx, "0"));
    if      (op == "rat_add") rr = ra + rb;
    else if (op == "rat_sub") rr = ra - rb;
    else if (op == "rat_mul") rr = ra * rb;
    else if (op == "rat_div") rr = ra / rb;
    else { std::cerr << "unknown rat op: " << op << "\n"; return 1; }
    std::cout << emit_result_json(op, vars, rr.to_string()) << "\n";
    return 0;
}

// ---------- Phase 2b: PoleDegree + RatResidue ----------

// Parse a rational-expression string as a Rat.  For now the input must
// be either a polynomial or a single "num/den" fraction (top-level `/`
// split, same convention as Rat::parse).
static hf::Rat parse_rat_expr(const hf::PolyCtx& ctx, const std::string& s) {
    return hf::Rat::parse(ctx, s);
}

int handle_pole_degree(const std::string& body) {
    // Request: {"op":"pole_degree", "f":<str>, "var":<str>, "vars":[...]}
    // Response: {"op":"pole_degree", "result":<int-or-"infinity">, "vars":[...]}
    std::string f   = json_str_field(body, "f");
    std::string var = json_str_field(body, "var");
    if (f.empty() || var.empty()) {
        std::cerr << "pole_degree: need \"f\" and \"var\"\n"; return 1;
    }
    auto vars = json_str_array(body, "vars");
    if (vars.empty()) vars = autoscan_vars({&f});
    bool present = false;
    for (const auto& v : vars) if (v == var) { present = true; break; }
    if (!present) vars.push_back(var);

    hf::PolyCtx ctx(vars);
    hf::Rat r = parse_rat_expr(ctx, f);
    size_t idx = 0;
    for (; idx < vars.size(); ++idx) if (vars[idx] == var) break;

    long d = r.pole_degree(idx);
    std::string val = (d == LONG_MAX) ? "\"infinity\"" : std::to_string(d);
    std::ostringstream o;
    o << "{\"op\":\"pole_degree\",\"result\":" << val
      << ",\"vars\":[";
    for (size_t i = 0; i < vars.size(); ++i) {
        if (i) o << ",";
        o << "\"" << json_escape(vars[i]) << "\"";
    }
    o << "]}\n";
    std::cout << o.str();
    return 0;
}

int handle_rat_residue(const std::string& body) {
    // Request: {"op":"rat_residue", "f":<str>, "var":<str>, "vars":[...]}
    // Response: {"op":"rat_residue", "result":<canonical rat string>, "vars":[...]}
    std::string f   = json_str_field(body, "f");
    std::string var = json_str_field(body, "var");
    if (f.empty() || var.empty()) {
        std::cerr << "rat_residue: need \"f\" and \"var\"\n"; return 1;
    }
    auto vars = json_str_array(body, "vars");
    if (vars.empty()) vars = autoscan_vars({&f});
    bool present = false;
    for (const auto& v : vars) if (v == var) { present = true; break; }
    if (!present) vars.push_back(var);

    hf::PolyCtx ctx(vars);
    hf::Rat r = parse_rat_expr(ctx, f);
    size_t idx = 0;
    for (; idx < vars.size(); ++idx) if (vars[idx] == var) break;

    hf::Rat residue = r.rat_residue(idx);
    std::cout << emit_result_json("rat_residue", vars, residue.to_string()) << "\n";
    return 0;
}

// ---------- Phase 3b: Shuffle algebra ----------

// Parse a letter-list (JSON array of strings) into a Word over `ctx`.
static hf::Word parse_word(const hf::PolyCtx& ctx,
                           const std::vector<std::string>& letters) {
    hf::Word w;
    w.letters.reserve(letters.size());
    for (const auto& s : letters) {
        w.letters.push_back(hf::Rat::parse(ctx, s));
    }
    return w;
}

// Serialize a Wordlist as JSON:
//   [ {"coef":"<rat>", "word":["<l1>",...]}, ... ]
static std::string emit_wordlist(const hf::Wordlist& wl) {
    std::ostringstream o;
    o << "[";
    for (size_t i = 0; i < wl.terms.size(); ++i) {
        if (i) o << ",";
        o << "{\"coef\":\"" << json_escape(wl.terms[i].coef.to_string()) << "\""
          << ",\"word\":[";
        for (size_t j = 0; j < wl.terms[i].word.size(); ++j) {
            if (j) o << ",";
            o << "\"" << json_escape(wl.terms[i].word[j].to_string()) << "\"";
        }
        o << "]}";
    }
    o << "]";
    return o.str();
}

int handle_shuffle_words(const std::string& body) {
    // Request: {"op":"shuffle_words","v":["l1",...],"w":["l1",...],"vars":[...]}
    auto v_letters = json_str_array(body, "v");
    auto w_letters = json_str_array(body, "w");
    auto vars      = json_str_array(body, "vars");
    if (vars.empty()) {
        // Pool all letters for autoscan
        std::vector<const std::string*> ptrs;
        for (auto& s : v_letters) ptrs.push_back(&s);
        for (auto& s : w_letters) ptrs.push_back(&s);
        std::set<std::string> seen;
        std::regex re("[A-Za-z][A-Za-z0-9_]*");
        for (auto* sp : ptrs) {
            for (auto it = std::sregex_iterator(sp->begin(), sp->end(), re);
                 it != std::sregex_iterator(); ++it) {
                seen.insert((*it)[0]);
            }
        }
        vars.assign(seen.begin(), seen.end());
        if (vars.empty()) vars.push_back("x");
    }
    hf::PolyCtx ctx(vars);
    hf::Word V = parse_word(ctx, v_letters);
    hf::Word W = parse_word(ctx, w_letters);
    hf::Wordlist out = hf::shuffle_words(V, W);

    std::ostringstream o;
    o << "{\"op\":\"shuffle_words\",\"result\":"
      << emit_wordlist(out) << ",\"vars\":[";
    for (size_t i = 0; i < vars.size(); ++i) {
        if (i) o << ",";
        o << "\"" << json_escape(vars[i]) << "\"";
    }
    o << "]}\n";
    std::cout << o.str();
    return 0;
}

// Extract the body of a top-level "<key>":[ ... ] array from a JSON
// string, tracking [] depth so nested arrays inside each element are
// handled correctly. Returns the inner (between the outer [ and ])
// text, or "" if not found. Does NOT honor quoted-string handling
// of [ / ] inside strings; flat-schema convention means our strings
// contain only the characters listed in json_escape (plus alnum) so
// this is safe.
static std::string extract_top_array(const std::string& body,
                                     const std::string& key) {
    // Search for the KEY with colon so that letters inside word-arrays
    // (which may spell identically) don't match as the key. e.g. key
    // "b" must not match the letter "b" in "word":["a","b"].
    std::regex re("\"" + key + "\"\\s*:\\s*\\[");
    std::smatch m;
    if (!std::regex_search(body, m, re)) return {};
    size_t open_bracket = m.position(0) + m.length(0) - 1;
    int depth = 0;
    size_t start = open_bracket + 1;
    for (size_t i = open_bracket; i < body.size(); ++i) {
        if (body[i] == '[') depth++;
        else if (body[i] == ']') {
            depth--;
            if (depth == 0) return body.substr(start, i - start);
        }
    }
    return {};
}

// Parse a wordlist JSON array: [ {"coef":"r","word":["l1",...]}, ... ]
static hf::Wordlist parse_wordlist(const hf::PolyCtx& ctx,
                                   const std::string& body,
                                   const std::string& key) {
    hf::Wordlist out;
    std::string inner = extract_top_array(body, key);
    if (inner.empty()) return out;
    // Split `inner` into top-level {...} objects (curly-brace-balanced).
    std::vector<std::string> terms;
    int depth = 0;
    size_t start = 0;
    for (size_t i = 0; i < inner.size(); ++i) {
        if (inner[i] == '{') {
            if (depth == 0) start = i;
            depth++;
        } else if (inner[i] == '}') {
            depth--;
            if (depth == 0) terms.push_back(inner.substr(start, i - start + 1));
        }
    }
    for (const auto& term : terms) {
        std::string coef_s = json_str_field(term, "coef");
        std::vector<std::string> letters = json_str_array(term, "word");
        // Iter-96 §T7 26th chunk Track-diagnostic-dump partial bridge
        // portion: VALUE-family macro relocation to bridge/env_flags.hpp
        // (second extension of that header; first extension since
        // iter-94 creation). The TRUTHY predicate semantics
        // `if (std::getenv(...))` are preserved verbatim. Scope-of-
        // effect: gates per-term shuffle_product input-parsing trace
        // emission only — see bridge/env_flags.hpp narrative for the
        // intentional narrow scope and the iter-96 review-note
        // Recommendation 3 caveat that future readers should not
        // assume cross-cutting "master verbose" semantics from the
        // generic name without explicit narrative caveat at any new
        // call site they add.
        if (HF_FLAG_DEBUG) {
            std::cerr << "[dbg] term='" << term << "'\n"
                      << "      coef='" << coef_s << "'"
                      << "  letters=[";
            for (auto& l : letters) std::cerr << "'" << l << "',";
            std::cerr << "]\n";
        }
        out.terms.push_back(
            hf::WordlistTerm{hf::Rat::parse(ctx, coef_s),
                              parse_word(ctx, letters)});
    }
    return out;
}

int handle_shuffle_product(const std::string& body) {
    // Request: {"op":"shuffle_product","a":[...], "b":[...], "vars":[...]}
    auto vars = json_str_array(body, "vars");
    if (vars.empty()) vars.push_back("x");
    hf::PolyCtx ctx(vars);
    hf::Wordlist A = parse_wordlist(ctx, body, "a");
    hf::Wordlist B = parse_wordlist(ctx, body, "b");
    hf::Wordlist out = hf::shuffle_product(A, B);
    std::ostringstream o;
    o << "{\"op\":\"shuffle_product\",\"result\":"
      << emit_wordlist(out) << ",\"vars\":[";
    for (size_t i = 0; i < vars.size(); ++i) {
        if (i) o << ",";
        o << "\"" << json_escape(vars[i]) << "\"";
    }
    o << "]}\n";
    std::cout << o.str();
    return 0;
}

int handle_concat_mul(const std::string& body) {
    auto vars = json_str_array(body, "vars");
    if (vars.empty()) vars.push_back("x");
    hf::PolyCtx ctx(vars);
    hf::Wordlist A = parse_wordlist(ctx, body, "a");
    hf::Wordlist B = parse_wordlist(ctx, body, "b");
    hf::Wordlist out = hf::concat_mul(A, B);
    std::ostringstream o;
    o << "{\"op\":\"concat_mul\",\"result\":"
      << emit_wordlist(out) << ",\"vars\":[";
    for (size_t i = 0; i < vars.size(); ++i) {
        if (i) o << ",";
        o << "\"" << json_escape(vars[i]) << "\"";
    }
    o << "]}\n";
    std::cout << o.str();
    return 0;
}

// ---------- Phase 3d: differentiation ----------

int handle_diff_hlog(const std::string& body) {
    // Request: {"op":"diff_hlog","z":<str>,"word":["l1",...],"var":<str>,"vars":[...]}
    // Response: {"op":"diff_hlog","result":[{"coef":"...","z":"...","word":[...]},...]}
    std::string z_s = json_str_field(body, "z");
    std::string var = json_str_field(body, "var");
    auto letters = json_str_array(body, "word");
    auto vars = json_str_array(body, "vars");
    if (vars.empty()) vars.push_back("x");
    bool present = false;
    for (auto& v : vars) if (v == var) { present = true; break; }
    if (!present) vars.push_back(var);

    hf::PolyCtx ctx(vars);
    hf::Rat z = hf::Rat::parse(ctx, z_s);
    hf::Word word = parse_word(ctx, letters);
    size_t idx = 0;
    for (; idx < vars.size(); ++idx) if (vars[idx] == var) break;
    auto dh = hf::diff_hlog(z, word, idx);

    // Build a sympy-friendly total string: each term becomes
    //   (coef) * H(z, w1, w2, ...)     (H is an inert function symbol)
    // Empty-word terms degrade to just `coef` (since Hlog[z,{}] = 1).
    std::ostringstream total;
    bool first = true;
    for (const auto& t : dh) {
        if (!first) total << " + ";
        first = false;
        total << "(" << t.coef.to_string() << ")";
        if (!t.hlog.word.empty()) {
            total << "*H(" << t.hlog.z.to_string();
            for (const auto& l : t.hlog.word.letters) {
                total << "," << l.to_string();
            }
            total << ")";
        }
    }
    std::string total_str = dh.empty() ? "0" : total.str();

    std::ostringstream o;
    o << "{\"op\":\"diff_hlog\",\"result\":[";
    for (size_t i = 0; i < dh.size(); ++i) {
        if (i) o << ",";
        o << "{\"coef\":\"" << json_escape(dh[i].coef.to_string()) << "\""
          << ",\"z\":\""    << json_escape(dh[i].hlog.z.to_string()) << "\""
          << ",\"word\":[";
        for (size_t j = 0; j < dh[i].hlog.word.size(); ++j) {
            if (j) o << ",";
            o << "\"" << json_escape(dh[i].hlog.word[j].to_string()) << "\"";
        }
        o << "]}";
    }
    o << "],\"total_str\":\"" << json_escape(total_str) << "\"";
    o << ",\"vars\":[";
    for (size_t i = 0; i < vars.size(); ++i) {
        if (i) o << ",";
        o << "\"" << json_escape(vars[i]) << "\"";
    }
    o << "]}\n";
    std::cout << o.str();
    return 0;
}

// Parse a JSON integer array: "ns":[1, 2, 3]
static std::vector<long> json_int_array(const std::string& body,
                                        const std::string& key) {
    std::vector<long> out;
    std::regex re("\"" + key + "\"\\s*:\\s*\\[([^\\]]*)\\]");
    std::smatch m;
    if (!std::regex_search(body, m, re)) return out;
    std::string inner = m[1];
    std::regex re_int("-?\\d+");
    for (auto it = std::sregex_iterator(inner.begin(), inner.end(), re_int);
         it != std::sregex_iterator(); ++it) {
        out.push_back(std::stol((*it)[0]));
    }
    return out;
}

int handle_diff_mpl(const std::string& body) {
    // Request: {"op":"diff_mpl","ns":[1,2,...],"zs":["..",...],"var":<str>,"vars":[...]}
    // Response: {"op":"diff_mpl","result":[{"coef":"...","ns":[...],"zs":[...]},...]}
    auto ns = json_int_array(body, "ns");
    auto zs_s = json_str_array(body, "zs");
    std::string var = json_str_field(body, "var");
    auto vars = json_str_array(body, "vars");
    if (vars.empty()) vars.push_back("x");
    bool present = false;
    for (auto& v : vars) if (v == var) { present = true; break; }
    if (!present) vars.push_back(var);

    hf::PolyCtx ctx(vars);
    std::vector<hf::Rat> zs;
    for (auto& s : zs_s) zs.push_back(hf::Rat::parse(ctx, s));
    size_t idx = 0;
    for (; idx < vars.size(); ++idx) if (vars[idx] == var) break;
    auto dm = hf::diff_mpl(ns, zs, idx);

    // Build a sympy-friendly total string: each term becomes
    //   (coef) * M_n1_n2_...(z1, z2, ...)
    // Empty-index terms degrade to `coef` (convention Mpl[{},{}] = 1).
    std::ostringstream total;
    bool first = true;
    for (const auto& t : dm) {
        if (!first) total << " + ";
        first = false;
        total << "(" << t.coef.to_string() << ")";
        if (!t.mpl.indices.empty()) {
            total << "*M";
            for (long n : t.mpl.indices) {
                total << "_" << (n < 0 ? "m" : "") <<
                    (n < 0 ? std::to_string(-n) : std::to_string(n));
            }
            total << "(";
            for (size_t j = 0; j < t.mpl.args.size(); ++j) {
                if (j) total << ",";
                total << t.mpl.args[j].to_string();
            }
            total << ")";
        }
    }
    std::string total_str = dm.empty() ? "0" : total.str();

    std::ostringstream o;
    o << "{\"op\":\"diff_mpl\",\"result\":[";
    for (size_t i = 0; i < dm.size(); ++i) {
        if (i) o << ",";
        o << "{\"coef\":\"" << json_escape(dm[i].coef.to_string()) << "\""
          << ",\"ns\":[";
        for (size_t j = 0; j < dm[i].mpl.indices.size(); ++j) {
            if (j) o << ",";
            o << dm[i].mpl.indices[j];
        }
        o << "],\"zs\":[";
        for (size_t j = 0; j < dm[i].mpl.args.size(); ++j) {
            if (j) o << ",";
            o << "\"" << json_escape(dm[i].mpl.args[j].to_string()) << "\"";
        }
        o << "]}";
    }
    o << "],\"total_str\":\"" << json_escape(total_str) << "\"";
    o << ",\"vars\":[";
    for (size_t i = 0; i < vars.size(); ++i) {
        if (i) o << ",";
        o << "\"" << json_escape(vars[i]) << "\"";
    }
    o << "]}\n";
    std::cout << o.str();
    return 0;
}

// ---------- Phase 3c: endpoint conversions ----------

int handle_convert_zero_one(const std::string& body) {
    auto vars = json_str_array(body, "vars");
    if (vars.empty()) vars.push_back("x");
    hf::PolyCtx ctx(vars);
    hf::Wordlist wl = parse_wordlist(ctx, body, "wl");
    hf::Wordlist out = hf::convert_zero_one(wl);
    std::ostringstream o;
    o << "{\"op\":\"convert_zero_one\",\"result\":"
      << emit_wordlist(out) << ",\"vars\":[";
    for (size_t i = 0; i < vars.size(); ++i) {
        if (i) o << ",";
        o << "\"" << json_escape(vars[i]) << "\"";
    }
    o << "]}\n";
    std::cout << o.str();
    return 0;
}

int handle_convert_ab_to_zero_inf(const std::string& body) {
    // Request: {"op":"convert_ab_to_zero_inf","wl":[...],"A":"<rat>","B":"<rat>","vars":[...]}
    auto vars = json_str_array(body, "vars");
    if (vars.empty()) vars.push_back("x");
    std::string A_s = json_str_field(body, "A");
    std::string B_s = json_str_field(body, "B");
    if (A_s.empty() || B_s.empty()) {
        std::cerr << "convert_ab_to_zero_inf: need \"A\" and \"B\"\n"; return 1;
    }
    hf::PolyCtx ctx(vars);
    hf::Wordlist wl = parse_wordlist(ctx, body, "wl");
    hf::Rat A = hf::Rat::parse(ctx, A_s);
    hf::Rat B = hf::Rat::parse(ctx, B_s);
    hf::Wordlist out = hf::convert_ab_to_zero_inf(ctx, wl, A, B);
    std::ostringstream o;
    o << "{\"op\":\"convert_ab_to_zero_inf\",\"result\":" << emit_wordlist(out)
      << ",\"vars\":[";
    for (size_t i = 0; i < vars.size(); ++i) {
        if (i) o << ",";
        o << "\"" << json_escape(vars[i]) << "\"";
    }
    o << "]}\n";
    std::cout << o.str();
    return 0;
}

int handle_convert_1inf_to_01(const std::string& body) {
    auto vars = json_str_array(body, "vars");
    if (vars.empty()) vars.push_back("x");
    hf::PolyCtx ctx(vars);
    hf::Wordlist wl = parse_wordlist(ctx, body, "wl");
    hf::Wordlist out = hf::convert_1inf_to_01(wl);
    std::ostringstream o;
    o << "{\"op\":\"convert_1inf_to_01\",\"result\":"
      << emit_wordlist(out) << ",\"vars\":[";
    for (size_t i = 0; i < vars.size(); ++i) {
        if (i) o << ",";
        o << "\"" << json_escape(vars[i]) << "\"";
    }
    o << "]}\n";
    std::cout << o.str();
    return 0;
}

int handle_collect_words(const std::string& body) {
    auto vars = json_str_array(body, "vars");
    if (vars.empty()) vars.push_back("x");
    hf::PolyCtx ctx(vars);
    hf::Wordlist A = parse_wordlist(ctx, body, "wl");
    hf::Wordlist out = hf::collect_words(A);
    std::ostringstream o;
    o << "{\"op\":\"collect_words\",\"result\":"
      << emit_wordlist(out) << ",\"vars\":[";
    for (size_t i = 0; i < vars.size(); ++i) {
        if (i) o << ",";
        o << "\"" << json_escape(vars[i]) << "\"";
    }
    o << "]}\n";
    std::cout << o.str();
    return 0;
}

// ---------- Phase 2c: PartialFractions ----------

int handle_partial_fractions(const std::string& body) {
    // Track 8.1b chunk-2b (iter-48): handler body factored into
    // hyperflint::handlers::partial_fractions per iter-44 plan §B
    // chunk-2 deferral.  CLI shim writes the (no-envelope) response
    // body to stdout with a trailing newline; the C ABI wrapper in
    // src/bridge/c_abi.cpp produces the envelope-stamped sibling.
    // Per chunk-2 disambiguation (iter-44 review note A2): CLI does NOT
    // add the schema_version/hf_version envelope so the legacy CLI
    // output stays byte-identical.  Errors are reported as JSON in
    // the handler return value (no std::cerr / non-zero exit code);
    // the byte-identical CLI snapshot ctest gates regressions.
    std::cout << hf::handlers::partial_fractions(body) << "\n";
    return 0;
}

// ---------- Phase 2a: LinearFactors ----------

int handle_linear_factors(const std::string& body) {
    // Track 8.1b chunk-3 (iter-50): handler body factored into
    // hyperflint::handlers::linear_factors per iter-44 plan §B chunk-3
    // deferral.  CLI shim writes the (no-envelope) response body to
    // stdout with a trailing newline; the C ABI wrapper in
    // src/bridge/c_abi.cpp produces the envelope-stamped sibling.
    // Per chunk-2 disambiguation (iter-44 review note A2): CLI does NOT
    // add the schema_version/hf_version envelope so the legacy CLI
    // output stays byte-identical.  Errors are reported as JSON in
    // the handler return value (no std::cerr / non-zero exit code);
    // the iter-51 chunk-3b byte-identical CLI snapshot ctest will gate
    // regressions.
    std::cout << hf::handlers::linear_factors(body) << "\n";
    return 0;
}

int handle_rat_sum(const std::string& body) {
    // Sum a list of rationals.  Each term is given as a pair of strings
    // [num, den]. The caller is responsible for pre-splitting --
    // HyperFLINT's current parser doesn't handle full expression
    // syntax (see PLAN.md: full-expression parsing arrives in Phase 3
    // with the Hlog word algebra).
    //
    // Request:
    //   {"op":"rat_sum",
    //    "terms":[["1","x"], ["1","y"], ["1","z"]],
    //    "vars":["x","y","z"]}
    // Response:
    //   {"op":"rat_sum", "result":"(y*z + x*z + x*y)/(x*y*z)", "vars":[...]}

    // Parse the "terms" array manually (flat JSON -- each term is a
    // 2-element string array).
    std::regex re("\"terms\"\\s*:\\s*\\[(.*)\\]");
    std::smatch m;
    std::vector<std::pair<std::string, std::string>> terms;
    if (std::regex_search(body, m, re)) {
        std::string inner = m[1];
        // Walk the comma-separated pair list.
        // Each pair looks like  ["num","den"] .
        std::regex re_pair("\\[\\s*\"((?:[^\"\\\\]|\\\\.)*)\"\\s*,"
                           "\\s*\"((?:[^\"\\\\]|\\\\.)*)\"\\s*\\]");
        for (auto it = std::sregex_iterator(inner.begin(), inner.end(), re_pair);
             it != std::sregex_iterator(); ++it) {
            terms.emplace_back((*it)[1], (*it)[2]);
        }
    }
    if (terms.empty()) {
        std::cerr << "rat_sum: need a non-empty \"terms\" array of [num,den] pairs\n";
        return 1;
    }
    auto vars = json_str_array(body, "vars");
    if (vars.empty()) {
        std::vector<const std::string*> ptrs;
        for (auto& t : terms) { ptrs.push_back(&t.first); ptrs.push_back(&t.second); }
        std::set<std::string> seen;
        std::regex re_id("[A-Za-z][A-Za-z0-9_]*");
        for (auto* sp : ptrs) {
            for (auto it = std::sregex_iterator(sp->begin(), sp->end(), re_id);
                 it != std::sregex_iterator(); ++it) {
                seen.insert((*it)[0]);
            }
        }
        vars.assign(seen.begin(), seen.end());
    }
    if (vars.empty()) vars.push_back("x");

    hf::PolyCtx ctx(vars);
    hf::Rat acc(hf::Poly(ctx, "0"));
    for (auto& t : terms) {
        acc = acc + hf::Rat(hf::Poly(ctx, t.first), hf::Poly(ctx, t.second));
    }
    std::cout << emit_result_json("rat_sum", vars, acc.to_string()) << "\n";
    return 0;
}

int handle_resultant(const std::string& body) {
    std::string a   = json_str_field(body, "a");
    std::string b   = json_str_field(body, "b");
    std::string var = json_str_field(body, "var");
    if (a.empty() || b.empty() || var.empty()) {
        std::cerr << "resultant: need \"a\", \"b\", \"var\"\n"; return 1;
    }
    auto vars = json_str_array(body, "vars");
    if (vars.empty()) vars = autoscan_vars({&a, &b});
    bool present = false;
    for (const auto& v : vars) if (v == var) { present = true; break; }
    if (!present) vars.push_back(var);

    hf::PolyCtx ctx(vars);
    hf::Poly pa(ctx, a), pb(ctx, b);
    size_t idx = 0;
    for (; idx < vars.size(); ++idx) if (vars[idx] == var) break;
    hf::Poly r = pa.resultant(pb, idx);
    std::cout << emit_result_json("resultant", vars, r.to_string()) << "\n";
    return 0;
}

int handle_subst(const std::string& body) {
    // Request: {"op":"subst","a":<expr>,"var":<str>,"value":<rational>,"vars":[...]}
    std::string a   = json_str_field(body, "a");
    std::string var = json_str_field(body, "var");
    std::string val = json_str_field(body, "value");
    if (a.empty() || var.empty() || val.empty()) {
        std::cerr << "subst: need \"a\", \"var\", \"value\"\n"; return 1;
    }
    auto vars = json_str_array(body, "vars");
    if (vars.empty()) vars = autoscan_vars({&a});
    bool present = false;
    for (const auto& v : vars) if (v == var) { present = true; break; }
    if (!present) vars.push_back(var);

    hf::PolyCtx ctx(vars);
    hf::Poly pa(ctx, a);
    size_t idx = 0;
    for (; idx < vars.size(); ++idx) if (vars[idx] == var) break;
    hf::Poly r = pa.substitute_one_rat(idx, val);
    std::cout << emit_result_json("subst", vars, r.to_string()) << "\n";
    return 0;
}

// ---------- Phase 4a: word-level series expansions ----------

// Emit a SeriesTable as nested JSON arrays of Rat strings:
//   [[ "c00","c01",... ], [ "c10","c11",... ], ... ]
static std::string emit_series_table(const hf::SeriesTable& t) {
    std::ostringstream o;
    o << "[";
    for (size_t i = 0; i < t.size(); ++i) {
        if (i) o << ",";
        o << "[";
        for (size_t j = 0; j < t[i].size(); ++j) {
            if (j) o << ",";
            o << "\"" << json_escape(t[i][j].to_string()) << "\"";
        }
        o << "]";
    }
    o << "]";
    return o.str();
}

// Shared body: parse {"word":[...], "min_order":<int>, "vars":[...]},
// dispatch to expand_zero_word or expand_inf_word based on `op_name`.
static int handle_expand(const std::string& body,
                         const char* op_name,
                         hf::SeriesTable (*fn)(const hf::PolyCtx&,
                                                const hf::Word&, long)) {
    auto letters = json_str_array(body, "word");
    auto vars = json_str_array(body, "vars");
    if (vars.empty()) vars.push_back("x");
    std::regex re("\"min_order\"\\s*:\\s*(-?\\d+)");
    std::smatch m;
    if (!std::regex_search(body, m, re)) {
        std::cerr << op_name << ": need \"min_order\"\n"; return 1;
    }
    long min_order = std::stol(m[1]);

    hf::PolyCtx ctx(vars);
    hf::Word word = parse_word(ctx, letters);
    // Use the (ctx, word, min_order) overload so the empty-word case
    // returns {{rat_one(ctx)}} instead of an empty table — see Phase
    // 6d-v-vi-0 dangling-stack-ctx fix.
    hf::SeriesTable table = fn(ctx, word, min_order);

    std::ostringstream o;
    o << "{\"op\":\"" << op_name << "\""
      << ",\"table\":" << emit_series_table(table)
      << ",\"vars\":[";
    for (size_t i = 0; i < vars.size(); ++i) {
        if (i) o << ",";
        o << "\"" << json_escape(vars[i]) << "\"";
    }
    o << "]}\n";
    std::cout << o.str();
    return 0;
}

int handle_expand_zero_word(const std::string& body) {
    return handle_expand(body, "expand_zero_word", &hf::expand_zero_word_in_ctx);
}

int handle_expand_inf_word(const std::string& body) {
    return handle_expand(body, "expand_inf_word", &hf::expand_inf_word_in_ctx);
}

// ---------- Phase 4b: MplSum ----------

int handle_mpl_sum(const std::string& body) {
    // Request: {"op":"mpl_sum","ns":[...],"zs":["..",...],"max_n":<int>,"vars":[...]}
    auto ns = json_int_array(body, "ns");
    auto zs_s = json_str_array(body, "zs");
    std::regex re("\"max_n\"\\s*:\\s*(-?\\d+)");
    std::smatch m;
    if (!std::regex_search(body, m, re)) {
        std::cerr << "mpl_sum: need \"max_n\"\n"; return 1;
    }
    long max_n = std::stol(m[1]);
    auto vars = json_str_array(body, "vars");
    if (vars.empty()) vars.push_back("x");

    hf::PolyCtx ctx(vars);
    std::vector<hf::Rat> zs;
    for (const auto& s : zs_s) zs.push_back(hf::Rat::parse(ctx, s));

    // For the zs.empty() edge case, library's mpl_sum fabricates an
    // empty PolyCtx which FLINT can't handle. Seed it with the CLI ctx
    // by providing a ghost element the library doesn't actually touch.
    hf::Rat result = zs.empty()
        ? (ns.empty() ? hf::Rat{hf::Poly(ctx, "1")}
                      : hf::Rat{hf::Poly(ctx, "0")})
        : hf::mpl_sum(ns, zs, max_n);

    std::ostringstream o;
    o << "{\"op\":\"mpl_sum\",\"result\":\""
      << json_escape(result.to_string()) << "\",\"vars\":[";
    for (size_t i = 0; i < vars.size(); ++i) {
        if (i) o << ",";
        o << "\"" << json_escape(vars[i]) << "\"";
    }
    o << "]}\n";
    std::cout << o.str();
    return 0;
}

// ---------- Phase 4c: HlogZeroExpand, HlogSeries, MplSeries ----------

// Emit an ExpansionSeries as a JSON array of term objects.
static std::string emit_expansion(const hf::ExpansionSeries& terms) {
    std::ostringstream o;
    o << "[";
    for (size_t i = 0; i < terms.size(); ++i) {
        if (i) o << ",";
        o << "{\"log_power\":" << terms[i].log_power
          << ",\"arg_power\":" << terms[i].arg_power
          << ",\"coef\":\"" << json_escape(terms[i].coef.to_string())
          << "\"}";
    }
    o << "]";
    return o.str();
}

// Emit the ExpansionSeries as a human-readable total-sum string:
//     coef_{k,j} * Log(arg)^k / k! * arg^j    (summed).
// Matches Mma's `HlogZeroExpand` InputForm which uses `Log[arg]`; we
// use `Log(arg)` so sympy can parse it directly (sympy treats `log` as
// a function; capitalized `Log` is an inert symbol there).  The
// comparator normalizes.
static std::string emit_expansion_total(const hf::ExpansionSeries& terms,
                                        const std::string& arg_str) {
    if (terms.empty()) return "0";
    std::ostringstream o;
    bool first = true;
    for (const auto& t : terms) {
        if (!first) o << " + ";
        first = false;
        o << "(" << t.coef.to_string() << ")";
        if (t.log_power > 0) {
            o << "*Log(" << arg_str << ")^" << t.log_power;
            // factorial divide
            long f = 1;
            for (long i = 1; i <= t.log_power; ++i) f *= i;
            o << "/" << f;
        }
        if (t.arg_power > 0) {
            o << "*(" << arg_str << ")^" << t.arg_power;
        }
    }
    return o.str();
}

int handle_hlog_zero_expand(const std::string& body) {
    // Request: {"op":"hlog_zero_expand","arg":<str>,"word":[...],"order":<int>,"vars":[...]}
    std::string arg_s = json_str_field(body, "arg");
    auto letters = json_str_array(body, "word");
    std::regex re("\"order\"\\s*:\\s*(-?\\d+)");
    std::smatch m;
    if (!std::regex_search(body, m, re)) {
        std::cerr << "hlog_zero_expand: need \"order\"\n"; return 1;
    }
    long order = std::stol(m[1]);
    auto vars = json_str_array(body, "vars");
    if (vars.empty()) vars.push_back("x");

    hf::PolyCtx ctx(vars);
    hf::Rat arg = hf::Rat::parse(ctx, arg_s);
    hf::Word word = parse_word(ctx, letters);

    hf::ExpansionSeries terms = hf::hlog_zero_expand(arg, word, order);
    std::string total_str = emit_expansion_total(terms, arg.to_string());

    std::ostringstream o;
    o << "{\"op\":\"hlog_zero_expand\""
      << ",\"arg\":\"" << json_escape(arg.to_string()) << "\""
      << ",\"total_str\":\"" << json_escape(total_str) << "\""
      << ",\"terms\":" << emit_expansion(terms)
      << ",\"vars\":[";
    for (size_t i = 0; i < vars.size(); ++i) {
        if (i) o << ",";
        o << "\"" << json_escape(vars[i]) << "\"";
    }
    o << "]}\n";
    std::cout << o.str();
    return 0;
}

static const char* hlog_branch_name(hf::HlogSeriesBranch b) {
    switch (b) {
        case hf::HlogSeriesBranch::kZeroLimit:      return "zero_limit";
        case hf::HlogSeriesBranch::kUnchanged:      return "unchanged";
        case hf::HlogSeriesBranch::kTaylorDeferred: return "taylor_deferred";
    }
    return "unknown";
}

int handle_hlog_series(const std::string& body) {
    // Request: {"op":"hlog_series","arg":<str>,"word":[...],"var":<str>,"order":<int>,"vars":[...]}
    std::string arg_s = json_str_field(body, "arg");
    std::string var   = json_str_field(body, "var");
    auto letters = json_str_array(body, "word");
    std::regex re("\"order\"\\s*:\\s*(-?\\d+)");
    std::smatch m;
    if (!std::regex_search(body, m, re)) {
        std::cerr << "hlog_series: need \"order\"\n"; return 1;
    }
    long order = std::stol(m[1]);
    auto vars = json_str_array(body, "vars");
    if (vars.empty()) vars.push_back("x");
    bool present = false;
    for (auto& v : vars) if (v == var) { present = true; break; }
    if (!present) vars.push_back(var);

    hf::PolyCtx ctx(vars);
    hf::Rat arg = hf::Rat::parse(ctx, arg_s);
    hf::Word word = parse_word(ctx, letters);
    size_t idx = 0;
    for (; idx < vars.size(); ++idx) if (vars[idx] == var) break;

    hf::HlogSeriesResult r = hf::hlog_series(arg, word, idx, order);

    // Build the total_str.  For kUnchanged we emit the inert Hlog(arg, w1, ...).
    // For kZeroLimit we emit the ExpansionSeries sum (no further var truncation).
    // For kTaylorDeferred we emit the empty string and signal the comparator.
    std::string total_str;
    if (r.branch == hf::HlogSeriesBranch::kUnchanged) {
        std::ostringstream tmp;
        tmp << "Hlog(" << arg.to_string();
        for (const auto& l : word.letters) tmp << "," << l.to_string();
        tmp << ")";
        total_str = tmp.str();
    } else if (r.branch == hf::HlogSeriesBranch::kZeroLimit) {
        total_str = emit_expansion_total(r.terms, arg.to_string());
    } else {
        total_str = "";
    }

    std::ostringstream o;
    o << "{\"op\":\"hlog_series\""
      << ",\"arg\":\"" << json_escape(arg.to_string()) << "\""
      << ",\"branch\":\"" << hlog_branch_name(r.branch) << "\""
      << ",\"total_str\":\"" << json_escape(total_str) << "\""
      << ",\"terms\":" << emit_expansion(r.terms)
      << ",\"vars\":[";
    for (size_t i = 0; i < vars.size(); ++i) {
        if (i) o << ",";
        o << "\"" << json_escape(vars[i]) << "\"";
    }
    o << "]}\n";
    std::cout << o.str();
    return 0;
}

static const char* mpl_branch_name(hf::MplSeriesBranch b) {
    switch (b) {
        case hf::MplSeriesBranch::kMplSum:          return "mpl_sum";
        case hf::MplSeriesBranch::kUnchanged:       return "unchanged";
        case hf::MplSeriesBranch::kPolyLogDeferred: return "polylog_deferred";
        case hf::MplSeriesBranch::kLogSingularity:  return "log_singularity";
        case hf::MplSeriesBranch::kTaylorDeferred:  return "taylor_deferred";
    }
    return "unknown";
}

int handle_mpl_series(const std::string& body) {
    // Request: {"op":"mpl_series","ns":[...],"zs":[...],"var":<str>,"order":<int>,"vars":[...]}
    auto ns   = json_int_array(body, "ns");
    auto zs_s = json_str_array(body, "zs");
    std::string var = json_str_field(body, "var");
    std::regex re("\"order\"\\s*:\\s*(-?\\d+)");
    std::smatch m;
    if (!std::regex_search(body, m, re)) {
        std::cerr << "mpl_series: need \"order\"\n"; return 1;
    }
    long order = std::stol(m[1]);
    auto vars = json_str_array(body, "vars");
    if (vars.empty()) vars.push_back("x");
    bool present = false;
    for (auto& v : vars) if (v == var) { present = true; break; }
    if (!present) vars.push_back(var);

    hf::PolyCtx ctx(vars);
    std::vector<hf::Rat> zs;
    for (const auto& s : zs_s) zs.push_back(hf::Rat::parse(ctx, s));
    size_t idx = 0;
    for (; idx < vars.size(); ++idx) if (vars[idx] == var) break;

    hf::MplSeriesResult r = hf::mpl_series(ns, zs, idx, order);

    std::ostringstream o;
    o << "{\"op\":\"mpl_series\""
      << ",\"branch\":\"" << mpl_branch_name(r.branch) << "\"";
    if (r.branch == hf::MplSeriesBranch::kMplSum) {
        o << ",\"scalar\":\"" << json_escape(r.scalar->to_string()) << "\"";
    }
    o << ",\"vars\":[";
    for (size_t i = 0; i < vars.size(); ++i) {
        if (i) o << ",";
        o << "\"" << json_escape(vars[i]) << "\"";
    }
    o << "]}\n";
    std::cout << o.str();
    return 0;
}

// ---------- Phase 5a: word regularization ----------

int handle_regzero_word(const std::string& body) {
    auto letters = json_str_array(body, "word");
    auto vars = json_str_array(body, "vars");
    if (vars.empty()) vars.push_back("x");
    hf::PolyCtx ctx(vars);
    hf::Word word = parse_word(ctx, letters);
    hf::Wordlist out = hf::regzero_word_in_ctx(ctx, word);
    std::ostringstream o;
    o << "{\"op\":\"regzero_word\",\"result\":" << emit_wordlist(out)
      << ",\"vars\":[";
    for (size_t i = 0; i < vars.size(); ++i) {
        if (i) o << ",";
        o << "\"" << json_escape(vars[i]) << "\"";
    }
    o << "]}\n";
    std::cout << o.str();
    return 0;
}

int handle_reg0(const std::string& body) {
    auto vars = json_str_array(body, "vars");
    if (vars.empty()) vars.push_back("x");
    hf::PolyCtx ctx(vars);
    hf::Wordlist wl = parse_wordlist(ctx, body, "wl");
    hf::Wordlist out = hf::reg0(wl);
    std::ostringstream o;
    o << "{\"op\":\"reg0\",\"result\":" << emit_wordlist(out)
      << ",\"vars\":[";
    for (size_t i = 0; i < vars.size(); ++i) {
        if (i) o << ",";
        o << "\"" << json_escape(vars[i]) << "\"";
    }
    o << "]}\n";
    std::cout << o.str();
    return 0;
}

// Helper: parse {"letter":<str>, "substitute":<str>} using the given ctx.
static hf::Letter parse_letter_field(const hf::PolyCtx& ctx,
                                     const std::string& body,
                                     const std::string& key,
                                     const std::string& dflt) {
    std::string s = json_str_field(body, key);
    if (s.empty()) s = dflt;
    return hf::Rat::parse(ctx, s);
}

int handle_reg_side(const std::string& body, const char* op_name,
                    bool head) {
    auto vars = json_str_array(body, "vars");
    if (vars.empty()) vars.push_back("x");
    hf::PolyCtx ctx(vars);
    hf::Wordlist wl = parse_wordlist(ctx, body, "wl");
    hf::Letter letter     = parse_letter_field(ctx, body, "letter", "0");
    hf::Letter substitute = parse_letter_field(ctx, body, "substitute", "0");
    hf::Wordlist out = head
        ? hf::reg_head(wl, letter, substitute)
        : hf::reg_tail(wl, letter, substitute);
    std::ostringstream o;
    o << "{\"op\":\"" << op_name << "\",\"result\":" << emit_wordlist(out)
      << ",\"vars\":[";
    for (size_t i = 0; i < vars.size(); ++i) {
        if (i) o << ",";
        o << "\"" << json_escape(vars[i]) << "\"";
    }
    o << "]}\n";
    std::cout << o.str();
    return 0;
}

int handle_reg_head(const std::string& body) {
    return handle_reg_side(body, "reg_head", /*head=*/true);
}
int handle_reg_tail(const std::string& body) {
    return handle_reg_side(body, "reg_tail", /*head=*/false);
}

// ---------- Phase 5c-i: Regulator I/O helpers ----------

// Serialize a Regulator as
//   [ {"coef":"<rat>", "key":[ ["<l1>",...], ... ]}, ... ]
// where each inner array is a Word (list of letter strings).
static std::string emit_regulator(const hf::Regulator& r) {
    std::ostringstream o;
    o << "[";
    for (size_t i = 0; i < r.size(); ++i) {
        if (i) o << ",";
        o << "{\"coef\":\"" << json_escape(r[i].coef.to_string()) << "\""
          << ",\"key\":[";
        for (size_t j = 0; j < r[i].key.size(); ++j) {
            if (j) o << ",";
            o << "[";
            for (size_t k = 0; k < r[i].key[j].size(); ++k) {
                if (k) o << ",";
                o << "\"" << json_escape(r[i].key[j][k].to_string()) << "\"";
            }
            o << "]";
        }
        o << "]}";
    }
    o << "]";
    return o.str();
}

// Extract a nested array "key":[ [l1,l2,...], [l1,l2,...], ... ] from
// a single {"coef":"...","key":[...]} object (as a JSON string body)
// and return the list of parsed Words. parse_word_array expects the
// JSON body ALREADY narrowed to the single regulator-term object.
static std::vector<hf::Word> parse_regkey_words(
    const hf::PolyCtx& ctx, const std::string& term) {
    std::vector<hf::Word> out;
    // Locate "key": [ ... ]  bracket-balanced top-level array.
    std::regex re("\"key\"\\s*:\\s*\\[");
    std::smatch m;
    if (!std::regex_search(term, m, re)) return out;
    size_t open = m.position(0) + m.length(0) - 1;
    int depth = 0;
    size_t start = open + 1, end = std::string::npos;
    for (size_t i = open; i < term.size(); ++i) {
        if (term[i] == '[') depth++;
        else if (term[i] == ']') {
            depth--;
            if (depth == 0) { end = i; break; }
        }
    }
    if (end == std::string::npos) return out;
    std::string inner = term.substr(start, end - start);
    // Split `inner` into top-level [ ... ] word-arrays.
    depth = 0;
    size_t wstart = 0;
    for (size_t i = 0; i < inner.size(); ++i) {
        if (inner[i] == '[') {
            if (depth == 0) wstart = i;
            depth++;
        } else if (inner[i] == ']') {
            depth--;
            if (depth == 0) {
                // inner[wstart..i] is "[\"l1\",\"l2\",...]"
                std::string word_json = inner.substr(wstart, i - wstart + 1);
                // Reuse json_str_array by wrapping as {"w":[...]}.
                std::string wrapped = "{\"w\":" + word_json + "}";
                auto letters = json_str_array(wrapped, "w");
                out.push_back(parse_word(ctx, letters));
            }
        }
    }
    return out;
}

// Parse a Regulator: [ {"coef":..,"key":[[...],...]}, ... ]
static hf::Regulator parse_regulator(const hf::PolyCtx& ctx,
                                     const std::string& body,
                                     const std::string& key) {
    hf::Regulator out;
    std::string inner = extract_top_array(body, key);
    if (inner.empty()) return out;
    // Split `inner` into top-level {...} objects (curly-balanced).
    std::vector<std::string> terms;
    int depth = 0;
    size_t start = 0;
    for (size_t i = 0; i < inner.size(); ++i) {
        if (inner[i] == '{') {
            if (depth == 0) start = i;
            depth++;
        } else if (inner[i] == '}') {
            depth--;
            if (depth == 0) terms.push_back(inner.substr(start, i - start + 1));
        }
    }
    for (const auto& t : terms) {
        std::string coef_s = json_str_field(t, "coef");
        auto words = parse_regkey_words(ctx, t);
        out.push_back(hf::RegTerm{hf::Rat::parse(ctx, coef_s),
                                    std::move(words)});
    }
    return out;
}

int handle_shuffle_symbolic(const std::string& body) {
    // Request: {"op":"shuffle_symbolic","a":<regulator>,"b":<regulator>,"vars":[...]}
    auto vars = json_str_array(body, "vars");
    if (vars.empty()) vars.push_back("x");
    hf::PolyCtx ctx(vars);
    hf::Regulator A = parse_regulator(ctx, body, "a");
    hf::Regulator B = parse_regulator(ctx, body, "b");
    hf::Regulator out = hf::shuffle_symbolic(A, B);
    std::ostringstream o;
    o << "{\"op\":\"shuffle_symbolic\",\"result\":" << emit_regulator(out)
      << ",\"vars\":[";
    for (size_t i = 0; i < vars.size(); ++i) {
        if (i) o << ",";
        o << "\"" << json_escape(vars[i]) << "\"";
    }
    o << "]}\n";
    std::cout << o.str();
    return 0;
}

// Serialize a TransformResult as
//   [ {"shuffle":[{"coef":..,"word":[...]},...],
//      "regulator":[{"coef":..,"key":[[...],...]},...]}, ... ]
static std::string emit_transform_result(const hf::TransformResult& tr) {
    std::ostringstream o;
    o << "[";
    for (size_t i = 0; i < tr.size(); ++i) {
        if (i) o << ",";
        o << "{\"shuffle\":" << emit_wordlist(tr[i].shuffle)
          << ",\"regulator\":" << emit_regulator(tr[i].regulator) << "}";
    }
    o << "]";
    return o.str();
}

// Forward decls (emit_regulator_sym is defined later alongside other
// sym helpers; emit_transform_result_sym wraps it).
static std::string emit_regulator_sym(const hf::RegulatorSym& r);
static std::string emit_transform_result_sym(const hf::TransformResultSym& tr) {
    std::ostringstream o;
    o << "[";
    for (size_t i = 0; i < tr.size(); ++i) {
        if (i) o << ",";
        o << "{\"shuffle\":" << emit_wordlist(tr[i].shuffle)
          << ",\"regulator\":" << emit_regulator_sym(tr[i].regulator) << "}";
    }
    o << "]";
    return o.str();
}

// Bug #6 lift: the per-step regulator chain (`reglim_word`,
// `transform_word`, `transform_shuffle`, `integration_step`) returns
// SymCoef-valued results now. Pre-existing CLI fixtures (that
// compared against Mma/Maple Regulator JSON) stay compatible because
// both `emit_regulator_sym` (below, for the sym coef string form)
// and `_sympy_canon_coef` in the comparator understand
// `I*Pi*delta[var]` residues. When the output is Rat-pure,
// `SymCoef::to_string` emits a plain rational string that's
// indistinguishable from the old Rat output. The narrow-ctx
// handlers below load the MZV reduction table and extend `vars`
// via `build_mzv_var_list` so mzv_* symbols emitted by Fragment
// P1/P2 can be parsed by Rat::parse.
static std::string resolve_mzv_data_path(const std::string& body);

int handle_transform_word(const std::string& body) {
    // Request: {"op":"transform_word","word":[...],"var":<str>,"vars":[...],
    //           "algebraic_letters":bool (optional)}
    auto letters = json_str_array(body, "word");
    std::string var = json_str_field(body, "var");
    auto vars = json_str_array(body, "vars");
    if (vars.empty()) vars.push_back(var.empty() ? "x" : var);
    if (!var.empty()) {
        bool present = false;
        for (const auto& v : vars) if (v == var) { present = true; break; }
        if (!present) vars.push_back(var);
    }
    // Bug #6: reglim_word's Fragment P1/P2 emits mzv_* symbols;
    // extend ctx with the MZV basis so Rat::parse succeeds on the
    // algebraic-string coefficients. Load the full reduction
    // table here (same path resolution as the hyperflint op).
    std::string data_path = resolve_mzv_data_path(body);
    hf::MzvReductionTable table = hf::load_mzv_reductions(data_path);
    vars = hf::build_mzv_var_list(table, vars);
    // Phase 7-vi-a: algebraic_letters flag — pre-allocate Wm/Wp pool.
    bool introduce_al = false;
    {
        std::regex re("\"algebraic_letters\"\\s*:\\s*(true|false)");
        std::smatch m;
        if (std::regex_search(body, m, re)) {
            introduce_al = (m[1] == "true");
        }
    }
    std::vector<std::string> ctx_vars = introduce_al
        ? hf::build_algebraic_letter_var_list(vars)
        : vars;
    hf::PolyCtx ctx(ctx_vars);
    hf::Word word = parse_word(ctx, letters);
    size_t idx = 0;
    if (!var.empty()) {
        for (; idx < ctx_vars.size(); ++idx) if (ctx_vars[idx] == var) break;
    }
    try {
        // Iter-52 C0c.1 Increment β: caller-side fresh transient ZWTable.
        auto _lf_zw = std::make_shared<hf::ZWTable>(ctx);
        hf::TransformResultSym tr = hf::transform_word(ctx, word, idx,
                                                         table,
                                                         _lf_zw,
                                                         introduce_al);
        std::ostringstream o;
        o << "{\"op\":\"transform_word\",\"result\":" << emit_transform_result_sym(tr)
          << ",\"vars\":[";
        for (size_t i = 0; i < vars.size(); ++i) {
            if (i) o << ",";
            o << "\"" << json_escape(vars[i]) << "\"";
        }
        o << "]}\n";
        std::cout << o.str();
    } catch (const hf::TransformFailed&) {
        std::cout << "{\"op\":\"transform_word\",\"failed\":true,\"vars\":[";
        for (size_t i = 0; i < vars.size(); ++i) {
            if (i) std::cout << ",";
            std::cout << "\"" << json_escape(vars[i]) << "\"";
        }
        std::cout << "]}\n";
    }
    return 0;
}

// Parse a ShuffleList from a top-level `"wordlist":[...]` array. Each
// entry is {"coef":"..", "shuffle":[[l1,..], ...]}. Shared with
// integration_step and hyperflint.
static hf::ShuffleList parse_shuffle_list(const hf::PolyCtx& ctx,
                                          const std::string& body,
                                          const std::string& key) {
    hf::ShuffleList input;
    std::string inner = extract_top_array(body, key);
    if (inner.empty()) return input;
    int depth = 0;
    size_t start = 0;
    std::vector<std::string> entries;
    for (size_t i = 0; i < inner.size(); ++i) {
        if (inner[i] == '{') {
            if (depth == 0) start = i;
            depth++;
        } else if (inner[i] == '}') {
            depth--;
            if (depth == 0) entries.push_back(inner.substr(start, i - start + 1));
        }
    }
    for (const auto& e : entries) {
        std::string coef_s = json_str_field(e, "coef");
        hf::ShuffleEntry ent{hf::Rat::parse(ctx, coef_s), {}};
        std::string sh_inner = extract_top_array(e, "shuffle");
        if (!sh_inner.empty()) {
            int wd = 0;
            size_t ws = 0;
            for (size_t i = 0; i < sh_inner.size(); ++i) {
                if (sh_inner[i] == '[') {
                    if (wd == 0) ws = i;
                    wd++;
                } else if (sh_inner[i] == ']') {
                    wd--;
                    if (wd == 0) {
                        std::string w_json = sh_inner.substr(ws, i - ws + 1);
                        std::string wrapped = "{\"w\":" + w_json + "}";
                        auto letters = json_str_array(wrapped, "w");
                        ent.shuffle.push_back(parse_word(ctx, letters));
                    }
                }
            }
        }
        input.push_back(std::move(ent));
    }
    return input;
}

// Forward declarations for helpers defined later in this file.
static std::string resolve_mzv_data_path(const std::string& body);
static std::string emit_regulator_sym(const hf::RegulatorSym& r);

int handle_hyperflint(const std::string& body) {
    // Phase γ.2: forwards to the transport-neutral handler so both the
    // CLI and the LibraryLink shared library exercise identical logic.
    // The handler returns a JSON response string; the CLI prints it to
    // stdout with a trailing newline.
    std::cout << hyperflint::handlers::hyperflint_sym(body) << "\n";
    return 0;
}

int handle_integration_step(const std::string& body) {
    // Request:
    //   {"op":"integration_step",
    //    "wordlist":[{"coef":"<rat>", "shuffle":[[l1,..], ...]}, ...],
    //    "var":<str>, "vars":[...],
    //    "algebraic_letters":bool (optional)}
    std::string var = json_str_field(body, "var");
    if (var.empty()) { std::cerr << "integration_step: need \"var\"\n"; return 1; }
    auto vars = json_str_array(body, "vars");
    if (vars.empty()) vars.push_back(var);
    bool present = false;
    for (const auto& v : vars) if (v == var) { present = true; break; }
    if (!present) vars.push_back(var);

    bool introduce_al = false;
    {
        std::regex re("\"algebraic_letters\"\\s*:\\s*(true|false)");
        std::smatch m;
        if (std::regex_search(body, m, re)) {
            introduce_al = (m[1] == "true");
        }
    }
    // Bug #6: extend ctx with MZV basis (see handle_transform_word).
    std::string data_path = resolve_mzv_data_path(body);
    hf::MzvReductionTable table = hf::load_mzv_reductions(data_path);
    vars = hf::build_mzv_var_list(table, vars);
    std::vector<std::string> ctx_vars = introduce_al
        ? hf::build_algebraic_letter_var_list(vars)
        : vars;
    hf::PolyCtx ctx(ctx_vars);
    // Parse ShuffleList: [{"coef":"..","shuffle":[[l1,..], ...]}, ...]
    hf::ShuffleList input;
    std::string inner = extract_top_array(body, "wordlist");
    if (!inner.empty()) {
        int depth = 0;
        size_t start = 0;
        std::vector<std::string> entries;
        for (size_t i = 0; i < inner.size(); ++i) {
            if (inner[i] == '{') {
                if (depth == 0) start = i;
                depth++;
            } else if (inner[i] == '}') {
                depth--;
                if (depth == 0) entries.push_back(inner.substr(start, i - start + 1));
            }
        }
        for (const auto& e : entries) {
            std::string coef_s = json_str_field(e, "coef");
            hf::ShuffleEntry ent{hf::Rat::parse(ctx, coef_s), {}};
            // "shuffle" is a top-level array of word-arrays.
            std::string sh_inner = extract_top_array(e, "shuffle");
            if (!sh_inner.empty()) {
                int wd = 0;
                size_t ws = 0;
                for (size_t i = 0; i < sh_inner.size(); ++i) {
                    if (sh_inner[i] == '[') {
                        if (wd == 0) ws = i;
                        wd++;
                    } else if (sh_inner[i] == ']') {
                        wd--;
                        if (wd == 0) {
                            std::string w_json = sh_inner.substr(ws, i - ws + 1);
                            std::string wrapped = "{\"w\":" + w_json + "}";
                            auto letters = json_str_array(wrapped, "w");
                            ent.shuffle.push_back(parse_word(ctx, letters));
                        }
                    }
                }
            }
            input.push_back(std::move(ent));
        }
    }
    size_t idx = 0;
    for (; idx < ctx_vars.size(); ++idx) if (ctx_vars[idx] == var) break;

    // Phase 5e-iii: optional "check_divergences" flag routes through the
    // boundary-divergence scan inside integration_step.
    //
    // DP.3 (divergence policy, 2026-06-03): the intended default for bare
    // requests is TRUE (independent single-integral usage; silent wrong
    // answers on divergent inputs are the worst failure mode), but the
    // flip is BLOCKED on HF-DIVCHECK-PARITY: the scan currently
    // false-positives on generic multi-pole CONVERGENT integrands (even
    // 1/((x+1)(x+2)) over x, value ln 2), because test_zero_function_sym
    // does not reduce cross-key bins to period VALUES the way the
    // reference TestZeroFunction (HyperIntica.wl:5050) does -- it tests
    // per-key after fibration over the supplied schedule vars only. See
    // notes/hf_divcheck_parity.md. Flip this default only after parity.
    bool check_div = false;
    {
        std::regex re("\"check_divergences\"\\s*:\\s*(true|false)");
        std::smatch m;
        if (std::regex_search(body, m, re)) check_div = (m[1] == "true");
    }

    // Caveat-3 plumbing: optional "remaining_vars" field lets callers
    // pass the outer integration schedule so the divergence check can
    // project residues onto the fibration basis over those vars (catches
    // shuffle-identity cancellations across distinct keys). Empty list
    // (absent field) falls back to base-case per-term is_zero.
    std::vector<size_t> remaining_idx;
    if (check_div) {
        auto remaining = json_str_array(body, "remaining_vars");
        for (const auto& rv : remaining) {
            for (size_t i = 0; i < ctx_vars.size(); ++i) {
                if (ctx_vars[i] == rv) { remaining_idx.push_back(i); break; }
            }
        }
    }

    try {
        // Iter-52 C0c.1 Increment β: caller-side fresh transient ZWTable.
        auto _lf_zw = std::make_shared<hf::ZWTable>(ctx);
        hf::RegulatorSym out =
            hf::integration_step(ctx, input, idx, table, _lf_zw,
                                  /*check_divergences=*/check_div,
                                  introduce_al,
                                  remaining_idx);
        std::ostringstream o;
        o << "{\"op\":\"integration_step\",\"result\":" << emit_regulator_sym(out)
          << ",\"vars\":[";
        for (size_t i = 0; i < vars.size(); ++i) {
            if (i) o << ",";
            o << "\"" << json_escape(vars[i]) << "\"";
        }
        o << "]}\n";
        std::cout << o.str();
    } catch (const hf::HyperFLINTDivergentIntegral& e) {
        std::cout << "{\"op\":\"integration_step\",\"divergent\":true"
                  << ",\"reason\":\"" << json_escape(e.what()) << "\""
                  << ",\"vars\":[";
        for (size_t i = 0; i < vars.size(); ++i) {
            if (i) std::cout << ",";
            std::cout << "\"" << json_escape(vars[i]) << "\"";
        }
        std::cout << "]}\n";
    } catch (const hf::IntegrationStepFailed&) {
        std::cout << "{\"op\":\"integration_step\",\"failed\":true,\"vars\":[";
        for (size_t i = 0; i < vars.size(); ++i) {
            if (i) std::cout << ",";
            std::cout << "\"" << json_escape(vars[i]) << "\"";
        }
        std::cout << "]}\n";
    }
    return 0;
}

static std::string resolve_mzv_data_path(const std::string& body) {
    std::string data_path = json_str_field(body, "mzv_data_path");
    if (!data_path.empty()) return data_path;
    const char* env = std::getenv("HYPERFLINT_DATA_DIR");
    if (env && *env) return std::string(env) + "/mzv_reductions.json";
    return "data/mzv_reductions.json";
}

int handle_test_zero_function(const std::string& body) {
    // Request: {"op":"test_zero_function","regulator":[...],"vars":[...]}
    auto user_vars = json_str_array(body, "vars");
    std::string data_path = resolve_mzv_data_path(body);
    hf::MzvReductionTable table = hf::load_mzv_reductions(data_path);
    std::vector<std::string> vars = hf::build_mzv_var_list(table, user_vars);
    hf::PolyCtx ctx(vars);
    hf::Regulator r = parse_regulator(ctx, body, "regulator");
    hf::Rat result = hf::test_zero_function(ctx, r, table);
    std::cout << emit_result_json("test_zero_function", vars,
                                   result.to_string()) << "\n";
    return 0;
}

int handle_evaluate_periods(const std::string& body) {
    // Request:
    //   {"op":"evaluate_periods",
    //    "regulator":[{"coef":..,"key":[[...],...]},...],
    //    "vars":[...]}
    auto user_vars = json_str_array(body, "vars");
    std::string data_path = resolve_mzv_data_path(body);
    hf::MzvReductionTable table = hf::load_mzv_reductions(data_path);
    std::vector<std::string> vars = hf::build_mzv_var_list(table, user_vars);
    hf::PolyCtx ctx(vars);
    hf::Regulator r = parse_regulator(ctx, body, "regulator");
    hf::Regulator out = hf::evaluate_periods(ctx, r, table);
    std::ostringstream o;
    o << "{\"op\":\"evaluate_periods\",\"result\":" << emit_regulator(out)
      << ",\"vars\":[";
    for (size_t i = 0; i < vars.size(); ++i) {
        if (i) o << ",";
        o << "\"" << json_escape(vars[i]) << "\"";
    }
    o << "]}\n";
    std::cout << o.str();
    return 0;
}

// Forward decl (definition lives below with the other SymCoef helpers).
static std::string sym_coef_to_mma_string(const hf::SymCoef& s);

// Phase 6d-v-ii: emit a RegulatorSym in the same JSON shape as a
// regular Regulator. SymCoef coef -> flat algebraic string the
// comparator's _sympy_canon_coef rewrites understand.
static std::string emit_regulator_sym(const hf::RegulatorSym& r) {
    std::ostringstream o;
    o << "[";
    for (size_t i = 0; i < r.size(); ++i) {
        if (i) o << ",";
        // Reuse sym_coef_to_mma_string by pretending it lives in `hf`.
        // (It's a static helper above, so call directly.)
        o << "{\"coef\":\"" << json_escape(sym_coef_to_mma_string(r[i].coef)) << "\""
          << ",\"key\":[";
        for (size_t j = 0; j < r[i].key.size(); ++j) {
            if (j) o << ",";
            o << "[";
            for (size_t k = 0; k < r[i].key[j].size(); ++k) {
                if (k) o << ",";
                o << "\"" << json_escape(r[i].key[j][k].to_string()) << "\"";
            }
            o << "]";
        }
        o << "]}";
    }
    o << "]";
    return o.str();
}

int handle_break_up_contour_sym(const std::string& body) {
    // Request: {"op":"break_up_contour_sym",
    //           "wl":[{"coef":"<rat>","word":["l1",...]}, ...],
    //           "on_axis":[{"letter":"<rat>","im_var":"<varname>"}, ...],
    //           "vars":[...]}
    auto user_vars = json_str_array(body, "vars");
    std::string data_path = resolve_mzv_data_path(body);
    hf::MzvReductionTable table = hf::load_mzv_reductions(data_path);
    std::vector<std::string> vars = hf::build_mzv_var_list(table, user_vars);

    hf::PolyCtx ctx(vars);
    hf::Wordlist wl = parse_wordlist(ctx, body, "wl");
    hf::WordlistSym wl_sym = hf::to_wordlist_sym(wl);

    // Parse on_axis: list of {letter, im_var} objects.
    std::vector<hf::OnAxisSymEntry> on_axis;
    std::string oa_inner = extract_top_array(body, "on_axis");
    if (!oa_inner.empty()) {
        int depth = 0;
        size_t start = 0;
        std::vector<std::string> entries;
        for (size_t i = 0; i < oa_inner.size(); ++i) {
            if (oa_inner[i] == '{') {
                if (depth == 0) start = i;
                depth++;
            } else if (oa_inner[i] == '}') {
                depth--;
                if (depth == 0) entries.push_back(oa_inner.substr(start, i - start + 1));
            }
        }
        for (const auto& e : entries) {
            std::string letter_s = json_str_field(e, "letter");
            std::string im_var   = json_str_field(e, "im_var");
            if (letter_s.empty() || im_var.empty()) {
                std::cerr << "break_up_contour_sym: each on_axis entry needs"
                             " \"letter\" and \"im_var\"\n";
                return 1;
            }
            on_axis.push_back(hf::OnAxisSymEntry{
                hf::Rat::parse(ctx, letter_s),
                hf::SymCoef::delta_factor(ctx, im_var)});
        }
    }

    // C0b.4 (iter-42): break_up_contour_sym now takes a mandatory
    // `std::shared_ptr<ZWTable>` parameter. The CLI bridge handler for
    // `op=break_up_contour_sym` is a standalone debug/probe path (not in
    // the production integration); allocate a per-call table when the
    // env-gate is on to mirror pre-iter-42 lifetime, leave null at
    // default-OFF (the body lambda short-circuits before any `zw_tab`
    // use when `runtime::scalar_rep_enabled()` is false).
    std::shared_ptr<hf::ZWTable> bcs_zw_tab_local;
    if (hf::runtime::scalar_rep_enabled()) {
        // Iter-44 (2026-05-09): HF_SCALAR_REP_REQUIRE_PERSISTENT=1
        // assertion (Concern-2 mitigation). See scalar_rep.hpp.
        if (hf::runtime::require_persistent_enabled()) {
            std::cerr << "[HF_SCALAR_REP_REQUIRE_PERSISTENT=1]"
                << " bridge/cli/main.cpp:op=break_up_contour_sym"
                << " (around line 2082): allocating per-call ZWTable"
                << " in CLI debug/probe path. This callsite is"
                << " standalone (not on the integration hot path)"
                << " and may be left transitional indefinitely;"
                << " when migrating, supply the table from a"
                << " driver-level allocation analogous to"
                << " hyperflint_sym." << std::endl;
            std::abort();
        }
        bcs_zw_tab_local = std::make_shared<hf::ZWTable>(ctx);
    }
    hf::RegulatorSym out =
        hf::break_up_contour_sym(ctx, wl_sym, on_axis, table,
                                   bcs_zw_tab_local);
    std::ostringstream o;
    o << "{\"op\":\"break_up_contour_sym\",\"result\":" << emit_regulator_sym(out)
      << ",\"vars\":[";
    for (size_t i = 0; i < vars.size(); ++i) {
        if (i) o << ",";
        o << "\"" << json_escape(vars[i]) << "\"";
    }
    o << "]}\n";
    std::cout << o.str();
    return 0;
}

int handle_break_up_contour(const std::string& body) {
    // Request: {"op":"break_up_contour","wl":[{"coef":..,"word":[...]}],
    //           "on_axis":[],    // Phase 6d-ii: empty only
    //           "vars":[...]}
    auto user_vars = json_str_array(body, "vars");
    std::string data_path = resolve_mzv_data_path(body);
    hf::MzvReductionTable table = hf::load_mzv_reductions(data_path);
    std::vector<std::string> vars = hf::build_mzv_var_list(table, user_vars);
    hf::PolyCtx ctx(vars);
    hf::Wordlist wl = parse_wordlist(ctx, body, "wl");
    // on_axis not yet supported beyond empty.
    std::vector<hf::OnAxisEntry> on_axis;
    hf::Regulator out = hf::break_up_contour(ctx, wl, on_axis, table);
    std::ostringstream o;
    o << "{\"op\":\"break_up_contour\",\"result\":" << emit_regulator(out)
      << ",\"vars\":[";
    for (size_t i = 0; i < vars.size(); ++i) {
        if (i) o << ",";
        o << "\"" << json_escape(vars[i]) << "\"";
    }
    o << "]}\n";
    std::cout << o.str();
    return 0;
}

int handle_zero_inf_period(const std::string& body) {
    auto letters = json_str_array(body, "word");
    auto user_vars = json_str_array(body, "vars");
    std::string data_path = resolve_mzv_data_path(body);

    hf::MzvReductionTable table = hf::load_mzv_reductions(data_path);
    std::vector<std::string> vars = hf::build_mzv_var_list(table, user_vars);

    hf::PolyCtx ctx(vars);
    hf::Word word = parse_word(ctx, letters);
    hf::Rat r = hf::zero_inf_period(ctx, word, table);
    std::cout << emit_result_json("zero_inf_period", vars, r.to_string()) << "\n";
    return 0;
}

int handle_zero_one_period(const std::string& body) {
    // Request: {"op":"zero_one_period","word":["0","-1",...], "vars":[...]}
    auto letters = json_str_array(body, "word");
    auto user_vars = json_str_array(body, "vars");
    std::string data_path = resolve_mzv_data_path(body);

    hf::MzvReductionTable table = hf::load_mzv_reductions(data_path);
    std::vector<std::string> vars = hf::build_mzv_var_list(table, user_vars);

    hf::PolyCtx ctx(vars);
    hf::Word word = parse_word(ctx, letters);
    hf::Rat r = hf::zero_one_period(ctx, word, table);
    std::cout << emit_result_json("zero_one_period", vars, r.to_string()) << "\n";
    return 0;
}

int handle_apply_mzv_reductions(const std::string& body) {
    // Request: {"op":"apply_mzv_reductions","f":"<rat>","vars":[...],
    //           "mzv_data_path": <optional path>}
    // Default mzv_data_path is <binary-dir>/../data/mzv_reductions.json
    // via the HYPERFLINT_DATA_DIR env var.
    std::string f = json_str_field(body, "f");
    if (f.empty()) { std::cerr << "apply_mzv_reductions: need \"f\"\n"; return 1; }
    auto vars = json_str_array(body, "vars");
    if (vars.empty()) vars = autoscan_vars({&f});

    std::string data_path = resolve_mzv_data_path(body);

    hf::PolyCtx ctx(vars);
    hf::Rat rf = hf::Rat::parse(ctx, f);
    hf::MzvReductionTable table = hf::load_mzv_reductions(data_path);
    hf::Rat r = hf::apply_mzv_reductions(ctx, table, rf);
    std::cout << emit_result_json("apply_mzv_reductions", vars, r.to_string()) << "\n";
    return 0;
}

// ---------- Phase 6d-v-i: SymCoef ----------

// Parse a single monomial JSON object {"prefactor":"<rat>", "pi":<int>,
// "i":<int>, "logs":[[n,k],...], "deltas":[["v",k],...]}.
static hf::SymCoef parse_sym_monomial(const hf::PolyCtx& ctx,
                                      const std::string& obj) {
    std::string pref = json_str_field(obj, "prefactor");
    if (pref.empty()) pref = "1";
    hf::SymCoef m = hf::SymCoef::from_rat(hf::Rat::parse(ctx, pref));

    std::regex pi_re("\"pi\"\\s*:\\s*(-?\\d+)");
    std::smatch sm;
    if (std::regex_search(obj, sm, pi_re)) {
        long n = std::stol(sm[1]);
        if (n < 0) throw std::runtime_error("sym_arith: pi power must be >= 0");
        for (long i = 0; i < n; ++i) m = m * hf::SymCoef::pi_factor(ctx);
    }

    std::regex i_re("\"i\"\\s*:\\s*(-?\\d+)");
    if (std::regex_search(obj, sm, i_re)) {
        long n = std::stol(sm[1]);
        if (n < 0) throw std::runtime_error("sym_arith: i power must be >= 0");
        for (long i = 0; i < n; ++i) m = m * hf::SymCoef::im_factor(ctx);
    }

    // logs: array of [n, k] pairs
    std::string logs_inner = extract_top_array(obj, "logs");
    if (!logs_inner.empty()) {
        int depth = 0;
        size_t start = 0;
        std::vector<std::string> entries;
        for (size_t i = 0; i < logs_inner.size(); ++i) {
            if (logs_inner[i] == '[') {
                if (depth == 0) start = i;
                depth++;
            } else if (logs_inner[i] == ']') {
                depth--;
                if (depth == 0) entries.push_back(logs_inner.substr(start, i - start + 1));
            }
        }
        std::regex pair_re("\\[\\s*(-?\\d+)\\s*,\\s*(-?\\d+)\\s*\\]");
        for (const auto& e : entries) {
            std::smatch pm;
            if (!std::regex_match(e, pm, pair_re)) {
                throw std::runtime_error("sym_arith: malformed logs entry: " + e);
            }
            long n = std::stol(pm[1]);
            long k = std::stol(pm[2]);
            if (k < 0) throw std::runtime_error("sym_arith: log power must be >= 0");
            for (long i = 0; i < k; ++i) m = m * hf::SymCoef::log_factor(ctx, n);
        }
    }

    // deltas: array of ["name", k] pairs
    std::string deltas_inner = extract_top_array(obj, "deltas");
    if (!deltas_inner.empty()) {
        int depth = 0;
        size_t start = 0;
        std::vector<std::string> entries;
        for (size_t i = 0; i < deltas_inner.size(); ++i) {
            if (deltas_inner[i] == '[') {
                if (depth == 0) start = i;
                depth++;
            } else if (deltas_inner[i] == ']') {
                depth--;
                if (depth == 0) entries.push_back(deltas_inner.substr(start, i - start + 1));
            }
        }
        std::regex pair_re("\\[\\s*\"([^\"]+)\"\\s*,\\s*(-?\\d+)\\s*\\]");
        for (const auto& e : entries) {
            std::smatch pm;
            if (!std::regex_match(e, pm, pair_re)) {
                throw std::runtime_error("sym_arith: malformed deltas entry: " + e);
            }
            std::string vname = pm[1];
            long k = std::stol(pm[2]);
            if (k < 0) throw std::runtime_error("sym_arith: delta power must be >= 0");
            for (long i = 0; i < k; ++i) m = m * hf::SymCoef::delta_factor(ctx, vname);
        }
    }
    return m;
}

// Parse a SymCoef as a list of monomials at body["<key>"]; sums them.
static hf::SymCoef parse_sym_coef(const hf::PolyCtx& ctx,
                                   const std::string& body,
                                   const std::string& key) {
    hf::SymCoef out(ctx);
    std::string inner = extract_top_array(body, key);
    if (inner.empty()) return out;
    int depth = 0;
    size_t start = 0;
    std::vector<std::string> entries;
    for (size_t i = 0; i < inner.size(); ++i) {
        if (inner[i] == '{') {
            if (depth == 0) start = i;
            depth++;
        } else if (inner[i] == '}') {
            depth--;
            if (depth == 0) entries.push_back(inner.substr(start, i - start + 1));
        }
    }
    for (const auto& e : entries) {
        out = out + parse_sym_monomial(ctx, e);
    }
    return out;
}

// Render a SymCoef as a flat algebraic-string Mma can parse (Pi, I,
// Log[n], delta[var]). Empty SymCoef -> "0".
static std::string sym_coef_to_mma_string(const hf::SymCoef& s) {
    if (s.is_zero()) return "0";
    std::ostringstream o;
    bool first = true;
    for (const auto& m : s.terms()) {
        std::string pre = m.prefactor.to_string();
        // Wrap negative or non-trivial prefactors in parens for safe parsing.
        bool needs_paren = pre.find('+') != std::string::npos
                        || pre.find('-') != std::string::npos
                        || pre.find('/') != std::string::npos
                        || pre.find('*') != std::string::npos;
        if (!first) o << " + ";
        if (needs_paren) o << "(" << pre << ")"; else o << pre;
        if (m.pi_power == 1)      o << "*Pi";
        else if (m.pi_power != 0) o << "*Pi^" << m.pi_power;
        if (m.i_power == 1)       o << "*I";
        else if (m.i_power != 0)  o << "*I^" << m.i_power;
        for (const auto& kv : m.log_powers) {
            if (kv.second == 1) o << "*Log[" << kv.first << "]";
            else                o << "*Log[" << kv.first << "]^" << kv.second;
        }
        for (const auto& kv : m.delta_powers) {
            if (kv.second == 1) o << "*delta[" << kv.first << "]";
            else                o << "*delta[" << kv.first << "]^" << kv.second;
        }
        first = false;
    }
    return o.str();
}

int handle_sym_arith(const std::string& body) {
    // Request: {"op":"sym_arith",
    //           "a":[{monomial}, ...],
    //           "b":[{monomial}, ...],
    //           "mode":"add"|"sub"|"mul",
    //           "vars":[...]}
    std::string mode = json_str_field(body, "mode");
    if (mode.empty()) { std::cerr << "sym_arith: need \"mode\"\n"; return 1; }
    auto user_vars = json_str_array(body, "vars");
    std::string data_path = resolve_mzv_data_path(body);
    hf::MzvReductionTable table = hf::load_mzv_reductions(data_path);
    std::vector<std::string> vars = hf::build_mzv_var_list(table, user_vars);

    hf::PolyCtx ctx(vars);
    hf::SymCoef A = parse_sym_coef(ctx, body, "a");
    hf::SymCoef B = parse_sym_coef(ctx, body, "b");
    hf::SymCoef R(ctx);
    if      (mode == "add") R = A + B;
    else if (mode == "sub") R = A - B;
    else if (mode == "mul") R = A * B;
    else { std::cerr << "sym_arith: unknown mode '" << mode << "'\n"; return 1; }

    std::string result_str = sym_coef_to_mma_string(R);
    std::cout << emit_result_json("sym_arith", vars, result_str) << "\n";
    return 0;
}

int handle_sym_reduce(const std::string& body) {
    // Request: {"op":"sym_reduce",
    //           "a":[{monomial}, ...],
    //           "vars":[...]}
    // Applies simplify_symcoef (Pi^(2k) -> (6 mzv_2)^k) and returns the
    // simplified SymCoef as a flat string.
    auto user_vars = json_str_array(body, "vars");
    std::string data_path = resolve_mzv_data_path(body);
    hf::MzvReductionTable table = hf::load_mzv_reductions(data_path);
    std::vector<std::string> vars = hf::build_mzv_var_list(table, user_vars);

    hf::PolyCtx ctx(vars);
    hf::SymCoef A = parse_sym_coef(ctx, body, "a");
    hf::SymCoef R = hf::simplify_symcoef(A, table);
    std::string result_str = sym_coef_to_mma_string(R);
    std::cout << emit_result_json("sym_reduce", vars, result_str) << "\n";
    return 0;
}

// ---------- Phase 7-i: AlgebraicLetterTable ----------

int handle_algebraic_letters_clear(const std::string& /*body*/) {
    hf::AlgebraicLetterTable::global().clear();
    std::cout << "{\"op\":\"algebraic_letters_clear\",\"size\":0}\n";
    return 0;
}

int handle_algebraic_letters_show(const std::string& /*body*/) {
    auto& T = hf::AlgebraicLetterTable::global();
    std::ostringstream o;
    o << "{\"op\":\"algebraic_letters_show\",\"entries\":[";
    bool first = true;
    for (long idx : T.indices()) {
        if (!first) o << ",";
        first = false;
        const auto& e = T.at(idx);
        o << "{\"idx\":" << e.idx
          << ",\"polynomial\":\"" << json_escape(e.polynomial.to_string()) << "\""
          << ",\"var_idx\":" << e.var_idx
          << ",\"sum\":\""        << json_escape(e.sum_value.to_string())     << "\""
          << ",\"product\":\""    << json_escape(e.product_value.to_string()) << "\""
          << ",\"discriminant\":\"" << json_escape(e.discriminant.to_string()) << "\""
          << ",\"wm\":\""    << json_escape(e.wm_name())    << "\""
          << ",\"wp\":\""    << json_escape(e.wp_name())    << "\""
          << ",\"wm_over_wp\":\"" << json_escape(e.wm_over_wp_name()) << "\""
          << "}";
    }
    o << "],\"size\":" << T.size() << "}\n";
    std::cout << o.str();
    return 0;
}

int handle_convert_to_hlog_reg_inf(const std::string& body) {
    // Request:  {"op":"convert_to_hlog_reg_inf",
    //            "expr":"<Mma-style expression>","vars":[...]}
    // Response: {"op":"convert_to_hlog_reg_inf",
    //            "result":[{"coef":"<rat>","key":[[l1,..],...]}, ...],
    //            "vars":[<augmented>]}
    //
    // Phase 3-b. Parses `expr`, dispatches through the 7-case Hlog
    // branch (plus minimal Plus/Times/Power/Leaf handling needed for
    // Case-4 re-entry). Emits a Regulator in HF's standard JSON shape.
    std::string expr_str = json_str_field(body, "expr");
    if (expr_str.empty()) {
        std::cerr << "convert_to_hlog_reg_inf: need \"expr\"\n"; return 1;
    }
    auto user_vars = json_str_array(body, "vars");
    try {
        auto parsed = hf::convert::parse_expression(expr_str, user_vars);
        hf::Regulator r = hf::convert::convert_to_hlog_reg_inf(
            parsed.expr, *parsed.ctx);
        std::ostringstream o;
        o << "{\"op\":\"convert_to_hlog_reg_inf\""
          << ",\"result\":" << emit_regulator(r)
          << ",\"vars\":[";
        for (size_t i = 0; i < parsed.augmented_vars.size(); ++i) {
            if (i) o << ",";
            o << "\"" << json_escape(parsed.augmented_vars[i]) << "\"";
        }
        o << "]}\n";
        std::cout << o.str();
    } catch (const hf::convert::ConvertFailed& e) {
        // The "error" field in a response triggers compare.py's
        // "backend errored" gate; use "reason" to preserve the
        // diagnostic string without short-circuiting the comparator.
        std::ostringstream o;
        o << "{\"op\":\"convert_to_hlog_reg_inf\""
          << ",\"failed\":true"
          << ",\"reason\":\"" << json_escape(e.what()) << "\""
          << ",\"vars\":[";
        for (size_t i = 0; i < user_vars.size(); ++i) {
            if (i) o << ",";
            o << "\"" << json_escape(user_vars[i]) << "\"";
        }
        o << "]}\n";
        std::cout << o.str();
    } catch (const hf::convert::ParseError& e) {
        std::ostringstream o;
        o << "{\"op\":\"convert_to_hlog_reg_inf\""
          << ",\"failed\":true"
          << ",\"reason\":\"" << json_escape(e.what()) << "\""
          << ",\"vars\":[";
        for (size_t i = 0; i < user_vars.size(); ++i) {
            if (i) o << ",";
            o << "\"" << json_escape(user_vars[i]) << "\"";
        }
        o << "]}\n";
        std::cout << o.str();
    }
    return 0;
}

int handle_fibration_basis(const std::string& body) {
    // Request:
    //   {"op":"fibration_basis",
    //    "wordlist":[{"coef":"<rat>","key":[[l1,..],...]}, ...],
    //    "vars_int":["x1","x2",...],        (reduction order; may be empty)
    //    "vars":[...]}                      (ambient PolyCtx)
    // Response:
    //   {"op":"fibration_basis",
    //    "vars_int":["x1",...],
    //    "terms":[{"key":[[w_1_letters], ...],"coef":"<rat>"}, ...],
    //    "vars":[<full ctx vars>]}
    //
    // Phase 6d-v-v. Each `terms[i].key` is a list of Words, one per
    // var in vars_int; the Word is the Hlog argument-word for that
    // var in the corresponding factor. Empty Word means "factor is 1".
    auto user_vars = json_str_array(body, "vars");
    auto vars_int  = json_str_array(body, "vars_int");
    if (user_vars.empty()) {
        for (const auto& v : vars_int) user_vars.push_back(v);
    }
    for (const auto& vi : vars_int) {
        bool present = false;
        for (const auto& v : user_vars) if (v == vi) { present = true; break; }
        if (!present) user_vars.push_back(vi);
    }
    if (user_vars.empty()) user_vars.push_back("x");

    std::string data_path = resolve_mzv_data_path(body);
    hf::MzvReductionTable table = hf::load_mzv_reductions(data_path);
    std::vector<std::string> vars = hf::build_mzv_var_list(table, user_vars);
    hf::PolyCtx ctx(vars);

    // Parse "wordlist":[{"coef":..,"key":[[..],..]},..] as a Regulator.
    hf::Regulator input;
    std::string inner = extract_top_array(body, "wordlist");
    if (!inner.empty()) {
        int depth = 0;
        size_t start = 0;
        std::vector<std::string> entries;
        for (size_t i = 0; i < inner.size(); ++i) {
            if (inner[i] == '{') {
                if (depth == 0) start = i;
                depth++;
            } else if (inner[i] == '}') {
                depth--;
                if (depth == 0) entries.push_back(
                    inner.substr(start, i - start + 1));
            }
        }
        for (const auto& e : entries) {
            std::string coef_s = json_str_field(e, "coef");
            hf::RegTerm t{hf::Rat::parse(ctx, coef_s),
                          parse_regkey_words(ctx, e)};
            input.push_back(std::move(t));
        }
    }

    std::vector<size_t> var_indices;
    for (const auto& vi : vars_int) {
        size_t idx = 0;
        for (; idx < vars.size(); ++idx) if (vars[idx] == vi) break;
        var_indices.push_back(idx);
    }

    // Phase 6d-v-v-ii: route through the SymCoef-valued core when the
    // request carries `"use_sym":true`. The Rat-only path throws on
    // Fragment-P2 residues; the sym path handles them natively, at a
    // slightly higher per-term cost in canonicalization.
    bool use_sym = false;
    {
        std::regex re("\"use_sym\"\\s*:\\s*(true|false)");
        std::smatch m;
        if (std::regex_search(body, m, re)) use_sym = (m[1] == "true");
    }

    if (use_sym) {
        try {
            hf::RegulatorSym input_sym;
            input_sym.reserve(input.size());
            for (const auto& t : input) {
                input_sym.push_back(
                    hf::RegTermSym{hf::SymCoef::from_rat(t.coef), t.key});
            }
            hf::FibrationBasisResultSym fbr =
                hf::fibration_basis_sym(ctx, input_sym, var_indices, table);
            std::ostringstream o;
            o << "{\"op\":\"fibration_basis\""
              << ",\"variant\":\"sym\""
              << ",\"vars_int\":[";
            for (size_t i = 0; i < vars_int.size(); ++i) {
                if (i) o << ",";
                o << "\"" << json_escape(vars_int[i]) << "\"";
            }
            o << "],\"terms\":[";
            for (size_t ti = 0; ti < fbr.terms.size(); ++ti) {
                if (ti) o << ",";
                const auto& [key, coef] = fbr.terms[ti];
                o << "{\"key\":[";
                for (size_t wi = 0; wi < key.size(); ++wi) {
                    if (wi) o << ",";
                    o << "[";
                    for (size_t li = 0; li < key[wi].size(); ++li) {
                        if (li) o << ",";
                        o << "\"" << json_escape(key[wi][li].to_string()) << "\"";
                    }
                    o << "]";
                }
                o << "],\"coef\":\"" << json_escape(coef.to_string()) << "\"}";
            }
            o << "],\"vars\":[";
            for (size_t i = 0; i < vars.size(); ++i) {
                if (i) o << ",";
                o << "\"" << json_escape(vars[i]) << "\"";
            }
            o << "]}\n";
            std::cout << o.str();
        } catch (const std::exception& e) {
            std::ostringstream o;
            o << "{\"op\":\"fibration_basis\""
              << ",\"variant\":\"sym\""
              << ",\"failed\":true,\"reason\":\""
              << json_escape(e.what()) << "\"}\n";
            std::cout << o.str();
        }
        return 0;
    }

    try {
        hf::FibrationBasisResult fbr =
            hf::fibration_basis(ctx, input, var_indices, table);

        std::ostringstream o;
        o << "{\"op\":\"fibration_basis\""
          << ",\"vars_int\":[";
        for (size_t i = 0; i < vars_int.size(); ++i) {
            if (i) o << ",";
            o << "\"" << json_escape(vars_int[i]) << "\"";
        }
        o << "],\"terms\":[";
        for (size_t ti = 0; ti < fbr.terms.size(); ++ti) {
            if (ti) o << ",";
            const auto& [key, coef] = fbr.terms[ti];
            o << "{\"key\":[";
            for (size_t wi = 0; wi < key.size(); ++wi) {
                if (wi) o << ",";
                o << "[";
                for (size_t li = 0; li < key[wi].size(); ++li) {
                    if (li) o << ",";
                    o << "\"" << json_escape(key[wi][li].to_string()) << "\"";
                }
                o << "]";
            }
            o << "],\"coef\":\"" << json_escape(coef.to_string()) << "\"}";
        }
        o << "],\"vars\":[";
        for (size_t i = 0; i < vars.size(); ++i) {
            if (i) o << ",";
            o << "\"" << json_escape(vars[i]) << "\"";
        }
        o << "]}\n";
        std::cout << o.str();
    } catch (const std::exception& e) {
        std::ostringstream o;
        o << "{\"op\":\"fibration_basis\""
          << ",\"failed\":true"
          << ",\"reason\":\"" << json_escape(e.what()) << "\""
          << ",\"vars\":[";
        for (size_t i = 0; i < vars.size(); ++i) {
            if (i) o << ",";
            o << "\"" << json_escape(vars[i]) << "\"";
        }
        o << "]}\n";
        std::cout << o.str();
    }
    return 0;
}

int handle_parse_expr(const std::string& body) {
    // Request: {"op":"parse_expr","expr":"<Mma-style string>","vars":[...]}
    // Response: {"op":"parse_expr","canonical":"<canonical form>","vars":[...]}
    //
    // Phase 3-a: round-trip diagnostic. Parses the Mma-style input
    // into the Expr AST and emits the canonical string form. The
    // Mma cross-validation runner emits the same canonical form by
    // walking its own ExpressionQ representation; the comparator
    // byte-compares the two.
    std::string expr_str = json_str_field(body, "expr");
    if (expr_str.empty()) {
        std::cerr << "parse_expr: need \"expr\"\n"; return 1;
    }
    auto user_vars = json_str_array(body, "vars");
    try {
        auto result = hf::convert::parse_expression(expr_str, user_vars);
        std::ostringstream o;
        o << "{\"op\":\"parse_expr\""
          << ",\"canonical\":\"" << json_escape(result.expr.to_canonical_string()) << "\""
          << ",\"vars\":[";
        for (size_t i = 0; i < result.augmented_vars.size(); ++i) {
            if (i) o << ",";
            o << "\"" << json_escape(result.augmented_vars[i]) << "\"";
        }
        o << "]}\n";
        std::cout << o.str();
    } catch (const hf::convert::ParseError& e) {
        std::ostringstream o;
        o << "{\"op\":\"parse_expr\""
          << ",\"error\":\"" << json_escape(e.what()) << "\""
          << ",\"vars\":[";
        for (size_t i = 0; i < user_vars.size(); ++i) {
            if (i) o << ",";
            o << "\"" << json_escape(user_vars[i]) << "\"";
        }
        o << "]}\n";
        std::cout << o.str();
    }
    return 0;
}

int handle_combine_wm_wp_ratios(const std::string& body) {
    // Request: {"op":"combine_wm_wp_ratios",
    //           "expr":"<rat>",
    //           "allocations":[{"polynomial":"...", "var":"<name>"}, ...],
    //           "vars":[...]}
    std::string expr_str = json_str_field(body, "expr");
    if (expr_str.empty()) {
        std::cerr << "combine_wm_wp_ratios: need \"expr\"\n"; return 1;
    }
    auto user_vars = json_str_array(body, "vars");
    if (user_vars.empty()) user_vars.push_back("x");
    std::vector<std::string> vars =
        hf::build_algebraic_letter_var_list(user_vars);
    hf::PolyCtx ctx(vars);

    hf::AlgebraicLetterTable::global().clear();
    std::string allocs_inner = extract_top_array(body, "allocations");
    if (!allocs_inner.empty()) {
        int depth = 0;
        size_t start = 0;
        std::vector<std::string> entries;
        for (size_t i = 0; i < allocs_inner.size(); ++i) {
            if (allocs_inner[i] == '{') {
                if (depth == 0) start = i;
                depth++;
            } else if (allocs_inner[i] == '}') {
                depth--;
                if (depth == 0)
                    entries.push_back(allocs_inner.substr(start, i - start + 1));
            }
        }
        for (const auto& e : entries) {
            std::string poly_s = json_str_field(e, "polynomial");
            std::string var_s  = json_str_field(e, "var");
            long var_idx = -1;
            for (size_t i = 0; i < vars.size(); ++i) {
                if (vars[i] == var_s) { var_idx = static_cast<long>(i); break; }
            }
            if (var_idx < 0) {
                std::cerr << "combine_wm_wp_ratios: allocation var '"
                          << var_s << "' not in PolyCtx\n";
                return 1;
            }
            hf::Poly p(ctx, poly_s);
            hf::AlgebraicLetterTable::global().allocate(p, var_idx);
        }
    }

    hf::Rat input = hf::Rat::parse(ctx, expr_str);
    hf::Rat reduced = hf::combine_wm_wp_ratios(input);
    std::cout << emit_result_json("combine_wm_wp_ratios", vars,
                                   reduced.to_string()) << "\n";
    return 0;
}

int handle_back_substitute(const std::string& body) {
    // Request: {"op":"back_substitute",
    //           "expr":"<rat over Wm_i, Wp_i, plus user vars>",
    //           "allocations":[{"polynomial":"...", "var":"<name>"}, ...],
    //           "vars":[...]}
    // Same allocation-replay pattern as simplify_with_vieta.
    std::string expr_str = json_str_field(body, "expr");
    if (expr_str.empty()) {
        std::cerr << "back_substitute: need \"expr\"\n"; return 1;
    }
    auto user_vars = json_str_array(body, "vars");
    if (user_vars.empty()) user_vars.push_back("x");
    std::vector<std::string> vars =
        hf::build_algebraic_letter_var_list(user_vars);
    hf::PolyCtx ctx(vars);

    hf::AlgebraicLetterTable::global().clear();
    std::string allocs_inner = extract_top_array(body, "allocations");
    if (!allocs_inner.empty()) {
        int depth = 0;
        size_t start = 0;
        std::vector<std::string> entries;
        for (size_t i = 0; i < allocs_inner.size(); ++i) {
            if (allocs_inner[i] == '{') {
                if (depth == 0) start = i;
                depth++;
            } else if (allocs_inner[i] == '}') {
                depth--;
                if (depth == 0)
                    entries.push_back(allocs_inner.substr(start, i - start + 1));
            }
        }
        for (const auto& e : entries) {
            std::string poly_s = json_str_field(e, "polynomial");
            std::string var_s  = json_str_field(e, "var");
            long var_idx = -1;
            for (size_t i = 0; i < vars.size(); ++i) {
                if (vars[i] == var_s) { var_idx = static_cast<long>(i); break; }
            }
            if (var_idx < 0) {
                std::cerr << "back_substitute: allocation var '" << var_s
                          << "' not in PolyCtx\n";
                return 1;
            }
            hf::Poly p(ctx, poly_s);
            hf::AlgebraicLetterTable::global().allocate(p, var_idx);
        }
    }

    hf::Rat input = hf::Rat::parse(ctx, expr_str);
    hf::Rat reduced = hf::back_substitute(input);
    std::cout << emit_result_json("back_substitute", vars,
                                   reduced.to_string()) << "\n";
    return 0;
}

int handle_simplify_with_vieta(const std::string& body) {
    // Request: {"op":"simplify_with_vieta",
    //           "expr":"<rat over Wm_i, Wp_i, plus user vars>",
    //           "allocations":[{"polynomial":"...", "var":"<name>"}, ...],
    //           "vars":[...]}
    //
    // The `allocations` array reproduces the AlgebraicLetterTable state
    // for this CLI invocation (since each invocation starts with a
    // fresh table). Each entry's index is the 1-based position in the
    // array; the user must reference Wm_<i> / Wp_<i> in `expr`
    // accordingly.
    std::string expr_str = json_str_field(body, "expr");
    if (expr_str.empty()) {
        std::cerr << "simplify_with_vieta: need \"expr\"\n"; return 1;
    }
    auto user_vars = json_str_array(body, "vars");
    if (user_vars.empty()) user_vars.push_back("x");
    std::vector<std::string> vars =
        hf::build_algebraic_letter_var_list(user_vars);
    hf::PolyCtx ctx(vars);

    // Reset the table and apply the allocations in order.
    hf::AlgebraicLetterTable::global().clear();
    std::string allocs_inner = extract_top_array(body, "allocations");
    if (!allocs_inner.empty()) {
        int depth = 0;
        size_t start = 0;
        std::vector<std::string> entries;
        for (size_t i = 0; i < allocs_inner.size(); ++i) {
            if (allocs_inner[i] == '{') {
                if (depth == 0) start = i;
                depth++;
            } else if (allocs_inner[i] == '}') {
                depth--;
                if (depth == 0)
                    entries.push_back(allocs_inner.substr(start, i - start + 1));
            }
        }
        for (const auto& e : entries) {
            std::string poly_s = json_str_field(e, "polynomial");
            std::string var_s  = json_str_field(e, "var");
            long var_idx = -1;
            for (size_t i = 0; i < vars.size(); ++i) {
                if (vars[i] == var_s) { var_idx = static_cast<long>(i); break; }
            }
            if (var_idx < 0) {
                std::cerr << "simplify_with_vieta: allocation var '"
                          << var_s << "' not in PolyCtx\n";
                return 1;
            }
            hf::Poly p(ctx, poly_s);
            hf::AlgebraicLetterTable::global().allocate(p, var_idx);
        }
    }

    hf::Rat input = hf::Rat::parse(ctx, expr_str);
    hf::Rat reduced = hf::simplify_with_vieta(input);
    std::cout << emit_result_json("simplify_with_vieta", vars,
                                   reduced.to_string()) << "\n";
    return 0;
}

int handle_algebraic_letters_allocate(const std::string& body) {
    // Request: {"op":"algebraic_letters_allocate",
    //           "polynomial":"<expr>", "var":"<name>", "vars":[...]}
    // Returns: {"idx":<allocated 1-based index>}.
    std::string poly_str = json_str_field(body, "polynomial");
    std::string var_name = json_str_field(body, "var");
    if (poly_str.empty() || var_name.empty()) {
        std::cerr << "algebraic_letters_allocate: need \"polynomial\" and \"var\"\n";
        return 1;
    }
    auto user_vars = json_str_array(body, "vars");
    bool present = false;
    for (const auto& v : user_vars) if (v == var_name) { present = true; break; }
    if (!present) user_vars.push_back(var_name);
    std::vector<std::string> vars =
        hf::build_algebraic_letter_var_list(user_vars);

    hf::PolyCtx ctx(vars);
    hf::Poly p(ctx, poly_str);
    long var_idx = -1;
    for (size_t i = 0; i < vars.size(); ++i) {
        if (vars[i] == var_name) { var_idx = static_cast<long>(i); break; }
    }
    if (var_idx < 0) {
        std::cerr << "algebraic_letters_allocate: var not in PolyCtx\n";
        return 1;
    }
    long idx = hf::AlgebraicLetterTable::global().allocate(p, var_idx);
    const auto& e = hf::AlgebraicLetterTable::global().at(idx);
    std::cout << "{\"op\":\"algebraic_letters_allocate\""
              << ",\"idx\":" << idx
              << ",\"sum\":\""          << json_escape(e.sum_value.to_string())     << "\""
              << ",\"product\":\""      << json_escape(e.product_value.to_string()) << "\""
              << ",\"discriminant\":\"" << json_escape(e.discriminant.to_string())  << "\""
              << ",\"wm\":\""           << json_escape(e.wm_name())                 << "\""
              << ",\"wp\":\""           << json_escape(e.wp_name())                 << "\""
              << "}\n";
    return 0;
}

int handle_series_expansion(const std::string& body) {
    // Request: {"op":"series_expansion","f":<str>,"var":<str>,"max_order":<int>,"vars":[...]}
    std::string f   = json_str_field(body, "f");
    std::string var = json_str_field(body, "var");
    if (f.empty() || var.empty()) {
        std::cerr << "series_expansion: need \"f\" and \"var\"\n"; return 1;
    }
    std::regex re("\"max_order\"\\s*:\\s*(-?\\d+)");
    std::smatch m;
    if (!std::regex_search(body, m, re)) {
        std::cerr << "series_expansion: need \"max_order\"\n"; return 1;
    }
    long max_order = std::stol(m[1]);
    auto vars = json_str_array(body, "vars");
    if (vars.empty()) vars = autoscan_vars({&f});
    bool present = false;
    for (const auto& v : vars) if (v == var) { present = true; break; }
    if (!present) vars.push_back(var);

    hf::PolyCtx ctx(vars);
    hf::Rat rf = hf::Rat::parse(ctx, f);
    size_t idx = 0;
    for (; idx < vars.size(); ++idx) if (vars[idx] == var) break;

    hf::Rat r = hf::series_expansion(ctx, rf, idx, max_order);
    std::cout << emit_result_json("series_expansion", vars, r.to_string()) << "\n";
    return 0;
}

int handle_integrate_ii(const std::string& body) {
    // Request: {"op":"integrate_ii","wl":[...],"var":<str>,"vars":[...],
    //           "algebraic_letters":bool (optional)}
    std::string var = json_str_field(body, "var");
    if (var.empty()) { std::cerr << "integrate_ii: need \"var\"\n"; return 1; }
    auto vars = json_str_array(body, "vars");
    if (vars.empty()) vars.push_back(var);
    bool present = false;
    for (const auto& v : vars) if (v == var) { present = true; break; }
    if (!present) vars.push_back(var);
    bool introduce_al = false;
    {
        std::regex re("\"algebraic_letters\"\\s*:\\s*(true|false)");
        std::smatch m;
        if (std::regex_search(body, m, re)) {
            introduce_al = (m[1] == "true");
        }
    }
    std::vector<std::string> ctx_vars = introduce_al
        ? hf::build_algebraic_letter_var_list(vars)
        : vars;
    hf::PolyCtx ctx(ctx_vars);
    hf::Wordlist wl_in = parse_wordlist(ctx, body, "wl");
    size_t idx = 0;
    for (; idx < ctx_vars.size(); ++idx) if (ctx_vars[idx] == var) break;
    try {
        // Iter-52 C0c.1 Increment β: caller-side fresh transient ZWTable.
        auto _lf_zw = std::make_shared<hf::ZWTable>(ctx);
        hf::Wordlist out = hf::integrate_ii(ctx, wl_in, idx, _lf_zw, introduce_al);
        std::ostringstream o;
        o << "{\"op\":\"integrate_ii\",\"result\":" << emit_wordlist(out)
          << ",\"vars\":[";
        for (size_t i = 0; i < vars.size(); ++i) {
            if (i) o << ",";
            o << "\"" << json_escape(vars[i]) << "\"";
        }
        o << "]}\n";
        std::cout << o.str();
    } catch (const hf::IntegrateIIFailed&) {
        std::cout << "{\"op\":\"integrate_ii\",\"failed\":true,\"vars\":[";
        for (size_t i = 0; i < vars.size(); ++i) {
            if (i) std::cout << ",";
            std::cout << "\"" << json_escape(vars[i]) << "\"";
        }
        std::cout << "]}\n";
    }
    return 0;
}

int handle_transform_shuffle(const std::string& body) {
    // Request: {"op":"transform_shuffle","wordlist":[[l1,..],...],"var":<str>,"vars":[...],
    //           "algebraic_letters":bool (optional)}
    std::string var = json_str_field(body, "var");
    auto vars = json_str_array(body, "vars");
    if (vars.empty()) vars.push_back(var.empty() ? "x" : var);
    if (!var.empty()) {
        bool present = false;
        for (const auto& v : vars) if (v == var) { present = true; break; }
        if (!present) vars.push_back(var);
    }
    // Bug #6: extend ctx with MZV basis (see handle_transform_word).
    std::string data_path = resolve_mzv_data_path(body);
    hf::MzvReductionTable table = hf::load_mzv_reductions(data_path);
    vars = hf::build_mzv_var_list(table, vars);
    bool introduce_al = false;
    {
        std::regex re("\"algebraic_letters\"\\s*:\\s*(true|false)");
        std::smatch m;
        if (std::regex_search(body, m, re)) {
            introduce_al = (m[1] == "true");
        }
    }
    std::vector<std::string> ctx_vars = introduce_al
        ? hf::build_algebraic_letter_var_list(vars)
        : vars;
    hf::PolyCtx ctx(ctx_vars);
    // Parse "wordlist": [[l1,...], [l1,...], ...]  top-level array of arrays.
    std::vector<hf::Word> words;
    {
        std::string inner = extract_top_array(body, "wordlist");
        if (!inner.empty()) {
            int depth = 0;
            size_t wstart = 0;
            for (size_t i = 0; i < inner.size(); ++i) {
                if (inner[i] == '[') {
                    if (depth == 0) wstart = i;
                    depth++;
                } else if (inner[i] == ']') {
                    depth--;
                    if (depth == 0) {
                        std::string word_json = inner.substr(wstart, i - wstart + 1);
                        std::string wrapped = "{\"w\":" + word_json + "}";
                        auto letters = json_str_array(wrapped, "w");
                        words.push_back(parse_word(ctx, letters));
                    }
                }
            }
        }
    }
    size_t idx = 0;
    if (!var.empty()) {
        for (; idx < ctx_vars.size(); ++idx) if (ctx_vars[idx] == var) break;
    }
    try {
        // Iter-52 C0c.1 Increment β: caller-side fresh transient ZWTable.
        auto _lf_zw = std::make_shared<hf::ZWTable>(ctx);
        hf::TransformResultSym tr =
            hf::transform_shuffle(ctx, words, idx, table, _lf_zw, introduce_al);
        std::ostringstream o;
        o << "{\"op\":\"transform_shuffle\",\"result\":" << emit_transform_result_sym(tr)
          << ",\"vars\":[";
        for (size_t i = 0; i < vars.size(); ++i) {
            if (i) o << ",";
            o << "\"" << json_escape(vars[i]) << "\"";
        }
        o << "]}\n";
        std::cout << o.str();
    } catch (const hf::TransformFailed&) {
        std::cout << "{\"op\":\"transform_shuffle\",\"failed\":true,\"vars\":[";
        for (size_t i = 0; i < vars.size(); ++i) {
            if (i) std::cout << ",";
            std::cout << "\"" << json_escape(vars[i]) << "\"";
        }
        std::cout << "]}\n";
    }
    return 0;
}

int handle_reglim_word(const std::string& body) {
    // Request: {"op":"reglim_word","word":[...],"var":<str>,"vars":[...]}
    auto letters = json_str_array(body, "word");
    std::string var = json_str_field(body, "var");
    auto vars = json_str_array(body, "vars");
    if (vars.empty()) vars.push_back(var.empty() ? "x" : var);
    if (!var.empty()) {
        bool present = false;
        for (const auto& v : vars) if (v == var) { present = true; break; }
        if (!present) vars.push_back(var);
    }
    // Bug #6: extend ctx with MZV basis (see handle_transform_word).
    std::string data_path = resolve_mzv_data_path(body);
    hf::MzvReductionTable table = hf::load_mzv_reductions(data_path);
    vars = hf::build_mzv_var_list(table, vars);
    hf::PolyCtx ctx(vars);
    hf::Word word = parse_word(ctx, letters);
    size_t idx = 0;
    if (!var.empty()) {
        for (; idx < vars.size(); ++idx) if (vars[idx] == var) break;
    }
    // Iter-66 C0c.1 path 1a sub-iter 3 (ABI cascade per iter-63 audit
    // MEMO §5.3): `reglim_word`'s signature widened to take a
    // `std::shared_ptr<ZWTable> zw_tab`. Allocate a fresh transient
    // here, mirroring the `_lf_zw` pattern in `handle_transform_word`
    // (main.cpp:1751-1755) and other CLI entrypoints.
    auto _lf_zw = std::make_shared<hf::ZWTable>(ctx);
    hf::RegulatorSym out = hf::reglim_word(ctx, word, idx, table, _lf_zw);
    std::ostringstream o;
    o << "{\"op\":\"reglim_word\",\"result\":" << emit_regulator_sym(out)
      << ",\"vars\":[";
    for (size_t i = 0; i < vars.size(); ++i) {
        if (i) o << ",";
        o << "\"" << json_escape(vars[i]) << "\"";
    }
    o << "]}\n";
    std::cout << o.str();
    return 0;
}

// ---------- Phase 5b: DifferentiateWordlist ----------

int handle_differentiate_wordlist(const std::string& body) {
    // Request: {"op":"differentiate_wordlist","wl":[...],"var":<str>,"vars":[...]}
    std::string var = json_str_field(body, "var");
    if (var.empty()) {
        std::cerr << "differentiate_wordlist: need \"var\"\n"; return 1;
    }
    auto vars = json_str_array(body, "vars");
    if (vars.empty()) vars.push_back(var);
    bool present = false;
    for (const auto& v : vars) if (v == var) { present = true; break; }
    if (!present) vars.push_back(var);

    hf::PolyCtx ctx(vars);
    hf::Wordlist wl = parse_wordlist(ctx, body, "wl");
    size_t idx = 0;
    for (; idx < vars.size(); ++idx) if (vars[idx] == var) break;

    // Fast path: an empty input wordlist needs no Rat construction.
    hf::Wordlist out = hf::differentiate_wordlist(wl, idx);

    std::ostringstream o;
    o << "{\"op\":\"differentiate_wordlist\",\"result\":" << emit_wordlist(out)
      << ",\"vars\":[";
    for (size_t i = 0; i < vars.size(); ++i) {
        if (i) o << ",";
        o << "\"" << json_escape(vars[i]) << "\"";
    }
    o << "]}\n";
    std::cout << o.str();
    return 0;
}

void usage() {
    std::cerr <<
        "HyperFLINT CLI (Phase 1a)\n"
        "\n"
        "  hyperflint eval-json < request.json\n"
        "     Read one JSON request from stdin, write JSON to stdout.\n"
        "\n"
        "Supported ops:\n"
        "  factor   {\"op\":\"factor\",\"expr\":<str>,\"vars\":[...]}\n"
        "  add/sub/mul  {\"op\":<op>,\"a\":<str>,\"b\":<str>,\"vars\":[...]}\n"
        "  neg      {\"op\":\"neg\",\"a\":<str>,\"vars\":[...]}\n"
        "  pow      {\"op\":\"pow\",\"a\":<str>,\"n\":<int>,\"vars\":[...]}\n"
        "\n"
        "Pretty-printing helpers:\n"
        "  hyperflint factor <expr>\n";
}

}  // namespace

// Phase-1.2 retrofit (2026-05-05): redirect GMP+FLINT allocations to
// mimalloc. Must be called before any GMP/FLINT memory is allocated;
// audit note: confirmed no namespace-scope statics
// touch GMP/FLINT, so the very first line of main() is safe. See
// bridge/cli/gmp_mimalloc_init.cpp for details.
extern "C" void hf_init_mimalloc_for_gmp_flint(void);

// Phase 4 §A.4 FOLD-I (iter-17): mimalloc weak symbols at file scope,
// mirroring HyperFLINT/src/integrator/hyper_int.cpp:35-39. The CLI links
// mimalloc statically via gmp_mimalloc_init.cpp; these declarations let
// the `mi-peak-test` subcommand probe peak_rss without forcing a hard
// dependency at link time.
extern "C" {
__attribute__((weak)) void mi_collect(int /*force*/);
__attribute__((weak)) void mi_process_info(size_t*, size_t*, size_t*,
                                             size_t*, size_t*,
                                             size_t*, size_t*, size_t*);
__attribute__((weak)) void* mi_malloc(size_t size);
__attribute__((weak)) void  mi_free(void* p);
}

int main(int argc, char** argv) {
    hf_init_mimalloc_for_gmp_flint();
    // HF FF Phase 5 §A.1 iter-49 REQ-1 BINDING: probe init MUST follow the
    // Phase 0.5 retrofit at gmp_mimalloc_init.cpp:145. The probe snapshots
    // the GMP function pointers via mp_get_memory_functions THEN re-registers
    // its wrapper layer which delegates through the saved pointers. Calling
    // this BEFORE hf_init_mimalloc_for_gmp_flint() would capture the libsystem
    // defaults (null/malloc) and break the M10 attribution lineage; calling
    // AFTER any FLINT/GMP allocation would leave a portion of the M10
    // denominator unlabelled. Default-OFF unless HF_DAG_HASHCONS_PROBE=1.
    hyperflint::hf_probe_init();
    // HF FF Phase 5 REC-1 (iter-83) REQ-2 BINDING chain-ordering: REC-1 init
    // MUST follow hf_probe_init() so the GMP-allocator chain is
    // REC-1 -> dag_hashcons_probe -> Phase 0.5 retrofit -> mimalloc. Reversed
    // ordering would double-classify GMP allocations against the iter-49
    // labelled-bytes counter. Default-OFF unless HF_REC1_TRACK_MPZ_POOL=1.
    hyperflint::hf_rec1_init();
    if (argc < 2) { usage(); return 1; }
    std::string cmd = argv[1];

    if (cmd == "--help" || cmd == "-h") { usage(); return 0; }

    // iter-29 Track 5.1: --version emits the build-variant stamp
    // compiled in by CMakePresets.json (HF_BUILD_VARIANT_STR). Allows
    // post-hoc attribution of a binary to the preset that built it.
    // iter-57 §T10: also emits HF_VERSION (MAJOR.MINOR.PATCH.BUILD)
    // from cmake's HF_VERSION cache var, so external tooling can
    // attribute a binary to a semver-tagged release without parsing
    // the eval-json envelope.
    if (cmd == "--version" || cmd == "-v") {
#ifndef HF_BUILD_VARIANT_STR
#define HF_BUILD_VARIANT_STR "unknown"
#endif
#ifndef HF_VERSION_STRING
#define HF_VERSION_STRING "unknown"
#endif
        std::cout << "HF_VERSION: " HF_VERSION_STRING "\n"
                  << "HF_BUILD_VARIANT: " HF_BUILD_VARIANT_STR "\n";
        return 0;
    }

    if (cmd == "eval-json") {
        std::string body, line;
        while (std::getline(std::cin, line)) body += line + "\n";
        std::string op = json_str_field(body, "op");
        if (op.empty()) {
            std::cerr << "eval-json: missing \"op\"\n";
            return 1;
        }
        try {
            if (op == "factor")                 return handle_factor(body);
            if (op == "discriminant")           return handle_discriminant(body);
            if (op == "find_lr_orders")         return handle_find_lr_orders(body);
            if (op == "find_lr_orders_scan")    return handle_find_lr_orders_scan(body);
            if (op == "factor_table")          return handle_factor_table(body);
            if (op == "add" || op == "sub" ||
                op == "mul")                    return handle_binary(op, body);
            if (op == "neg")                    return handle_neg(body);
            if (op == "pow")                    return handle_pow(body);
            if (op == "derivative")             return handle_derivative(body);
            if (op == "eval")                   return handle_eval(body);
            if (op == "subst")                  return handle_subst(body);
            if (op == "divexact")               return handle_divexact(body);
            if (op == "gcd")                    return handle_gcd(body);
            if (op == "resultant")              return handle_resultant(body);
            if (op == "rat_add" || op == "rat_sub" ||
                op == "rat_mul" || op == "rat_div")
                                                return handle_rat_binary(op, body);
            if (op == "rat_sum")                return handle_rat_sum(body);
            if (op == "linear_factors")         return handle_linear_factors(body);
            if (op == "pole_degree")            return handle_pole_degree(body);
            if (op == "rat_residue")            return handle_rat_residue(body);
            if (op == "partial_fractions")      return handle_partial_fractions(body);
            if (op == "shuffle_words")          return handle_shuffle_words(body);
            if (op == "shuffle_product")        return handle_shuffle_product(body);
            if (op == "concat_mul")             return handle_concat_mul(body);
            if (op == "collect_words")          return handle_collect_words(body);
            if (op == "convert_zero_one")       return handle_convert_zero_one(body);
            if (op == "convert_1inf_to_01")     return handle_convert_1inf_to_01(body);
            if (op == "convert_ab_to_zero_inf") return handle_convert_ab_to_zero_inf(body);
            if (op == "diff_hlog")              return handle_diff_hlog(body);
            if (op == "diff_mpl")               return handle_diff_mpl(body);
            if (op == "expand_zero_word")       return handle_expand_zero_word(body);
            if (op == "expand_inf_word")        return handle_expand_inf_word(body);
            if (op == "mpl_sum")                return handle_mpl_sum(body);
            if (op == "hlog_zero_expand")       return handle_hlog_zero_expand(body);
            if (op == "hlog_series")            return handle_hlog_series(body);
            if (op == "mpl_series")             return handle_mpl_series(body);
            if (op == "regzero_word")           return handle_regzero_word(body);
            if (op == "reg0")                   return handle_reg0(body);
            if (op == "reg_head")               return handle_reg_head(body);
            if (op == "reg_tail")               return handle_reg_tail(body);
            if (op == "differentiate_wordlist") return handle_differentiate_wordlist(body);
            if (op == "shuffle_symbolic")       return handle_shuffle_symbolic(body);
            if (op == "reglim_word")            return handle_reglim_word(body);
            if (op == "transform_word")         return handle_transform_word(body);
            if (op == "transform_shuffle")      return handle_transform_shuffle(body);
            if (op == "integrate_ii")           return handle_integrate_ii(body);
            if (op == "series_expansion")       return handle_series_expansion(body);
            if (op == "integration_step")       return handle_integration_step(body);
            if (op == "hyperflint")             return handle_hyperflint(body);
            if (op == "apply_mzv_reductions")   return handle_apply_mzv_reductions(body);
            if (op == "zero_one_period")        return handle_zero_one_period(body);
            if (op == "zero_inf_period")        return handle_zero_inf_period(body);
            if (op == "break_up_contour")       return handle_break_up_contour(body);
            if (op == "break_up_contour_sym")   return handle_break_up_contour_sym(body);
            if (op == "evaluate_periods")       return handle_evaluate_periods(body);
            if (op == "test_zero_function")     return handle_test_zero_function(body);
            if (op == "sym_arith")              return handle_sym_arith(body);
            if (op == "sym_reduce")             return handle_sym_reduce(body);
            if (op == "algebraic_letters_clear")
                return handle_algebraic_letters_clear(body);
            if (op == "algebraic_letters_show")
                return handle_algebraic_letters_show(body);
            if (op == "algebraic_letters_allocate")
                return handle_algebraic_letters_allocate(body);
            if (op == "simplify_with_vieta")
                return handle_simplify_with_vieta(body);
            if (op == "back_substitute")
                return handle_back_substitute(body);
            if (op == "combine_wm_wp_ratios")
                return handle_combine_wm_wp_ratios(body);
            if (op == "parse_expr")             return handle_parse_expr(body);
            if (op == "convert_to_hlog_reg_inf")
                return handle_convert_to_hlog_reg_inf(body);
            if (op == "fibration_basis")       return handle_fibration_basis(body);
            std::cout << "{\"op\":\"" << op
                      << "\",\"error\":\"unknown op\"}\n";
            return 2;
        } catch (const std::exception& e) {
            std::cout << "{\"op\":\"" << op
                      << "\",\"error\":\"" << json_escape(e.what())
                      << "\"}\n";
            return 2;
        }
    }

    // Pretty-print compat with Phase 0:  hyperflint factor <expr>
    if (cmd == "factor") {
        if (argc < 3) { std::cerr << "factor: expected <expr>\n"; return 1; }
        std::string expr = argv[2];
        std::vector<std::string> vars;
        for (int i = 3; i + 1 < argc; i += 2) {
            if (std::string(argv[i]) == "--vars") {
                std::string v = argv[i + 1];
                std::string cur;
                for (char c : v) {
                    if (c == ',') { if (!cur.empty()) vars.push_back(cur); cur.clear(); }
                    else cur += c;
                }
                if (!cur.empty()) vars.push_back(cur);
            }
        }
        if (vars.empty()) vars = autoscan_vars({&expr});
        if (vars.empty()) vars.push_back("x");
        try {
            hf::PolyCtx ctx(vars);
            hf::Poly p(ctx, expr);
            auto fac = hf::factor(p);
            std::cout << fac.constant;
            for (const auto& [base, exp] : fac.factors) {
                std::cout << " * (" << base << ")";
                if (exp != 1) std::cout << "^" << exp;
            }
            std::cout << "\n";
            return 0;
        } catch (const std::exception& e) {
            std::cerr << "error: " << e.what() << "\n";
            return 2;
        }
    }

    // Phase 4 §A.4 FOLD-I (iter-17): mimalloc peak_rss reliability smoke
    // test. Allocates a large block via malloc (mimalloc-overridden), frees
    // it, calls mi_collect(force=1), then reads mi_process_info to check
    // whether peak_rss reflects the historical peak ("lifetime peak") or
    // the post-free residency ("since-last-reset"). The Phase 4 §A.4
    // per-step probe only emits at post_collect, so we need to know which
    // mode mimalloc uses to interpret the trace correctly.
    //
    // Usage: hyperflint mi-peak-test [size_mb]   (default 512 MB)
    if (cmd == "mi-peak-test") {
        if (!mi_collect || !mi_process_info || !mi_malloc || !mi_free) {
            std::cerr << "mi-peak-test: mimalloc weak symbols unresolved "
                      << "(mi_collect="     << (void*)mi_collect
                      << " mi_process_info=" << (void*)mi_process_info
                      << " mi_malloc="       << (void*)mi_malloc
                      << " mi_free="         << (void*)mi_free
                      << ")\n";
            return 2;
        }
        size_t mb = 512;
        if (argc >= 3) {
            try { mb = static_cast<size_t>(std::stoul(argv[2])); }
            catch (...) { std::cerr << "mi-peak-test: bad size_mb\n"; return 1; }
        }
        const size_t N = mb * size_t(1) << 20;  // bytes

        auto sample = [&](const char* phase) {
            size_t elapsed = 0, user = 0, system = 0;
            size_t cur_rss = 0, peak_rss = 0;
            size_t cur_commit = 0, peak_commit = 0, page_faults = 0;
            mi_process_info(&elapsed, &user, &system,
                            &cur_rss, &peak_rss,
                            &cur_commit, &peak_commit, &page_faults);
            std::fprintf(stderr,
                "[mi-peak-test] phase=%s rss=%.1fMB peak_rss=%.1fMB "
                "commit=%.1fMB peak_commit=%.1fMB\n",
                phase,
                cur_rss / 1048576.0, peak_rss / 1048576.0,
                cur_commit / 1048576.0, peak_commit / 1048576.0);
            std::fflush(stderr);
            return std::pair<size_t, size_t>{cur_rss, peak_rss};
        };

        auto [rss_pre,   peak_pre]   = sample("pre_alloc");
        char* buf = static_cast<char*>(mi_malloc(N));
        if (!buf) { std::cerr << "mi-peak-test: mi_malloc failed\n"; return 2; }
        // Touch every page so RSS actually grows (lazy commit otherwise).
        for (size_t i = 0; i < N; i += 4096) buf[i] = static_cast<char>(i);
        auto [rss_after_touch, peak_after_touch] = sample("post_alloc_touch");
        mi_free(buf);
        auto [rss_after_free, peak_after_free] = sample("post_free");
        mi_collect(/*force=*/1);
        auto [rss_after_collect, peak_after_collect] = sample("post_collect");

        // Verdict:
        //   "lifetime"      peak_rss after free+collect ≈ peak after touch
        //   "since-reset"   peak_rss after free+collect << peak after touch
        const double ratio = (peak_after_touch > 0)
            ? double(peak_after_collect) / double(peak_after_touch) : 0.0;
        const bool is_lifetime = ratio > 0.8;
        std::fprintf(stderr,
            "[mi-peak-test] verdict=%s ratio=%.3f size_mb_requested=%zu\n",
            is_lifetime ? "LIFETIME_PEAK" : "SINCE_LAST_RESET",
            ratio, mb);
        std::cout << "{\"verdict\":\""
                  << (is_lifetime ? "LIFETIME_PEAK" : "SINCE_LAST_RESET")
                  << "\",\"size_mb_requested\":" << mb
                  << ",\"peak_rss_pre\":" << peak_pre
                  << ",\"peak_rss_post_touch\":" << peak_after_touch
                  << ",\"peak_rss_post_free\":" << peak_after_free
                  << ",\"peak_rss_post_collect\":" << peak_after_collect
                  << ",\"ratio\":" << ratio << "}\n";
        return 0;
    }

    std::cerr << "unknown command: " << cmd << "\n";
    usage();
    return 1;
}
