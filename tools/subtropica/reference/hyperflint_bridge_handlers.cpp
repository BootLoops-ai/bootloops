// Phase γ.1: transport-neutral JSON handlers for HyperFLINT ops.
//
// Currently exposed:
//   - find_lr_orders
//
// Both the CLI (`hyperflint eval-json`) and the LibraryLink shared
// library call into these.  The CLI wrapper prints the returned string
// to stdout + newline; the LibraryLink wrapper returns it to Mma via
// MArgument_setUTF8String.  Errors travel in the returned JSON (`"error"`
// field), never stderr; neither transport aborts.

#include "hyperflint/bridge/handlers.hpp"
#include "hyperflint/bridge/env_flags.hpp"  // iter-94 Track-OMP bridge portion: HF_FLAG_MAX_THREADS_PER_CALL (NEW first bridge-domain env_flags header; §5.1 rule-1 BINDING placement)

#include "hyperflint/c_abi.h"  // HF_SCHEMA_VERSION SSOT (Track 8.1b chunk-1, iter-46)
#include "hyperflint/algebra/algebraic_letters.hpp"
#include "hyperflint/algebra/linear_factors.hpp"  // clear_linear_factors_cache
#include "hyperflint/algebra/partial_fractions.hpp"  // Track 8.1b chunk-2b iter-48
#include "hyperflint/algebra/shuffle.hpp"
#include "hyperflint/core/zw_table.hpp"  // Track 8.1b chunk-2b iter-48: handlers::partial_fractions transient
#include "hyperflint/convert/convert_hlog.hpp"
#include "hyperflint/convert/parse.hpp"
#include "hyperflint/reduce/mzv_expansion.hpp"   // HF basis-ctx campaign (PHASE_2 iter 10)
#include "hyperflint/core/poly.hpp"
#include "hyperflint/core/rat.hpp"
#include "hyperflint/core/factored_rat.hpp"  // A3: deferred-denominator side-channel
#include "hyperflint/integrator/lr_scan.hpp"  // find_lr_orders_scan op (Doppio-port bridge, 2026-06-06)
#include "hyperflint/integrator/factor_table.hpp"  // factor_table op (spec 2026-06-11)
#include "hyperflint/core/addpf_probe.hpp"  // period-tuples Phase 0 census
#include "hyperflint/core/symcoef.hpp"
#include "hyperflint/integrator/ctx_probe.hpp"
#include "hyperflint/integrator/env_flags.hpp"  // iter-77 Track-probe-ctx (cross-domain rule-3)
#include "hyperflint/integrator/hyper_int.hpp"
#include "hyperflint/integrator/integration_step.hpp"  // NarrowCtxTooNarrow
#include "hyperflint/integrator/lr_search.hpp"
#include "hyperflint/integrator/regularize.hpp"
#include "hyperflint/integrator/step_strategy.hpp"  // Track 6.5 wire-in iter-40
#include "hyperflint/reduce/mzv_reduce.hpp"             // clear_rhs_cache
#include "hyperflint/core/period_table.hpp"             // period_powers emission (Phase 2)
#include "hyperflint/reduce/period_scratch.hpp"          // period_tuples_enabled (Phase 2)
#include "hyperflint/runtime/narrow_ctx_flag.hpp"       // reset_narrow_ctx_flag
#include "hyperflint/runtime/env_flags.hpp"             // HF_FLAG_NARROW_CTX (iter-69)

#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <flint/flint.h>
#include <iostream>
#include <memory>
#include <mutex>      // std::lock_guard/std::mutex (not transitively included on libstdc++)
#include <regex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>

// Track 8.1 (iter-43): defensive fallback so handlers.cpp still
// compiles if some downstream target somehow links it without picking
// up the HF_VERSION_STRING compile-def from the hyperflint /
// hyperflint_nomp / hyperflint_nomp_sd libraries.  Production builds
// always supply HF_VERSION_STRING via CMakeLists.txt top-of-file
// HF_VERSION cache var.  "unknown" is a load-bearing sentinel: the
// Mma-side gate in SubTropica.wl warns on hf_version == "unknown" as
// "binary was not built with HF_VERSION_STRING — verify your build."
#ifndef HF_VERSION_STRING
#define HF_VERSION_STRING "unknown"
#endif

namespace hyperflint {
namespace handlers {

// Track 8.1 (iter-43): eval-json response schema version.
//
// Semantics (request side):
//   - Optional request field `"schema_version_min": <int>`.  If the
//     request asserts a minimum schema version greater than what this
//     binary serves, the handler returns an error_json with a
//     structured message — the caller (Mma / future C ABI) MUST gate
//     subsequent parses on the response NOT being an error before
//     reading the regular fields.
//
// Semantics (response side):
//   - Every successful response now carries:
//       "schema_version": <int>   (currently 1)
//       "hf_version": "<string>"  (from HF_VERSION_STRING, default
//                                  "unknown" if the build is older
//                                  than this define)
//   - Error responses also carry these fields so the caller can
//     diagnose mismatches against an *error* without first having to
//     parse a successful payload.
//
// Bump policy:
//   - Backwards-compatible field additions: keep schema_version=1.
//   - Renames, removals, or semantic changes to existing fields:
//     bump schema_version and update SubTropica.wl's
//     $SubTropicaHFSchemaVersionExpected in lockstep.
//
// Track 8.1 scope (iter-43): find_lr_orders body + error_json +
// error_json_op all stamp the envelope.  error_json_op was folded
// IN-ITER per iter-43 review note Q1/Q8 advisory
// (the one-line shim is mechanical and avoids leaving an
// inconsistent-error-surface gap for hyperflint_sym).  Per-op
// success-path emission for hyperflint_sym remains a future-iter
// item if the op's response shape ever stabilises behind a documented
// schema; today it returns opaque payloads and adding the envelope
// would require an explicit response-builder refactor.  No silent
// inconsistency between error paths.
// §K.1 SSOT retrofit (iter-46 Track 8.1b chunk-1): kSchemaVersion sources
// HF_SCHEMA_VERSION from c_abi.h so the eval-json envelope (Track 8.1) and
// the forthcoming stable C ABI (Track 8.1b chunks 2-4) share a single bump
// point.  Replaces the literal `1` per iter-44 review note B3 BINDING fold.
constexpr int kSchemaVersion = HF_SCHEMA_VERSION;
inline const char* kHFVersion() { return HF_VERSION_STRING; }

namespace {

// ---------- minimal JSON helpers (duplicated from bridge/cli/main.cpp
// for now; promote to shared bridge/ utility if a third caller appears).

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

bool json_bool_field(const std::string& body, const std::string& key) {
    // Track 8.1b chunk-2b (iter-48): mirror of the c_abi.cpp anon-namespace
    // helper; promotion to shared utility deferred to chunk-4 per iter-44
    // plan §F.
    std::regex re("\"" + key + "\"\\s*:\\s*(true|false)");
    std::smatch m;
    if (!std::regex_search(body, m, re)) return false;
    return m[1] == "true";
}

std::vector<std::string> autoscan_vars_single(const std::string& expr) {
    // Single-string autoscan: identifier-pattern scan used when the
    // caller doesn't supply an explicit "vars" array.  CLI's main.cpp
    // exposes an initializer-list variant; we keep this minimal form
    // local because partial_fractions only scans one expression.
    std::set<std::string> seen;
    std::regex re("[A-Za-z][A-Za-z0-9_]*");
    for (auto it = std::sregex_iterator(expr.begin(), expr.end(), re);
         it != std::sregex_iterator(); ++it) {
        seen.insert((*it)[0]);
    }
    return std::vector<std::string>(seen.begin(), seen.end());
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

std::string extract_top_array(const std::string& body,
                               const std::string& key) {
    std::string search = "\"" + key + "\"";
    size_t k = body.find(search);
    if (k == std::string::npos) return {};
    size_t colon = body.find(':', k + search.size());
    if (colon == std::string::npos) return {};
    size_t bracket = body.find('[', colon);
    if (bracket == std::string::npos) return {};
    int depth = 1;
    size_t i = bracket + 1;
    while (i < body.size() && depth > 0) {
        if (body[i] == '[') depth++;
        else if (body[i] == ']') depth--;
        ++i;
    }
    if (depth != 0) return {};
    return body.substr(bracket + 1, i - bracket - 2);
}

std::string error_json(const std::string& msg) {
    // Track 8.1 (iter-43): error responses also carry schema_version
    // and hf_version so the caller can route on mismatch vs other
    // failure classes without first having to parse a successful
    // payload.  Field ordering matches the success path emitted by
    // find_lr_orders below: op, schema_version, hf_version, error.
    std::ostringstream o;
    o << "{\"op\":\"find_lr_orders\""
      << ",\"schema_version\":" << kSchemaVersion
      << ",\"hf_version\":\"" << json_escape(kHFVersion()) << "\""
      << ",\"error\":\"" << json_escape(msg) << "\"}";
    return o.str();
}

// Budget-exceeded response (2026-06-20).  DISTINCT from a normal NOLR
// result: it carries "error" (so the Mma wrapper's
// `KeyExistsQ[resp,"error"]` -> $Failed path fires, exactly like any
// other failure), PLUS "budget_exceeded":true and "reason" so a caller
// that wants to distinguish a clean budget bail from a genuine engine
// error can do so.  CRITICAL: this NEVER emits best_order / score /
// nolr, so a budget bail can never be mistaken for a NOLR verdict.
// Field ordering matches error_json (op, schema_version, hf_version,
// then the failure fields) so a streaming parser sees the envelope head
// identically.
std::string budget_exceeded_json(const std::string& reason) {
    std::ostringstream o;
    o << "{\"op\":\"find_lr_orders\""
      << ",\"schema_version\":" << kSchemaVersion
      << ",\"hf_version\":\"" << json_escape(kHFVersion()) << "\""
      << ",\"budget_exceeded\":true"
      << ",\"reason\":\"" << json_escape(reason) << "\""
      << ",\"error\":\"" << json_escape(reason) << "\"}";
    return o.str();
}

std::string error_json_op(const std::string& op, const std::string& msg) {
    // Track 8.1 (iter-43, in-iter fold per review note
    // rec Q1/Q8): error_json_op now ALSO carries schema_version +
    // hf_version so the eval-json error envelope is uniform across
    // every op handler (find_lr_orders, hyperflint_sym, and any
    // future op).  Without this, a hyperflint_sym error path returns
    // a bare {"op":"hyperflint","error":"..."} and any C-ABI consumer
    // cannot diagnose whether the failure is schema / version /
    // anything else.  One-line shim; matches error_json above.
    std::ostringstream o;
    o << "{\"op\":\"" << json_escape(op) << "\""
      << ",\"schema_version\":" << kSchemaVersion
      << ",\"hf_version\":\"" << json_escape(kHFVersion()) << "\""
      << ",\"error\":\"" << json_escape(msg) << "\"}";
    return o.str();
}

// HF MZV-rewrite C-prep.4 (iter-32): F2 parse-boundary safety rail.
//
// Iter-31 post-commit advisory review note F4
// flagged the JSON-bridge `Rat::parse` entry points
// at handlers.cpp:186/211/694 as the production entry vector for
// non-canonical-content Rat instances. T4c/T5c at iter-32 (test_rat_
// content_invariance.cpp) PASSED: even at TRUE F2 falsifier inputs
// where each input Rat carries shared num/den Z-content, the rep-swap
// and legacy backends produce byte-identical Rat::to_string. The
// canonical-emission writer (sym_coef_canonical_string in this file)
// is therefore not load-bearing on parse-level integer-content
// coprime-ness within each input Rat.
//
// This rail is defense-in-depth: in debug builds (NDEBUG unset), it
// asserts that every JSON-bridge-parsed Rat satisfies the to_string
// idempotency invariant -- i.e., that re-parsing the canonical
// to_string output produces a Rat with the same to_string. T1 of
// test_rat_content_invariance establishes this on synthetic inputs;
// this rail extends the check to every JSON-bridge production input.
//
// In release builds (NDEBUG set), this is a zero-cost no-op.
//
// See c_prep_4_content_audit_memo.md §5.6 (iter-32 closure rail).
inline void debug_check_parse_idempotent(
    [[maybe_unused]] const hyperflint::PolyCtx& ctx,
    [[maybe_unused]] const hyperflint::Rat& parsed,
    [[maybe_unused]] const char* site_tag) {
#ifndef NDEBUG
    const std::string& s1 = parsed.to_string();
    hyperflint::Rat r2 = hyperflint::Rat::parse(ctx, s1);
    if (s1 != r2.to_string()) {
        std::cerr << "[hyperflint][F2-parse-rail] to_string idempotency "
                     "broken at " << site_tag
                  << ": parsed=\"" << s1
                  << "\", reparsed=\"" << r2.to_string() << "\"\n";
        std::abort();
    }
#endif
}

// ---------- hyperflint_sym supporting helpers (migrated from
// bridge/cli/main.cpp anon namespace so the library + LibraryLink can
// both reach them).  Keep `static` / anon-namespace linkage; they're
// implementation detail of handlers.cpp.

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

std::string resolve_mzv_data_path(const std::string& body) {
    std::string data_path = json_str_field(body, "mzv_data_path");
    if (!data_path.empty()) return data_path;
    const char* env = std::getenv("HYPERFLINT_DATA_DIR");
    if (env && *env) return std::string(env) + "/mzv_reductions.json";
    return "data/mzv_reductions.json";
}

// Phase γ.2: process-global MZV-table cache.  The CLI reloads the table
// every process; LibraryLink keeps the same process alive across Mma
// calls, so a per-process cache turns ~15 ms/call of JSON parsing into
// ~0.05 ms/call.  Keyed on path + mtime so a stale .dylib with a fresh
// reductions file still rebuilds.  Single-threaded CLI never races; the
// LibraryLink transport is called from the main Mma kernel thread and
// subkernel kernels get their own process with their own cache, so no
// mutex is needed.
const hyperflint::MzvReductionTable&
get_cached_mzv_table(const std::string& path) {
    static std::string cached_path;
    static std::time_t cached_mtime = 0;
    static std::unique_ptr<hyperflint::MzvReductionTable> cached;
    std::time_t mtime = 0;
    struct stat st;
    if (::stat(path.c_str(), &st) == 0) mtime = st.st_mtime;
    if (!cached || path != cached_path || mtime != cached_mtime) {
        cached = std::make_unique<hyperflint::MzvReductionTable>(
            hyperflint::load_mzv_reductions(path));
        cached_path  = path;
        cached_mtime = mtime;
    }
    return *cached;
}

hyperflint::Word parse_word(const hyperflint::PolyCtx& ctx,
                             const std::vector<std::string>& letters) {
    hyperflint::Word w;
    w.letters.reserve(letters.size());
    for (const auto& s : letters) {
        // HF MZV-rewrite C-prep.4 iter-32 F2 parse-boundary safety rail:
        // debug-only to_string idempotency check. See helper comment.
        hyperflint::Rat r = hyperflint::Rat::parse(ctx, s);
        debug_check_parse_idempotent(ctx, r, "parse_word/letter");
        w.letters.push_back(std::move(r));
    }
    return w;
}

hyperflint::ShuffleList parse_shuffle_list(const hyperflint::PolyCtx& ctx,
                                           const std::string& body,
                                           const std::string& key) {
    hyperflint::ShuffleList input;
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
        // HF MZV-rewrite C-prep.4 iter-32 F2 parse-boundary safety rail.
        hyperflint::Rat coef_r = hyperflint::Rat::parse(ctx, coef_s);
        debug_check_parse_idempotent(ctx, coef_r,
                                      "parse_shuffle_list/coef");
        hyperflint::ShuffleEntry ent{std::move(coef_r), {}};
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

std::string sym_coef_to_mma_string(const hyperflint::SymCoef& s) {
    if (s.is_zero()) return "0";
    std::ostringstream o;
    bool first = true;
    for (const auto& m : s.terms()) {
        std::string pre = m.prefactor.to_string();
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
        // Period-tuples Phase 2: render period generators by their atom
        // names (same tokens the wide-ctx path embeds in the Rat string,
        // so the Mathematica-side parser reuses its existing handling).
        for (const auto& kv : m.period_powers) {
            const std::string& nm =
                hyperflint::PeriodTable::instance().key_for(kv.first);
            if (kv.second == 1) o << "*" << nm;
            else                o << "*" << nm << "^" << kv.second;
        }
        first = false;
    }
    return o.str();
}

// Phase 7-vi-b: serialize the process-global AlgebraicLetterTable so
// the Mma caller can register matching entries in
// HyperIntica`$HyperAlgebraicLetterTable (index-remapped onto
// $HyperAlgebraicLetterCounter) and keep SimplifyWithVieta /
// GetAlgebraicBackSubRules / stEchoAlgebraicLettersSummary working
// unchanged on Wm[i]/Wp[i] symbols that came out of HF.  Each entry
// emits the raw degree-2 polynomial data; Mma derives WmValue/WpValue
// locally as (-b ∓ √disc)/(2·lc).
//
// Output shape: JSON array of objects with keys
//   "idx":      1-based HF-local index  (remapped on Mma side)
//   "poly":     the degree-2 polynomial as a parseable string
//   "var":      name of the special variable inside `poly`
//   "lc":       leading coefficient string
//   "sum":      Vieta sum (-b/lc) string
//   "product":  Vieta product (c/lc) string
//   "disc":     discriminant (b² - 4·lc·c) string
std::string emit_algebraic_letter_table() {
    const auto& tbl = hyperflint::AlgebraicLetterTable::global();
    std::ostringstream o;
    o << "[";
    bool first = true;
    for (long i : tbl.indices()) {
        const auto& e = tbl.at(i);
        if (!first) o << ",";
        first = false;
        o << "{\"idx\":"     << e.idx
          << ",\"poly\":\""  << json_escape(e.polynomial.to_string())     << "\""
          << ",\"var\":\""   << json_escape(e.polynomial.ctx().vars()[e.var_idx]) << "\""
          << ",\"lc\":\""    << json_escape(e.lc.to_string())              << "\""
          << ",\"sum\":\""   << json_escape(e.sum_value.to_string())       << "\""
          << ",\"product\":\"" << json_escape(e.product_value.to_string()) << "\""
          << ",\"disc\":\""  << json_escape(e.discriminant.to_string())    << "\""
          << "}";
    }
    o << "]";
    return o.str();
}

std::string emit_regulator_sym(const hyperflint::RegulatorSym& r) {
    std::ostringstream o;
    o << "[";
    for (size_t i = 0; i < r.size(); ++i) {
        if (i) o << ",";
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

// HF MZV-rewrite C-prep.4 (iter-27) -- canonical-emission writer.
//
// Path (c) per c_prep_4_scoping_memo.md: enforce structural-canonical
// input by calling SymCoef::canonicalize() before emitting. On already-
// canonical input the byte sequence matches sym_coef_to_mma_string
// exactly; on non-canonical input (post-OMP-merge cross-thread reorder
// at C0a) it absorbs the structural-permutation drift into a single
// canonical byte sequence.
//
// REFINEMENT, not a weakening: the comparator stays sha256 cross-cell
// match on the canonical-emission output; on inputs that already pass
// bit-identity today this writer produces the same bytes as
// sym_coef_to_mma_string. On structurally-distinct-but-algebraically-
// equivalent inputs the comparator still rejects (e.g. Pi^(2k) vs zeta
// absorption is NOT performed by canonicalize()). See iter-27 drift
// check M1 + memo §5 case 4.
//
// Iter-30 (HF MZV-rewrite C-prep.4 content audit close,
// c_prep_4_content_audit_memo.md §5.4): T4+T5 verified that the two
// production hot paths (`Rat::add_repswap` == `add_via_q_underscore`
// and `Rat::add_legacy` == cross-mult+gcd_cofactors) produce
// byte-identical `Rat::to_string` on algebraically-equal operands at
// both nvars=60 (rep-swap dispatch regime) and nvars=12 (Smirnov tst2
// legacy dispatch regime).  The Rat instances reaching this writer in
// production come exclusively through these arithmetic paths plus the
// `Rat::Rat(Poly,Poly)` ctor on arithmetic results, all of which
// preserve content normalization across paths.
//
// **Parse-level invariant (T3 caveat)**: `fmpq_mpoly_get_str_pretty`
// is integer-content-faithful for Rats produced by the production
// arithmetic paths, but NOT for Rats produced by `Rat::parse` of an
// arbitrary non-coprime-content input string (e.g.
// `Rat::parse("(2*x)/(2*x+2)")` lands on `"2*x/(2*x + 2)"` while
// `Rat::parse("x/(x+1)")` lands on `"x/(x + 1)"`).  If a future caller
// starts feeding `Rat::parse(non_coprime_content_string)` directly
// into `sym_coef_canonical_string`, add `from_canonical_normalize_content`
// (~30 LOC: pull integer content out of num and den, divide both by
// the common gcd_Z, re-emit) as a preprocessor at the writer entry
// before doing so. See c_prep_4_content_audit_memo.md §5 for the
// audit verdict and §5.2 for the remediation sketch.
std::string sym_coef_canonical_string(const hyperflint::SymCoef& s) {
    return sym_coef_to_mma_string(s.canonicalize());
}

std::string emit_regulator_sym_canonical(const hyperflint::RegulatorSym& r) {
    // canonicalize_regulator_sym (break_up_contour.cpp:222) sorts by
    // regkey_content_key, collects duplicate keys, drops zero-coef
    // entries; the inner key entries are canonicalize_regkey-sorted.
    hyperflint::RegulatorSym canon = hyperflint::canonicalize_regulator_sym(r);
    std::ostringstream o;
    o << "[";
    for (size_t i = 0; i < canon.size(); ++i) {
        if (i) o << ",";
        o << "{\"coef\":\""
          << json_escape(sym_coef_canonical_string(canon[i].coef)) << "\""
          << ",\"key\":[";
        for (size_t j = 0; j < canon[i].key.size(); ++j) {
            if (j) o << ",";
            o << "[";
            for (size_t k = 0; k < canon[i].key[j].size(); ++k) {
                if (k) o << ",";
                o << "\"" << json_escape(canon[i].key[j][k].to_string()) << "\"";
            }
            o << "]";
        }
        o << "]}";
    }
    o << "]";
    return o.str();
}

}  // namespace

// ---------- find_lr_orders ----------

std::string find_lr_orders(const std::string& body) {
    try {
        // Track 8.1 (iter-43): optional request-side gate.  If the
        // caller asserts a minimum schema version greater than what
        // this binary supports, fail fast before any algebraic work.
        // The error_json path also stamps schema_version + hf_version
        // so the caller can diagnose the mismatch without ambiguity.
        // Negative or non-numeric values fall through to std::stoi
        // throwing — caught by the outer try/catch as a parse error.
        {
            std::regex re_sv(
                R"~("schema_version_min"\s*:\s*([0-9]+))~");
            std::smatch m_sv;
            if (std::regex_search(body, m_sv, re_sv)) {
                const int requested = std::stoi(m_sv[1]);
                if (requested > kSchemaVersion) {
                    return error_json(
                        "schema_version_min=" + std::to_string(requested)
                        + " exceeds supported schema_version="
                        + std::to_string(kSchemaVersion)
                        + " (hf_version=" + kHFVersion() + ")");
                }
            }
        }

        auto xvars = json_str_array(body, "xvars");
        if (xvars.empty()) return error_json("need \"xvars\"");
        auto coeff_vars = json_str_array(body, "coeff_vars");

        // Group shape: either "groups":[[...],[...]] (multi-group: one
        // group per ADDEND; the LR search finds an order reducible for
        // EVERY group = the INTERSECTION) or "polys":[...] (single group
        // = the UNION of all letters = the FUSED face). For a sum of
        // denominator-disjoint counterterm addends the union is NOLR by
        // construction even when each addend is reducible in a common
        // order, so use "groups" for tbox/lazy-sum faces -- "polys" is a
        // convenience for genuinely single-integrand inputs only.
        // See docs/cross-subsystem-invariants.md (INV-LAZYSUM-GROUPS) and
        // HyperFLINT/docs/op-contracts.md#find_lr_orders.
        std::vector<std::vector<std::string>> group_strs;
        std::string groups_inner = extract_top_array(body, "groups");
        if (!groups_inner.empty()) {
            int depth = 0;
            size_t start = 0;
            for (size_t i = 0; i < groups_inner.size(); ++i) {
                if (groups_inner[i] == '[') {
                    if (depth == 0) start = i;
                    depth++;
                } else if (groups_inner[i] == ']') {
                    depth--;
                    if (depth == 0) {
                        std::string sub_obj = "{\"xs\":" +
                            groups_inner.substr(start, i - start + 1) + "}";
                        group_strs.push_back(json_str_array(sub_obj, "xs"));
                    }
                }
            }
        } else {
            auto polys_str = json_str_array(body, "polys");
            if (polys_str.empty()) return error_json("need \"polys\" or \"groups\"");
            group_strs.push_back(std::move(polys_str));
        }
        if (group_strs.empty()) return error_json("empty group list");

        std::vector<std::string> all_vars = xvars;
        for (const auto& cv : coeff_vars) all_vars.push_back(cv);

        hyperflint::PolyCtx ctx(all_vars);
        std::vector<size_t> xvar_indices;
        xvar_indices.reserve(xvars.size());
        for (const auto& v : xvars) xvar_indices.push_back(ctx.index_of(v));

        std::vector<std::vector<hyperflint::Poly>> group_polys;
        group_polys.reserve(group_strs.size());
        for (size_t g = 0; g < group_strs.size(); ++g) {
            std::vector<hyperflint::Poly> parsed;
            parsed.reserve(group_strs[g].size());
            for (const auto& s : group_strs[g]) {
                try {
                    parsed.emplace_back(ctx, s);
                } catch (const std::exception& e) {
                    return error_json(std::string("poly parse failed in group ")
                        + std::to_string(g) + ": " + e.what());
                }
            }
            group_polys.push_back(std::move(parsed));
        }

        // Phase 7-vii: optional algebraic_letters flag.
        bool allow_al = false;
        {
            std::regex re("\"algebraic_letters\"\\s*:\\s*(true|false)");
            std::smatch m;
            if (std::regex_search(body, m, re)) {
                allow_al = (m[1] == "true");
            }
        }

        // Carry-discharge (Doppio FindRoots) keep rule (2026-06-07),
        // DEFAULT OFF: absent => false (spec 2026-06-10-carry-option-
        // design.md 4a.1 — the 2026-06-07 WIP default-ON escaped into
        // the v1.2.3 release with no WL sender; absent => false
        // restores released 1.2.2.x semantics for all field-less
        // callers).  Only active when algebraic_letters is on; an
        // explicit "carry_discharge":true selects the carry DFS.
        // Mirrors fr_judge exactly via the shared
        // lr_scan::step_fr_judge primitive — see find_lr_orders.
        bool carry_discharge = false;
        {
            std::regex re("\"carry_discharge\"\\s*:\\s*(true|false)");
            std::smatch m;
            if (std::regex_search(body, m, re)) {
                carry_discharge = (m[1] == "true");
            }
        }

        // ScorePruneFactor (2026-06-17): optional "score_prune_factor":<num>
        // relative branch-and-bound cutoff over the subset DP.  Absent =>
        // +inf (no pruning); the Mma wrapper omits the field for Infinity.
        double score_prune_factor =
            std::numeric_limits<double>::infinity();
        {
            std::regex re("\"score_prune_factor\"\\s*:\\s*"
                          "([0-9]+\\.?[0-9]*([eE][+-]?[0-9]+)?)");
            std::smatch m;
            if (std::regex_search(body, m, re))
                score_prune_factor = std::atof(m[1].str().c_str());
        }

        // Order-resolved singularities (2026-06-07): optional emit_sings
        // flag.  When set, a SingCollector is threaded through the LR
        // walk and accumulates every IRREDUCIBLE kinematic divisor (free
        // of all integration variables, ANY degree) encountered, in the
        // engine's canonical proportionality form.  Absent => collector
        // is nullptr, the walk takes its byte-identical path, and the
        // response carries no "sings" field (gate #1).  The collector
        // observes factors BEFORE/ASIDE from the deg-2 letter cap that
        // bounds the VERDICT; it never alters control flow or the order.
        bool emit_sings = false;
        {
            std::regex re("\"emit_sings\"\\s*:\\s*(true|false)");
            std::smatch m;
            if (std::regex_search(body, m, re)) {
                emit_sings = (m[1] == "true");
            }
        }
        hyperflint::lr_search::SingCollector sings_collector;
        hyperflint::lr_search::SingCollector* sings_ptr =
            emit_sings ? &sings_collector : nullptr;

        // VERIFY-ORDER mode (2026-06-13): optional "verify_order":[v1,v2,...]
        // field.  When present, we VERIFY that this one specific order is
        // linearly reducible (verify_order_is_lr; O(n) st_fubini_lr calls, no
        // O(2^n) search) and SKIP find_lr_orders entirely.  The response adds
        // "order_is_lr" (+ blocking info); the standard envelope fields
        // (best_order empty, nolr, score null) are inert in this mode.  The
        // carry executor's order-pinning guard uses this to certify the
        // PINNED shared order directly (cheap, exact) instead of a free
        // search + best-order comparison.
        std::vector<std::string> verify_order_names;
        bool verify_requested = false;
        {
            std::regex re("\"verify_order\"\\s*:\\s*\\[([^\\]]*)\\]");
            std::smatch m;
            if (std::regex_search(body, m, re)) {
                verify_requested = true;
                std::string inner = m[1].str();
                std::regex name_re("\"([^\"]+)\"");
                for (std::sregex_iterator it(inner.begin(), inner.end(), name_re),
                         end; it != end; ++it) {
                    verify_order_names.push_back((*it)[1].str());
                }
            }
        }

        // Value-initialize (2026-06-20): in VERIFY-ORDER mode
        // find_lr_orders is NOT called, so `result` is never assigned.
        // LrResult's only scalar WITHOUT a default member initializer is
        // `score` (double); a plain `LrResult result;` leaves it
        // indeterminate, and the envelope below reads result.score /
        // result.nolr() unconditionally (the verify short-circuit lives on
        // the Mma side, AFTER the JSON is built).  That produced a garbage
        // "score" (e.g. 2.1e-314) and a bogus "nolr":false in the
        // verify-mode payload.  `LrResult result{};` zero-inits score (and
        // runs the NSDMIs) so verify-mode emits score:null, nolr:false.
        hyperflint::lr_search::LrResult result{};
        hyperflint::lr_search::OrderVerifyResult verify_res;
        double compute_s = 0.0;
        auto t0 = std::chrono::steady_clock::now();
        if (verify_requested) {
            std::vector<size_t> order_idx;
            order_idx.reserve(verify_order_names.size());
            for (const auto& v : verify_order_names)
                order_idx.push_back(ctx.index_of(v));
            verify_res = hyperflint::lr_search::verify_order_is_lr(
                group_polys, xvar_indices, order_idx, allow_al);
        } else {
            result = hyperflint::lr_search::find_lr_orders(
                group_polys, xvar_indices, allow_al, sings_ptr,
                carry_discharge, score_prune_factor);
        }
        auto t1 = std::chrono::steady_clock::now();
        compute_s = std::chrono::duration<double>(t1 - t0).count();

        /* @ARCH:step-strategy v=1 */
        // Track 6.5 / 6.3-followup wire-in (iter-40, satisfies iter-35 rec 3
        // BINDING by iter-45 contract).  Construct StepInputs from the
        // request body + LR-search result, classify via pick_step_strategy
        // (single source of truth in
        // include/hyperflint/integrator/step_strategy.hpp), and emit the
        // resulting StepStrategy enum name in the response JSON.  Mma
        // callers consume this field to dispatch LR_NoOpt / LR_OptOrdered
        // (HF-native) vs Fubini_Lungo / Fubini_Espresso (Mma-delegated
        // fallback) without re-deriving the rule.  See ARCHITECTURE.md §vi.
        std::string method_lr_hint_str = "Lungo";  // Mma default
        {
            std::regex re_hint(
                R"~("method_lr_hint"\s*:\s*"([A-Za-z_]+)")~");
            std::smatch m_hint;
            if (std::regex_search(body, m_hint, re_hint)) {
                method_lr_hint_str = m_hint[1].str();
            }
        }
        std::size_t n_factors_total = 0;
        for (const auto& gp : group_polys) n_factors_total += gp.size();
        hyperflint::integrator::StepInputs strategy_inputs;
        strategy_inputs.degree_budget = allow_al ? 2 : 1;
        strategy_inputs.n_factors     = n_factors_total;
        strategy_inputs.n_letters     = xvar_indices.size();
        strategy_inputs.lr_found      = !result.nolr();
        strategy_inputs.method_lr_hint =
            (method_lr_hint_str == "Espresso")
                ? hyperflint::integrator::StepInputs::MethodLR::Espresso
                : hyperflint::integrator::StepInputs::MethodLR::Lungo;
        const hyperflint::integrator::StepStrategy strategy_choice =
            hyperflint::integrator::pick_step_strategy(strategy_inputs);
        const char* strategy_name = "UNKNOWN";
        switch (strategy_choice) {
            case hyperflint::integrator::StepStrategy::LR_OptOrdered:
                strategy_name = "LR_OptOrdered"; break;
            case hyperflint::integrator::StepStrategy::LR_NoOpt:
                strategy_name = "LR_NoOpt"; break;
            case hyperflint::integrator::StepStrategy::Fubini_Lungo:
                strategy_name = "Fubini_Lungo"; break;
            case hyperflint::integrator::StepStrategy::Fubini_Espresso:
                strategy_name = "Fubini_Espresso"; break;
        }
        /* @ARCH:end step-strategy */

        std::ostringstream o;
        // Track 8.1 (iter-43): schema_version + hf_version emitted
        // immediately after "op" so Mma callers and any future C ABI
        // consumer can gate on the field before doing the full parse.
        // Field-ordering rationale: keeps the version envelope at the
        // head of the JSON so a streaming parser can short-circuit on
        // mismatch without scanning the rest of the payload.
        o << "{\"op\":\"find_lr_orders\""
          << ",\"schema_version\":" << kSchemaVersion
          << ",\"hf_version\":\"" << json_escape(kHFVersion()) << "\""
          << ",\"best_order\":[";
        for (size_t i = 0; i < result.order.size(); ++i) {
            if (i) o << ",";
            const size_t ctx_idx = result.order[i];
            o << "\"" << json_escape(all_vars[ctx_idx]) << "\"";
        }
        o << "],\"score\":";
        // In VERIFY-ORDER mode there is no search result, so the score /
        // nolr fields are inert: emit score:null, nolr:false explicitly
        // (the order-verdict lives in order_is_lr below).  `result` is
        // value-initialized for this mode (score==0), so this guard is a
        // defensive belt-and-braces over the value-init: it guarantees the
        // verify-mode envelope never carries a search-shaped score even if
        // the default ever changes.  The non-verify path is unchanged.
        if (verify_requested || result.nolr() || !std::isfinite(result.score)) {
            o << "null";
        } else {
            o << result.score;
        }
        o << ",\"nolr\":"
          << ((!verify_requested && result.nolr()) ? "true" : "false")
          << ",\"strategy\":\"" << strategy_name << "\""
          << ",\"timing_compute_s\":" << compute_s
          << ",\"nXVars\":" << xvars.size()
          << ",\"nGroups\":" << group_polys.size()
          << ",\"nPolys\":[";
        for (size_t g = 0; g < group_polys.size(); ++g) {
            if (g) o << ",";
            o << group_polys[g].size();
        }
        o << "]";
        // Phase 7-vii: emit the deg-2 polys collected during the LR
        // walk so the SubTropica side knows which polys to introduce
        // as Wm/Wp at integration time.  Mma's STFasterFubini2 returns
        // these as `result[[2]]` under FindRoots=True.
        if (allow_al) {
            o << ",\"root_polys\":[";
            for (size_t i = 0; i < result.root_polys.size(); ++i) {
                if (i) o << ",";
                o << "\"" << json_escape(result.root_polys[i].to_string()) << "\"";
            }
            o << "]";
        }
        // Carry-discharge (Doppio FindRoots) profile of the chosen order.
        // Emitted ONLY when the carry tier actually ran (algebraic_letters
        // on AND carry_discharge on), so the Strict (carry_discharge:false)
        // response stays byte-identical to the pre-change envelope (the
        // regression gate).  These mirror lr_scan's ScanOrder counters:
        // carried_sqrts = deferred sqrt-obligations folded along the path,
        // kin_sqrts = pure-kinematic sqrt letters, terminal_quads =
        // quadratics whose roots are letters of the final answer.  WL
        // consumers read them via Lookup, so adding them to the default
        // response is additive (no existing reader breaks).
        if (allow_al && carry_discharge) {
            o << ",\"carried_sqrts\":" << result.carried_sqrts
              << ",\"kin_sqrts\":" << result.kin_sqrts
              << ",\"terminal_quads\":" << result.terminal_quads;
            // Phase 2 (spec 2026-06-11-carry-phase2 §3.2): the distinct
            // carried obligations along best_order (leaf-replay,
            // proportionality-deduped; size <= carried_sqrts).  Same
            // string serialization as root_polys above.  Additive and
            // carry-gated, so the Strict envelope stays byte-identical;
            // kSchemaVersion stays 2 (bump policy reserves bumps for
            // renames/removals/semantic changes).
            o << ",\"carried_polys\":[";
            for (size_t i = 0; i < result.obligation_polys.size(); ++i) {
                if (i) o << ",";
                o << "\""
                  << json_escape(result.obligation_polys[i].to_string())
                  << "\"";
            }
            o << "]";
        }
        // Order-resolved singularities: emit the collected canonical
        // irreducible kinematic divisors and their count.  Only present
        // when the caller set emit_sings (so the default response is
        // byte-identical).  Order is first-encounter order in the DP
        // walk (SingCollector::ordered), which is deterministic for a
        // fixed input + memo configuration.
        if (emit_sings) {
            o << ",\"sings\":[";
            for (size_t i = 0; i < sings_collector.ordered.size(); ++i) {
                if (i) o << ",";
                o << "\"" << json_escape(sings_collector.ordered[i]) << "\"";
            }
            o << "],\"sings_total\":" << sings_collector.ordered.size();
        }
        // VERIFY-ORDER response (2026-06-13): present only when the request
        // carried "verify_order".  order_is_lr is the single binding bit; the
        // blocking fields explain a false (which step/letter, and whether the
        // failure was a deg-2 forbidden-pending-variable dependence vs a plain
        // degree-too-high).  Additive + request-gated, so the default
        // envelope is byte-identical.
        if (verify_requested) {
            o << ",\"order_is_lr\":"
              << (verify_res.is_lr ? "true" : "false")
              << ",\"verify_malformed\":"
              << (verify_res.malformed ? "true" : "false")
              << ",\"verify_blocking_step\":" << verify_res.blocking_step
              << ",\"verify_blocking_degree\":" << verify_res.blocking_degree
              << ",\"verify_forbidden_dep\":"
              << (verify_res.forbidden_dep ? "true" : "false")
              << ",\"verify_blocking_letter\":\""
              << json_escape(verify_res.blocking_letter) << "\"";
        }
        o << "}";
        return o.str();
    } catch (const hyperflint::lr_search::LrBudgetExceeded& e) {
        // Budget safety net (2026-06-20): the exhaustive LR search (or the
        // verify walk) hit its time / operand-size budget and bailed
        // CLEANLY mid-flight.  Serialize a DISTINCT failed response
        // (budget_exceeded + reason + error) so the Mma wrapper sees
        // $Failed, NOT a NOLR verdict.  MUST come before the generic
        // std::exception catch (LrBudgetExceeded derives from
        // std::runtime_error, so the base handler would otherwise swallow
        // it and mislabel it a plain "error").  RAII (the Poly wrappers)
        // releases any in-flight FLINT memory as the exception unwinds.
        return budget_exceeded_json(e.what());
    } catch (const std::exception& e) {
        return error_json(e.what());
    } catch (...) {
        return error_json("unknown exception in find_lr_orders");
    }
}

// ---------- factor_table ----------
//
// Factor-prediction table for the Fubini reduction (spec
// docs/superpowers/specs/2026-06-11-stfactorpredictor-design.md).
// Single-chain replay along the supplied LR order; tabulates monic
// pair differences (deg-1) and per-letter coefficient/disc factor
// lists, trial-divided against the stage pool with exact fallback.
//
// Request:
//   {"op":"factor_table", "groups":[[...]] | "polys":[...],
//    "xvars":[...], "coeff_vars":[...], "order":[...],
//    "algebraic_letters":bool?, "max_pairs":N?, "max_singletons":N?,
//    "max_response_mb":N?}
// Response: spec section 4.4 (interned "polys" array + "stages" +
// "pairs" + "singletons" + "stats", schema/version envelope).
// All guards are loud errors, never truncation.

namespace {

size_t json_size_field(const std::string& body, const std::string& key,
                       size_t dflt) {
    std::regex re("\"" + key + "\"\\s*:\\s*([0-9]+)");
    std::smatch m;
    if (std::regex_search(body, m, re)) {
        try {
            return std::stoull(m[1]);
        } catch (const std::exception&) {
            throw std::runtime_error("bad value for " + key
                                     + " (out of range)");
        }
    }
    // Key present but not a nonnegative integer: loud error, never a
    // silent fall-back to the default (review fold).
    std::regex re_key("\"" + key + "\"\\s*:");
    if (std::regex_search(body, re_key))
        throw std::runtime_error("bad value for " + key
                                 + " (expected a nonnegative integer)");
    return dflt;
}

void emit_factored_object(std::ostringstream& o,
                          const hyperflint::factor_table::FactoredObject& fo) {
    o << "\"c\":\"" << json_escape(fo.c) << "\",\"factors\":[";
    for (size_t j = 0; j < fo.factors.size(); ++j)
        o << (j ? "," : "") << "[" << fo.factors[j].first << ","
          << fo.factors[j].second << "]";
    o << "],\"oop\":" << (fo.oop ? "true" : "false");
}

}  // namespace

std::string factor_table(const std::string& body) {
    static const char* kOp = "factor_table";
    try {
        {
            std::regex re_sv(R"~("schema_version_min"\s*:\s*([0-9]+))~");
            std::smatch m_sv;
            if (std::regex_search(body, m_sv, re_sv)) {
                const int requested = std::stoi(m_sv[1]);
                if (requested > kSchemaVersion) {
                    return error_json_op(kOp,
                        "schema_version_min=" + std::to_string(requested)
                        + " exceeds supported schema_version="
                        + std::to_string(kSchemaVersion)
                        + " (hf_version=" + kHFVersion() + ")");
                }
            }
        }

        auto xvars = json_str_array(body, "xvars");
        if (xvars.empty()) return error_json_op(kOp, "need \"xvars\"");
        auto coeff_vars = json_str_array(body, "coeff_vars");

        std::vector<std::vector<std::string>> group_strs;
        std::string groups_inner = extract_top_array(body, "groups");
        if (!groups_inner.empty()) {
            int depth = 0;
            size_t start = 0;
            for (size_t i = 0; i < groups_inner.size(); ++i) {
                if (groups_inner[i] == '[') {
                    if (depth == 0) start = i;
                    depth++;
                } else if (groups_inner[i] == ']') {
                    depth--;
                    if (depth == 0) {
                        std::string sub_obj = "{\"xs\":" +
                            groups_inner.substr(start, i - start + 1) + "}";
                        group_strs.push_back(json_str_array(sub_obj, "xs"));
                    }
                }
            }
        } else {
            auto polys_str = json_str_array(body, "polys");
            if (polys_str.empty())
                return error_json_op(kOp, "need \"polys\" or \"groups\"");
            group_strs.push_back(std::move(polys_str));
        }
        if (group_strs.empty()) return error_json_op(kOp, "empty group list");

        auto order_strs = json_str_array(body, "order");
        if (order_strs.empty()) return error_json_op(kOp, "need \"order\"");
        {
            // Permutation validation: string multisets must agree
            // (spec 4.2; exact-name matching, no aliasing).
            std::vector<std::string> a = order_strs, b = xvars;
            std::sort(a.begin(), a.end());
            std::sort(b.begin(), b.end());
            if (a != b)
                return error_json_op(kOp,
                    "order is not a permutation of xvars");
        }

        std::vector<std::string> all_vars = xvars;
        for (const auto& cv : coeff_vars) all_vars.push_back(cv);

        hyperflint::PolyCtx ctx(all_vars);
        std::vector<size_t> order_indices;
        order_indices.reserve(order_strs.size());
        for (const auto& v : order_strs) {
            const size_t i = ctx.index_of(v);
            if (i == SIZE_MAX)
                return error_json_op(kOp,
                    "unknown variable in order: " + v);
            order_indices.push_back(i);
        }

        std::vector<std::vector<hyperflint::Poly>> group_polys;
        group_polys.reserve(group_strs.size());
        for (size_t g = 0; g < group_strs.size(); ++g) {
            std::vector<hyperflint::Poly> parsed;
            parsed.reserve(group_strs[g].size());
            for (const auto& s : group_strs[g]) {
                try {
                    parsed.emplace_back(ctx, s);
                } catch (const std::exception& e) {
                    return error_json_op(kOp,
                        std::string("poly parse failed in group ")
                        + std::to_string(g) + ": " + e.what());
                }
            }
            group_polys.push_back(std::move(parsed));
        }

        bool allow_al = false;
        {
            std::regex re("\"algebraic_letters\"\\s*:\\s*(true|false)");
            std::smatch m;
            if (std::regex_search(body, m, re)) allow_al = (m[1] == "true");
        }

        hyperflint::factor_table::Limits lim;
        lim.max_pairs = json_size_field(body, "max_pairs", lim.max_pairs);
        lim.max_singletons =
            json_size_field(body, "max_singletons", lim.max_singletons);
        lim.max_response_mb =
            json_size_field(body, "max_response_mb", lim.max_response_mb);

        // Per-request memo hygiene: cold memos, bounded memory, clean
        // trace stats (spec 4.1; warm sharing deliberately deferred).
        // reset_lr_budget starts the steady_clock deadline fresh and
        // clears any stale (already-past) deadline from a prior budgeted
        // call in the same process; no-op when the budget env vars are
        // unset.
        hyperflint::lr_search::reset_lr_memos();
        hyperflint::lr_search::reset_lr_trace();
        hyperflint::lr_search::reset_lr_budget();

        hyperflint::factor_table::FactorTable t =
            hyperflint::factor_table::build(ctx, group_polys, order_indices,
                                            allow_al, lim);

        std::ostringstream o;
        o << "{\"op\":\"factor_table\""
          << ",\"schema_version\":" << kSchemaVersion
          << ",\"hf_version\":\"" << json_escape(kHFVersion()) << "\""
          << ",\"order\":[";
        for (size_t i = 0; i < order_strs.size(); ++i)
            o << (i ? "," : "") << "\"" << json_escape(order_strs[i]) << "\"";
        o << "],\"polys\":[";
        for (size_t i = 0; i < t.intern_strs.size(); ++i)
            o << (i ? "," : "") << "\"" << json_escape(t.intern_strs[i])
              << "\"";
        o << "],\"stages\":[";
        for (size_t s = 0; s < t.stages.size(); ++s) {
            const auto& st = t.stages[s];
            o << (s ? "," : "") << "{\"var\":\""
              << json_escape(ctx.vars()[st.var_idx])
              << "\",\"admissible\":[";
            for (size_t g = 0; g < st.admissible.size(); ++g) {
                o << (g ? "," : "") << "[";
                for (size_t i = 0; i < st.admissible[g].size(); ++i)
                    o << (i ? "," : "") << st.admissible[g][i];
                o << "]";
            }
            o << "],\"pool\":[";
            for (size_t i = 0; i < st.pool.size(); ++i)
                o << (i ? "," : "") << st.pool[i];
            o << "],\"n_pairs\":" << st.n_pairs
              << ",\"n_singletons\":" << st.n_singletons
              << ",\"n_inadmissible\":" << st.n_inadmissible
              << ",\"t_build_s\":" << st.t_build_s << "}";
        }
        o << "],\"pairs\":[";
        for (size_t i = 0; i < t.pairs.size(); ++i) {
            const auto& pe = t.pairs[i];
            o << (i ? "," : "") << "{\"var\":\""
              << json_escape(ctx.vars()[pe.var_idx])
              << "\",\"f\":" << pe.f_id << ",\"g\":" << pe.g_id << ",";
            emit_factored_object(o, pe.diff);
            o << "}";
        }
        o << "],\"singletons\":[";
        for (size_t i = 0; i < t.singletons.size(); ++i) {
            const auto& se = t.singletons[i];
            o << (i ? "," : "") << "{\"var\":\""
              << json_escape(ctx.vars()[se.var_idx])
              << "\",\"id\":" << se.id << ",\"deg\":" << se.deg
              << ",\"coeffs\":[";
            for (size_t j = 0; j < se.coeffs.size(); ++j) {
                o << (j ? "," : "") << "{\"power\":" << se.coeffs[j].power
                  << ",";
                emit_factored_object(o, se.coeffs[j].fo);
                o << "}";
            }
            o << "]";
            if (se.has_disc) {
                o << ",\"disc\":{";
                emit_factored_object(o, se.disc);
                o << "}";
            }
            o << "}";
        }
        o << "],\"stats\":{\"pairs_total\":" << t.stats.pairs_total
          << ",\"singletons_total\":" << t.stats.singletons_total
          << ",\"oop\":" << t.stats.oop
          << ",\"pair_fallbacks\":" << t.stats.pair_fallbacks
          << ",\"trial_s\":" << t.stats.trial_s
          << ",\"fallback_s\":" << t.stats.fallback_s << "}}";

        // Post-hoc bound: the response is fully materialized before
        // this check; the in-build max_pairs / max_singletons guards
        // bound the size in practice (review fold: documented, not a
        // streaming cap).
        std::string out = o.str();
        if (out.size() > lim.max_response_mb * 1024ull * 1024ull)
            return error_json_op(kOp,
                "max_response_mb exceeded (response "
                + std::to_string(out.size() / (1024ull * 1024ull))
                + " MB, cap " + std::to_string(lim.max_response_mb)
                + " MB)");
        return out;
    } catch (const std::exception& e) {
        return error_json_op(kOp, e.what());
    } catch (...) {
        return error_json_op(kOp, "unknown exception in factor_table");
    }
}

// ---------- find_lr_orders_scan ----------

// Doppio-port phase 3 bridge op (2026-06-06): expose
// lr_scan::find_lr_orders_scan — the projective Cheng-Wu GAUGE SCAN
// with the Doppio keep rules (Strict / FindRoots carried-sqrt tiers).
// This is the PRE-GAUGE engine: groups are the RAW factor lists of the
// ungauged n-variable system (boundary monomials added internally),
// and per-poly twist exponents a + b*eps are REQUIRED for the
// projectivity gate (sum a_i d_i == -n, sum b_i d_i == 0).
//
// Request:
//   {"op":"find_lr_orders_scan",
//    "schema_version_min": 1,                  (optional fast-fail gate)
//    "groups":[["<p1>", ...], ...],            (RAW groups)
//    "xvars":["x1", ...],
//    "coeff_vars":["s", ...],                  (optional)
//    "exps":[[[a,b], ...], ...],               (REQUIRED; shape == groups)
//    "keep_rule":"Strict"|"FindRoots",         (default "Strict")
//    "euler_filter":true|false,                (default false; needs msolve)
//    "max_orders":N}                           (default 8192)
// Response:
//   {"op":"find_lr_orders_scan", schema/version envelope,
//    "projective":bool, "truncated":bool,
//    "orders":[{"order":[...names...],"gauge":"<name>","score":F,
//               "carried_sqrts":n,"kin_sqrts":n,"terminal_quads":n},...],
//    "timing_compute_s":F, "nXVars":K, "nGroups":G}
// orders is score-ascending; empty + projective=false means the input
// failed the projectivity gate (NOT a NOLR verdict); empty +
// projective=true means no admissible (gauge, order) pair exists under
// the requested keep rule.
std::string find_lr_orders_scan(const std::string& body) {
    static const char* kOp = "find_lr_orders_scan";
    try {
        {
            std::regex re_sv(R"~("schema_version_min"\s*:\s*([0-9]+))~");
            std::smatch m_sv;
            if (std::regex_search(body, m_sv, re_sv)) {
                const int requested = std::stoi(m_sv[1]);
                if (requested > kSchemaVersion) {
                    return error_json_op(kOp,
                        "schema_version_min=" + std::to_string(requested)
                        + " exceeds supported schema_version="
                        + std::to_string(kSchemaVersion)
                        + " (hf_version=" + kHFVersion() + ")");
                }
            }
        }

        auto xvars = json_str_array(body, "xvars");
        if (xvars.empty()) return error_json_op(kOp, "need \"xvars\"");
        auto coeff_vars = json_str_array(body, "coeff_vars");

        // groups: list of lists of poly strings (same parser shape as
        // find_lr_orders, multi-group form only — the scan is defined
        // on per-term groups).
        std::vector<std::vector<std::string>> group_strs;
        {
            std::string groups_inner = extract_top_array(body, "groups");
            if (groups_inner.empty())
                return error_json_op(kOp, "need \"groups\"");
            int depth = 0;
            size_t start = 0;
            for (size_t i = 0; i < groups_inner.size(); ++i) {
                if (groups_inner[i] == '[') {
                    if (depth == 0) start = i;
                    depth++;
                } else if (groups_inner[i] == ']') {
                    depth--;
                    if (depth == 0) {
                        std::string sub_obj = "{\"xs\":" +
                            groups_inner.substr(start, i - start + 1) + "}";
                        group_strs.push_back(json_str_array(sub_obj, "xs"));
                    }
                }
            }
        }
        if (group_strs.empty()) return error_json_op(kOp, "empty group list");

        // exps: [[[a,b],...],...], shape matching groups.  Integer pairs.
        std::vector<std::vector<hyperflint::lr_scan::ScanExponent>> exps;
        {
            std::string exps_inner = extract_top_array(body, "exps");
            if (exps_inner.empty())
                return error_json_op(kOp,
                    "need \"exps\" (per-poly twist exponents [a,b]; "
                    "the projectivity gate is load-bearing)");
            int depth = 0;
            size_t gstart = 0;
            for (size_t i = 0; i < exps_inner.size(); ++i) {
                if (exps_inner[i] == '[') {
                    if (depth == 0) gstart = i;
                    depth++;
                } else if (exps_inner[i] == ']') {
                    depth--;
                    if (depth == 0) {
                        // one group's [[a,b],[a,b],...]
                        const std::string g =
                            exps_inner.substr(gstart, i - gstart + 1);
                        std::vector<hyperflint::lr_scan::ScanExponent> ge;
                        int d2 = 0;
                        size_t pstart = 0;
                        for (size_t j = 1; j + 1 < g.size(); ++j) {
                            if (g[j] == '[') {
                                if (d2 == 0) pstart = j;
                                d2++;
                            } else if (g[j] == ']') {
                                d2--;
                                if (d2 == 0) {
                                    const std::string pair =
                                        g.substr(pstart + 1, j - pstart - 1);
                                    const size_t comma = pair.find(',');
                                    if (comma == std::string::npos)
                                        return error_json_op(kOp,
                                            "malformed exps pair: " + pair);
                                    hyperflint::lr_scan::ScanExponent e;
                                    e.a = std::stol(pair.substr(0, comma));
                                    e.b = std::stol(pair.substr(comma + 1));
                                    ge.push_back(e);
                                }
                            }
                        }
                        exps.push_back(std::move(ge));
                    }
                }
            }
        }
        if (exps.size() != group_strs.size())
            return error_json_op(kOp,
                "exps group count (" + std::to_string(exps.size())
                + ") != groups (" + std::to_string(group_strs.size()) + ")");
        for (size_t g = 0; g < exps.size(); ++g)
            if (exps[g].size() != group_strs[g].size())
                return error_json_op(kOp,
                    "exps[" + std::to_string(g) + "] length "
                    + std::to_string(exps[g].size()) + " != group size "
                    + std::to_string(group_strs[g].size()));

        // keep_rule / euler_filter / max_orders
        hyperflint::lr_scan::KeepRule rule =
            hyperflint::lr_scan::KeepRule::Strict;
        {
            std::string kr = json_str_field(body, "keep_rule");
            if (kr == "FindRoots")
                rule = hyperflint::lr_scan::KeepRule::FindRoots;
            else if (!kr.empty() && kr != "Strict")
                return error_json_op(kOp, "unknown keep_rule: " + kr);
        }
        bool euler_filter = false;
        {
            std::regex re("\"euler_filter\"\\s*:\\s*(true|false)");
            std::smatch m;
            if (std::regex_search(body, m, re))
                euler_filter = (m[1] == "true");
        }
        size_t max_orders = 8192;
        {
            std::regex re(R"~("max_orders"\s*:\s*([0-9]+))~");
            std::smatch m;
            if (std::regex_search(body, m, re))
                max_orders = static_cast<size_t>(std::stoul(m[1]));
        }

        std::vector<std::string> all_vars = xvars;
        for (const auto& cv : coeff_vars) all_vars.push_back(cv);
        hyperflint::PolyCtx ctx(all_vars);
        std::vector<size_t> xvar_indices;
        xvar_indices.reserve(xvars.size());
        for (const auto& v : xvars) xvar_indices.push_back(ctx.index_of(v));

        std::vector<std::vector<hyperflint::Poly>> group_polys;
        group_polys.reserve(group_strs.size());
        for (size_t g = 0; g < group_strs.size(); ++g) {
            std::vector<hyperflint::Poly> parsed;
            parsed.reserve(group_strs[g].size());
            for (const auto& s : group_strs[g]) {
                try {
                    parsed.emplace_back(ctx, s);
                } catch (const std::exception& e) {
                    return error_json_op(kOp,
                        std::string("poly parse failed in group ")
                        + std::to_string(g) + ": " + e.what());
                }
            }
            group_polys.push_back(std::move(parsed));
        }

        auto t0 = std::chrono::steady_clock::now();
        hyperflint::lr_scan::ScanResult result =
            hyperflint::lr_scan::find_lr_orders_scan(
                group_polys, xvar_indices, exps, rule, euler_filter,
                max_orders);
        auto t1 = std::chrono::steady_clock::now();
        const double compute_s =
            std::chrono::duration<double>(t1 - t0).count();

        std::ostringstream o;
        o << "{\"op\":\"" << kOp << "\""
          << ",\"schema_version\":" << kSchemaVersion
          << ",\"hf_version\":\"" << json_escape(kHFVersion()) << "\""
          << ",\"projective\":" << (result.projective ? "true" : "false")
          << ",\"truncated\":" << (result.truncated ? "true" : "false")
          << ",\"orders\":[";
        for (size_t i = 0; i < result.orders.size(); ++i) {
            const auto& so = result.orders[i];
            if (i) o << ",";
            o << "{\"order\":[";
            for (size_t j = 0; j < so.order.size(); ++j) {
                if (j) o << ",";
                o << "\"" << json_escape(all_vars[so.order[j]]) << "\"";
            }
            o << "],\"gauge\":\"" << json_escape(all_vars[so.gauge]) << "\""
              << ",\"score\":" << so.score
              << ",\"carried_sqrts\":" << so.carried_sqrts
              << ",\"kin_sqrts\":" << so.kin_sqrts
              << ",\"terminal_quads\":" << so.terminal_quads
              << "}";
        }
        o << "],\"timing_compute_s\":" << compute_s
          << ",\"nXVars\":" << xvars.size()
          << ",\"nGroups\":" << group_polys.size()
          << "}";
        return o.str();
    } catch (const std::exception& e) {
        return error_json_op(kOp, e.what());
    } catch (...) {
        return error_json_op(kOp, "unknown exception in find_lr_orders_scan");
    }
}

// ---------- partial_fractions ----------

// Track 8.1b chunk-2b (iter-48): transport-neutral partial_fractions
// handler.  Factored out of bridge/cli/main.cpp:916 per iter-44 plan §B
// chunk-2.  Returns the CLI-form response (no envelope, no trailing
// newline) so the CLI shim can print + newline; the C ABI wrapper in
// src/bridge/c_abi.cpp splices the schema_version + hf_version envelope
// in front of the payload (or, in a future chunk, calls a sibling
// envelope-emitting variant directly).
//
// Field ordering and content match the legacy main.cpp output verbatim
// so the iter-48 byte-identical CLI snapshot ctest gates regressions.
// Error reporting: on missing `f` or `var`, returns
// `{"op":"partial_fractions","error":"<msg>"}` (no envelope —
// matches the legacy `std::cerr` + non-zero exit-code path semantically
// but moves the message into JSON so non-CLI callers can route on it).
std::string partial_fractions(const std::string& body) {
    try {
        std::string f   = json_str_field(body, "f");
        std::string var = json_str_field(body, "var");
        if (f.empty() || var.empty()) {
            return std::string(
                "{\"op\":\"partial_fractions\","
                "\"error\":\"need \\\"f\\\" and \\\"var\\\" fields\"}");
        }
        auto vars = json_str_array(body, "vars");
        if (vars.empty()) vars = autoscan_vars_single(f);
        bool present = false;
        for (const auto& v : vars) if (v == var) { present = true; break; }
        if (!present) vars.push_back(var);

        // Phase 7-vi-a: same `algebraic_letters` flag as the other
        // pipeline stages; when true, deg-2 irreducible denominators
        // split into Wm/Wp pairs.
        bool introduce_al = json_bool_field(body, "algebraic_letters");
        std::vector<std::string> ctx_vars = introduce_al
            ? hyperflint::build_algebraic_letter_var_list(vars)
            : vars;
        hyperflint::PolyCtx ctx(ctx_vars);
        hyperflint::Rat r = hyperflint::Rat::parse(ctx, f);
        size_t idx = 0;
        for (; idx < ctx_vars.size(); ++idx) {
            if (ctx_vars[idx] == var) break;
        }
        // Iter-52 C0c.1 Increment β convention: caller-side fresh
        // transient ZWTable.  Bridge/CLI handlers don't share state
        // across requests, so a per-call transient is appropriate.
        auto _lf_zw = std::make_shared<hyperflint::ZWTable>(ctx);
        // Fully qualified call: the enclosing namespace is
        // hyperflint::handlers, where `partial_fractions` would name
        // *this* function via the unqualified-lookup rule.  Use the
        // global `::hyperflint::partial_fractions` to reach the
        // algebra-level algorithm in algebra/partial_fractions.hpp.
        auto pf = ::hyperflint::partial_fractions(r, idx, _lf_zw,
                                                  introduce_al);

        std::ostringstream o;
        o << "{\"op\":\"partial_fractions\""
          << ",\"var\":\"" << json_escape(var) << "\""
          << ",\"polynomial_part\":\""
          << json_escape(pf.polynomial_part.to_string()) << "\""
          << ",\"poles\":[";
        for (size_t i = 0; i < pf.poles.size(); ++i) {
            if (i) o << ",";
            const auto& P = pf.poles[i];
            o << "{\"pole\":\"" << json_escape(P.pole.to_string()) << "\""
              << ",\"multiplicity\":" << P.multiplicity
              << ",\"coefs\":[";
            for (size_t j = 0; j < P.coefs.size(); ++j) {
                if (j) o << ",";
                o << "\"" << json_escape(P.coefs[j].to_string()) << "\"";
            }
            o << "]}";
        }
        o << "],\"vars\":[";
        for (size_t i = 0; i < vars.size(); ++i) {
            if (i) o << ",";
            o << "\"" << json_escape(vars[i]) << "\"";
        }
        o << "]}";
        return o.str();
    } catch (const std::exception& e) {
        return error_json_op("partial_fractions", e.what());
    } catch (...) {
        return error_json_op("partial_fractions",
                             "unknown exception in partial_fractions");
    }
}

// ---------- linear_factors ----------

// Track 8.1b chunk-3 (iter-50): transport-neutral linear_factors
// handler.  Factored out of bridge/cli/main.cpp:933 per iter-44 plan §B
// chunk-3.  Mirror of partial_fractions above (chunk-2b iter-48):
// returns the CLI-form response (no envelope, no trailing newline), the
// CLI shim prints + newline, and the C ABI wrapper in
// src/bridge/c_abi.cpp splices the schema_version + hf_version envelope
// for external consumers.
//
// Field ordering and content match the legacy main.cpp output verbatim
// (constant, linear, nonlinear, vars) so the iter-51 chunk-3b byte-
// identical CLI snapshot ctest will gate regressions.  Error reporting:
// on missing `poly` or `var`, returns
// `{"op":"linear_factors","error":"<msg>"}` (no envelope — matches the
// chunk-2b partial_fractions convention; replaces the legacy
// std::cerr + non-zero exit-code path with a JSON error response so
// non-CLI callers can route on it).
std::string linear_factors(const std::string& body) {
    try {
        std::string poly_s = json_str_field(body, "poly");
        std::string var    = json_str_field(body, "var");
        if (poly_s.empty() || var.empty()) {
            return std::string(
                "{\"op\":\"linear_factors\","
                "\"error\":\"need \\\"poly\\\" and \\\"var\\\" fields\"}");
        }
        auto user_vars = json_str_array(body, "vars");
        if (user_vars.empty()) user_vars = autoscan_vars_single(poly_s);
        bool present = false;
        for (const auto& v : user_vars) if (v == var) { present = true; break; }
        if (!present) user_vars.push_back(var);

        // Phase 7-ii: when introduce_algebraic_letters=true, the deg-2
        // branch needs Wm_<i> / Wp_<i> atoms in the PolyCtx.
        bool introduce_al = json_bool_field(body, "introduce_algebraic_letters");
        std::vector<std::string> vars = introduce_al
            ? hyperflint::build_algebraic_letter_var_list(user_vars)
            : user_vars;

        hyperflint::PolyCtx ctx(vars);
        hyperflint::Poly p(ctx, poly_s);
        size_t idx = 0;
        for (; idx < vars.size(); ++idx) if (vars[idx] == var) break;
        // Iter-52 C0c.1: caller-side fresh transient ZWTable for the
        // mandatory `zw_tab` parameter (Option A).  Bridge/CLI handlers
        // don't share state across requests, so a per-call transient is
        // appropriate, mirroring the pre-iter-52 lambda-internal
        // allocation semantics.
        auto _lf_zw = std::make_shared<hyperflint::ZWTable>(p.ctx());
        // axis-C-lf-constant-defer: standalone bridge-CLI op needs the
        // serialised constant for compare.py cross-tests, so pass
        // compute_constant=true.  Integration callers default to false.
        // Fully qualified call: the enclosing namespace is
        // hyperflint::handlers, where `linear_factors` would name *this*
        // function via unqualified lookup; use ::hyperflint::linear_factors
        // to reach the algebra-level algorithm in
        // algebra/linear_factors.hpp.
        auto lf = ::hyperflint::linear_factors(p, idx, _lf_zw, introduce_al,
                                                /*compute_constant=*/true);

        std::ostringstream o;
        o << "{\"op\":\"linear_factors\""
          << ",\"constant\":\"" << json_escape(lf.constant) << "\""
          << ",\"linear\":[";
        for (size_t i = 0; i < lf.linear.size(); ++i) {
            if (i) o << ",";
            o << "[" << lf.linear[i].multiplicity
              << ",\"" << json_escape(lf.linear[i].pole.num().to_string()) << "\""
              << ",\"" << json_escape(lf.linear[i].pole.den().to_string()) << "\"]";
        }
        o << "],\"nonlinear\":[";
        for (size_t i = 0; i < lf.nonlinear.size(); ++i) {
            if (i) o << ",";
            o << "[" << lf.nonlinear[i].multiplicity
              << ",\"" << json_escape(lf.nonlinear[i].polynomial.to_string()) << "\""
              << "," << lf.nonlinear[i].degree_in_var << "]";
        }
        o << "],\"vars\":[";
        for (size_t i = 0; i < vars.size(); ++i) {
            if (i) o << ",";
            o << "\"" << json_escape(vars[i]) << "\"";
        }
        o << "]}";
        return o.str();
    } catch (const std::exception& e) {
        return error_json_op("linear_factors", e.what());
    } catch (...) {
        return error_json_op("linear_factors",
                             "unknown exception in linear_factors");
    }
}

// ---------- hyperflint_sym ----------

std::string hyperflint_sym(const std::string& body) {
    // Memory operational lever: per-call FLINT thread count from env
    // HF_MAX_THREADS_PER_CALL.  When set, calls
    // flint_set_num_threads(N) before any FLINT work in this call.
    // Use case: Mma master-kernel STIntegrate dispatch needs single-
    // threaded HF to keep per-call peak RSS at ~970 MB (vs ~2.5 GB
    // at OMP=13).  Env-driven so it's overridable per-process; the
    // SubTropica.wl bridge sets it explicitly when calling HF in
    // memory-constrained contexts (master kernel, parallel
    // subkernels).  Not gated by introduce_al / narrow ctx — works
    // on every code path.
    // iter-94 Track-OMP macro-layer LAND: env-var literal relocated to
    // bridge/env_flags.hpp under §5.1 rule-1 BINDING from the
    // iter-94 Q-19 B1 review note (this is the
    // FIRST bridge-domain env_flags header). The POSITIVE_INTEGER value-
    // family semantics (atoi-then-(n>=1) guard) are preserved verbatim.
    if (const char* mt = HF_FLAG_MAX_THREADS_PER_CALL) {
        if (mt && *mt) {
            int n = std::atoi(mt);
            if (n >= 1) flint_set_num_threads(n);
        }
    }
    // R20 step 2 instrumentation: env-gated request-body dump.
    // When HF_PROBE_DUMP_DIR is set to a directory, every entry of
    // hyperflint_sym writes its JSON body to <dir>/face_<NNN>.json
    // where NNN is an atomic per-process counter. Used for the
    // wide-ctx campaign's per-face variance probe.
    if (const char* dir = HF_FLAG_PROBE_DUMP_DIR) {
        if (dir && *dir) {
            static std::atomic<int> g_probe_dump_counter{0};
            int idx = g_probe_dump_counter.fetch_add(1);
            char fname[2048];
            std::snprintf(fname, sizeof(fname),
                "%s/face_%03d_pid%d.json",
                dir, idx, getpid());
            FILE* f = std::fopen(fname, "w");
            if (f) {
                std::fwrite(body.data(), 1, body.size(), f);
                std::fclose(f);
            }
        }
    }
    try {
        auto user_vars = json_str_array(body, "vars");
        auto vars_int = json_str_array(body, "vars_int");
        auto vars_int_from = json_str_array(body, "vars_int_from");
        auto vars_int_to   = json_str_array(body, "vars_int_to");
        if (user_vars.empty()) {
            for (const auto& v : vars_int) user_vars.push_back(v);
        }
        for (const auto& vi : vars_int) {
            bool present = false;
            for (const auto& v : user_vars) if (v == vi) { present = true; break; }
            if (!present) user_vars.push_back(vi);
        }
        if (user_vars.empty()) user_vars.push_back("x");

        const bool have_ranges =
            !vars_int_from.empty() || !vars_int_to.empty();
        if (have_ranges) {
            if (vars_int_from.size() != vars_int.size() ||
                vars_int_to.size()   != vars_int.size()) {
                return error_json_op("hyperflint",
                    "vars_int_from/vars_int_to length mismatch with vars_int");
            }
        }

        bool introduce_al = false;
        {
            std::regex re("\"algebraic_letters\"\\s*:\\s*(true|false)");
            std::smatch m;
            if (std::regex_search(body, m, re)) {
                introduce_al = (m[1] == "true");
            }
        }
        // Always clear the singleton at entry — required for correctness
        // under in-process transport (LibraryLink).  Cheap when unused.
        hyperflint::AlgebraicLetterTable::global().clear();

        // R24 rev 2 / chain 17 — defuse pointer-reuse latent UB on the
        // (ctx*, rhs_str*) caches.  When a `PolyCtx` is destroyed at
        // call boundary, its raw pointer becomes dangling but cache
        // entries keyed by it survive.  A future `PolyCtx` placed at
        // the same heap address by the allocator would produce a key
        // collision and return a `Rat` valued in the OLD ctx (UB,
        // crashes possible, silent wrong answers possible).  Clearing
        // at handler entry costs one ms-scale wide-ctx parse pass on
        // first miss per call but eliminates the hazard entirely.  Also
        // reset the narrow-ctx flag so a flag set by a previous failed
        // call (in-process LibraryLink) cannot poison this call.
        hyperflint::clear_rhs_cache();
        hyperflint::clear_linear_factors_cache();
        hyperflint::reset_narrow_ctx_flag();
        // Lazy-sum risk (b), 2026-06-13: same in-process poison hazard for
        // the nonlinear-denominator flag (integration_step also resets it
        // per call, but a parse-fail path may never reach integration_step).
        hyperflint::reset_nonlinear_den_flag();
        // PHASE_4 round-3 BLOCKER fix (2026-05-28): same hazard as
        // clear_rhs_cache. The (ctx*, table*) cache populated by
        // apply_mzv_reductions's no-op guard must be cleared at
        // bridge entry so stale entries keyed by a now-destroyed
        // PolyCtx address cannot collide with a freshly-allocated
        // PolyCtx at the same heap address.
        hyperflint::clear_ctx_has_no_lhs_cache();

        std::string data_path = resolve_mzv_data_path(body);
        // Phase γ.2: cached MZV reduction table.  First call loads;
        // subsequent calls in the same process reuse.  ~15 ms → ~0.05
        // ms on the Phase-γ micro-bench.  CLI invocation (one process
        // per call) sees no benefit but also no regression.
        const hyperflint::MzvReductionTable& table = get_cached_mzv_table(data_path);

        // HF basis-ctx campaign (PHASE_2 iter 10): env-gated opt-in to
        // the slim-ctx path. When HF_USE_BASIS_CTX=1, build a slim ctx
        // (basis + user_vars, no LHS) AND install the active expansion
        // table so to_mzv_one_word's arm-2 fires at every mint. Mutex
        // with HF_NARROW_CTX (slim is a strictly tighter narrow); mutex
        // with introduce_al (algebraic letters retain wide ctx per
        // design §5.2 A-2 carve-out).
        const char* slim_env = std::getenv("HF_USE_BASIS_CTX");
        const bool use_slim = slim_env && *slim_env && slim_env[0] != '0'
            && !introduce_al;
        const hyperflint::MzvExpansionTable* slim_exp = nullptr;
        if (use_slim) {
            // Cache the expansion table as a process singleton so
            // repeated bridge calls reuse the same instance (the
            // basis_ctx inside is shared_ptr-held; safe across calls).
            // C1 round-3 advisory fold (2026-05-28): static-init is
            // thread-safe (C++11) but the rebuild-on-data_path-change
            // block is not atomic. Mutex-protect for defence in depth
            // even though LibraryLink documents main-thread-only entry.
            static std::mutex                                  cached_exp_mu;
            static std::unique_ptr<hyperflint::MzvExpansionTable> cached_exp;
            static std::string                                 cached_exp_path;
            {
                std::lock_guard<std::mutex> lk(cached_exp_mu);
                if (!cached_exp || cached_exp_path != data_path) {
                    cached_exp = std::make_unique<hyperflint::MzvExpansionTable>(
                        hyperflint::load_mzv_expansion(data_path));
                    cached_exp_path = data_path;
                }
                slim_exp = cached_exp.get();
            }
        }

        // R20 Route (i): narrow per-call PolyCtx, env-gated for opt-in.
        // Discovers actually-touched MZV symbols by scanning integrand
        // string + transitive closure on reduction graph.  ~3-7x per-term
        // RSS reduction at the cost of ~ms-scale discovery.  Only safe
        // when introduce_al is False (algebraic letters need full pool).
        std::vector<std::string> base_vars;
        std::string expr_str_for_discovery = json_str_field(body, "expr");
        const char* narrow_env = HF_FLAG_NARROW_CTX;
        const bool use_narrow = narrow_env && *narrow_env && narrow_env[0] != '0'
            && !introduce_al && !expr_str_for_discovery.empty()
            && !use_slim;  // slim takes precedence
        if (hyperflint::period_tuples_enabled() && !introduce_al) {
            // Period-tuples Phase 2 (spec 2026-06-04 §2.1): SLIM ctx --
            // kinematic vars only; period content lives structurally in
            // SymMonomial::period_powers via the scratch-ring mint.
            // introduce_al (Wm/Wp) requests fall back to the legacy wide
            // path below (pilot exclusion, spec §2.5).
            base_vars = user_vars;
        } else if (use_slim) {
            base_vars = hyperflint::build_basis_var_list(*slim_exp, user_vars);
        } else if (use_narrow) {
            base_vars = hyperflint::build_narrow_var_list(
                table, user_vars, expr_str_for_discovery);
        } else {
            base_vars = introduce_al
                ? hyperflint::build_full_var_list(table, user_vars)
                : hyperflint::build_mzv_var_list(table, user_vars);
        }
        // Period-tuples Phase 0 Task 0.2: kinematic prefix length for the
        // density census (atoms occupy indices >= user_vars.size()).
        hyperflint::addpf_probe::set_user_var_count(user_vars.size());
        // HF_CTX_PAD_VARS (period-tuples Phase 0 falsifier, 2026-06-04):
        // append N synthetic never-referenced variables to measure the
        // width-attributable share of RSS/wall directly (same fixture,
        // PAD=0 vs PAD=N, diff peak RSS). Values must be byte-identical
        // since the pad vars never appear in any polynomial. Default OFF.
        if (const char* pad_env = std::getenv("HF_CTX_PAD_VARS")) {
            const long npad = std::strtol(pad_env, nullptr, 10);
            for (long i = 0; i < npad; ++i)
                base_vars.push_back("hfpad_" + std::to_string(i));
        }
        // RAII guard: install active expansion for the lifetime of this
        // bridge call. Mint site reads it via get_active_mzv_expansion()
        // as a fallback when no explicit expansion parameter is provided.
        std::unique_ptr<hyperflint::ActiveMzvExpansionScope> exp_scope;
        if (use_slim) {
            exp_scope = std::make_unique<
                hyperflint::ActiveMzvExpansionScope>(slim_exp);
            // PHASE_3 / MF-3: scan the full body for any LHS / out-of-
            // table MZV token before parsing begins. Cheap (~µs against
            // a 3KB JSON). Throws with clear site identification.
            hyperflint::assert_no_lhs_tokens(body, slim_exp,
                                              "hyperflint_sym");
        }

        std::unique_ptr<hyperflint::PolyCtx> ctx_holder;
        std::vector<std::string> vars;
        hyperflint::ShuffleList input;
        // LAZY-SUM (HF_LAZY_SUM, default off): when on AND the expr parses
        // to a top-level Plus, integrate each addend separately and sum the
        // result tables -- the R-class parse-fusion cure. lazy_inputs holds
        // one ShuffleList per top-level addend; empty => the single fused
        // `input` path (historical behaviour).
        std::vector<hyperflint::ShuffleList> lazy_inputs;
        bool lazy_sum = false;
        {
            const char* e = std::getenv("HF_LAZY_SUM");
            lazy_sum = (e && e[0] && std::string(e) != "0");
        }

        std::string f_str    = json_str_field(body, "f");
        std::string expr_str = json_str_field(body, "expr");
        {
            int given = (!expr_str.empty() ? 1 : 0) + (!f_str.empty() ? 1 : 0)
                + (!extract_top_array(body, "wordlist").empty() ? 1 : 0);
            if (given > 1) {
                std::cerr << "hyperflint: warning: multiple input forms "
                             "provided (expr/f/wordlist); expr > f > wordlist "
                             "priority applied\n";
            }
        }

        if (!expr_str.empty()) {
            try {
                auto parsed = hyperflint::convert::parse_expression(
                    expr_str, base_vars, lazy_sum);
                vars = std::move(parsed.augmented_vars);
                const bool lazy_split =
                    lazy_sum
                    && parsed.expr.kind() == hyperflint::convert::ExprKind::Plus
                    && parsed.expr.num_children() > 1;
                if (std::getenv("HF_LAZY_SUM_DEBUG")) {
                    std::cerr << "[lazy-sum] lazy_sum=" << lazy_sum
                              << " expr_kind=" << (int)parsed.expr.kind()
                              << " num_children=" << parsed.expr.num_children()
                              << " lazy_split=" << lazy_split << "\n";
                }
                if (lazy_split) {
                    // Convert EACH top-level addend into its own ShuffleList
                    // in the shared ctx. No fused Rat::add, no merge across
                    // addends -- each addend is integrated on its own and the
                    // result tables are summed downstream.
                    for (size_t c = 0; c < parsed.expr.num_children(); ++c) {
                        hyperflint::Regulator reg =
                            hyperflint::convert::convert_to_hlog_reg_inf(
                                parsed.expr.child(c), *parsed.ctx);
                        hyperflint::ShuffleList si;
                        si.reserve(reg.size());
                        for (auto& t : reg) {
                            si.push_back(hyperflint::ShuffleEntry{
                                std::move(t.coef), std::move(t.key)});
                        }
                        lazy_inputs.push_back(std::move(si));
                    }
                } else {
                    hyperflint::Regulator reg =
                        hyperflint::convert::convert_to_hlog_reg_inf(
                            parsed.expr, *parsed.ctx);
                    input.reserve(reg.size());
                    for (auto& t : reg) {
                        input.push_back(hyperflint::ShuffleEntry{
                            std::move(t.coef), std::move(t.key)});
                    }
                }
                ctx_holder = std::move(parsed.ctx);
            } catch (const hyperflint::convert::ConvertFailed& e) {
                std::ostringstream o;
                o << "{\"op\":\"hyperflint\",\"failed\":true"
                  << ",\"reason\":\"" << json_escape(e.what()) << "\""
                  << ",\"vars\":[";
                for (size_t i = 0; i < base_vars.size(); ++i) {
                    if (i) o << ",";
                    o << "\"" << json_escape(base_vars[i]) << "\"";
                }
                o << "]}";
                return o.str();
            } catch (const hyperflint::convert::ParseError& e) {
                std::ostringstream o;
                o << "{\"op\":\"hyperflint\",\"failed\":true"
                  << ",\"reason\":\"" << json_escape(e.what()) << "\""
                  << ",\"vars\":[";
                for (size_t i = 0; i < base_vars.size(); ++i) {
                    if (i) o << ",";
                    o << "\"" << json_escape(base_vars[i]) << "\"";
                }
                o << "]}";
                return o.str();
            }
        } else {
            vars = base_vars;
            ctx_holder = std::make_unique<hyperflint::PolyCtx>(vars);
            hyperflint::PolyCtx& ctx_ref = *ctx_holder;
            if (!f_str.empty()) {
                // A3 (HF perf campaign — stay-factored lever). For a BARE
                // single-integration-variable integrand of shape
                // (NUM)^p/(DEN)^q whose DEN is linear in that variable (the
                // LR / qbox face shape), parse via FactoredRat::parse so the
                // expanded DEN^q (the parse-time pow_fps wall) is NEVER formed,
                // and carry the factored denominator as a side-channel; the
                // integrand `coef` then holds the NUMERATOR only. integrate_ii
                // dispatches the first partial-fractions call to
                // partial_fractions_factored_den (which itself safe-degrades
                // for any non-matching shape, so the result is always correct).
                // Conservative gating: single integration variable, no variable
                // rescaling (have_ranges false), and the factored parse yields
                // exactly one denominator factor linear in the integration var.
                // Conservative: skip A3 when a divergence check is requested
                // (check_divergences_pass would integrate the numerator-only
                // entry without threading the side-channel). All A3 targets run
                // check_divergences = false.
                bool cd_requested = false;
                {
                    std::regex re_cd(
                        "\"check_divergences\"\\s*:\\s*(true|false)");
                    std::smatch m_cd;
                    if (std::regex_search(body, m_cd, re_cd))
                        cd_requested = (m_cd[1] == "true");
                }
                // A3 kill-switch (HF_DISABLE_A3=1): force the ordinary
                // Rat::parse path, for A/B measurement and as a safety opt-out.
                static const bool a3_disabled = [] {
                    const char* s = std::getenv("HF_DISABLE_A3");
                    return s && s[0] && s[0] != '0';
                }();
                bool routed_factored = false;
                if (!a3_disabled && vars_int.size() == 1 && !have_ranges &&
                    !cd_requested) {
                    size_t vidx = 0;
                    for (; vidx < vars.size(); ++vidx)
                        if (vars[vidx] == vars_int[0]) break;
                    if (vidx < vars.size()) {
                        try {
                            hyperflint::FactoredRat fr =
                                hyperflint::FactoredRat::parse(ctx_ref, f_str);
                            const auto& dfs = fr.den_factors();
                            if (dfs.size() == 1 &&
                                dfs[0].base.degree_in_var(vidx) == 1) {
                                hyperflint::Rat num_only{
                                    hyperflint::Poly(fr.numerator())};
                                hyperflint::ShuffleEntry se{
                                    std::move(num_only), {}, fr};
                                input.push_back(std::move(se));
                                routed_factored = true;
                            }
                        } catch (const std::exception&) {
                            routed_factored = false;  // fall through to Rat::parse
                        }
                    }
                }
                if (!routed_factored) {
                    // HF MZV-rewrite C-prep.4 iter-32 F2 parse-boundary
                    // safety rail (top-level f_str entry).
                    hyperflint::Rat f_r = hyperflint::Rat::parse(ctx_ref, f_str);
                    debug_check_parse_idempotent(ctx_ref, f_r,
                                                  "main_handler/f_str");
                    input.push_back(
                        hyperflint::ShuffleEntry{std::move(f_r), {}});
                }
            } else {
                input = parse_shuffle_list(ctx_ref, body, "wordlist");
            }
        }
        hyperflint::PolyCtx& ctx = *ctx_holder;

        std::vector<size_t> var_indices;
        for (const auto& vi : vars_int) {
            size_t idx = 0;
            for (; idx < vars.size(); ++idx) if (vars[idx] == vi) break;
            var_indices.push_back(idx);
        }

        // Unify the single-input and lazy-split paths: integrate every
        // entry of `addend_inputs` and sum the result tables. Non-lazy =>
        // one entry (the fused `input`).
        std::vector<hyperflint::ShuffleList> addend_inputs;
        if (lazy_inputs.empty())
            addend_inputs.push_back(std::move(input));
        else
            addend_inputs = std::move(lazy_inputs);

        if (have_ranges) {
            for (auto& ai : addend_inputs) {
                for (size_t k = 0; k < vars_int.size(); ++k) {
                    const std::string& from = vars_int_from[k];
                    const std::string& to   = vars_int_to[k];
                    if (from == "0" && (to == "Infinity" || to == "+Infinity" ||
                                         to == "oo"))
                        continue;
                    ai = hyperflint::rescale_interval(ctx, ai, var_indices[k],
                                                       from, to);
                    if (ai.empty()) break;
                }
            }
        }

        // DP.3 (divergence policy, 2026-06-03): the intended bare-request
        // default is TRUE, but the flip is BLOCKED on HF-DIVCHECK-PARITY
        // (the scan false-positives on generic multi-pole convergent
        // integrands) -- see bridge/cli/main.cpp and
        // notes/hf_divcheck_parity.md. Flip only after parity.
        bool check_div = false;
        {
            std::regex re("\"check_divergences\"\\s*:\\s*(true|false)");
            std::smatch m;
            if (std::regex_search(body, m, re)) check_div = (m[1] == "true");
        }

        // LAZY-SUM x check_divergences guard (adversarial review 2026-06-13,
        // fix B): the lazy split integrates each top-level addend separately,
        // and R-class addends are INDIVIDUALLY divergent by construction
        // (only their sum is finite). A per-addend divergence check would
        // therefore spuriously report the whole face divergent. Refuse the
        // combination loudly rather than mis-report. (addend_inputs.size()>1
        // == the lazy split path fired; non-lazy is always a single input.
        // Note lazy_inputs has already been move()d into addend_inputs.)
        if (addend_inputs.size() > 1 && check_div) {
            std::ostringstream o;
            o << "{\"op\":\"hyperflint\",\"failed\":true,\"reason\":\""
              << "HF_LAZY_SUM is incompatible with check_divergences: "
              << "top-level addends are individually divergent by "
              << "construction; the per-addend check cannot certify the "
              << "finite sum. Disable one of the two.\",\"vars\":[";
            for (size_t i = 0; i < vars.size(); ++i) {
                if (i) o << ",";
                o << "\"" << json_escape(vars[i]) << "\"";
            }
            o << "]}";
            return o.str();
        }

        // HF MZV-rewrite C-prep.4 (iter-27): opt-in canonical-emission
        // writer. Default false preserves existing JSON-byte output.
        bool canonical_emission = false;
        {
            std::regex re("\"canonical_emission\"\\s*:\\s*(true|false)");
            std::smatch m;
            if (std::regex_search(body, m, re))
                canonical_emission = (m[1] == "true");
        }

        // DP.3 spectator projection (2026-06-03): user variables that are
        // never integrated (free kinematic parameters) feed the
        // divergence check's fibration-basis zero test. Built from
        // base_vars (the user-supplied "vars" list, BEFORE atom-pool
        // augmentation) so MZV constants / Log2 / minted letters are
        // never projected over. Cheap: only computed when the check is
        // armed; empty when every user var is integrated.
        std::vector<size_t> spectator_idx;
        if (check_div) {
            for (const auto& bv : base_vars) {
                bool is_int = false;
                for (const auto& vi : vars_int)
                    if (vi == bv) { is_int = true; break; }
                if (is_int) continue;
                for (size_t i = 0; i < vars.size(); ++i)
                    if (vars[i] == bv) { spectator_idx.push_back(i); break; }
            }
        }

        try {
            auto _bench_t0 = std::chrono::steady_clock::now();
            hyperflint::RegulatorSym out;
            if (addend_inputs.size() == 1) {
                out = hyperflint::hyperflint_sym(
                    ctx, addend_inputs[0], var_indices, table, introduce_al,
                    check_div, spectator_idx);
            } else {
                // LAZY-SUM: integrate each top-level addend, then sum the
                // result tables. Each addend is its own integration; the
                // heavy parse/PF fusion of the denominator-disjoint addends
                // never forms. canonicalize_regulator_sym merges same-key
                // terms (SymCoef add) AND sorts -- byte-matching the fused
                // path's terminal canonicalize, so the default emitter output
                // is order-canonical, not addend-concatenation order
                // (adversarial review 2026-06-13, fix A).
                //
                // LAZY-SUM risk (b), 2026-06-13: if an individual addend is
                // not integrable in the recorded order (degree>=3 nonlinear
                // denominator factor), hyperflint_sym throws
                // NonlinearDenominatorUnsupported (a typed IntegrationStepFailed
                // raised post-barrier by integration_step -- NOT a terminate;
                // see the Phase-1 flag-and-rethrow there). It propagates to
                // the outer handler catch below, which emits {"failed": true}
                // -- so the sweep RECORDS the face and CONTINUES (the runner
                // only halts on signal/crash), instead of the SIGABRT that
                // halted v8 on ord_0_face_68. A fused-fallback that RECOVERS
                // the face value (re-integrating the whole face fused, where
                // the cancellation makes the denominator tractable) is the
                // intended Phase-2 enhancement; it is NOT wired here yet
                // because a first attempt hit a >300 GB re-integration blowup
                // in the in-handler re-parse path that is not yet root-caused
                // (the same fused expr integrates in ~32 MB when run lazy-off
                // from a clean process). Tracked in
                // notes/hf_tree_merge/PLAN_lazy_sum_and_blinding.md.
                for (auto& ai : addend_inputs) {
                    hyperflint::RegulatorSym oi = hyperflint::hyperflint_sym(
                        ctx, ai, var_indices, table, introduce_al, check_div,
                        spectator_idx);
                    for (auto& term : oi) out.push_back(std::move(term));
                }
                out = hyperflint::canonicalize_regulator_sym(out);
            }
            auto _bench_t1 = std::chrono::steady_clock::now();
            double compute_s =
                std::chrono::duration<double>(_bench_t1 - _bench_t0).count();
            // Round-19 wide-ctx probe dump.
            if (hyperflint::ctx_probe_enabled()) {
                std::cerr << hyperflint::ctx_probe_dump_and_clear(ctx)
                          << "\n";
            }
            std::ostringstream o;
            o << "{\"op\":\"hyperflint\",\"result\":"
              << (canonical_emission
                      ? emit_regulator_sym_canonical(out)
                      : emit_regulator_sym(out))
              << ",\"timing_compute_s\":" << compute_s
              << ",\"vars\":[";
            for (size_t i = 0; i < vars.size(); ++i) {
                if (i) o << ",";
                o << "\"" << json_escape(vars[i]) << "\"";
            }
            o << "]";
            // Phase 7-vi-b: emit the algebraic-letter table when the
            // caller turned it on.  Always emit (even if empty) so the
            // Mma-side parser has a stable response shape.
            if (introduce_al) {
                o << ",\"algebraic_letters\":" << emit_algebraic_letter_table();
            }
            o << "}";
            return o.str();
        } catch (const hyperflint::NarrowCtxTooNarrow&) {
            // R24 rev 2 / chain 17 — narrow ctx is missing some MZV
            // variable referenced at runtime.  Emit a structured-error
            // JSON so the Mma-side `STHyperFlint` can detect it and
            // re-issue `RunProcess` with `HF_NARROW_CTX=0`
            // (out-of-process retry — no in-process recursion, no
            // setenv MT-unsafety).  R1: cleanup FLINT pools so the
            // failed call doesn't leak per-thread arenas across
            // LibraryLink retries.
            flint_cleanup_master();
            std::ostringstream o;
            o << "{\"op\":\"hyperflint\",\"narrow_ctx_insufficient\":true"
              << ",\"vars\":[";
            for (size_t i = 0; i < vars.size(); ++i) {
                if (i) o << ",";
                o << "\"" << json_escape(vars[i]) << "\"";
            }
            o << "]}";
            return o.str();
        } catch (const hyperflint::HyperFLINTDivergentIntegral& e) {
            // R26 R1 -- release FLINT thread-pool arenas on the failure
            // path so subsequent LibraryLink calls don't see leaked RSS.
            flint_cleanup_master();
            std::ostringstream o;
            o << "{\"op\":\"hyperflint\",\"divergent\":true"
              << ",\"reason\":\"" << json_escape(e.what()) << "\""
              << ",\"vars\":[";
            for (size_t i = 0; i < vars.size(); ++i) {
                if (i) o << ",";
                o << "\"" << json_escape(vars[i]) << "\"";
            }
            o << "]}";
            return o.str();
        } catch (const hyperflint::IntegrationStepFailed&) {
            flint_cleanup_master();  // R26 R1 -- see above.
            std::ostringstream o;
            o << "{\"op\":\"hyperflint\",\"failed\":true,\"vars\":[";
            for (size_t i = 0; i < vars.size(); ++i) {
                if (i) o << ",";
                o << "\"" << json_escape(vars[i]) << "\"";
            }
            o << "]}";
            return o.str();
        }
    } catch (const std::exception& e) {
        flint_cleanup_master();  // R26 R1 -- catch-all failure path.
        return error_json_op("hyperflint", e.what());
    } catch (...) {
        flint_cleanup_master();  // R26 R1 -- unknown-exception failure path.
        return error_json_op("hyperflint", "unknown exception");
    }
}

}  // namespace handlers
}  // namespace hyperflint
