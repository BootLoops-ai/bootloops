# verify_cabi.jl — post-build acceptance for deps/libhf_cabi_static_f2.so.
# Checks, all in ONE process (that is the
# point — the F2 bug class only bites on in-process reuse):
#   1. hf_version_string() == "1.2.8" (strict; source builds mis-stamp 1.2.0)
#   2. response envelope carries "hf_version":"1.2.8" (double check)
#   3. unary-minus canary through the in-process parser (pfrac of -x^2)
#   4. F2 pair-repro sequence (1-var pfrac -> 2-var pfrac -> identical 2-var
#      pfrac, one process) — the unpatched lib dies or corrupts here;
#      repeat responses must be byte-identical and error-free.
# Run: ulimit -v 32505856; julia deps/verify_cabi.jl
using Libdl

const LIB = joinpath(@__DIR__, "libhf_cabi_static_f2.so")
h = Libdl.dlopen(LIB)   # RTLD_LOCAL default: hf internals stay private
pfree = Libdl.dlsym(h, :hf_free_string)

function abi(op::Symbol, req::String)
    p = ccall(Libdl.dlsym(h, op), Ptr{UInt8}, (Cstring,), req)
    p == C_NULL && error("$op: NULL return (catastrophic alloc failure)")
    s = unsafe_string(p)
    ccall(pfree, Cvoid, (Ptr{UInt8},), p)
    return s
end

fails = String[]
chk(name, ok) = ok ? println("PASS  $name") : (push!(fails, name); println("FAIL  $name"))

# 1. static version string — strict equality (source builds stamp 1.2.0 by default)
ver = unsafe_string(ccall(Libdl.dlsym(h, :hf_version_string), Cstring, ()))
chk("version_string == 1.2.8 (got $ver)", ver == "1.2.8")

# 2+3. unary-minus canary through the in-process parser + envelope stamp
r = abi(:hf_partial_fractions,
        """{"op":"partial_fractions","f":"-x^2","vars":["x"],"var":"x"}""")
chk("envelope hf_version 1.2.8", occursin("\"hf_version\":\"1.2.8\"", r))
chk("unary-minus canary (-x^2 stays -x^2)",
    occursin("-x^2", r) && !occursin("\"error\"", r))

# 4. F2 pair-repro sequence, one process
r1 = abi(:hf_partial_fractions,
         """{"op":"partial_fractions","f":"1/(x*(1 + x))","vars":["x"],"var":"x"}""")
r2 = abi(:hf_partial_fractions,
         """{"op":"partial_fractions","f":"1/(x*(1 + x)*(x - y))","vars":["x","y"],"var":"x"}""")
r3 = abi(:hf_partial_fractions,
         """{"op":"partial_fractions","f":"1/(x*(1 + x)*(x - y))","vars":["x","y"],"var":"x"}""")
chk("F2 pair sequence error-free",
    all(!occursin("\"error\"", x) for x in (r1, r2, r3)))
chk("F2 repeat byte-identical", r2 == r3)

# linear_factors face of F2 (segv face on unpatched lib), same process
l1 = abi(:hf_linear_factors,
         """{"op":"linear_factors","poly":"2*x^2 + 3*x*y + y^2","var":"x","vars":["x","y"]}""")
l2 = abi(:hf_linear_factors,
         """{"op":"linear_factors","poly":"2*x^2 + 3*x*y + y^2","var":"x","vars":["x","y"]}""")
chk("lf pair error-free + byte-identical",
    !occursin("\"error\"", l1) && l1 == l2)

isempty(fails) ? println("VERIFY_CABI: ALL PASS") :
    (println("VERIFY_CABI: FAILURES: ", fails); exit(1))
