# GiNaC G(a;z) oracle bridge — shells out to a tiny C++ CLI (ginac_gpl, built
# on first use from ginac_gpl.cpp; needs g++ and the GiNaC/CLN libraries).
# The binary links GPL libraries and is never distributed; this package only
# runs it as a subprocess.  This is the field-standard arbitrary-precision
# reference; used in tests.

const _GINAC_DIR = @__DIR__
const _GINAC_BIN = joinpath(_GINAC_DIR, "ginac_gpl")

function _ensure_ginac_built()
    isfile(_GINAC_BIN) && return _GINAC_BIN
    src = joinpath(_GINAC_DIR, "ginac_gpl.cpp")
    flags = try strip(read(`pkg-config --cflags --libs ginac`, String))
            catch; "-lginac -lcln" end
    cmd = `bash -c "g++ -O2 -std=c++17 -o $(_GINAC_BIN) $(src) $flags"`
    run(cmd)
    return _GINAC_BIN
end

_ratstr(x::Rational) = string(numerator(x)) * "/" * string(denominator(x))
_ratstr(x::Integer)  = string(x)
function _ratstr(x::Real)
    bx = big(x)
    return Base.MPFR.string_mpfr(bx, "%.60Re")
end
_ratstr(x::BigFloat) = iszero(x) ? "0" : Base.MPFR.string_mpfr(x, "%.$(precision(x) ÷ 3)Re")

# Real inputs MUST hand GiNaC an exact-integer 0 imaginary part: a CLN float
# zero (e.g. "0.0e+00") leaves the letter as `re + 0.0*I`, and GiNaC's G then
# refuses to pick an i0 side for an on-path letter (evalf returns symbolic G).
_re_im(x) = x isa Complex ? (real(x), imag(x)) : (x, 0)

"""
    gpl_ginac(a, z; digits) -> Complex{BigFloat}

Reference value of G(a;z) from GiNaC at `digits` decimal digits. Letters and z
should be exact (`Integer`, `Rational`, or `Complex` thereof) for best fidelity.
"""
function gpl_ginac(a::AbstractVector, z; digits::Int = 60)
    bin = _ensure_ginac_built()
    args = String[string(digits)]
    zr, zi = _re_im(z); push!(args, _ratstr(zr), _ratstr(zi))
    for ai in a
        ar, aim = _re_im(ai); push!(args, _ratstr(ar), _ratstr(aim))
    end
    out = read(`$bin $args`, String)
    lines = split(strip(out), '\n')
    setprecision(BigFloat, ceil(Int, digits * 3.33) + 64) do
        Complex(parse(BigFloat, replace(lines[1], "E" => "e")),
                parse(BigFloat, replace(lines[2], "E" => "e")))
    end
end
