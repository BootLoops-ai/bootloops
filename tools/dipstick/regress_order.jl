# DIPSTICK member — `order` regression (7/7: six jet controls + the 3L-banana
# validate case).
# Usage: julia --project=<your Oscar/Singular project> regress_order.jl <pf_rank_file>
# Run it against BOTH routes: tools/dipstick/pf_rank.jl (member) and
# tools/pf_rank.jl (shim) — results must be identical.
# Runs the 6 jet controls (verbatim targets) against
# pf_probe/pf_probe_gen/pf_probe_esc — the drain helpers below are
# reference copies — plus c7 = the pf_rank_validate.jl 3L-banana case
# (linear-path preservation: order 3, holonomic_rank 15).
pfr = abspath(ARGS[1])
println("# tool under test: ", basename(pfr))
include(pfr)

# ---------- drain helpers (reference copies; needed by c5/c6) -----------------
_degv(P, x) = iszero(P) ? 0 : maximum(e[Oscar.var_index(x)]
                                      for e in Nemo.exponent_vectors(P); init=0)

function _oddpart(P, R, znames; log)
    fac = Oscar.factor(P)
    keep = one(R)
    for (f, e) in fac
        zdep = any(_degv(f, g) > 0 for g in Oscar.gens(R) if string(g) in znames)
        if iseven(e)
            push!(log, Dict("drop"=>"even^$(e)", "factor"=>string(f)))
        elseif !zdep
            push!(log, Dict("drop"=>"odd-pure-kin (algebraic prefactor sqrt)",
                            "factor"=>string(f)))
        else
            keep *= f^(isone(e) ? 1 : e)
        end
        if isodd(e) && e > 1 && zdep
            keep = Oscar.divexact(keep, f^(e-1))
            push!(log, Dict("drop"=>"square-part of odd^$(e)", "factor"=>string(f)))
        end
    end
    keep
end

function drain(B, R, znames)
    log = Any[]; chain = String[]
    P = _oddpart(B, R, znames; log=log)
    while true
        live = [g for g in Oscar.gens(R) if string(g) in znames && _degv(P, g) > 0]
        length(live) <= 1 && break
        cand = [z for z in live if _degv(P, z) == 2]
        if isempty(cand)
            push!(log, Dict("stall"=>"no deg-2 candidate",
                            "live"=>[string(z) for z in live],
                            "degs"=>[_degv(P, z) for z in live]))
            break
        end
        best = nothing
        for z in cand
            a = Oscar.divexact(Oscar.derivative(Oscar.derivative(P, z), z), R(2))
            b = Oscar.derivative(P, z) - 2*a*z
            c = P - a*z^2 - b*z
            D = b^2 - 4*a*c
            nxt = iszero(D) ? a : D
            sublog = Any[]
            nxt = iszero(nxt) ? nxt : _oddpart(nxt, R, znames; log=sublog)
            sc = length(Nemo.exponent_vectors(nxt) |> collect)
            if best === nothing || sc < best[1]
                best = (sc, z, nxt, sublog, iszero(D))
            end
        end
        sc, z, nxt, sublog, wasdeg = best
        push!(chain, string(z))
        push!(log, Dict("drain"=>string(z), "disc_zero(perfect square)"=>wasdeg,
                        "n_terms_after"=>sc))
        append!(log, sublog)
        if iszero(nxt) || sc == 0
            push!(log, Dict("note"=>"residual vanished/constant — fully drained"))
            P = nxt
            break
        end
        P = nxt
    end
    live = [string(g) for g in Oscar.gens(R) if string(g) in znames && _degv(P, g) > 0]
    P, chain, live, log
end

# ---------- controls (c1..c6 jet controls; c7 = validate case) ---------------
npass = 0; nfail = 0
function tally(name, cond)
    global npass, nfail
    cond ? (npass += 1) : (nfail += 1)
    println(name, " ", cond ? "PASS" : "FAIL"); flush(stdout)
end

# c1: regression — pf_probe on sunrise LP, linear-in-s v1 path (expect 2, hr 7)
G_sun = "x1*x2+x1*x3+x2*x3 + s*x1*x2*x3 - (x1+x2+x3)*(x1*x2+x1*x3+x2*x3)"
t0=time(); r1 = pf_probe(G_sun, ["s"], ["x1","x2","x3"], "s"; Lmax=3)
println("c1 sunrise v1-probe: order=", r1.pf_order, " 2p=", r1.pf_order_2p,
        " hr=", r1.holonomic_rank, " (expect 2, hr 7) ", round(time()-t0,digits=1), "s")
tally("c1", r1.pf_order==2 && r1.stable && r1.holonomic_rank==7)

# c2: pf_probe_gen on the same G (linear s through the general jet)
R2, (x1,x2,x3,s) = Oscar.polynomial_ring(Oscar.QQ, ["x1","x2","x3","s"])
G2 = x1*x2+x1*x3+x2*x3 + s*x1*x2*x3 - (x1+x2+x3)*(x1*x2+x1*x3+x2*x3)
t0=time(); r2 = pf_probe_gen(G2, R2, ["x1","x2","x3"], ["s"], "s"; Lmax=3)
println("c2 sunrise gen-jet:  order=", r2.pf_order, " 2p=", r2.pf_order_2p,
        " (expect 2) ", round(time()-t0,digits=1), "s")
tally("c2", r2.pf_order==2 && r2.stable)

# c3: Legendre pullback w^2 = z(z-1)(z-s^2) — deg_s=2, expect 2
R3, (w,z,s3) = Oscar.polynomial_ring(Oscar.QQ, ["w","z","s"])
G3 = w^2 - z*(z-1)*(z-s3^2)
t0=time(); r3 = pf_probe_gen(G3, R3, ["w","z"], ["s"], "s"; Lmax=3)
println("c3 Legendre(s^2):    order=", r3.pf_order, " 2p=", r3.pf_order_2p,
        " (expect 2) ", round(time()-t0,digits=1), "s")
tally("c3", r3.pf_order==2 && r3.stable)

# c4: w^2 = (z-s^2)(z-s^2-1) — constant period, deg_s=4, expect 1
G4 = w^2 - (z-s3^2)*(z-s3^2-1)
t0=time(); r4 = pf_probe_gen(G4, R3, ["w","z"], ["s"], "s"; Lmax=3)
println("c4 const-period:     order=", r4.pf_order, " 2p=", r4.pf_order_2p,
        " (expect 1) ", round(time()-t0,digits=1), "s")
tally("c4", r4.pf_order==1 && r4.stable)

# c5: drain mechanics — nested quadric, known chain + order 1
R5, (z1,z2,z3,s5) = Oscar.polynomial_ring(Oscar.QQ, ["z1","z2","z3","s"])
B5 = z1^2 - (z2^2 - (z3^2 - s5*(s5+1)))
P5, ch5, zr5, _ = drain(B5, R5, ["z1","z2","z3"])
R5f, (w5,z3f,s5f) = Oscar.polynomial_ring(Oscar.QQ, ["w__","z3","s"])
img5 = elem_type(R5f)[R5f(0), R5f(0), z3f, s5f]
G5 = w5^2 - Oscar.evaluate(P5, img5)
r5 = pf_probe_gen(G5, R5f, ["w__","z3"], ["s"], "s"; Lmax=3)
println("c5 drain+probe:      chain=", ch5, " residual=", zr5, " order=",
        r5.pf_order, " (expect [z1,z2]->z3, 1)")
tally("c5", zr5 == ["z3"] && r5.pf_order == 1)

# c6: elliptic through the FULL drain path (disc quartic -> expect order 2)
R6, (za,zb,s6) = Oscar.polynomial_ring(Oscar.QQ, ["z1","z2","s"])
B6 = za^2 - zb*(zb-1)*(zb-s6)*(zb-2)
P6, ch6, zr6, _ = drain(B6, R6, ["z1","z2"])
R6f, (w6,zbf,s6f) = Oscar.polynomial_ring(Oscar.QQ, ["w__","z2","s"])
img6 = elem_type(R6f)[R6f(0), zbf, s6f]
G6 = w6^2 - Oscar.evaluate(P6, img6)
r6 = pf_probe_gen(G6, R6f, ["w__","z2"], ["s"], "s"; Lmax=3)
println("c6 drain->elliptic:  chain=", ch6, " order=", r6.pf_order, " (expect 2)")
tally("c6", r6.pf_order == 2 && r6.stable)

# c7: pre-existing pf_rank_validate.jl case — 3L equal-mass banana (K3),
# v1 linear path preservation: order 3, holonomic_rank 15.
G_ban3 = "(x1*x2*x3+x1*x2*x4+x1*x3*x4+x2*x3*x4) + " *
         "(x1+x2+x3+x4)*(x1*x2*x3+x1*x2*x4+x1*x3*x4+x2*x3*x4) - s*x1*x2*x3*x4"
t0=time(); r7 = pf_probe(G_ban3, ["s"], ["x1","x2","x3","x4"], "s"; Lmax=4)
println("c7 3L-banana v1:     order=", r7.pf_order, " 2p=", r7.pf_order_2p,
        " hr=", r7.holonomic_rank, " (expect 3, hr 15) ", round(time()-t0,digits=1), "s")
tally("c7", r7.pf_order==3 && r7.stable && r7.holonomic_rank==15)

println("\nRESULT: ", npass, "/", npass+nfail, " PASS")
nfail == 0 || exit(1)
