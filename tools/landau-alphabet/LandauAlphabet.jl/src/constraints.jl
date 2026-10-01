# Generic exact linear-constraint assembly over Q.
#
# A ConstraintSystem collects linear conditions  Σ_j A_{ij} c_j = b_i  on the
# rational coefficient vector c of an ansatz; `reduce_space` returns the
# affine solution space (particular solution + nullspace basis) via Nemo.

export ConstraintSystem, add_constraint!, reduce_space, solution_dim

mutable struct ConstraintSystem
    ncoeff::Int
    rows::Vector{Vector{QQFieldElem}}
    rhs::Vector{QQFieldElem}
    labels::Vector{String}
end

ConstraintSystem(n::Int) = ConstraintSystem(n, Vector{Vector{QQFieldElem}}(), QQFieldElem[], String[])

function add_constraint!(cs::ConstraintSystem, row::Vector{QQFieldElem},
                         rhs::QQFieldElem = QQ(0); label::String = "")
    length(row) == cs.ncoeff || error("row length mismatch")
    push!(cs.rows, row)
    push!(cs.rhs, rhs)
    push!(cs.labels, label)
    cs
end

"""
    reduce_space(cs) -> (particular::Vector{QQFieldElem}, null::Vector{Vector{QQFieldElem}})

Solve the system exactly. Throws if inconsistent.
"""
function reduce_space(cs::ConstraintSystem)
    m = length(cs.rows)
    A = zero_matrix(QQ, m, cs.ncoeff)
    b = zero_matrix(QQ, m, 1)
    for i in 1:m
        for j in 1:cs.ncoeff
            A[i, j] = cs.rows[i][j]
        end
        b[i, 1] = cs.rhs[i]
    end
    fl, x = can_solve_with_solution(A, b, side = :right)
    fl || error("constraint system inconsistent")
    nd, ns = nullspace(A)
    part = [x[j, 1] for j in 1:cs.ncoeff]
    nulls = [[ns[j, k] for j in 1:cs.ncoeff] for k in 1:nd]
    (part, nulls)
end

solution_dim(cs::ConstraintSystem) = length(reduce_space(cs)[2])
