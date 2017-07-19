export symmlq_iterable

import Base.LinAlg.BLAS.axpy!
import Base: start, next, done

type SymmLQIterable
    A
    
    # Vectors
    b
    x
    v_next
    v_curr
    v_prev
    w

    β₁
    β₂

    # Given's rotations
    cos
    sin
    g1
    g2

    resnorm
    reltol
    maxiter
end

converged(s::SymmLQIterable) = s.resnorm ≤ s.reltol

done(s::SymmLQIterable, iteration::Int) = iteration > s.maxiter || converged(s)

start(::SymmLQIterable) = 0

function next(s::SymmLQIterable, iteration::Int)
    A_mul_B!(s.v_next, s.A, s.v_curr)
    axpy!(-s.β₁, s.v_prev, s.v_next)

    # Orthogonalize (once)
    # v_next := (I - v_curr v_curr')v_next
    α = dot(s.v_curr, s.v_next)
    axpy!(-α, s.v_curr, s.v_next)
    s.β₁ = norm(s.v_next)

    # Move the vecs around.
    copy!(s.v_prev, s.v_curr)
    copy!(s.v_curr, s.v_next)

    # Normalize
    scale!(s.v_curr, inv(s.β₁))

    # Implicitly decompose the tridiagonal lanczos matrix
    l1 = s.sin * α - s.cos * s.β₂
    l2 = s.sin * s.β₁
    α̃ = -s.sin * s.β₂ - s.cos * α
    s.β₂ = s.cos * s.β₁
    l0 = √(α̃ ^ 2 + s.β₁ ^ 2)
    s.cos = α̃ / l0
    s.sin = s.β₁ / l0

    g = s.g2 - l1 * s.g1
    s.g2 = -l2 * s.g1
    s.g1 = g / l0

    # Update x
    axpy!(s.g1 * s.cos, s.w, s.x)
    axpy!(s.g1 * s.sin, s.v_curr, s.x)

    # Update w
    scale!(s.w, s.sin)
    axpy!(-s.cos, s.v_curr, s.w)

    # Implicitly computed residual norm
    s.resnorm = √(g ^ 2 + s.g2 ^ 2)

    s.resnorm, iteration + 1
end

function symmlq_iterable(A, b; maxiter = min(20, size(A, 1)), tol = sqrt(eps(real(eltype(b)))))
    T = eltype(b)
    x = zeros(b)
    resnorm = norm(b)

    v_curr = copy(b)
    scale!(v_curr, inv(resnorm))
    v_next = zeros(b)
    v_prev = zeros(b)
    w = copy(v_curr)

    givens_c = -one(T)
    givens_s = zero(T)
    g1 = zero(T)
    g2 = resnorm

    reltol = resnorm * tol

    SymmLQIterable(
        A, b, x, v_next, v_curr, v_prev, w,
        zero(T), zero(T),
        givens_c, givens_s,
        g1, g2,
        resnorm, reltol, maxiter
    )
end