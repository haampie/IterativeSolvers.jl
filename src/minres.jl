export symmlq_iterable

import Base.LinAlg.BLAS.axpy!
import Base: start, next, done

type MinResIterable
    A
    
    # Vectors
    b
    x
    v_next
    v_curr
    v_prev
    
    w₀
    w₁
    w₂

    β₁
    β₂

    # Given's rotations
    cos
    sin

    ρ
    reltol
    maxiter
end

converged(s::MinResIterable) = abs(s.ρ) ≤ s.reltol

done(s::MinResIterable, iteration::Int) = iteration > s.maxiter || converged(s)

start(::MinResIterable) = 0

function next(s::MinResIterable, iteration::Int)
    # v_next = A * v_next
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

    # Search vector
    copy!(s.w₁, s.w₂)
    axpy!(-l1, s.w₀, s.w₁)

    copy!(s.w₂, s.v_curr)
    axpy!(-l2, s.w₀, s.w₂)
    copy!(s.w₀, s.w₁)
    scale!(inv(l0), s.w₀)

    # Update x
    axpy!(s.ρ * s.cos, s.w₀, s.x)

    s.ρ *= s.sin

    # Implicitly computed residual & iteration
    abs(s.ρ), iteration + 1
end

function minres_iterable(A, b; maxiter = min(20, size(A, 1)), tol = sqrt(eps(real(eltype(b)))))
    T = eltype(b)
    x = zeros(b)
    resnorm = norm(b)

    v_curr = copy(b)
    scale!(v_curr, inv(resnorm))
    v_next = zeros(b)
    v_prev = zeros(b)

    w₀ = zeros(b)
    w₁ = zeros(b)
    w₂ = copy(v_curr)

    givens_c = -one(T)
    givens_s = zero(T)

    reltol = resnorm * tol

    MinResIterable(
        A, b, x, v_next, v_curr, v_prev, 
        w₀, w₁, w₂,
        zero(T), zero(T),
        givens_c, givens_s,
        resnorm, reltol, maxiter
    )
end