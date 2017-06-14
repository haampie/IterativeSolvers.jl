export improved_cg

function improved_cg(A, b; log = true, maxiter = 1500, tol = sqrt(eps(real(eltype(b)))))
    T = eltype(b)
    n = size(A, 1)

    x = zeros(T, n) # Initial "guess"
    c = zeros(T, n)
    u = zeros(T, n)
    r = zeros(T, n)
    copy!(r, b) # Initial residual vector
    ρ = one(T)
    
    iter = 0

    last_residual = norm(r)
    reltol = last_residual * tol

    if log
        history::Vector{real(T)} = []
    end

    while last_residual > reltol && iter < maxiter
        # Preconditioner K = I
        copy!(c, r)

        σ = -ρ
        ρ = dot(c, r)
        β = ρ / σ

        # u ← r - βu
        scale!(u, -β)
        axpy!(one(T), r, u)

        # c = A * u
        A_mul_B!(c, A, u)

        α = ρ / dot(c, u)
    
        axpy!(α, u, x)
        axpy!(-α, c, r)

        iter += 1
        last_residual = norm(r)
        if log
            push!(history, last_residual)
        end
    end

    log ? (x, history) : x
end