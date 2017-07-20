export minres_iterable, minres

import Base.LinAlg.BLAS.axpy!
import Base: start, next, done

function rotation(a, b)
    r = √(a^2 + b^2)
    a / r, b / r, r
end

function minres2(A, b; maxiter = 10)
    n = length(b)

    # This'll hold the active part of the Lanczos matrix
    H = zeros(maxiter + 1, maxiter)
    T = zeros(maxiter + 1, maxiter)

    x = zeros(b)
    v_prev = similar(b)
    v_curr = copy(b)
    v_next = similar(b)
    w_prev = similar(b)
    w_curr = similar(b)
    w_next = similar(b)

    # Create the right-hand side
    rhs = zeros(maxiter + 1)
    rhs[1] = norm(v_curr)

    # Normalize the first Krylov basis vector
    v_curr /= rhs[1]

    # The normalization constant of v_prev (initially zero)
    prev_norm = 0.0

    # Givens rotations
    c_prev, s_prev = 1.0, 0.0
    c_curr, s_curr = 1.0, 0.0

    for i = 1 : maxiter
        # v_next = A * v_curr - prev_norm * v_prev
        A_mul_B!(v_next, A, v_curr)
        axpy!(-prev_norm, v_prev, v_next)
        
        # Orthogonalize w.r.t. v_curr
        T[i, i] = dot(v_curr, v_next)
        axpy!(-T[i, i], v_curr, v_next)

        # Normalize
        T[i + 1, i] = prev_norm = norm(v_next)
        
        if i < maxiter
            T[i, i + 1] = prev_norm
        end

        v_next /= prev_norm

        # Previous 2 rotations
        if i > 2
            # Use the fact that T[iter - 2, iter] = 0
            T[i - 2, i] = s_prev * T[i - 1, i]
            T[i - 1, i] = c_prev * T[i - 1, i]
        end

        if i > 1
            tmp = -s_curr * T[i - 1, i] + c_curr * T[i, i]
            T[i - 1, i] = c_curr * T[i - 1, i] + s_curr * T[i, i]
            T[i, i] = tmp
        end

        # New rotation
        c, s, r = rotation(T[i, i], T[i + 1, i])
        T[i, i] = r
        T[i + 1, i] = 0.0

        # Applied to the rhs:
        rhs[i + 1] = -s * rhs[i]
        rhs[i] = c * rhs[i]

        # Update W = V * inv(R)
        copy!(w_next, v_curr)
        
        if i > 1
            axpy!(-T[i - 1, i], w_curr, w_next)
        end

        if i > 2
            axpy!(-T[i - 2, i], w_prev, w_next)
        end
        
        scale!(w_next, inv(T[i, i]))

        # Update solution x
        axpy!(rhs[i], w_next, x)

        # Move on: next -> curr, curr -> prev
        v_prev, v_curr, v_next = v_curr, v_next, v_prev
        w_prev, w_curr, w_next = w_curr, w_next, w_prev
        c_prev, s_prev, c_curr, s_curr = c_curr, s_curr, c, s

        println("iter = ", i, " res = ", norm(b - A * x), " approx = ", abs(rhs[i + 1]))
    end

    return T
end