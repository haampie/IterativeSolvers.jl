export minres_iterable, minres

import Base.LinAlg.BLAS.axpy!
import Base: start, next, done

function rotation(a, b)
    r = √(a^2 + b^2)
    a / r, b / r, r
end

function minres2(A, b; maxiter = 10)
    n = length(b)

    x = zeros(b)
    v_prev = similar(b)
    v_curr = copy(b)
    v_next = similar(b)
    w_prev = similar(b)
    w_curr = similar(b)
    w_next = similar(b)

    previous_norm = 0.0

    # Last column of the R matrix: QR = H where H is 
    # the tridiagonal Lanczos matrix.
    R = zeros(4)

    # Last two entries of the right-hand side
    rhs = [norm(v_curr); 0.0]

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

        if i > 1
            axpy!(-prev_norm, v_prev, v_next)
        end
        
        # Orthogonalize w.r.t. v_curr
        R[3] = dot(v_curr, v_next)
        axpy!(-R[3], v_curr, v_next)

        # Normalize
        R[4] = prev_norm = norm(v_next)

        v_next /= prev_norm

        # Previous 2 rotations
        if i > 2
            # Use the fact that T[iter - 2, iter] = 0
            R[1] = s_prev * R[2]
            R[2] = c_prev * R[2]
        end

        if i > 1
            tmp = -s_curr * R[2] + c_curr * R[3]
            R[2] = c_curr * R[2] + s_curr * R[3]
            R[3] = tmp
        end

        # New rotation
        c, s, R[3] = rotation(R[3], R[4])

        # Applied to the rhs:
        rhs[2] = -s * rhs[1]
        rhs[1] = c * rhs[1]

        # Update W = V * inv(R)
        copy!(w_next, v_curr)
        
        # Could be BLAS-2 rather than BLAS-1
        if i > 1
            axpy!(-R[2], w_curr, w_next)
        end

        if i > 2
            axpy!(-R[1], w_prev, w_next)
        end
        
        scale!(w_next, inv(R[3]))

        # Update solution x
        axpy!(rhs[1], w_next, x)

        # Move on: next -> curr, curr -> prev
        v_prev, v_curr, v_next = v_curr, v_next, v_prev
        w_prev, w_curr, w_next = w_curr, w_next, w_prev
        c_prev, s_prev, c_curr, s_curr = c_curr, s_curr, c, s

        R[2] = prev_norm
        rhs[1] = rhs[2]
    end
end