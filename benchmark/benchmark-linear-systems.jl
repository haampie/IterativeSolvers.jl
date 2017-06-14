module LinearSystemsBench

using BenchmarkTools
using IterativeSolvers

function indefinite(n)
    # Generate an indefinite "hard" matrix
    srand(1)
    A = 3 * speye(n) + sprand(n, n, 12.0 / n)
    A = (A + A') / 2
    x = ones(n)
    b = A * x

    A, b
end

function posdef(n)
    A, b = indefinite(n)

    # Shift the spectrum a bit to make A positive definite

    return A, b
end

function gmres(; n = 10000, tol = 1e-5, restart::Int = 15, maxiter::Int = 1500)
    A, b = indefinite(n)
    outer = div(maxiter, restart)

    println("Matrix of size ", n, " with ~", nnz(A) / n, " nonzeros per row")
    println("Tolerance = ", tol, "; restart = ", restart, "; max #iterations = ", maxiter)
    
    impr = @benchmark improved_gmres($A, $b, tol = $tol, restart = $restart, outer = $outer)
    old = @benchmark gmres($A, $b, tol = $tol, restart = $restart, maxiter = $maxiter, log = true)

    impr, old
end

function cg(; n = 100000, tol = 1e-10, maxiter::Int = 1500)
    A, b = posdef(n)

    println("Matrix of size ", n, " with ~", nnz(A) / n, " nonzeros per row")
    println("Tolerance = ", tol, "; max #iterations = ", maxiter)
    
    x1, his1 = IterativeSolvers.improved_cg(A, b, maxiter = 3, tol = tol, log = false)
    x2, his2 = IterativeSolvers.cg(A, b, maxiter = 3, tol = tol, log = false)

    @show norm(x1 - x2)

    impr = @benchmark IterativeSolvers.improved_cg($A, $b, maxiter = $maxiter, tol = $tol, log = false)
    old = @benchmark IterativeSolvers.cg($A, $b, maxiter = $maxiter, tol = $tol, log = false)

    impr, old
end

function gmres_allocations(; n = 10000, tol = 1e-5, restart::Int = 15, maxiter::Int = 1500)
    A, b = indefinite(n)
    outer = div(maxiter, restart)

    # Dry run
    impr = gmres(A, b, tol = tol, restart = restart, maxiter = maxiter, log = true)
    Profile.clear_malloc_data()

    # Run the test
    impr = gmres(A, b, tol = tol, restart = restart, maxiter = maxiter, log = true)
end

function cg_allocations(; n = 10000, tol = 1e-10, maxiter::Int = 1500)
    A, b = posdef(n)
    
    # Dry run
    impr = improved_cg(A, b, maxiter = maxiter, tol = tol, log = true)
    Profile.clear_malloc_data()

    # Run the test
    impr = improved_cg(A, b, maxiter = maxiter, tol = tol, log = true)
end
end
