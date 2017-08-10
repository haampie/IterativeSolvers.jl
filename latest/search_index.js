var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#IterativeSolvers.jl-1",
    "page": "Home",
    "title": "IterativeSolvers.jl",
    "category": "section",
    "text": "IterativeSolvers.jl is a Julia package that provides iterative algorithms for solving linear systems, eigenproblems, and singular value problems. The purpose of this package is to provide efficient Julia implementations for iterative methods. The package aims to accept a wide variety of input types and that's why most arguments don't specify a specific type.For bug reports, feature requests and questions please submit an issue. If you're interested in contributing, please see the Contributing guide.For more information on future methods have a look at the package roadmap on deterministic methods, for randomized algorithms check here."
},

{
    "location": "index.html#Linear-Solvers-1",
    "page": "Home",
    "title": "Linear Solvers",
    "category": "section",
    "text": "Stationary methodsJacobi\nGauss-Seidel\nSuccessive over-relaxation (SOR)\nSymmetric successive over-relaxation (SSOR)Non stationary methodsConjugate Gradients (CG)\nMINRES\nBiCGStabl(l)\nIDR(s)\nRestarted GMRES\nChebyshev iteration\nLSMR\nLSQR"
},

{
    "location": "index.html#Eigenproblem-Solvers-1",
    "page": "Home",
    "title": "Eigenproblem Solvers",
    "category": "section",
    "text": "(Inverse) power iteration"
},

{
    "location": "index.html#Singular-Value-Decomposition-1",
    "page": "Home",
    "title": "Singular Value Decomposition",
    "category": "section",
    "text": "Golub-Kahan-Lanczos\nRandomized singular value decomposition"
},

{
    "location": "index.html#Randomized-1",
    "page": "Home",
    "title": "Randomized",
    "category": "section",
    "text": "Condition number estimate\nExtremal eigenvalue estimates\nNorm estimate"
},

{
    "location": "index.html#Documentation-Outline-1",
    "page": "Home",
    "title": "Documentation Outline",
    "category": "section",
    "text": ""
},

{
    "location": "index.html#Manual-1",
    "page": "Home",
    "title": "Manual",
    "category": "section",
    "text": "Pages = [\n    \"user_manual.md\",\n]\nDepth = 2"
},

{
    "location": "index.html#Library-1",
    "page": "Home",
    "title": "Library",
    "category": "section",
    "text": "Pages = [\"library/public.md\", \"library/internal.md\"]\nDepth = 2"
},

{
    "location": "index.html#main-index-1",
    "page": "Home",
    "title": "Index",
    "category": "section",
    "text": ""
},

{
    "location": "index.html#Functions-1",
    "page": "Home",
    "title": "Functions",
    "category": "section",
    "text": "Pages = [\"library/public.md\", \"library/internals.md\"]\nOrder = [:function]"
},

{
    "location": "index.html#Types-1",
    "page": "Home",
    "title": "Types",
    "category": "section",
    "text": "Pages = [\"library/public.md\", \"library/internals.md\"]\nOrder = [:type]"
},

{
    "location": "user_manual.html#",
    "page": "Manual",
    "title": "Manual",
    "category": "page",
    "text": ""
},

{
    "location": "user_manual.html#Manual-1",
    "page": "Manual",
    "title": "Manual",
    "category": "section",
    "text": ""
},

{
    "location": "user_manual.html#Installation-1",
    "page": "Manual",
    "title": "Installation",
    "category": "section",
    "text": "The package can be installed with a simple instruction.julia> Pkg.add(\"IterativeSolvers\")After installing the package, if you wish to use the latest features of the package you must switch to the master branch with Pkg.clone.julia> Pkg.checkout(\"IterativeSolvers\")"
},

{
    "location": "user_manual.html#Interface-1",
    "page": "Manual",
    "title": "Interface",
    "category": "section",
    "text": "All linear-algebraic routines will take as input a linear operator A that maps vectors to vectors. Typically A is a Matrix or a SparseMatrixCSC, but since A is not explicitly typed, any linear operator that supports matrix operations can be used as well. This makes it possible to apply solvers matrix-free. In  IterativeSolvers.jl we strongly recommend LinearMaps.jl  for non-matrix types of A.For matrix-free types of A the following interface is expected to be defined:A*v computes the matrix-vector product on a v::AbstractVector;\nA_mul_B!(y, A, v) computes the matrix-vector product on a v::AbstractVector in-place;\neltype(A) returns the element type implicit in the equivalent matrix representation of A;\nsize(A, d) returns the nominal dimensions along the dth axis in the equivalent matrix representation of A."
},

{
    "location": "user_manual.html#Solvers-1",
    "page": "Manual",
    "title": "Solvers",
    "category": "section",
    "text": "All linear solvers have a common function declaration (with a few exceptions).solver(A, b::Vector; kwargs...)\nsolver!(x, A, b::Vector; kwargs...)In the case of eigenproblems or singular value decompositions:eigsolver(A; kwargs...)\neigsolver!(x, A; kwargs...)A is a linear operator as described above.b is the vector to be solved.x is a vector for the initial guess. In the case of a mutating call this parameter will be overwritten.Output will be the solution to the system."
},

{
    "location": "user_manual.html#Additional-arguments-1",
    "page": "Manual",
    "title": "Additional arguments",
    "category": "section",
    "text": "Keyword names will vary depending on the method, however some of them will always have the same spelling:tol: (relative) stopping tolerance of the method;\nverbose: print information during the iterations;\nmaxiter: maximum number of allowed iterations;\nPl and Pr: left and right preconditioner. See Preconditioning;\nlog::Bool = false: output an extra element of type ConvergenceHistory containing the convergence history."
},

{
    "location": "user_manual.html#log-keyword-1",
    "page": "Manual",
    "title": "log keyword",
    "category": "section",
    "text": "Most solvers contain the log keyword. This is to be used when obtaining more information is required, to use it place the set log to true.x, ch = cg(Master, rand(10, 10), rand(10) log=true)\nsvd, L, ch = svdl(Master, rand(100, 100), log=true)The function will now return one more parameter of type ConvergenceHistory."
},

{
    "location": "user_manual.html#ConvergenceHistory-1",
    "page": "Manual",
    "title": "ConvergenceHistory",
    "category": "section",
    "text": "A ConvergenceHistory instance stores information of a solver.Number of iterations.ch.itersConvergence status.ch.isconvergedStopping tolerances. (A Symbol key is needed to access)ch[:tol]Maximum number of iterations per restart. (Only on restarted methods)nrests(ch)Number of matrix-vectors and matrix-transposed-vector products.nprods(ch)Data stored on each iteration, accessed information can be either a vector or matrix. This data can be a lot of things, most commonly residual. (A Symbol key is needed to access)ch[:resnorm] #Vector or Matrix\nch[:resnorm, x] #Vector or Matrix element\nch[:resnorm, x, y] #Matrix elementThe available keys of each method is described in the Public Documentation."
},

{
    "location": "user_manual.html#Plotting-1",
    "page": "Manual",
    "title": "Plotting",
    "category": "section",
    "text": "ConvergeHistory provides a recipe to use with the package Plots.jl, this makes it really easy to plot on different plot backends. There are two recipes provided:One for the whole ConvergenceHistory.plot(ch)The other one to plot data binded to a key._, ch = gmres(rand(10,10), rand(10), maxiter = 100, log=true)\nplot(ch, :resnorm, sep = :blue)Plot additional keywordssep::Symbol = :white: color of the line separator in restarted methods."
},

{
    "location": "library/cg.html#",
    "page": "Conjugate Gradients",
    "title": "Conjugate Gradients",
    "category": "page",
    "text": ""
},

{
    "location": "library/cg.html#Conjugate-Gradients-(CG)-1",
    "page": "Conjugate Gradients",
    "title": "Conjugate Gradients (CG)",
    "category": "section",
    "text": "Conjugate Gradients solves Ax = b approximately for x where A is a symmetric, positive-definite linear operator and b the right-hand side vector. The method uses short recurrences and therefore has fixed memory costs and fixed computational costs per iteration."
},

{
    "location": "library/cg.html#IterativeSolvers.cg",
    "page": "Conjugate Gradients",
    "title": "IterativeSolvers.cg",
    "category": "Function",
    "text": "cg(A, b; kwargs...) -> x, [history]\n\nSame as cg!, but allocates a solution vector x initialized with zeros.\n\n\n\n"
},

{
    "location": "library/cg.html#IterativeSolvers.cg!",
    "page": "Conjugate Gradients",
    "title": "IterativeSolvers.cg!",
    "category": "Function",
    "text": "cg!(x, A, b; kwargs...) -> x, [history]\n\nArguments\n\nx: Initial guess, will be updated in-place;\nA: linear operator;\nb: right-hand side.\n\nKeywords\n\ninitially_zero::Bool: If true assumes that iszero(x) so that one  matrix-vector product can be saved when computing the initial  residual vector;\nPl = Identity(): left preconditioner of the method. Should be symmetric,  positive-definite like A.\ntol::Real = sqrt(eps(real(eltype(b)))): tolerance for stopping condition |r_k| / |r_0| ≤ tol;\nmaxiter::Integer = size(A,2): maximum number of iterations;\nverbose::Bool = false: print method information;\nlog::Bool = false: keep track of the residual norm in each iteration;\n\nOutput\n\nif log is false\n\nx: approximated solution.\n\nif log is true\n\nx: approximated solution.\nch: convergence history.\n\nConvergenceHistory keys\n\n:tol => ::Real: stopping tolerance.\n:resnom => ::Vector: residual norm at each iteration.\n\n\n\n"
},

{
    "location": "library/cg.html#Usage-1",
    "page": "Conjugate Gradients",
    "title": "Usage",
    "category": "section",
    "text": "cg\ncg!"
},

{
    "location": "library/cg.html#Implementation-details-1",
    "page": "Conjugate Gradients",
    "title": "Implementation details",
    "category": "section",
    "text": "The current implementation follows a rather standard approach. Note that preconditioned CG (or PCG) is slightly different from ordinary CG, because the former must compute the residual explicitly, while it is available as byproduct in the latter. Our implementation of CG ensures the minimal number of vector operations."
},

{
    "location": "library/chebyshev.html#",
    "page": "Chebyshev iteration",
    "title": "Chebyshev iteration",
    "category": "page",
    "text": ""
},

{
    "location": "library/chebyshev.html#IterativeSolvers.chebyshev",
    "page": "Chebyshev iteration",
    "title": "IterativeSolvers.chebyshev",
    "category": "Function",
    "text": "chebyshev(A, b, λmin::Real, λmax::Real; kwargs...) -> x, [history]\n\nSame as chebyshev!, but allocates a solution vector x initialized with zeros.\n\n\n\n"
},

{
    "location": "library/chebyshev.html#IterativeSolvers.chebyshev!",
    "page": "Chebyshev iteration",
    "title": "IterativeSolvers.chebyshev!",
    "category": "Function",
    "text": "chebyshev!(x, A, b, λmin::Real, λmax::Real; kwargs...) -> x, [history]\n\nSolve Ax = b for symmetric, definite matrices A using Chebyshev iteration.\n\nArguments\n\nx: initial guess, will be updated in-place;\nA: linear operator;\nb: right-hand side;\nλmin::Real: lower bound for the real eigenvalues\nλmax::Real: upper bound for the real eigenvalues\n\nKeywords\n\ninitially_zero::Bool = false: if true assumes that iszero(x) so that one  matrix-vector product can be saved when computing the initial  residual vector;\ntol: tolerance for stopping condition |r_k| / |r_0| ≤ tol.\nmaxiter::Int: maximum number of inner iterations of GMRES;\nPl = Identity(): left preconditioner;\nlog::Bool = false: keep track of the residual norm in each iteration;\nverbose::Bool = false: print convergence information during the iterations.\n\nReturn values\n\nif log is false\n\nx: approximate solution.\n\nif log is true\n\nx: approximate solution;\nhistory: convergence history.\n\n\n\n"
},

{
    "location": "library/chebyshev.html#Chebyshev-1",
    "page": "Chebyshev iteration",
    "title": "Chebyshev",
    "category": "section",
    "text": "Chebyshev iteration solves the problem Ax=b approximately for x where A is a symmetric, definite linear operator and b the right-hand side vector. The methods assumes the interval lambda_min lambda_max containing all eigenvalues of A is known, so that x can be iteratively constructed via a Chebyshev polynomial with zeros in this interval. This polynomial ultimately acts as a filter that removes components in the direction of the eigenvectors from the initial residual.The main advantage with respect to Conjugate Gradients is that BLAS1 operations such as inner products are avoided.chebyshev\nchebyshev!note: BLAS1 operations\nAlthough the method is often used to avoid computation of inner products, the stopping criterion is still based on the residual norm. Hence the current implementation is not free of BLAS1 operations."
},

{
    "location": "library/minres.html#",
    "page": "MINRES",
    "title": "MINRES",
    "category": "page",
    "text": ""
},

{
    "location": "library/minres.html#MINRES-1",
    "page": "MINRES",
    "title": "MINRES",
    "category": "section",
    "text": "MINRES is a short-recurrence version of GMRES for solving Ax = b approximately for x where A is a symmetric, Hermitian, skew-symmetric or skew-Hermitian linear operator and b the right-hand side vector."
},

{
    "location": "library/minres.html#IterativeSolvers.minres",
    "page": "MINRES",
    "title": "IterativeSolvers.minres",
    "category": "Function",
    "text": "minres(A, b; kwargs...) -> x, [history]\n\nSame as minres!, but allocates a solution vector x initialized with zeros.\n\n\n\n"
},

{
    "location": "library/minres.html#IterativeSolvers.minres!",
    "page": "MINRES",
    "title": "IterativeSolvers.minres!",
    "category": "Function",
    "text": "minres!(x, A, b; kwargs...) -> x, [history]\n\nSolve Ax = b for (skew-)Hermitian matrices A using MINRES.\n\nArguments\n\nx: initial guess, will be updated in-place;\nA: linear operator;\nb: right-hand side.\n\nKeywords\n\ninitially_zero::Bool = false: if true assumes that iszero(x) so that one  matrix-vector product can be saved when computing the initial  residual vector;\nskew_hermitian::Bool = false: if true assumes that A is skew-symmetric or skew-Hermitian;\ntol: tolerance for stopping condition |r_k| / |r_0| ≤ tol. Note that the residual is computed only approximately;\nmaxiter::Int: maximum number of inner iterations of GMRES;\nPl: left preconditioner;\nPr: right preconditioner;\nlog::Bool = false: keep track of the residual norm in each iteration;\nverbose::Bool = false: print convergence information during the iterations.\n\nReturn values\n\nif log is false\n\nx: approximate solution.\n\nif log is true\n\nx: approximate solution;\nhistory: convergence history.\n\n\n\n"
},

{
    "location": "library/minres.html#Usage-1",
    "page": "MINRES",
    "title": "Usage",
    "category": "section",
    "text": "minres\nminres!"
},

{
    "location": "library/minres.html#Implementation-details-1",
    "page": "MINRES",
    "title": "Implementation details",
    "category": "section",
    "text": "MINRES exploits the tridiagonal structure of the Hessenberg matrix. Although MINRES is mathematically equivalent to GMRES, it might not be equivalent in finite precision. MINRES updates the solution asx = x_0 + (V R^-1) (Q^*r_0e_1)where V is the orthonormal basis for the Krylov subspace and QR is the QR-decomposition of the Hessenberg matrix. Note that the brackets are placed slightly differently from how GMRES would update the residual.MINRES computes V and W = VR^-1 via a three-term recurrence, using only the last column of R Therefore we pre-allocate only six vectors, save only the last two entries of Q^*r_0e_1 and part of the last column of the Hessenberg matrix.note: Real and complex arithmetic\nIf A is Hermitian, then the Hessenberg matrix will be real. This is exploited in the current implementation.If A is skew-Hermitian, the diagonal of the Hessenberg matrix will be imaginary, and hence we use complex arithmetic in that case."
},

{
    "location": "library/bicgstabl.html#",
    "page": "BiCGStab(l)",
    "title": "BiCGStab(l)",
    "category": "page",
    "text": ""
},

{
    "location": "library/bicgstabl.html#BiCGStab(l)-1",
    "page": "BiCGStab(l)",
    "title": "BiCGStab(l)",
    "category": "section",
    "text": "BiCGStab(l) solves the problem Ax = b approximately for x where A is a general, linear operator and b the right-hand side vector. The methods combines BiCG with l GMRES iterations, resulting in a short-reccurence iteration. As a result the memory is fixed as well as the computational costs per iteration."
},

{
    "location": "library/bicgstabl.html#IterativeSolvers.bicgstabl",
    "page": "BiCGStab(l)",
    "title": "IterativeSolvers.bicgstabl",
    "category": "Function",
    "text": "bicgstabl(A, b, l; kwargs...) -> x, [history]\n\nSame as bicgstabl!, but allocates a solution vector x initialized with zeros.\n\n\n\n"
},

{
    "location": "library/bicgstabl.html#IterativeSolvers.bicgstabl!",
    "page": "BiCGStab(l)",
    "title": "IterativeSolvers.bicgstabl!",
    "category": "Function",
    "text": "bicgstabl!(x, A, b, l; kwargs...) -> x, [history]\n\nArguments\n\nA: linear operator;\nb: right hand side (vector);\nl::Int = 2: Number of GMRES steps.\n\nKeywords\n\nmax_mv_products::Int = min(30, size(A, 1)): maximum number of matrix vector products.\n\nFor BiCGStab(l) this is a less dubious term than \"number of iterations\";\n\nPl = Identity(): left preconditioner of the method;\ntol::Real = sqrt(eps(real(eltype(b)))): tolerance for stopping condition |r_k| / |r_0| ≤ tol.   Note that (1) the true residual norm is never computed during the iterations,   only an approximation; and (2) if a preconditioner is given, the stopping condition is based on the   preconditioned residual.\n\nReturn values\n\nif log is false\n\nx: approximate solution.\n\nif log is true\n\nx: approximate solution;\nhistory: convergence history.\n\n\n\n"
},

{
    "location": "library/bicgstabl.html#Usage-1",
    "page": "BiCGStab(l)",
    "title": "Usage",
    "category": "section",
    "text": "bicgstabl\nbicgstabl!"
},

{
    "location": "library/bicgstabl.html#Implementation-details-1",
    "page": "BiCGStab(l)",
    "title": "Implementation details",
    "category": "section",
    "text": "The method is based on the original article [Sleijpen1993], but does not implement later improvements. The normal equations arising from the GMRES steps are solved without orthogonalization. Hence the method should only be reliable for relatively small values of l.The r and u factors are pre-allocated as matrices of size n times (l + 1), so that BLAS2 methods can be used. Also the random shadow residual is pre-allocated as a vector. Hence the storage costs are approximately 2l + 3 vectors.[Sleijpen1993]: Sleijpen, Gerard LG, and Diederik R. Fokkema. \"BiCGstab(l) for  linear equations involving unsymmetric matrices with complex spectrum.\"  Electronic Transactions on Numerical Analysis 1.11 (1993): 2000."
},

{
    "location": "library/gmres.html#",
    "page": "Restarted GMRES",
    "title": "Restarted GMRES",
    "category": "page",
    "text": ""
},

{
    "location": "library/gmres.html#Restarted-GMRES-1",
    "page": "Restarted GMRES",
    "title": "Restarted GMRES",
    "category": "section",
    "text": "GMRES solves the problem Ax = b approximately for x where A is a general, linear operator and b the right-hand side vector. The method is optimal in the sense that it selects the solution with minimal residual from a Krylov subspace, but the price of optimality is increasing storage and computation effort per iteration. Restarts are necessary to fix these costs."
},

{
    "location": "library/gmres.html#IterativeSolvers.gmres",
    "page": "Restarted GMRES",
    "title": "IterativeSolvers.gmres",
    "category": "Function",
    "text": "gmres(A, b; kwargs...) -> x, [history]\n\nSame as gmres!, but allocates a solution vector x initialized with zeros.\n\n\n\n"
},

{
    "location": "library/gmres.html#IterativeSolvers.gmres!",
    "page": "Restarted GMRES",
    "title": "IterativeSolvers.gmres!",
    "category": "Function",
    "text": "gmres!(x, A, b; kwargs...) -> x, [history]\n\nSolves the problem Ax = b with restarted GMRES.\n\nArguments\n\nx: Initial guess, will be updated in-place;\nA: linear operator;\nb: right-hand side.\n\nKeywords\n\ninitially_zero::Bool: If true assumes that iszero(x) so that one  matrix-vector product can be saved when computing the initial  residual vector;\ntol: relative tolerance;\nrestart::Int: restarts GMRES after specified number of iterations;\nmaxiter::Int: maximum number of inner iterations of GMRES;\nPl: left preconditioner;\nPr: right preconditioner;\nlog::Bool: keep track of the residual norm in each iteration;\nverbose::Bool: print convergence information during the iterations.\n\nReturn values\n\nif log is false\n\nx: approximate solution.\n\nif log is true\n\nx: approximate solution;\nhistory: convergence history.\n\n\n\n"
},

{
    "location": "library/gmres.html#Usage-1",
    "page": "Restarted GMRES",
    "title": "Usage",
    "category": "section",
    "text": "gmres\ngmres!"
},

{
    "location": "library/gmres.html#Implementation-details-1",
    "page": "Restarted GMRES",
    "title": "Implementation details",
    "category": "section",
    "text": "The implementation pre-allocates a matrix V of size n by restart whose columns form an orthonormal basis for the Krylov subspace. This allows BLAS2 operations when updating the solution vector x. The Hessenberg matrix is also pre-allocated.Modified Gram-Schmidt is used to orthogonalize the columns of V.The computation of the residual norm is implemented in a non-standard way, namely keeping track of a vector gamma in the null-space of H_k^*, which is the adjoint of the (k + 1) times k Hessenberg matrix H_k at the kth iteration. Only when x needs to be updated is the Hessenberg matrix mutated with Givens rotations. Advanced users can therefore use the underlying GMRESIterable to access the Hessenberg matrix during the iterations.note: Note\nGMRES is designed for general operators A. Consider MINRES for indefinite, (skew-)symmetric or (skew-)Hermitian operators, and CG for positive-definite, symmetric or Hermitian matrices."
},

{
    "location": "library/stationary.html#",
    "page": "Stationary methods",
    "title": "Stationary methods",
    "category": "page",
    "text": ""
},

{
    "location": "library/stationary.html#Stationary-methods-1",
    "page": "Stationary methods",
    "title": "Stationary methods",
    "category": "section",
    "text": "Stationary methods are typically used as smoothers in multigrid methods, where only very few iterations are applied to get rid of high-frequency components in the error. The implementations of stationary methods have this goal in mind, which means there is no other stopping criterion besides the maximum number of iterations.note: CSC versus CSR\nJulia stores matrices column-major. In order to avoid cache misses, the implementations of our stationary methods traverse the matrices column-major. This deviates from classical textbook implementations. Also the (S)SOR methods cannot be computed fully in-place, but require a temporary vector.When it comes to SparseMatrixCSC, we precompute an integer array of the indices of the diagonal as well to avoid expensive searches in each iteration."
},

{
    "location": "library/stationary.html#IterativeSolvers.jacobi",
    "page": "Stationary methods",
    "title": "IterativeSolvers.jacobi",
    "category": "Function",
    "text": "jacobi(A, b)\n\nSolve A*x=b with the Jacobi method.\n\nch is a ConvergenceHistory object. Otherwise it will only return x. If log is set to true is given, method will output a tuple x, ch.\n\nArguments\n\nA: linear operator.\n\nKeywords\n\ntol::Real = size(A,2)^3*eps(): stopping tolerance.\n\nmaxiter::Integer = size(A,2)^2: maximum number of iterations.\n\nverbose::Bool = false: verbose flag.\n\nlog::Bool = false: output an extra element of type ConvergenceHistory containing extra information of the method execution.\n\nOutput\n\nif log is false\n\nx: approximated solution.\n\nif log is true\n\nx: approximated solution.\n\nch: convergence history.\n\nConvergenceHistory keys\n\n:tol => ::Real: stopping tolerance.\n\n:resnom => ::Vector: residual norm at each iteration.\n\n\n\n"
},

{
    "location": "library/stationary.html#IterativeSolvers.jacobi!",
    "page": "Stationary methods",
    "title": "IterativeSolvers.jacobi!",
    "category": "Function",
    "text": "jacobi!(x, A, b)\n\nOverwrite x.\n\nSolve A*x=b with the Jacobi method.\n\nch is a ConvergenceHistory object. Otherwise it will only return x. If log is set to true is given, method will output a tuple x, ch.\n\nArguments\n\nx: initial guess, overwrite final estimation.\n\nA: linear operator.\n\nKeywords\n\ntol::Real = size(A,2)^3*eps(): stopping tolerance.\n\nmaxiter::Integer = size(A,2)^2: maximum number of iterations.\n\nverbose::Bool = false: verbose flag.\n\nlog::Bool = false: output an extra element of type ConvergenceHistory containing extra information of the method execution.\n\nOutput\n\nif log is false\n\nx: approximated solution.\n\nif log is true\n\nx: approximated solution.\n\nch: convergence history.\n\nConvergenceHistory keys\n\n:tol => ::Real: stopping tolerance.\n\n:resnom => ::Vector: residual norm at each iteration.\n\n\n\n"
},

{
    "location": "library/stationary.html#Jacobi-1",
    "page": "Stationary methods",
    "title": "Jacobi",
    "category": "section",
    "text": "jacobi\njacobi!"
},

{
    "location": "library/stationary.html#IterativeSolvers.gauss_seidel",
    "page": "Stationary methods",
    "title": "IterativeSolvers.gauss_seidel",
    "category": "Function",
    "text": "gauss_seidel(A, b)\n\nSolve A*x=b with the Gauss-Seidel method.\n\nch is a ConvergenceHistory object. Otherwise it will only return x. If log is set to true is given, method will output a tuple x, ch.\n\nArguments\n\nA: linear operator.\n\nKeywords\n\ntol::Real = size(A,2)^3*eps(): stopping tolerance.\n\nmaxiter::Integer = size(A,2)^2: maximum number of iterations.\n\nverbose::Bool = false: verbose flag.\n\nlog::Bool = false: output an extra element of type ConvergenceHistory containing extra information of the method execution.\n\nOutput\n\nif log is false\n\nx: approximated solution.\n\nif log is true\n\nx: approximated solution.\n\nch: convergence history.\n\nConvergenceHistory keys\n\n:tol => ::Real: stopping tolerance.\n\n:resnom => ::Vector: residual norm at each iteration.\n\n\n\n"
},

{
    "location": "library/stationary.html#IterativeSolvers.gauss_seidel!",
    "page": "Stationary methods",
    "title": "IterativeSolvers.gauss_seidel!",
    "category": "Function",
    "text": "gauss_seidel!(x, A, b)\n\nOverwrite x.\n\nSolve A*x=b with the Gauss-Seidel method.\n\nch is a ConvergenceHistory object. Otherwise it will only return x. If log is set to true is given, method will output a tuple x, ch.\n\nArguments\n\nx: initial guess, overwrite final estimation.\n\nA: linear operator.\n\nKeywords\n\ntol::Real = size(A,2)^3*eps(): stopping tolerance.\n\nmaxiter::Integer = size(A,2)^2: maximum number of iterations.\n\nverbose::Bool = false: verbose flag.\n\nlog::Bool = false: output an extra element of type ConvergenceHistory containing extra information of the method execution.\n\nOutput\n\nif log is false\n\nx: approximated solution.\n\nif log is true\n\nx: approximated solution.\n\nch: convergence history.\n\nConvergenceHistory keys\n\n:tol => ::Real: stopping tolerance.\n\n:resnom => ::Vector: residual norm at each iteration.\n\n\n\n"
},

{
    "location": "library/stationary.html#Gauss-Seidel-1",
    "page": "Stationary methods",
    "title": "Gauss-Seidel",
    "category": "section",
    "text": "gauss_seidel\ngauss_seidel!"
},

{
    "location": "library/stationary.html#IterativeSolvers.sor",
    "page": "Stationary methods",
    "title": "IterativeSolvers.sor",
    "category": "Function",
    "text": "sor(A, b, ω)\n\nSolve A*x=b with the successive overrelaxation method.\n\nch is a ConvergenceHistory object. Otherwise it will only return x. If log is set to true is given, method will output a tuple x, ch.\n\nArguments\n\nshift::Number=0: shift to be applied to matrix A.\n\nA: linear operator.\n\nKeywords\n\ntol::Real = size(A,2)^3*eps(): stopping tolerance.\n\nmaxiter::Integer = size(A,2)^2: maximum number of iterations.\n\nverbose::Bool = false: verbose flag.\n\nlog::Bool = false: output an extra element of type ConvergenceHistory containing extra information of the method execution.\n\nOutput\n\nif log is false\n\nx: approximated solution.\n\nif log is true\n\nx: approximated solution.\n\nch: convergence history.\n\nConvergenceHistory keys\n\n:tol => ::Real: stopping tolerance.\n\n:resnom => ::Vector: residual norm at each iteration.\n\n\n\n"
},

{
    "location": "library/stationary.html#IterativeSolvers.sor!",
    "page": "Stationary methods",
    "title": "IterativeSolvers.sor!",
    "category": "Function",
    "text": "sor!(x, A, b, ω)\n\nOverwrite x.\n\nSolve A*x=b with the successive overrelaxation method.\n\nch is a ConvergenceHistory object. Otherwise it will only return x. If log is set to true is given, method will output a tuple x, ch.\n\nArguments\n\nx: initial guess, overwrite final estimation.\n\nshift::Number=0: shift to be applied to matrix A.\n\nA: linear operator.\n\nKeywords\n\ntol::Real = size(A,2)^3*eps(): stopping tolerance.\n\nmaxiter::Integer = size(A,2)^2: maximum number of iterations.\n\nverbose::Bool = false: verbose flag.\n\nlog::Bool = false: output an extra element of type ConvergenceHistory containing extra information of the method execution.\n\nOutput\n\nif log is false\n\nx: approximated solution.\n\nif log is true\n\nx: approximated solution.\n\nch: convergence history.\n\nConvergenceHistory keys\n\n:tol => ::Real: stopping tolerance.\n\n:resnom => ::Vector: residual norm at each iteration.\n\n\n\n"
},

{
    "location": "library/stationary.html#Successive-over-relaxation-(SOR)-1",
    "page": "Stationary methods",
    "title": "Successive over-relaxation (SOR)",
    "category": "section",
    "text": "sor\nsor!"
},

{
    "location": "library/stationary.html#IterativeSolvers.ssor",
    "page": "Stationary methods",
    "title": "IterativeSolvers.ssor",
    "category": "Function",
    "text": "ssor(A, b, ω)\n\nSolve A*x=b with the symmetric successive overrelaxation method.\n\nch is a ConvergenceHistory object. Otherwise it will only return x. If log is set to true is given, method will output a tuple x, ch.\n\nArguments\n\nshift::Number=0: shift to be applied to matrix A.\n\nA: linear operator.\n\nKeywords\n\ntol::Real = size(A,2)^3*eps(): stopping tolerance.\n\nmaxiter::Integer = size(A,2)^2: maximum number of iterations.\n\nverbose::Bool = false: verbose flag.\n\nlog::Bool = false: output an extra element of type ConvergenceHistory containing extra information of the method execution.\n\nOutput\n\nif log is false\n\nx: approximated solution.\n\nif log is true\n\nx: approximated solution.\n\nch: convergence history.\n\nConvergenceHistory keys\n\n:tol => ::Real: stopping tolerance.\n\n:resnom => ::Vector: residual norm at each iteration.\n\n\n\n"
},

{
    "location": "library/stationary.html#IterativeSolvers.ssor!",
    "page": "Stationary methods",
    "title": "IterativeSolvers.ssor!",
    "category": "Function",
    "text": "ssor!(x, A, b, ω)\n\nOverwrite x.\n\nSolve A*x=b with the symmetric successive overrelaxation method.\n\nch is a ConvergenceHistory object. Otherwise it will only return x. If log is set to true is given, method will output a tuple x, ch.\n\nArguments\n\nx: initial guess, overwrite final estimation.\n\nshift::Number=0: shift to be applied to matrix A.\n\nA: linear operator.\n\nKeywords\n\ntol::Real = size(A,2)^3*eps(): stopping tolerance.\n\nmaxiter::Integer = size(A,2)^2: maximum number of iterations.\n\nverbose::Bool = false: verbose flag.\n\nlog::Bool = false: output an extra element of type ConvergenceHistory containing extra information of the method execution.\n\nOutput\n\nif log is false\n\nx: approximated solution.\n\nif log is true\n\nx: approximated solution.\n\nch: convergence history.\n\nConvergenceHistory keys\n\n:tol => ::Real: stopping tolerance.\n\n:resnom => ::Vector: residual norm at each iteration.\n\n\n\n"
},

{
    "location": "library/stationary.html#Symmetric-successive-over-relaxation-(SSOR)-1",
    "page": "Stationary methods",
    "title": "Symmetric successive over-relaxation (SSOR)",
    "category": "section",
    "text": "ssor\nssor!"
},

{
    "location": "library/power_method.html#",
    "page": "Power method",
    "title": "Power method",
    "category": "page",
    "text": ""
},

{
    "location": "library/power_method.html#IterativeSolvers.powm",
    "page": "Power method",
    "title": "IterativeSolvers.powm",
    "category": "Function",
    "text": "powm(B; kwargs...) -> λ, x, [history]\n\nSee powm!. Calls powm!(B, x0; kwargs...) with  x0 initialized as a random, complex unit vector.\n\n\n\n"
},

{
    "location": "library/power_method.html#IterativeSolvers.powm!",
    "page": "Power method",
    "title": "IterativeSolvers.powm!",
    "category": "Function",
    "text": "powm!(B, x; shift = zero(eltype(B)), inverse::Bool = false, kwargs...) -> λ, x, [history]\n\nBy default finds the approximate eigenpair (λ, x) of B where |λ| is largest.\n\nArguments\n\nB: linear map, see the note below.\nx: normalized initial guess. Don't forget to use complex arithmetic when necessary.\n\nKeywords\n\ntol::Real = eps(real(eltype(B))) * size(B, 2) ^ 3: stopping tolerance for the residual norm;\nmaxiter::Integer = size(B,2): maximum number of iterations;\nlog::Bool: keep track of the residual norm in each iteration;\nverbose::Bool: print convergence information during the iterations.\n\nnote: Shift-and-invert\nWhen applying shift-and-invert to Ax = x with invert = true and shift = ..., note  that the role of B * b becomes computing inv(A - shift I) * b. So rather than  passing the linear map A itself, pass a linear map B that has the action of  shift-and-invert. The eigenvalue is transformed back to an eigenvalue of the actual  matrix A.\n\nReturn values\n\nif log is false\n\nλ::Number approximate eigenvalue computed as the Rayleigh quotient;\nx::Vector approximate eigenvector.\n\nif log is true\n\nλ::Number: approximate eigenvalue computed as the Rayleigh quotient;\nx::Vector: approximate eigenvector;\nhistory: convergence history.\n\nConvergenceHistory keys\n\n:tol => ::Real: stopping tolerance;\n:resnom => ::Vector: residual norm at each iteration.\n\nExamples\n\nusing LinearMaps\nσ = 1.0 + 1.3im\nA = rand(Complex128, 50, 50)\nF = lufact(A - σ * I)\nFmap = LinearMap{Complex128}((y, x) -> A_ldiv_B!(y, F, x), 50, ismutating = true)\nλ, x = powm(Fmap, inverse = true, shift = σ, tol = 1e-4, maxiter = 200)\n\n\n\n"
},

{
    "location": "library/power_method.html#IterativeSolvers.invpowm",
    "page": "Power method",
    "title": "IterativeSolvers.invpowm",
    "category": "Function",
    "text": "invpowm(B; shift = σ, kwargs...) -> λ, x, [history]\n\nFind the approximate eigenpair (λ, x) of A near shift, where B is a linear map that has the effect B * v = inv(A - σI) * v.\n\nThe method calls powm!(B, x0; inverse = true, shift = σ) with x0 a random, complex unit vector. See powm!\n\nExamples\n\nusing LinearMaps\nσ = 1.0 + 1.3im\nA = rand(Complex128, 50, 50)\nF = lufact(A - σ * I)\nFmap = LinearMap{Complex128}((y, x) -> A_ldiv_B!(y, F, x), 50, ismutating = true)\nλ, x = invpowm(Fmap, shift = σ, tol = 1e-4, maxiter = 200)\n\n\n\n"
},

{
    "location": "library/power_method.html#IterativeSolvers.invpowm!",
    "page": "Power method",
    "title": "IterativeSolvers.invpowm!",
    "category": "Function",
    "text": "invpowm!(B, x0; shift = σ, kwargs...) -> λ, x, [history]\n\nFind the approximate eigenpair (λ, x) of A near shift, where B is a linear map that has the effect B * v = inv(A - σI) * v.\n\nThe method calls powm!(B, x0; inverse = true, shift = σ). See powm!.\n\n\n\n"
},

{
    "location": "library/power_method.html#(Inverse)-power-method-1",
    "page": "Power method",
    "title": "(Inverse) power method",
    "category": "section",
    "text": "Solves the eigenproblem Ax = x approximately where A is a general linear map. By default converges towards the dominant eigenpair ( x) such that  is largest. Shift-and-invert can be applied to target a specific eigenvalue near shift in the complex plane.powm\npowm!\ninvpowm\ninvpowm!"
},

{
    "location": "library/power_method.html#Implementation-details-1",
    "page": "Power method",
    "title": "Implementation details",
    "category": "section",
    "text": "Storage requirements are 3 vectors: the approximate eigenvector x, the residual vector r and a temporary. The residual norm lags behind one iteration, as it is computed when Ax is performed. Therefore the final resdiual norm is even smaller."
},

{
    "location": "preconditioning.html#",
    "page": "Preconditioning",
    "title": "Preconditioning",
    "category": "page",
    "text": ""
},

{
    "location": "preconditioning.html#Preconditioning-1",
    "page": "Preconditioning",
    "title": "Preconditioning",
    "category": "section",
    "text": "Many iterative solvers have the option to provide left and right preconditioners (Pl and Pr resp.) in order to speed up convergence or prevent stagnation. They transform a problem Ax = b into a better conditioned system (P_l^-1AP_r^-1)y = P_l^-1b, where x = P_r^-1y.These preconditioners should support the operations A_ldiv_B!(y, P, x) computes P \\ x in-place of y;\nA_ldiv_B!(P, x) computes P \\ x in-place of x;\nand P \\ x.If no preconditioners are passed to the solver, the method will default to Pl = Pr = IterativeSolvers.Identity().IterativeSolvers.jl itself does not provide any other preconditioners besides Identity(), but recommends the following external packages:ILU.jl for incomplete LU decompositions (using drop tolerance);\nIncompleteSelectedInversion.jl for incomplete LDLt decompositions."
},

{
    "location": "library/public.html#",
    "page": "Public",
    "title": "Public",
    "category": "page",
    "text": ""
},

{
    "location": "library/public.html#Iterative-methods-1",
    "page": "Public",
    "title": "Iterative methods",
    "category": "section",
    "text": "Documentation for IterativeSolvers.jl's public interface.Pages = [\"public.md\"]\nDepth = 4"
},

{
    "location": "library/public.html#Index-1",
    "page": "Public",
    "title": "Index",
    "category": "section",
    "text": "Pages = [\"public.md\"]"
},

{
    "location": "library/public.html#Linear-Solvers-1",
    "page": "Public",
    "title": "Linear Solvers",
    "category": "section",
    "text": ""
},

{
    "location": "library/public.html#IterativeSolvers.idrs",
    "page": "Public",
    "title": "IterativeSolvers.idrs",
    "category": "Function",
    "text": "idrs(A, b)\n\nSolve A*x=b using the induced dimension reduction method.\n\nIf log is set to true is given, method will output a tuple x, ch. Where ch is a ConvergenceHistory object. Otherwise it will only return x.\n\nArguments\n\nA: linear operator.\n\nb: right hand side.\n\nKeywords\n\nPl = 1: left preconditioner of the method.\n\nPr = 1: right preconditioner of the method.\n\ntol::Real = sqrt(eps()): stopping tolerance.\n\nrestart::Integer = min(20,length(b)): maximum number of iterations per restart.\n\nmaxiter::Integer = min(20,length(b)): maximum number of iterations.\n\nverbose::Bool = false: print method information.\n\nlog::Bool = false: output an extra element of type ConvergenceHistory containing extra information of the method execution.\n\nOutput\n\nif log is false\n\nx: approximated solution.\n\nif log is true\n\nx: approximated solution.\n\nch: convergence history.\n\nConvergenceHistory keys\n\n:tol => ::Real: stopping tolerance.\n\n:resnom => ::Vector: residual norm at each iteration.\n\n\n\n"
},

{
    "location": "library/public.html#IterativeSolvers.idrs!",
    "page": "Public",
    "title": "IterativeSolvers.idrs!",
    "category": "Function",
    "text": "idrs!(x, A, b)\n\nOverwrite x.\n\nSolve A*x=b using the induced dimension reduction method.\n\nIf log is set to true is given, method will output a tuple x, ch. Where ch is a ConvergenceHistory object. Otherwise it will only return x.\n\nArguments\n\nx: initial guess, overwrite final estimation.\n\nA: linear operator.\n\nb: right hand side.\n\nKeywords\n\nPl = 1: left preconditioner of the method.\n\nPr = 1: right preconditioner of the method.\n\ntol::Real = sqrt(eps()): stopping tolerance.\n\nrestart::Integer = min(20,length(b)): maximum number of iterations per restart.\n\nmaxiter::Integer = min(20,length(b)): maximum number of iterations.\n\nverbose::Bool = false: print method information.\n\nlog::Bool = false: output an extra element of type ConvergenceHistory containing extra information of the method execution.\n\nOutput\n\nif log is false\n\nx: approximated solution.\n\nif log is true\n\nx: approximated solution.\n\nch: convergence history.\n\nConvergenceHistory keys\n\n:tol => ::Real: stopping tolerance.\n\n:resnom => ::Vector: residual norm at each iteration.\n\n\n\n"
},

{
    "location": "library/public.html#IDR(s)-1",
    "page": "Public",
    "title": "IDR(s)",
    "category": "section",
    "text": "idrs\nidrs!References[1] IDR(s): a family of simple and fast algorithms for solving large\n    nonsymmetric linear systems. P. Sonneveld and M. B. van Gijzen\n    SIAM J. Sci. Comput. Vol. 31, No. 2, pp. 1035--1062, 2008\n[2] Algorithm 913: An Elegant IDR(s) Variant that Efficiently Exploits\n    Bi-orthogonality Properties. M. B. van Gijzen and P. Sonneveld\n    ACM Trans. Math. Software,, Vol. 38, No. 1, pp. 5:1-5:19, 2011\n[3] This file is a translation of the following MATLAB implementation:\n    http://ta.twi.tudelft.nl/nw/users/gijzen/idrs.m\n[4] IDR(s)' webpage http://ta.twi.tudelft.nl/nw/users/gijzen/IDR.html"
},

{
    "location": "library/public.html#IterativeSolvers.lsmr",
    "page": "Public",
    "title": "IterativeSolvers.lsmr",
    "category": "Function",
    "text": "lsmr(A, b)\n\nMinimize ||Ax-b||^2 + λ^2 ||x||^2 for A*x=b.\n\nThe method is based on the Golub-Kahan bidiagonalization process. It is algebraically equivalent to applying MINRES to the normal equation (ATA+λ2I)x=ATb, but has better numerical properties, especially if A is ill-conditioned.\n\nIf log is set to true is given, method will output a tuple x, ch. Where ch is a ConvergenceHistory object. Otherwise it will only return x.\n\nArguments\n\nA: linear operator.\n\nb: right hand side.\n\nKeywords\n\nλ::Number = 0: lambda.\n\natol::Number = 1e-6, btol::Number = 1e-6: stopping tolerances. If both are 1.0e-9 (say), the final residual norm should be accurate to about 9 digits. (The final x will usually have fewer correct digits, depending on cond(A) and the size of damp).\n\nconlim::Number = 1e8: stopping tolerance.  lsmr terminates if an estimate of cond(A) exceeds conlim.  For compatible systems Ax = b, conlim could be as large as 1.0e+12 (say).  For least-squares problems, conlim should be less than 1.0e+8. Maximum precision can be obtained by setting atol = btol = conlim = zero, but the number of iterations may then be excessive.\n\nmaxiter::Integer = min(20,length(b)): maximum number of iterations.\n\nverbose::Bool = false: print method information.\n\nlog::Bool = false: output an extra element of type ConvergenceHistory containing extra information of the method execution.\n\nOutput\n\nif log is false\n\nx: approximated solution.\n\nif log is true\n\nx: approximated solution.\n\nch: convergence history.\n\nConvergenceHistory keys\n\n:atol => ::Real: atol stopping tolerance.\n\n:btol => ::Real: btol stopping tolerance.\n\n:ctol => ::Real: ctol stopping tolerance.\n\n:anorm => ::Real: anorm.\n\n:rnorm => ::Real: rnorm.\n\n:cnorm => ::Real: cnorm.\n\n:resnom => ::Vector: residual norm at each iteration.\n\n\n\n"
},

{
    "location": "library/public.html#IterativeSolvers.lsmr!",
    "page": "Public",
    "title": "IterativeSolvers.lsmr!",
    "category": "Function",
    "text": "lsmr!(x, A, b)\n\nOverwrite x.\n\nMinimize ||Ax-b||^2 + λ^2 ||x||^2 for A*x=b.\n\nThe method is based on the Golub-Kahan bidiagonalization process. It is algebraically equivalent to applying MINRES to the normal equation (ATA+λ2I)x=ATb, but has better numerical properties, especially if A is ill-conditioned.\n\nIf log is set to true is given, method will output a tuple x, ch. Where ch is a ConvergenceHistory object. Otherwise it will only return x.\n\nArguments\n\nx: initial guess, overwrite final estimation.\n\nA: linear operator.\n\nb: right hand side.\n\nKeywords\n\nλ::Number = 0: lambda.\n\natol::Number = 1e-6, btol::Number = 1e-6: stopping tolerances. If both are 1.0e-9 (say), the final residual norm should be accurate to about 9 digits. (The final x will usually have fewer correct digits, depending on cond(A) and the size of damp).\n\nconlim::Number = 1e8: stopping tolerance.  lsmr terminates if an estimate of cond(A) exceeds conlim.  For compatible systems Ax = b, conlim could be as large as 1.0e+12 (say).  For least-squares problems, conlim should be less than 1.0e+8. Maximum precision can be obtained by setting atol = btol = conlim = zero, but the number of iterations may then be excessive.\n\nmaxiter::Integer = min(20,length(b)): maximum number of iterations.\n\nverbose::Bool = false: print method information.\n\nlog::Bool = false: output an extra element of type ConvergenceHistory containing extra information of the method execution.\n\nOutput\n\nif log is false\n\nx: approximated solution.\n\nif log is true\n\nx: approximated solution.\n\nch: convergence history.\n\nConvergenceHistory keys\n\n:atol => ::Real: atol stopping tolerance.\n\n:btol => ::Real: btol stopping tolerance.\n\n:ctol => ::Real: ctol stopping tolerance.\n\n:anorm => ::Real: anorm.\n\n:rnorm => ::Real: rnorm.\n\n:cnorm => ::Real: cnorm.\n\n:resnom => ::Vector: residual norm at each iteration.\n\n\n\n"
},

{
    "location": "library/public.html#LSMR-1",
    "page": "Public",
    "title": "LSMR",
    "category": "section",
    "text": "lsmr\nlsmr!ReferencesAdapted from: http://web.stanford.edu/group/SOL/software/lsmr/"
},

{
    "location": "library/public.html#IterativeSolvers.lsqr",
    "page": "Public",
    "title": "IterativeSolvers.lsqr",
    "category": "Function",
    "text": "lsqr(A, b)\n\nLSQR solves Ax = b or min ||b - Ax||^2 if damp = 0, or   min ||(b) - (  A   )x||   otherwise.          ||(0)   (damp*I) ||^2.\n\nThe method is based on the Golub-Kahan bidiagonalization process. It is algebraically equivalent to applying CG to the normal equation (ATA+λ2I)x=ATb, but has better numerical properties, especially if A is ill-conditioned.\n\nIf log is set to true is given, method will output a tuple x, ch. Where ch is a ConvergenceHistory object. Otherwise it will only return x.\n\nArguments\n\nA: linear operator.\n\nb: right hand side.\n\nKeywords\n\ndamp::Number = 0: damping parameter.\n\natol::Number = 1e-6, btol::Number = 1e-6: stopping tolerances. If both are 1.0e-9 (say), the final residual norm should be accurate to about 9 digits. (The final x will usually have fewer correct digits, depending on cond(A) and the size of damp).\n\nconlim::Number = 1e8: stopping tolerance.  lsmr terminates if an estimate of cond(A) exceeds conlim.  For compatible systems Ax = b, conlim could be as large as 1.0e+12 (say).  For least-squares problems, conlim should be less than 1.0e+8. Maximum precision can be obtained by setting atol = btol = conlim = zero, but the number of iterations may then be excessive.\n\nmaxiter::Integer = min(20,length(b)): maximum number of iterations.\n\nverbose::Bool = false: print method information.\n\nlog::Bool = false: output an extra element of type ConvergenceHistory containing extra information of the method execution.\n\nOutput\n\nif log is false\n\nx: approximated solution.\n\nif log is true\n\nx: approximated solution.\n\nch: convergence history.\n\nConvergenceHistory keys\n\n:atol => ::Real: atol stopping tolerance.\n\n:btol => ::Real: btol stopping tolerance.\n\n:ctol => ::Real: ctol stopping tolerance.\n\n:anorm => ::Real: anorm.\n\n:rnorm => ::Real: rnorm.\n\n:cnorm => ::Real: cnorm.\n\n:resnom => ::Vector: residual norm at each iteration.\n\n\n\n"
},

{
    "location": "library/public.html#IterativeSolvers.lsqr!",
    "page": "Public",
    "title": "IterativeSolvers.lsqr!",
    "category": "Function",
    "text": "lsqr!(x, A, b)\n\nOverwrite x.\n\nLSQR solves Ax = b or min ||b - Ax||^2 if damp = 0, or   min ||(b) - (  A   )x||   otherwise.          ||(0)   (damp*I) ||^2.\n\nThe method is based on the Golub-Kahan bidiagonalization process. It is algebraically equivalent to applying CG to the normal equation (ATA+λ2I)x=ATb, but has better numerical properties, especially if A is ill-conditioned.\n\nIf log is set to true is given, method will output a tuple x, ch. Where ch is a ConvergenceHistory object. Otherwise it will only return x.\n\nArguments\n\nx: initial guess, overwrite final estimation.\n\nA: linear operator.\n\nb: right hand side.\n\nKeywords\n\ndamp::Number = 0: damping parameter.\n\natol::Number = 1e-6, btol::Number = 1e-6: stopping tolerances. If both are 1.0e-9 (say), the final residual norm should be accurate to about 9 digits. (The final x will usually have fewer correct digits, depending on cond(A) and the size of damp).\n\nconlim::Number = 1e8: stopping tolerance.  lsmr terminates if an estimate of cond(A) exceeds conlim.  For compatible systems Ax = b, conlim could be as large as 1.0e+12 (say).  For least-squares problems, conlim should be less than 1.0e+8. Maximum precision can be obtained by setting atol = btol = conlim = zero, but the number of iterations may then be excessive.\n\nmaxiter::Integer = min(20,length(b)): maximum number of iterations.\n\nverbose::Bool = false: print method information.\n\nlog::Bool = false: output an extra element of type ConvergenceHistory containing extra information of the method execution.\n\nOutput\n\nif log is false\n\nx: approximated solution.\n\nif log is true\n\nx: approximated solution.\n\nch: convergence history.\n\nConvergenceHistory keys\n\n:atol => ::Real: atol stopping tolerance.\n\n:btol => ::Real: btol stopping tolerance.\n\n:ctol => ::Real: ctol stopping tolerance.\n\n:anorm => ::Real: anorm.\n\n:rnorm => ::Real: rnorm.\n\n:cnorm => ::Real: cnorm.\n\n:resnom => ::Vector: residual norm at each iteration.\n\n\n\n"
},

{
    "location": "library/public.html#LSQR-1",
    "page": "Public",
    "title": "LSQR",
    "category": "section",
    "text": "lsqr\nlsqr!ReferencesAdapted from: http://web.stanford.edu/group/SOL/software/lsqr/\n\n1. C. C. Paige and M. A. Saunders (1982a).\n    LSQR: An algorithm for sparse linear equations and sparse least squares,\n    ACM TOMS 8(1), 43-71.\n\n2. C. C. Paige and M. A. Saunders (1982b).\n    Algorithm 583.  LSQR: Sparse linear equations and least squares problems,\n    ACM TOMS 8(2), 195-209.\n\n3. M. A. Saunders (1995).  Solution of sparse rectangular systems using\n    LSQR and CRAIG, BIT 35, 588-604."
},

{
    "location": "library/public.html#Eigen-Solvers-1",
    "page": "Public",
    "title": "Eigen Solvers",
    "category": "section",
    "text": ""
},

{
    "location": "library/public.html#IterativeSolvers.svdl",
    "page": "Public",
    "title": "IterativeSolvers.svdl",
    "category": "Function",
    "text": "svdl(A)\n\nCompute some singular values (and optionally vectors) using Golub-Kahan-Lanczos bidiagonalization \\cite{Golub1965} with thick restarting \\cite{Wu2000}.\n\nIf log is set to true is given, method will output a tuple X, L, ch. Where ch is a ConvergenceHistory object. Otherwise it will only return X, L.\n\nArguments\n\nA : The matrix or matrix-like object whose singular values are desired.\n\nKeywords\n\nnsv::Int = 6: number of singular values requested.\n\nv0 = random unit vector: starting guess vector in the domain of A. The length of q should be the number of columns in A.\n\nk::Int = 2nsv: maximum number of Lanczos vectors to compute before restarting.\n\nj::Int = nsv: number of vectors to keep at the end of the restart. We don't recommend j < nsv.\n\nmaxiter::Int = minimum(size(A)): maximum number of iterations to run.\n\nverbose::Bool = false: print information at each iteration.\n\ntol::Real = √eps(): maximum absolute error in each desired singular value.\n\nreltol::Real=√eps(): maximum error in each desired singular value relative to the estimated norm of the input matrix.\n\nmethod::Symbol=:ritz: restarting algorithm to use. Valid choices are:\n\n:ritz: Thick restart with Ritz values [Wu2000].\n:harmonic: Restart with harmonic Ritz values [Baglama2005].\n\nvecs::Symbol = :none: singular vectors to return.\n\n:both: Both left and right singular vectors are returned.\n:left: Only the left singular vectors are returned.\n:right: Only the right singular vectors are returned.\n:none: No singular vectors are returned.\n\ndolock::Bool=false: If true, locks converged Ritz values, removing them from the Krylov subspace being searched in the next macroiteration.\n\nlog::Bool = false: output an extra element of type ConvergenceHistory containing extra information of the method execution.\n\nOutput\n\nif log is false\n\nΣ: list of the desired singular values if vecs == :none (the default), otherwise returns an SVD object with the desired singular vectors filled in.\n\nL: computed partial factorizations of A.\n\nif log is true\n\nΣ: list of the desired singular values if vecs == :none (the default), otherwise returns an SVD object with the desired singular vectors filled in.\n\nL: computed partial factorizations of A.\n\nch::ConvergenceHistory: convergence history.\n\nConvergenceHistory keys\n\n:betas => betas: The history of the computed betas.\n\n:Bs => Bs: The history of the computed projected matrices.\n\n:ritz => ritzvalhist: Ritz values computed at each iteration.\n\n:conv => convhist: Convergence data.\n\n\n\n"
},

{
    "location": "library/public.html#Golub-Kahan-Lanczos-1",
    "page": "Public",
    "title": "Golub-Kahan-Lanczos",
    "category": "section",
    "text": "svdlImplementation notesThe implementation of thick restarting follows closely that of SLEPc as described in [Hernandez2008]. Thick restarting can be turned off by setting k = maxiter, but most of the time this is not desirable.The singular vectors are computed directly by forming the Ritz vectors from the product of the Lanczos vectors L.P/L.Q and the singular vectors of L.B. Additional accuracy in the singular triples can be obtained using inverse iteration.References@article{Golub1965,\n    author = {Golub, G. and Kahan, W.},\n    doi = {10.1137/0702016},\n    journal = {Journal of the Society for Industrial and Applied Mathematics\n        Series B Numerical Analysis},\n    volume = 2,\n    number = 2,\n    pages = {205--224},\n    title = {Calculating the Singular Values and Pseudo-Inverse of a Matrix},\n    year = 1965\n}\n\n@article{Wu2000,\n    author = {Wu, Kesheng and Simon, Horst},\n    journal = {SIAM Journal on Matrix Analysis and Applications},\n    number = 2,\n    pages = {602--616},\n    title = {Thick-Restart {L}anczos Method for Large Symmetric Eigenvalue Problems},\n    volume = 22,\n    year = 2000\n}\n\n@article{Baglama2005,\n    author = {Baglama, James and Reichel, Lothar},\n    doi = {10.1137/04060593X},\n    journal = {SIAM Journal on Scientific Computing},\n    number = 1,\n    pages = {19--42},\n    title = {Augmented Implicitly Restarted {L}anczos Bidiagonalization Methods},\n    volume = 27,\n    year = 2005\n}\n\n@article{Hernandez2008,\n    author = {Hern\\'{a}ndez, Vicente and Rom\\'{a}n, Jos\\'{e} E and Tom\\'{a}s,\n    Andr\\'{e}s},\n    journal = {Electronic Transactions on Numerical Analysis},\n    pages = {68--85},\n    title = {A Robust and Efficient Parallel {SVD} Solver based on Restarted\n        {L}anczos Bidiagonalization},\n    url = {http://etna.mcs.kent.edu/volumes/2001-2010/vol31/abstract.php?vol=31\\&pages=68-85},\n    volume = 31,\n    year = 2008\n}"
},

{
    "location": "library/public.html#Randomized-1",
    "page": "Public",
    "title": "Randomized",
    "category": "section",
    "text": ""
},

{
    "location": "library/public.html#IterativeSolvers.rcond",
    "page": "Public",
    "title": "IterativeSolvers.rcond",
    "category": "Function",
    "text": "rcond(A, iters=1)\n\nEstimate matrix condition number randomly.\n\nArguments\n\nA: matrix whose condition number to estimate. Must be square and support premultiply (A*⋅) and solve (A\\⋅).\n\niters::Int = 1: number of power iterations to run.\n\nKeywords\n\np::Real = 0.05: probability that estimate fails to hold as an upper bound.\n\nOutput\n\nInterval (x, y) which contains κ(A) with probability 1 - p.\n\nImplementation note\n\n\\cite{Dixon1983} originally describes this as a computation that can be done by computing the necessary number of power iterations given p and the desired accuracy parameter θ=y/x. However, these bounds were only derived under the assumptions of exact arithmetic. Empirically, iters≥4 has been seen to result in incorrect results in that the computed interval does not contain the true condition number. This implemention therefore makes iters an explicitly user-controllable parameter from which to infer the accuracy parameter and hence the interval containing κ(A).\n\nReferences\n\n\\cite[Theorem 2]{Dixon1983}\n\n@article{Dixon1983,\n    author = {Dixon, John D},\n    doi = {10.1137/0720053},\n    journal = {SIAM Journal on Numerical Analysis},\n    number = {4},\n    pages = {812--814},\n    title = {Estimating Extremal Eigenvalues and Condition Numbers of\n	Matrices},\n    volume = {20},\n    year = {1983}\n}\n\n\n\n"
},

{
    "location": "library/public.html#Condition-number-estimate-1",
    "page": "Public",
    "title": "Condition number estimate",
    "category": "section",
    "text": "rcond"
},

{
    "location": "library/public.html#IterativeSolvers.reigmin",
    "page": "Public",
    "title": "IterativeSolvers.reigmin",
    "category": "Function",
    "text": "reigmin(A, iters=1)\n\nEstimate minimal eigenvalue randomly.\n\nArguments\n\nA: Matrix whose maximal eigenvalue to estimate. Must be square and support premultiply (A*⋅).\n\niters::Int=1: Number of power iterations to run. (Recommended: iters ≤ 3)\n\nKeywords\n\np::Real=0.05: Probability that estimate fails to hold as an upper bound.\n\nOutput\n\nInterval (x, y) which contains the maximal eigenvalue of A with probability 1 - p.\n\nReferences\n\n\\cite[Corollary of Theorem 1]{Dixon1983}.\n\n\n\n"
},

{
    "location": "library/public.html#IterativeSolvers.reigmax",
    "page": "Public",
    "title": "IterativeSolvers.reigmax",
    "category": "Function",
    "text": "reigmax(A, iters=1)\n\nEstimate maximal eigenvalue randomly.\n\nArguments\n\nA: Matrix whose maximal eigenvalue to estimate. Must be square and support premultiply (A*⋅).\n\niters::Int=1: Number of power iterations to run. (Recommended: iters ≤ 3)\n\nKeywords\n\np::Real=0.05: Probability that estimate fails to hold as an upper bound.\n\nOutput\n\nInterval (x, y) which contains the maximal eigenvalue of A with probability 1 - p.\n\nReferences\n\n\\cite[Corollary of Theorem 1]{Dixon1983}.\n\n\n\n"
},

{
    "location": "library/public.html#Extremal-eigenvalue-estimates-1",
    "page": "Public",
    "title": "Extremal eigenvalue estimates",
    "category": "section",
    "text": "reigmin\nreigmax"
},

{
    "location": "library/public.html#IterativeSolvers.rnorm",
    "page": "Public",
    "title": "IterativeSolvers.rnorm",
    "category": "Function",
    "text": "rnorm(A, mvps)\n\nCompute a probabilistic upper bound on the norm of a matrix A. ‖A‖ ≤ α √(2/π) maxᵢ ‖Aωᵢ‖ with probability p=α^(-mvps).\n\nArguments\n\nA: matrix whose norm to estimate.\n\nmvps::Int: number of matrix-vector products to compute.\n\nKeywords\n\np::Real=0.05: probability of upper bound failing.\n\nOutput\n\nEstimate of ‖A‖.\n\nSee also\n\nsee rnorms for a different estimator that uses premultiplying by both A and A'.\n\nReferences\n\n\\cite[Lemma 4.1]{Halko2011}\n\n\n\n"
},

{
    "location": "library/public.html#IterativeSolvers.rnorms",
    "page": "Public",
    "title": "IterativeSolvers.rnorms",
    "category": "Function",
    "text": "rnorms(A, iters=1)\n\nEstimate matrix norm randomly using A'A.\n\nCompute a probabilistic upper bound on the norm of a matrix A.\n\nρ = √(‖(A'A)ʲω‖/‖(A'A)ʲ⁻¹ω‖)\n\nwhich is an estimate of the spectral norm of A produced by iters steps of the power method starting with normalized ω, is a lower bound on the true norm by a factor\n\nρ ≤ α ‖A‖\n\nwith probability greater than 1 - p, where p = 4\\sqrt(n/(iters-1)) α^(-2iters).\n\nArguments\n\nA: matrix whose norm to estimate.\n\niters::Int = 1: mumber of power iterations to perform.\n\nKeywords\n\np::Real = 0.05: probability of upper bound failing.\n\nAt = A': Transpose of A.\n\nOutput\n\nEstimate of ‖A‖.\n\nSee also\n\nsee rnorm for a different estimator that does not require premultiplying by A'\n\nReferences\n\nAppendix of \\cite{Liberty2007}.\n\n@article{Liberty2007,\n    authors = {Edo Liberty and Franco Woolfe and Per-Gunnar Martinsson\n    and Vladimir Rokhlin and Mark Tygert},\n    title = {Randomized algorithms for the low-rank approximation of matrices},\n    journal = {Proceedings of the National Academy of Sciences},\n    volume = {104},\n    issue = {51},\n    year = {2007},\n    pages = {20167--20172},\n    doi  = {10.1073/pnas.0709640104}\n}\n\n\n\n"
},

{
    "location": "library/public.html#Norm-estimate-1",
    "page": "Public",
    "title": "Norm estimate",
    "category": "section",
    "text": "rnorm\nrnorms"
},

{
    "location": "library/public.html#IterativeSolvers.reig",
    "page": "Public",
    "title": "IterativeSolvers.reig",
    "category": "Function",
    "text": "reig(A, l)\n\nCompute the spectral (Eigen) decomposition of A using a randomized algorithm.\n\nArguments\n\nA: input matrix.\n\nl::Int: number of eigenpairs to find.\n\nOutput\n\n::Base.LinAlg.Eigen: eigen decomposition.\n\nImplementation note\n\nThis is a wrapper around eigfact_onepass() which uses the randomized samples found using srft(l).\n\nReferences\n\n@article{Halko2011,\n    author = {Halko, N and Martinsson, P G and Tropp, J A},\n    doi = {10.1137/090771806},\n    journal = {SIAM Review},\n    month = jan,\n    number = {2},\n    pages = {217--288},\n    title = {Finding Structure with Randomness: Probabilistic Algorithms for Constructing Approximate Matrix Decompositions},\n    volume = {53},\n    year = {2011}\n}\n\n\n\n"
},

{
    "location": "library/public.html#IterativeSolvers.rsvdfact",
    "page": "Public",
    "title": "IterativeSolvers.rsvdfact",
    "category": "Function",
    "text": "rsvdfact(A, n, p=0)\n\nCompute partial singular value decomposition of A using a randomized algorithm.\n\nArguments\n\nA: input matrix.\n\nn::Int: number of singular value/vector pairs to find.\n\np::Int=0: number of extra vectors to include in computation.\n\nOutput\n\n::SVD: singular value decomposition.\n\nWarning\n\nThis variant of the randomized singular value decomposition is the most commonly found implementation but is not recommended for accurate computations, as it often has trouble finding the n largest singular pairs, but rather finds n large singular pairs which may not necessarily be the largest.\n\nImplementation note\n\nThis function calls rrange, which uses naive randomized rangefinding to compute a basis for a subspace of dimension n (Algorithm 4.1 of \\cite{Halko2011}), followed by svdfact_restricted(), which computes the exact SVD factorization on the restriction of A to this randomly selected subspace (Algorithm 5.1 of \\cite{Halko2011}).\n\nAlternatively, you can mix and match your own randomized algorithm using any of the randomized range finding algorithms to find a suitable subspace and feeding the result to one of the routines that computes the SVD restricted to that subspace.\n\nReferences\n\n@article{Halko2011,\n    author = {Halko, N and Martinsson, P G and Tropp, J A},\n    doi = {10.1137/090771806},\n    journal = {SIAM Review},\n    month = jan,\n    number = {2},\n    pages = {217--288},\n    title = {Finding Structure with Randomness: Probabilistic Algorithms for Constructing Approximate Matrix Decompositions},\n    volume = {53},\n    year = {2011}\n}\n\n\n\n"
},

{
    "location": "library/public.html#IterativeSolvers.rsvd_fnkz",
    "page": "Public",
    "title": "IterativeSolvers.rsvd_fnkz",
    "category": "Function",
    "text": "rsvd_fnkz(A, k)\n\nCompute the randomized SVD by iterative refinement from randomly selected columns/rows.\n\nArguments\n\nA: matrix whose SVD is desired.\n\nk::Int: desired rank of approximation (k ≤ min(m, n)).\n\nKeywords\n\nl::Int = k: number of columns/rows to sample at each iteration (1 ≤ l ≤ k).\n\nN::Int = minimum(size(A)): maximum number of iterations.\n\nϵ::Real = prod(size(A))*eps(): relative threshold for convergence, as measured by growth of the spectral norm.\n\nmethod::Symbol = :eig: problem to solve.\n\n:eig: eigenproblem.\n:svd: singular problem.\n\nverbose::Bool = false: print convergence information at each iteration.\n\nOutput\n\nSVD object of rank ≤ k.\n\nReferences\n\n@inproceedings{,\n    author={Friedland, S. and Niknejad, A. and Kaveh, Mostafa and Zare, H.},\n    booktitle={System of Systems Engineering, 2006 IEEE/SMC International Conference on},\n    title={Fast Monte-Carlo low rank approximations for matrices},\n    year={2006},\n    month={April},\n    pages={218--223},\n    doi={10.1109/SYSOSE.2006.1652299}\n}\n\n\n\n"
},

{
    "location": "library/public.html#Randomized-singular-value-decomposition-1",
    "page": "Public",
    "title": "Randomized singular value decomposition",
    "category": "section",
    "text": "reig\nrsvdfact\nrsvd_fnkz"
},

{
    "location": "library/public.html#Types-1",
    "page": "Public",
    "title": "Types",
    "category": "section",
    "text": ""
},

{
    "location": "library/public.html#IterativeSolvers.ConvergenceHistory",
    "page": "Public",
    "title": "IterativeSolvers.ConvergenceHistory",
    "category": "Type",
    "text": "Store general and in-depth information about an iterative method.\n\nFields\n\nmvps::Int: number of matrix vector products.\n\nmtvps::Int: number of transposed matrix-vector products\n\niters::Int: iterations taken by the method.\n\nrestart::T: restart relevant information.\n\nT == Int: iterations per restart.\nT == Void: methods without restarts.\n\nisconverged::Bool: convergence of the method.\n\ndata::Dict{Symbol,Any}: Stores all the information stored during the method execution. It stores tolerances, residuals and other information, e.g. ritz values in svdl.\n\nConstructors\n\nConvergenceHistory()\nConvergenceHistory(restart)\n\nCreate ConvergenceHistory with empty fields.\n\nArguments\n\nrestart: number of iterations per restart.\n\nPlots\n\nSupports plots using the Plots.jl package via a type recipe. Vectors are ploted as series and matrices as scatterplots.\n\nImplements\n\nBase: getindex, setindex!, push!\n\n\n\n"
},

{
    "location": "library/public.html#ConvergenceHistory-1",
    "page": "Public",
    "title": "ConvergenceHistory",
    "category": "section",
    "text": "ConvergenceHistory"
},

{
    "location": "library/internal.html#",
    "page": "Internal",
    "title": "Internal",
    "category": "page",
    "text": ""
},

{
    "location": "library/internal.html#Internal-Documentation-1",
    "page": "Internal",
    "title": "Internal Documentation",
    "category": "section",
    "text": "Documentation for IterativeSolvers.jl's internals.Pages = [\"internal.md\"]\nDepth = 4"
},

{
    "location": "library/internal.html#Index-1",
    "page": "Internal",
    "title": "Index",
    "category": "section",
    "text": "Pages = [\"internal.md\"]"
},

{
    "location": "library/internal.html#ConvergenceHistory-Internals-1",
    "page": "Internal",
    "title": "ConvergenceHistory Internals",
    "category": "section",
    "text": "TypealiasesIterativeSolvers.PlainHistory\nIterativeSolvers.RestartedHistoryFunctionsIterativeSolvers.nextiter!\nIterativeSolvers.reserve!\nIterativeSolvers.shrink!\nIterativeSolvers.setmvps\nIterativeSolvers.setmtvps\nIterativeSolvers.setconv\nIterativeSolvers.showplot"
},

{
    "location": "library/internal.html#Other-Functions-1",
    "page": "Internal",
    "title": "Other Functions",
    "category": "section",
    "text": "IterativeSolvers.idfact\nIterativeSolvers.isconverged\nIterativeSolvers.thickrestart!\nIterativeSolvers.harmonicrestart!\nIterativeSolvers.plot_collection\nIterativeSolvers.plotable\nIterativeSolvers.Adivtype\nIterativeSolvers.Amultype\nIterativeSolvers.randx\nIterativeSolvers.zerox\nIterativeSolvers.update!\nIterativeSolvers.initrand!(::Vector)"
},

{
    "location": "about/CONTRIBUTING.html#",
    "page": "Contributing",
    "title": "Contributing",
    "category": "page",
    "text": ""
},

{
    "location": "about/CONTRIBUTING.html#Contributing-1",
    "page": "Contributing",
    "title": "Contributing",
    "category": "section",
    "text": "Contributions are always welcome, as are feature requests and suggestions. Feel free to open issues and pull requests at any time. If you aren't familiar with git or Github please start now.It is important to note that almost every method in the package has documentation, to know what it does simply use ?<method> in the terminal.julia> using IterativeSolvers\n\nhelp?> IterativeSolvers.Adivtype\n  Adivtype(A, b)\n\n  Determine type of the division of an element of b against an element of A:\n\n  typeof(one(eltype(b))/one(eltype(A)))"
},

{
    "location": "about/CONTRIBUTING.html#Setting-workspace-up-1",
    "page": "Contributing",
    "title": "Setting workspace up",
    "category": "section",
    "text": "Julia's internal package manager makes it easy to install and modify packages from Github. Any package hosted on Github can be installed via Pkg.clone by providing the repository's URL, so installing a fork on your system is a simple task.Pkg.clone(\"https://github.com/johndoe/IterativeSolvers.jl\")It is to note here if you have the original package installed the fork will replace it, this is not a problem.Now find your fork's location.Pkg.dir(\"IterativeSolvers\")Once there you will notice you are on the master branch, whenever a package is imported Julia will use the code in the current branch, this means checking out other git branches will let you use/test whatever there is."
},

{
    "location": "about/CONTRIBUTING.html#Adding-or-modifying-iterative-methods-1",
    "page": "Contributing",
    "title": "Adding or modifying iterative methods",
    "category": "section",
    "text": "Each iterative method method must log information using the inner ConvergenceHistory type. When information is not necessary to be stored (plot is set to false) then instead of ConvergenceHistory create a DummyHistory, this type has the same calls ConvergenceHistory does but without storing anything.There are two types of ConvergenceHistory: plain and restarted. The only real difference between the two is how they are plotted and how the number of restarts is calculated, everything else is the same.Before logging information space must always be reserved.log = ConvergenceHistory()\nlog[:tol] = tol\nreserve!(log,:betas, maxiter) # Vector of length maxiter\nreserve!(log,:conv, maxiter, T=BitArray) # Vector of length maxiter\nreserve!(log,:ritz, maxiter, k) # Matrix of size (maxiter, k)To store information at each iteration use push!.push!(log, :conv, conv)\npush!(log, :ritz, F[:S][1:k])\npush!(log, :betas, L.β)To advance the log index to the next iteration use nextiter!.nextiter!(log)A more detailed explanation of all the functions is in both the public and internal documentation of ConvergenceHistory.The most rich example of the usage of ConvergenceHistory is in svdl."
},

{
    "location": "about/license.html#",
    "page": "License",
    "title": "License",
    "category": "page",
    "text": ""
},

{
    "location": "about/license.html#License-(MIT)-1",
    "page": "License",
    "title": "License (MIT)",
    "category": "section",
    "text": "Copyright (c) 2013--2016 The Julia Language\n\nPermission is hereby granted, free of charge, to any person obtaining a copy of\nthis software and associated documentation files (the \"Software\"), to deal in\nthe Software without restriction, including without limitation the rights to\nuse, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of\nthe Software, and to permit persons to whom the Software is furnished to do so,\nsubject to the following conditions:\n\nThe above copyright notice and this permission notice shall be included in all\ncopies or substantial portions of the Software.\n\nTHE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR\nIMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS\nFOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR\nCOPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER\nIN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN\nCONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE."
},

]}
