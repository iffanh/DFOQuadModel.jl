using TestItemRunner
using TestItems
using Test 

@run_package_tests verbose = true

@testitem "Subproblem 1" begin
    
    using DFOQuadModel
    using DynamicPolynomials

    f = (a...) -> a[1]^2 + a[2]^2
    ci1 = (a...) -> a[1] + a[2] - 0.5
    ci2 = (a...) -> a[1] - a[2] - 1

    X = [0.0 0.0; 
        1.0 0.0;
        0.0 1.0; 
        2.0 0.0;
        1.0 1.0; 
        0.0 2.0];

    DynamicPolynomials.@polyvar x[1:size(X)[2]]

    IS = DFOQuadModel.evaluate(X, f, [ci1, ci2], [])
    M = DFOQuadModel.build_models(IS, x)

    x_opt = DFOQuadModel.solve_subproblem(M, 1.0)

    tol = 1e-5
    @test abs(x_opt[1] - 0.75) < tol 
    @test abs(x_opt[2] - -0.25) < tol 

end