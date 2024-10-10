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

@testitem "Filter" begin
    
    using DFOQuadModel

    ℱ = []
    f = 2.0
    ϑ = 10.0
    ℱ, accepted = DFOQuadModel._update_filter!(ℱ, f, ϑ)

    @test ℱ == [(2.0, 10.0)]
    @test accepted == true

    function test_filter_update()
        ℱ = []
        test_points = [(2.0, 10.0, true),
                    (11.0, 2.0, true), 
                    (4.0, 3.0, true), 
                    (4.0, 3.0, false), 
                    (2.0, 2.0, true)]

        for (f, ϑ, boolean) in test_points
            ℱ, accepted = DFOQuadModel._update_filter!(ℱ, f, ϑ)
            @test accepted == boolean

            if accepted
                ℱ = DFOQuadModel.add_to_filter!(ℱ, f, ϑ)
            end

        end
    end

    test_filter_update()

end

@testitem "Rosenbrock" begin
    using DFOQuadModel

    function rosenbrock(x::Vector{Float64}; a::Float64=1.0, b::Float64=100.0)
        return (a - x[1])^2 + b * (x[2] - x[1]^2)^2
    end

    x0 = [-0.5, -0.5]
    dftr = DFOQuadModel.DFTR(rosenbrock, 
                            x0, 
                            [], 
                            []; 
                            max_iteration=200,
                            radius_factor_0=0.7,
                            radius_factor_1=0.9,
                            radius_factor_2=1.5,
                            ratio_threshold_1=1e-12,
                            ratio_threshold_2=1e-11,
                            stopping_radius=1e-16,
                            poisedness_threshold=100.0,
                            max_points="maximum")
    res = DFOQuadModel.optimize!(dftr)

    tol = 1e-2
    @test abs(res.x[1] - 1.0) < tol
    @test abs(res.x[2] - 1.0) < tol
    @test abs(res.of - 0.0) < tol
end