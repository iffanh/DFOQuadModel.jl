struct DFTR
    f::Function
    x0::Vector{Float64}
    ci::Vector{Function}
    ce::Vector{Function}
    hyperparameters::Dict{Symbol, Any}
end

function DFTR(
                f,
                x0,
                ci,
                ce
              ;
                max_iteration::Int=100,
                poisedness_threshold::Float64=100.0,
                initial_radius::Float64=1.0,
                stopping_radius::Float64=1e-5,
                radius_factor_0::Float64=0.6,
                radius_factor_1::Float64=0.8,
                radius_factor_2::Float64=1.2,
                ratio_threshold_1::Float64=0.25,
                ratio_threshold_2::Float64=0.75
            )
    hyperparameters = Dict(
        :max_iteration => max_iteration,
        :poisedness_threshold => poisedness_threshold, 
        :initial_radius => initial_radius,
        :stopping_radius => stopping_radius,
        :radius_factor_0 => radius_factor_0,
        :radius_factor_1 => radius_factor_1,
        :radius_factor_2 => radius_factor_2,
        :ratio_threshold_1 => ratio_threshold_1,
        :ratio_threshold_2 => ratio_threshold_2 
    )
    return DFTR(f, x0, ci, ce, hyperparameters)
end


struct InterpolationSet
    X::Matrix{Float64}
    F
    Ci
    Ce
end

struct InterpolationModels
    IS::InterpolationSet
    mf
    mci
    mce
end

function evaluate(X::Matrix{Float64}, f, ci, ce)
    F = Vector{Float64}()
    Ci = []
    Ce = []
    for i in eachindex(X[:,1])
        f_, ci_, ce_ = evaluate_single_point(X[i,:], f, ci, ce)
        push!(F, f_)
        push!(Ci, ci_)
        push!(Ce, ce_)
    end

    Ci = reduce(hcat, Ci)
    Ce = reduce(hcat, Ce)

    return InterpolationSet(X, F, Ci, Ce)
end

function evaluate_single_point(x::Vector{Float64}, f, ci, ce)

    f_ = f(x...)
    ci_ = [c(x...) for c in ci]
    ce_ = [c(x...) for c in ce]
    return f_, ci_, ce_
end 

function build_models(IS::InterpolationSet, x)

    X = IS.X
    F = IS.F
    Ci = IS.Ci 
    Ce = IS.Ce 

    m = 2
    n = size(X)[2]

    mf = LagrangePoly.generate_lagrange_poly(n, m, X, F, x)

    mci = []
    if size(Ci, 1) == 0
    else
        for i in range(1,size(Ci,1))
            mci_ = LagrangePoly.generate_lagrange_poly(n, m, X, Ci[i,:], x)
            push!(mci, mci_)
        end
    end

    mce = []
    if size(Ce, 1) == 0
    else
        for i in range(1,size(Ce, 1))
            mce_ = LagrangePoly.generate_lagrange_poly(n, m, X, Ce[i,:], x)
            push!(mce, mce_)
        end 
    end

    return InterpolationModels(IS, mf, mci, mce)

end

function preparation(x0::Vector{Float64}, f, ci, ce, r0::Float64)

    # STEP 0: Preparation

    DynamicPolynomials.@polyvar x[1:size(x0)[1]]

    # Step 0.1 : Check all the hyperparameters

    # Step 0.2 : Add one more point to initialize the model with two points 
    
    δx = zeros(size(x0)[1])
    δx[1] = x0[1] + r0
    x1 = x0 + δx 
    X = vcat(x0', x1')

    # Step 0.3 : Evaluate the points
    IS = evaluate(X, f, ci, ce)

    # Step 0.4 : Build the initial models
    M = build_models(IS, x)

    return M
end

function improve_model!(M::InterpolationModels, r::Float64, Λth::Float64)

    X = M.IS.X
    n = size(X)[2]
    m = 2

    lpolys = LagrangePoly.generate_lagrange_bases(n, m, M.IS.X)
    lpolys, X = LagrangePoly.model_improvement!(lpolys, M.IS.X, r, Λth)
    return X
end

function optimize!(optimizer::DFTR)


    f = optimizer.f
    x0 = optimizer.x0
    ci = optimizer.ci
    ce = optimizer.ce

    max_iter = optimizer.hyperparameters[:max_iteration]
    Λth = optimizer.hyperparameters[:poisedness_threshold]
    r0 = optimizer.hyperparameters[:initial_radius]
    r_min = optimizer.hyperparameters[:stopping_radius]


    γ0 = optimizer.hyperparameters[:radius_factor_0]
    γ1 = optimizer.hyperparameters[:radius_factor_1]
    γ2 = optimizer.hyperparameters[:radius_factor_2]
    η1 = optimizer.hyperparameters[:ratio_threshold_1]
    η2 = optimizer.hyperparameters[:ratio_threshold_2]


    M = preparation(x0, f, ci, ce, r0)

    ## Main loop
    
    r = r0
    for k in range(1, max_iter)

        # STEP 1: Improve model
        println("STEP1")
        X = improve_model!(M, r, Λth)
u
        # STEP 2: Evaluate the points
        println("STEP2")
        IS = evaluate(X, f, ci, ce)
        M = build_models(IS, x)

        # STEP 3: Solve the subproblem
        println("STEP3")
        xinc = solve_subproblem(M, r)


    end

     

end

# function main()
# end

