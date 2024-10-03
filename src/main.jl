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

function evaluate(x::Vector{Float64}, f, ci, ce)

    f_ = f(x...)
    ci_ = [c(x...) for c in ci]
    ce_ = [c(x...) for c in ce]
    return f_, ci_, ce_
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


    # STEP 0: Preparation
    # Step 0.1 : Check all the hyperparameters
    # Step 0.2 : Add one more point to initialize the model with two points 
    
    δx = zeros(size(x0)[1])
    δx[1] = x0[1] + r0
    x1 = x0 + δx 
    X = vcat(x0', x1')

    # Step 0.3 : Evaluate the points

    F = []
    Ci = []
    Ce = []
    for i in eachindex(X[:,1])
        f_, ci_, ce_ = evaluate(X[i,:], f, ci, ce)
        println(f_, ci_, ce_)
    end

    

    # Step 0.4 : Build the initial models


    ## Main loop
    # STEP 1: improve model 


end

function main()
end

objective = (x, y) -> (x-1)^2 + (y-2)^2
ineq_constraints = [(x,y) -> x + y - 2]
eq_constraints = [(x,y) -> x - y]
x0 = [0.0, 0.0]

dftr = DFTR(objective, x0, ineq_constraints, Vector{Float64}();)

optimize!(dftr)