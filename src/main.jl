struct DFTR
    f::Function
    x0::Vector{Float64}
    ci::Vector{Any}
    ce::Vector{Any}
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
                ratio_threshold_2::Float64=0.75,
                reduction_constant::Float64=0.5,
                max_points::String="default",
            )

    n = size(x0)[1]
    if max_points == "default"
        max_points = 2*n+1
    elseif max_points == "minimum"
        max_points = n+1
    elseif max_points == "maximum"
        max_points = round(Int, 0.5*((n+1)*(n+2)))
    end

    hyperparameters = Dict(
        :max_iteration => max_iteration,
        :poisedness_threshold => poisedness_threshold, 
        :initial_radius => initial_radius,
        :stopping_radius => stopping_radius,
        :radius_factor_0 => radius_factor_0,
        :radius_factor_1 => radius_factor_1,
        :radius_factor_2 => radius_factor_2,
        :ratio_threshold_1 => ratio_threshold_1,
        :ratio_threshold_2 => ratio_threshold_2,
        :reduction_constant => reduction_constant,
        :max_points => max_points, 
    )
    return DFTR(f, x0, ci, ce, hyperparameters)
end


struct InterpolationSet
    X::Matrix{Float64}
    F
    Ci
    Ce
    ϑs
end

struct InterpolationModels
    IS::InterpolationSet
    mf
    mci
    mce
end

function add_to_filter!(ℱ::Vector{Any}, f::Float64, ϑ::Float64)
    push!(ℱ, (f, ϑ))
    return ℱ
end

function _update_filter!(ℱ::Vector{Any}, f::Float64, ϑ::Float64)

    if length(ℱ) == 0
        push!(ℱ, (f, ϑ))
        accepted = true
        return ℱ, accepted
    else
        accepted = false
        γᵥ = 0.5
        conds = []
        for (fⱼ, ϑⱼ) in ℱ
            if ϑ < (1 - γᵥ)*ϑⱼ || f < fⱼ - γᵥ*ϑ
                cond = true
            else
                cond = false
            end
            append!(conds, cond)
        end 
        accepted = all(conds)
        # if accepted
        #     push!(ℱ, (f, ϑ))
        # end

        return ℱ, accepted
    end

end

function update_filter!(ℱ::Vector{Any}, IS::InterpolationSet)
    F = IS.F
    ϑs = IS.ϑs

    for (f, ϑ) in zip(F, ϑs)
        ℱ, accepted = _update_filter!(ℱ, f, ϑ)
    end

    return ℱ
end

function calculate_violation(ci::Vector{Any}, ce::Vector{Any})
    ϑ = max(0.0, max(-ci..., 0.0), max(abs.(ce)..., 0.0))
    return ϑ

end

function evaluate(X::Matrix{Float64}, f, ci, ce)
    F = Vector{Float64}()
    Ci = []
    Ce = []
    ϑs = Vector{Float64}()
    for i in eachindex(X[:,1])
        f_, ci_, ce_, ϑ = evaluate_single_point(X[i,:], f, ci, ce)
        push!(F, f_)
        push!(Ci, ci_)
        push!(Ce, ce_)
        push!(ϑs, ϑ)
    end

    Ci = reduce(hcat, Ci)
    Ce = reduce(hcat, Ce)

    sorted_indices = reorder_samples!(F, ϑs)
    X = X[sorted_indices,:]
    Ci = Ci[:,sorted_indices]
    Ce = Ce[:,sorted_indices]
    ϑs = ϑs[sorted_indices]

    return InterpolationSet(X, F, Ci, Ce, ϑs)
end

function evaluate_single_point(x::Vector{Float64}, f, ci, ce)

    f_ = f(x)
    ci_ = [c(x) for c in ci]
    ce_ = [c(x) for c in ce]
    ϑ = calculate_violation(ci_, ce_)
    return f_, ci_, ce_, ϑ
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

function preparation(x0::Vector{Float64}, f, ci, ce, r0::Float64, x)

    # STEP 0: Preparation

    # Step 0.1 : Check all the hyperparameters

    # Step 0.2 : Add one more point to initialize the model with two points 
    
    δx = zeros(size(x0)[1])
    δx[1] = x0[1] + r0
    x1 = x0 + δx 
    X = vcat(x0', x1')

    # Step 0.3 : Evaluate the points
    IS = evaluate(X, f, ci, ce)

    # Step 0.4 : Initialize Filter
    ℱ = []
    ℱ = update_filter!(ℱ, IS)

    # Step 0.5 : Build the initial models
    M = build_models(IS, x)

    return M, ℱ
end

function improve_model!(X, r::Float64, Λth::Float64)

    n = size(X)[2]
    m = 2

    lpolys = LagrangePoly.generate_lagrange_bases(n, m, X)
    lpolys, X = LagrangePoly.model_improvement!(lpolys, X, r, Λth)
    return X
end

function update_trust_region_radius!(Δ, ρ, η₁, η₂, γ₀, γ₁, γ₂)
    if ρ < η₁
        return γ₀*Δ
    elseif ρ < η₂
        return γ₁*Δ
    else
        return γ₂*Δ 
    end
end

function reorder_samples!(F, ϑs)

    indices = sortperm(1:length(F), by = i -> (ϑs[i], F[i]))

    return indices
end

function delete_point!(X::Matrix{Float64}, Δ::Float64, pₘ::Int64)
    #= 
    When there's too many points, we remove those that contribute to 
    the largest poisedness value.
    =#

    n = size(X)[2]
    m = 2

    if size(X)[1] > pₘ
        lpolys = LagrangePoly.generate_lagrange_bases(n, m, X)
        Λs = LagrangePoly.compute_poisedness(X[1,:], Δ, lpolys)
        index_to_delete = max(sortperm(Λs)...)
        X = X[1:end .!= index_to_delete, :]
        return X
    end
    
    return X
end

function optimize!(optimizer::DFTR)


    f = optimizer.f
    x0 = optimizer.x0
    ci = optimizer.ci
    ce = optimizer.ce

    max_iter = optimizer.hyperparameters[:max_iteration]
    Λth = optimizer.hyperparameters[:poisedness_threshold]
    Δ₀ = optimizer.hyperparameters[:initial_radius]
    r_min = optimizer.hyperparameters[:stopping_radius]


    γ₀ = optimizer.hyperparameters[:radius_factor_0]
    γ₁ = optimizer.hyperparameters[:radius_factor_1]
    γ₂ = optimizer.hyperparameters[:radius_factor_2]
    η₁ = optimizer.hyperparameters[:ratio_threshold_1]
    η₂ = optimizer.hyperparameters[:ratio_threshold_2]
    κᵥ = optimizer.hyperparameters[:reduction_constant]
    pₘ = optimizer.hyperparameters[:max_points]

    DynamicPolynomials.@polyvar x[1:size(x0)[1]]

    M, ℱ = preparation(x0, f, ci, ce, Δ₀, x)
    
    ## Main loop
    X = M.IS.X
    Δ = Δ₀
    for k in range(1, max_iter)

        # STEP 1: Improve model
        X = improve_model!(X, Δ, Λth)

        # println(X, size(X))

        # STEP 2: Evaluate the points
        IS = evaluate(X, f, ci, ce)
        M = build_models(IS, x)
        ℱ = update_filter!(ℱ, IS)

        xbest = IS.X[1,:]
        fbest = IS.F[1,:]
        ϑbest = IS.ϑs[1,:]
        npoint = size(IS.X)[1]
        println("It. $k, xbest = $xbest, fbest = $fbest, ϑbest = $ϑbest, Δ = $Δ, #points = $npoint")

        # STEP 3: Solve the subproblem and evaluate
        xinc = solve_subproblem(M, Δ)
        f_, ci_, ce_, ϑ_ = evaluate_single_point(xinc, f, ci, ce)
        ℱ, accepted = _update_filter!(ℱ, f_, ϑ_)

        println(xinc, f_, ϑ_)

        
        if !accepted
            Δ = γ₀*Δ

            X = vcat(X, xinc')
            X = delete_point!(X, Δ, pₘ)

        else
            if M.mf(X[1,:]) - M.mf(xinc) ≥ κᵥ*IS.ϑs[1]
                ℱ = add_to_filter!(ℱ, f_, ϑ_)
            end
            X = vcat(X, xinc')
            X = delete_point!(X, Δ, pₘ)

            # STEP 4: Update trust-region radius
            ρₖ = (IS.F[1] - f_)/(M.mf(X[1,:]) - M.mf(xinc))
            Δ = update_trust_region_radius!(Δ, ρₖ, η₁, η₂, γ₀, γ₁, γ₂)
        end

    end

end

# function main()
# end

