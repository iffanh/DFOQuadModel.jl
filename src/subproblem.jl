function solve_subproblem(M, r)

    X = M.IS.X
    center = X[1,:]

    solver = optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true)
    model = SumOfSquares.SOSModel(solver)

    JuMP.@variable(model, α)
    JuMP.@objective(model, Max, α)

    xf = DynamicPolynomials.variables(M.mf)
    x = DynamicPolynomials.variables(M.mci[1])

    the_sum = sum((x[j] - center[j])^2 for j in eachindex(center))
    tr_expr = r^2 - the_sum

    # Build the domain
    n_ineq = length([i for i in eachindex(M.mci)])
    n_eq = length([i for i in eachindex(M.mce)]) 

    if n_ineq > 0 && n_eq > 0
        eq = [ce for ce in M.mce]
        ineq = [ci for ci in M.mci]
        append!(ineq, [tr_expr])
        S = SemialgebraicSets.basic_semialgebraic_set(SemialgebraicSets.algebraic_set(eq), ineq)  
    elseif n_eq > 0
        eq = [ce for ce in M.mce]
        S = SemialgebraicSets.basic_semialgebraic_set(SemialgebraicSets.algebraic_set(eq), [tr_expr])  
    elseif n_ineq > 0
        ineq = [ci for ci in M.mci]
        append!(ineq, [tr_expr])
        S = SemialgebraicSets.basic_semialgebraic_set(SemialgebraicSets.FullSpace(), ineq)
    else
        S = SemialgebraicSets.basic_semialgebraic_set(SemialgebraicSets.FullSpace(), [tr_expr])  
    end

    JuMP.@constraint(model, cref, M.mf >= α, domain = S, maxdegree = 4)
    JuMP.optimize!(model)

    θ = moment_matrix(cref);
    η = atomic_measure(θ, 1e-3);    
    x_opt = η.atoms[1].center

    return x_opt

end


    
