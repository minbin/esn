using Plots
using DiffEqProblemLibrary
using DifferentialEquations
using QuadGK
using DelimitedFiles

function sdde(α,γr,C₀,τ,ß,R₀,tspan,st)
  function f(du, u, h, p, t)
    r = h(p, t)[1]
    rτ = h(p, t-τ)[1]
    du[1] = ((α * C₀^2)/(C₀ + rτ)^2) - ((γr * r)/(R₀ + r)) - (ß * r)
  end

  function h(p, t)
    α, γr, C₀, τ, ß, R₀ = p
    [-γr * t]
  end

  p = (α, γr, C₀, τ, ß, R₀)
  lags = [τ]

  prob = DDEProblem(f, h, tspan, p; constant_lags=lags)
  alg = MethodOfSteps(Rosenbrock23(autodiff=false))
  sol = solve(prob, saveat=st, alg; alg_hints=[:stiff])

  plot(sol, title = "Mass-action kinetics of delayed production and degradation", xlabel = "time", ylabel = "repressor", legend=false)

  savefig("base.png")
  writedlm("base.txt", sol.u)
end

sdde(300, 80, 10, 1, 0.1, 1, (0.0, 10000.0), 0.1)
