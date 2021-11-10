using FSimBase
using DifferentialEquations
using ComponentArrays
using Test
using LinearAlgebra
using DataFrames


function main()
    state0 = [1.0, 2.0]
    tf = 10
    Δt = 1
    @Loggable function dynamics!(dx, x, p, t)
        @log t
        @log x
        dx .= 0.99*x
    end
    simulator = Simulator(
                          state0, dynamics!;
                          tf=tf,
                          Problem=DiscreteProblem,
                          solver=FunctionMap(),
                         )
    # solve approach (automatically reinitialised)
    @time _df = solve(simulator; savestep=Δt)
    @show _df
    # interactive simulation
    ## step_until!
    df = DataFrame()
    reinit!(simulator)
    ts = 0:Δt:tf
    @time for t in ts
        step_until!(simulator, t)
        push!(simulator, df)  # flag=true is default
        # safer way:
        # flag = step_until!(simulator, t)
        # push!(simulator, df, flag)
        # or, equivalently,
        # push!(simulator, df, step_until!(simulator, t))  # compact form
    end
    @test df == _df
end
