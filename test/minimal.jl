using FSimBase
using DifferentialEquations
using ComponentArrays
using Test
using LinearAlgebra
using DataFrames


function main()
    state0 = [1, 2]
    p = 1
    tf = 1.0
    Δt = 0.01
    @Loggable function dynamics!(dx, x, p, t)
        @log t
        @log x
        dx .= -p.*x
    end
    simulator = Simulator(
                          state0, dynamics!, p;
                          tf=tf, solver=Tsit5(),
                         )
    # solve approach (automatically reinitialised)
    @time _df = solve(simulator;
                     savestep=Δt,
                    )
    @time _df = solve(simulator;
                     savestep=Δt,
                    )
    # interactive simulation
    ## step!
    df = DataFrame()
    reinit!(simulator)
    step!(simulator, Δt)
    @test simulator.integrator.t ≈ Δt
    ## step_until!
    reinit!(simulator)
    ts_weird = 0:Δt:tf+Δt
    @time for t in ts_weird
        step_until!(simulator, t; df=df)
    end
    @show df
    @test norm(_df.sol[end].x - df.sol[end].x) < 1e-6
    @test simulator.integrator.t ≈ tf
end

@testset "minimal" begin
    main()
end
