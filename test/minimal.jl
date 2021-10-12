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
                          solver=Tsit5(),
                          tf=tf,
                         )
    # solve approach (automatically reinitialised)
    @time _df = solve(simulator; savestep=Δt)
    # interactive simulation
    ## step!
    reinit!(simulator)
    step!(simulator, Δt)
    @test simulator.integrator.t ≈ Δt
    ## step_until! (callback-like)
    ts_weird = 0:Δt:tf+Δt
    df_ = DataFrame()
    reinit!(simulator)
    @time for t in ts_weird
        flag = step_until!(simulator, t)  # flag == false if step is inappropriate
        if simulator.integrator.u[1] < 5e-1
            break
        else
            push!(simulator, df_, flag)  # push data only when flag == true
        end
    end
    println(df_[end-5:end, :])
    ## step_until!
    df = DataFrame()
    reinit!(simulator)
    @time for t in ts_weird
        step_until!(simulator, t)
        push!(simulator, df)  # flag=true is default
        # safer way:
        # flag = step_until!(simulator, t)
        # push!(simulator, df, flag)
        # or, equivalently,
        # push!(simulator, df, step_until!(simulator, t))  # compact form
    end
    println(df[end-5:end, :])
    @test norm(_df.sol[end].x - df.sol[end].x) < 1e-6
    @test simulator.integrator.t ≈ tf
end

@testset "minimal" begin
    main()
end
