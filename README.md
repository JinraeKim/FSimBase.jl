# FSimBase
[FSimBase.jl](https://github.com/JinraeKim/FSimBase.jl) is
the lightweight base library for numerical simulation supporting nested dynamical systems and macro-based data logger,
compatible with [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl).

## Notes
- [FSimBase.jl](https://github.com/JinraeKim/FSimBase.jl) **works alone**!
For more functionality, see [FlightSims.jl](https://github.com/JinraeKim/FlightSims.jl).
- In [FSimBase.jl](https://github.com/JinraeKim/FSimBase.jl),
you must specify [the differential equation problem](https://diffeq.sciml.ai/stable/solvers/discrete_solve/#DiscreteProblems) and [the differential equation solver](https://diffeq.sciml.ai/stable/#Solver-Algorithms).
Only `ODEProblem` and `DiscreteProblem` are tested.

## APIs
Main APIs are provided in `src/APIs`.

### Make an environment
- `AbstractEnv`: an abstract type for user-defined and predefined environments.
In general, environments is a sub-type of `AbstractEnv`.
    ```julia
    struct LinearSystemEnv <: AbstractEnv
        A
        B
    end
    ```
- `State(env::AbstractEnv)`: return a function that produces structured states.
    ```julia
    function State(env::LinearSystemEnv)
        @unpack B = env
        n = size(B)[1]
        return function (x)
            @assert length(x) == n
            x
        end
    end
    ```
- `Dynamics!(env::AbstractEnv)`: return a function that maps in-place dynamics,
compatible with [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl).
User can extend these methods or simply define other methods.
    ```julia
    function Dynamics!(env::LinearSystemEnv)
        @unpack A, B = env
        @Loggable function dynamics!(dx, x, p, t; u)  # data would not be saved without @Loggable. Follow this form!
            @log state = x  # syntax sugar; macro-based logging
            @log input = u
            dx .= A*x + B*u
        end
    end
    ```
- (Optional) `Params(env::AbstractEnv)`: returns structured parameters of given environment `env`.

### Simulation, logging, and data saving & loading
**Simulator**
- `Simulator(state0, dyn, p; Problem, solver)` is a simulator struct that will be simulated by `solve` (non-interactive) or `step!` and `step_until!` (interactive).
`Problem = :ODE` and `Problem = :Discrete` imply [`ODEProblem`](https://diffeq.sciml.ai/stable/types/ode_types/) and [`DiscreteProblem`](https://diffeq.sciml.ai/stable/types/discrete_types/#Discrete-Problems), respectively.
For more details, see `src/APIs/simulation.jl`.

**Non-interactive interface (e.g., compatible with callbacks from DifferentialEquations.jl)**
- `solve(simulation::Simulator)` will solve (O)DE and provide `df::DataFrame`.
    - For now, only [**in-place** method (iip)](https://diffeq.sciml.ai/stable/basics/problem/#In-place-vs-Out-of-Place-Function-Definition-Forms) is supported.

**Interactive interface (you should be aware of how to use [`integrator` interface in DifferentialEquations.jl](https://diffeq.sciml.ai/stable/basics/integrator/#integrator))**
- `reinit!(simulator::Simulator)` will reinitialise `simulator::Simulator`.
- `step!(simulator::Simulator, Δt; stop_at_tdt=true)` will step the `simulator::Simulator` as `Δt`.
- `step_until!(simulator::Simulator, tf)` will step the `simulator::Simulator` until `tf`.
- `push!(simulator::Simulator, df::DataFrame)` will push a datum from `simulator` to `df`.

**Utilities**
- `apply_inputs(func; kwargs...)`
    - By using this, user can easily apply external inputs into environments. It is borrowed from [an MRAC example of ComponentArrays.jl](https://jonniedie.github.io/ComponentArrays.jl/stable/examples/adaptive_control/) and extended to be compatible with [SimulationLogger.jl](https://github.com/JinraeKim/SimulationLogger.jl).
    - (Limitations) for now, dynamical equations wrapped by `apply_inputs` will automatically generate logging function (even without `@Loggable`). In this case, all data will be an array of empty `NamedTuple`.
- Macros for logging data: `@Loggable`, `@log`, `@onlylog`, `@nested_log`, `@nested_onlylog`.
    - For more details, see [SimulationLogger.jl](https://github.com/JinraeKim/SimulationLogger.jl).

## Examples

### Custom environments, discrete problems, etc.
See directory `./test`.

### (TL; DR) Minimal example
```julia
using FSimBase
using DifferentialEquations
using ComponentArrays
using Test
using LinearAlgebra
using DataFrames


function main()
    state0 = [1.0, 2.0]
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
                          Problem=ODEProblem,
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
```
![ex_screenshot](./figures/custom_example.png)



## Related packages
### Dependencies
- [SimulationLogger.jl](https://github.com/JinraeKim/SimulationLogger.jl): A convenient logging tools compatible with [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl).
