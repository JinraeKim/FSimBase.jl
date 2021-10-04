# FSimBase
[FSimBase.jl](https://github.com/JinraeKim/FSimBase.jl) is
the lightweight base library for numerical simulation supporting nested dynamical systems and macro-based data logger,
compatible with [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl).

## Notes
- [FSimBase.jl](https://github.com/JinraeKim/FSimBase.jl) **works alone**!
For more functionality, see [FlightSims.jl](https://github.com/JinraeKim/FlightSims.jl).
- In [FSimBase.jl](https://github.com/JinraeKim/FSimBase.jl),
you must specify [the differential equation solver](https://diffeq.sciml.ai/stable/#Solver-Algorithms).

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
- `Simulator(state0, dyn, p)` is a simulator struct that will be simulated by `solve` (non-interactive) or `step!` and `step_until!` (interactive).
For more details, see `src/APIs/simulation.jl`.

**Non-interactive interface (e.g., compatible with callbacks from DifferentialEquations.jl)**
- `solve(simulation::Simulator)` will solve (O)DE and provide `df::DataFrame`.
    - For now, only [**in-place** method (iip)](https://diffeq.sciml.ai/stable/basics/problem/#In-place-vs-Out-of-Place-Function-Definition-Forms) is supported.

**Interactive interface (similar to `integrator` interface in DifferentialEquations.jl)**
- `reinit!(simulator::Simulator)` will reinitialise `simulator::Simulator`.
- `step!(simulator::Simulator, Δt, df=nothing; stop_at_tdt=true)` will step the `simulator::Simulator` as `Δt`.
**Data will be pushed into `df::DataFrame` if the argument `df` is provided.**
- `step_until!(simulator::Simulator, tf, df=nothing)` will step the `simulator::Simulator` until `tf`.
**Data will be pushed into `df::DataFrame` if the argument `df` is provided.**

**Utilities**
- `apply_inputs(func; kwargs...)`
    - By using this, user can easily apply external inputs into environments. It is borrowed from [an MRAC example of ComponentArrays.jl](https://jonniedie.github.io/ComponentArrays.jl/stable/examples/adaptive_control/) and extended to be compatible with [SimulationLogger.jl](https://github.com/JinraeKim/SimulationLogger.jl).
    - (Limitations) for now, dynamical equations wrapped by `apply_inputs` will automatically generate logging function (even without `@Loggable`). In this case, all data will be an array of empty `NamedTuple`.
- Macros for logging data: `@Loggable`, `@log`, `@onlylog`, `@nested_log`, `@nested_onlylog`.
    - For more details, see [SimulationLogger.jl](https://github.com/JinraeKim/SimulationLogger.jl).

## Examples

### (TL; DR) Minimal example
```julia
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
                          state0, dynamics!, Tsit5(), p;
                          tf=tf,
                         )
    # solve approach (automatically reinitialised)
    @time _df = solve(simulator; savestep=Δt)
    # interactive simulation
    ## step!
    reinit!(simulator)
    step!(simulator, Δt)
    @test simulator.integrator.t ≈ Δt
    ## step_until!
    ts_weird = 0:Δt:tf+Δt
    df = DataFrame()
    reinit!(simulator)
    @time for t in ts_weird
        step_until!(simulator, t, df)  # log data at the third element
    end
    @show df
    @test norm(_df.sol[end].x - df.sol[end].x) < 1e-6
    @test simulator.integrator.t ≈ tf
end

@testset "minimal" begin
    main()
end
```

### (TL; DR) An example with custom environments
```julia
using FSimBase

using LinearAlgebra  # for I, e.g., Matrix(I, n, n)
using ComponentArrays
using UnPack
using Transducers
using Plots
using DifferentialEquations
using Test


struct MyEnv <: AbstractEnv  # AbstractEnv exported from FSimBase
    a
    b
end

"""
FlightSims recommends you to use closures for State and Dynamics!. For more details, see https://docs.julialang.org/en/v1/devdocs/functions/.
"""
function State(env::MyEnv)
    return function (x1::Number, x2::Number)
        ComponentArray(x1=x1, x2=x2)
    end
end

function Dynamics!(env::MyEnv)
    @unpack a, b = env  # @unpack is very useful!
    @Loggable function dynamics!(dx, x, p, t; u)  # `Loggable` makes it loggable via SimulationLogger.jl (imported in FSimBase)
        @unpack x1, x2 = x
        @log x1  # to log x1
        @log x2  # to log x2
        dx.x1 = a*x2
        dx.x2 = b*u
    end
end

function my_controller(x, p, t)
    @unpack x1, x2 = x
    -(x1+x2)
end


function main()
    n = 2
    m = 1
    a, b = 1, 1
    A = -Matrix(I, n, n)
    B = Matrix(I, m, m)
    env = MyEnv(a, b)
    tf = 10.0
    Δt = 0.01
    x10, x20 = 10.0, 0.0
    x0 = State(env)(x10, x20)
    # simulator
    simulator = Simulator(
                          x0, apply_inputs(Dynamics!(env); u=my_controller), Tsit5();
                          tf=tf,
                         )
    @time df = solve(simulator; savestep=Δt)
    ts = df.time
    x1s = df.sol |> Map(datum -> datum.x1) |> collect
    x2s = df.sol |> Map(datum -> datum.x2) |> collect
    # plot
    p_x1 = plot(ts, x1s; label="x1")
    p_x2 = plot(ts, x2s; label="x2")
    p_x = plot(p_x1, p_x2, layout=(2, 1))
    # save
    dir_log = "figures"
    mkpath(dir_log)
    savefig(p_x, joinpath(dir_log, "custom_example.png"))
    display(p_x)
end

@testset "custom_example" begin
    main()
end
```
![ex_screenshot](./figures/custom_example.png)



## Related packages
### Dependencies
- [SimulationLogger.jl](https://github.com/JinraeKim/SimulationLogger.jl): A convenient logging tools compatible with [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl).
