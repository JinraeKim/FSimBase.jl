"""
Initialise simulator.
"""
function _initialise(state0, __dyn, p, t0, tf)
    tspan = (t0, tf)
    iip = true
    prob = ODEProblem{iip}(__dyn, state0, tspan, p)  # true: isinplace
    prob
end

"""
Extended to deal with empty Dict.
"""
function _namedtuple(x::Dict)
    if x == Dict()
        return NamedTuple()  # empty NamedTuple
    else
        return namedtuple(x)
    end
end

function recursive_namedtuple(x::Any)
    x
end

function recursive_namedtuple(d::Dict)
    _namedtuple(
                Dict(k => recursive_namedtuple(v) for (k, v) in d)
               )
end


"""
Simulator struct
"""
struct Simulator
    integrator::DEIntegrator
    log_func::Union{Function, Nothing}
    # df::DataFrame
    function Simulator(state0, dyn, p=nothing;
            t0=0.0, tf=1.0, solver=nothing, kwargs...,
        )
        # DEProblem
        __dyn = (dx, x, p, t) -> dyn(dx, x, p, t)
        prob = _initialise(state0, __dyn, p, t0, tf)
        log_func = nothing
        # save func (logging func)
        if hasmethod(dyn, Tuple{Any, Any, Any, Any, __LOG_INDICATOR__})
            __log_indicator__ = __LOG_INDICATOR__()  # just an indicator for logging
            log_func = function (x, t, integrator::DEIntegrator; kwargs...)
                x = copy(x)  # `x` merely denotes a "view"
                dyn(zero.(x), x, integrator.p, t, __log_indicator__; kwargs...)
            end
        end
        integrator = init(prob, solver; kwargs...)
        integrator
        new(integrator, log_func)
    end
end

"""
Reinitialise simulator.
"""
function DiffEqBase.reinit!(simulator::Simulator, args...; kwargs...)
    DiffEqBase.reinit!(simulator.integrator, args...; kwargs...)
    simulator
end

"""
Step `dt` time.
"""
function DiffEqBase.step!(simulator::Simulator, Δt, df=nothing; stop_at_tdt=true)
    DiffEqBase.step!(simulator.integrator, Δt, stop_at_tdt)
    if df !=nothing
        log!(df, simulator)
    end
end

"""
Step until `tf`.
"""
function step_until!(simulator::Simulator, tf, df=nothing;
        suppress_termination_warn=true,
        suppress_truncation_warn=false,
    )
    integrator = simulator.integrator
    t = integrator.t
    _tf = integrator.sol.prob.tspan[2]
    if (simulator.integrator.tdir > 0 && _tf < tf) || (simulator.integrator.tdir < 0 && _tf > tf)
        if !suppress_truncation_warn
            @warn("step! truncated up to tf = $(_tf)")
        end
        tf = _tf
    end
    if tf ≈ t
        if !suppress_termination_warn
            @warn("step! ignored; simulator seems already terminated")
        end
    else
        DiffEqBase.step!(simulator, tf - t)
        if df !=nothing
            log!(df, simulator)
        end
    end
end

function _push!(df::DataFrame, integrator::DEIntegrator, log_func)
    x = integrator.u
    t = integrator.t
    __log_dict__ = log_func(x, t, integrator)
    __log_nt__ = recursive_namedtuple(__log_dict__)
    push!(df, (; time=t, sol=__log_nt__))
end

function log!(df::DataFrame, simulator::Simulator)
    integrator = simulator.integrator
    log_func = simulator.log_func
    _push!(df, integrator, log_func)
end


"""
    solve(simulator; ...)

`solve` for `simulator::Simulator` (similar to that of DifferentialEquations.jl).
This method will automatically reinitialise `simulator`.
"""
function DiffEqBase.solve(simulator::Simulator;
        saveat=nothing,
        savestep=0.01,
        callback::DECallback=CallbackSet(),
        kwargs...,
    )
    DiffEqBase.reinit!(simulator)
    prob = simulator.integrator.sol.prob
    solver = simulator.integrator.alg
    log_func = simulator.log_func
    t0, tf = prob.tspan[1], prob.tspan[2]
    if saveat == nothing && savestep == nothing
        saveat = []  # default
    elseif saveat != nothing && savestep == nothing
        # nothing
    elseif saveat == nothing && savestep != nothing
        saveat = t0:savestep:tf
    elseif saveat != nothing && savestep != nothing
        error("Assign values of either `saveat` or `savestep`")
    end
    saved_values = SavedValues(Float64, Dict)
    if log_func != nothing
        cb_save = SavingCallback(log_func, saved_values;
                                 saveat=saveat, tdir=Int(sign(tf-t0)))
        callback = CallbackSet(callback, cb_save)  # save values "after all other callbacks"
    end
    _ = solve(prob, solver;
                callback=callback,
                kwargs...)
    df = DataFrame(
                   time = saved_values.t,
                   sol = saved_values.saveval |> Map(recursive_namedtuple) |> collect,
                  )
    df
end


"""
# Notes
- The basic concept is borrowed from [an MRAC example](https://jonniedie.github.io/ComponentArrays.jl/stable/examples/adaptive_control/).
- It is modified to be compatible with [SimulationLogger.jl](https://github.com/JinraeKim/SimulationLogger.jl).
# Limitations
- Conditional method definition is troublesome; see [#16](https://github.com/JinraeKim/FlightSims.jl/issues/16).
Instead of it, I decided to merely define a new method with argument `__log_indicator__::__LOG_INDICATOR__`,
which will provide "empty Dict" in the case of no logging.
"""
maybe_apply(f::Function, x, p, t) = f(x, p, t)
maybe_apply(f, x, p, t) = f
function apply_inputs(func; kwargs...)
    function simfunc(dx, x, p, t)
        func(dx, x, p, t; map(f -> maybe_apply(f, x, p, t), (; kwargs...))...)
    end
    function simfunc(dx, x, p, t, __log_indicator__::__LOG_INDICATOR__)
        if hasmethod(func, Tuple{Any, Any, Any, Any, __LOG_INDICATOR__})
            return func(dx, x, p, t, __log_indicator__; map(f -> maybe_apply(f, x, p, t), (; kwargs...))...)
        else
            return Dict()  # see the above method `NamedTupleTools.namedtuple`
        end
    end
    simfunc
end
