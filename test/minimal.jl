using FSimBase
const FSBase = FSimBase
using DifferentialEquations
using SimulationLogger: @log

using ComponentArrays


function main()
    state0 = [1, 2]
    p = 1
    @Loggable function dynamics!(dx, x, p, t)
        @log t
        @log x
        dx .= -p.*x
    end
    prob, df = FSimBase.sim(state0, dynamics!, p;
                            solver=Tsit5(),
                            tf=10.0,
                           )
    df
end