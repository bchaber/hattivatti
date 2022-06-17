using BenchmarkTools
include(first(ARGS))

t = @belapsed step()
MLUps = 1e-6 * H̃ * W̃ / t
@show MLUps