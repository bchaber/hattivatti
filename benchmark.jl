using BenchmarkTools

function benchmark()
  t = @belapsed step()
  MLUps = 1e-6 * H̃ * W̃ / t
  @show MLUps
  return nothing
end
