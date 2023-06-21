using ProgressMeter

function simulate(timesteps)
  init()
  @showprogress for t=1:timesteps
    step()
  end
  evaluate()
  return nothing
end

""" Useful relations

cs = Δx̃^2 / √(3.0)/ Δt̃^2   # speed of sound [l. u.]
ν  = cs^2 * Δt̃ * (τ̃ - 0.5) # kinematic viscosity [l. u.]
"""
