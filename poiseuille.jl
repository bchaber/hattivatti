# physical input parameters
const ν = 0.1 # kinematic viscosity [m^2/s]
const H = 1.0 # channel height [m]
const ρ = 1.0 # density [kg/m^3]
const a = 0.8 # gravity [m/s^2]
const F = a*ρ # force per unit volume [N/m^3] 
# derived parameters
const um = a * H^2 / (8.0 * ν)    # maximum velocity [m/s]
const Re = a * H^3 / (8.0 * ν^2)  # Reynolds number [1]

# simulation parameters
const H̃  = 50  # height of the channel [l. u.]
const W̃  = 100 # width of the channel [l. u.]
const ρ̃  = 1.0 # average density [1]
const τ̃  = 0.8 # relaxation time [1]
const Δx̃ = 1.0 # assumed (not used)
const Δt̃ = 1.0 # assumed (not used)
const λ̃  = τ̃ * Δt̃

# conversion factors
CH = H / H̃                      # [m]
Cρ = ρ / ρ̃                      # [kg / m^3] 
Ct = (τ̃ - 0.5) / 3.0 * CH^2 / ν # [s]
Cu = CH / Ct                    # [m / s]
Cν = Cu * CH                    # [m^2 / s]
CF = Cρ * CH / Ct^2             # [kg / m^3 * m / s^2 = N / m^3]

# normalized parameters
const F̃  = F / CF
const ν̃  = ν / Cν
const ũm = um / Cu
const R̃e = ũm * H̃ / ν̃ # normalized Reynolds number matches the physical one
@assert R̃e ≈ Re  "Reynolds number should be preserved in lattice units"
@assert ũm < 1.0 "Population moves MORE than one cell in a timestep"
@assert ũm < 0.1 "The results might not be accurate"

using Einsum
# constants
const D        = 2 # spatial dimensions
const Q        = 9 # discrete velocities 
const ex       = [0, 1, 1, 0,-1,-1,-1, 0, 1]
const ey       = [0, 0, 1, 1, 1, 0,-1,-1,-1]
const west     = [2, 3, 9] # W, NW, SW
const east     = [5, 6, 7] # NE, E, SE
const north    = [3, 4, 5] # NW, N, NE
const south    = [7, 8, 9] # SE, S, SW
const weights  = [4/9,1/9,1/36,1/9,1/36,1/9,1/36,1/9,1/36]
const opposite = [1, 6, 7, 8, 9, 2, 3, 4, 5]

## initial conditions
# populations on lattice (PDFs)
const f  = ones(H̃, W̃, Q) * ρ / Q # normalized distribution
const f′ = similar(f)        # post-collision distribution
# derived, physical quantities on lattice
const rho = ones(H̃, W̃) * ρ
const ux  = zeros(H̃, W̃)
const uy  = zeros(H̃, W̃)

const fi  = zeros(H̃, W̃)
function stream()
    @inbounds for i in 1:Q
        copyto!(fi,view(f′, :, :, i))
        circshift!(view(f,  :, :, i), fi, (ex[i], ey[i]))
    end

    return nothing
end

const ω̃ = 1.0 / τ̃
function collide() # [Luo 1997]
    @inbounds for n in 1:H̃, m in 1:W̃
        uu = ux[n, m]^2 + uy[n, m]^2
        @inbounds for i in 1:Q
            eu  = ex[i] * ux[n,m] + ey[i] * uy[n,m]
            feq = weights[i] * rho[n, m] * (1. + 3. * eu + 9/2 * eu^2 - 3/2 * uu)
            fi  = f[n, m, i]
            Si  = 3. * weights[i] * ey[i] * F̃
            f′[n, m, i] = fi + ω̃ * (feq - fi) + Δt̃ * Si
        end
    end

    return nothing
end
function bounceback()
    @inbounds for n=1, m in 1:W̃, i in west
        f[n, m, i] = f′[n, m, opposite[i]]
    end

    @inbounds for n=H̃, m in 1:W̃, i in east
        f[n, m, i] = f′[n, m, opposite[i]]
    end

    return nothing
end

function moments()
    @einsum rho[n, m] = f[n, m, i]
    @einsum ux[n, m]  = f[n, m, i] * ex[i] / rho[n, m]
    @einsum uy[n, m]  = f[n, m, i] * ey[i] / rho[n, m]
end

function init()
    @inbounds for n in 1:H̃, m in 1:W̃
        @inbounds for i in 1:Q
            f[n, m, i] = weights[i] * ρ̃ 
        end
    end
end

function step()
    moments()
    collide()
    stream()
    bounceback()
    
    return nothing
end

# simulation and post-processing
using BenchmarkTools
using ProgressMeter
using UnicodePlots
# simulate
init()
@showprogress for t=1:50_000
    step()
end
@show maximum(Cu * uy)
# benchmark
t = @belapsed step()
MLUps = H̃ * W̃ / t * 1e-6
@show MLUps
