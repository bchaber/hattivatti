# physical input parameters
const ν = 1e-1 # kinematic viscosity [m^2/s]
const H = 1.00 # channel height [m]
const ρ = 1e+0 # density [kg/m^3]
const g = 0.8  # gravity [m/s^2]

# derived parameters
const um = g * H^2 / (8.0 * ν)    # maximum velocity [m/s]
const Re = g * H^3 / (8.0 * ν^2)  # Reynolds number [1]

# simulation parameters
const H̃ = 50   # height of the channel [l. u.]
const W̃ = 100  # width of the channel [l. u.]
const ρ̃ = 1.0  # average density [1]
const τ = 0.8  # relaxation time [1]

# conversion factors
CH = H / H̃                      # [m]
Cρ = ρ / ρ̃                      # [kg / m^3] 
Ct = (τ - 0.5) / 3.0 * CH^2 / ν # [s]
Cu = CH / Ct                    # [m / s]
Cν = Cu * CH                    # [?]
Cg = Cρ * CH / Ct^2             # [?]

# normalized parameters
const g̃  = g / Cg
const ν̃  = ν / Cν
const ũm = um / Cu
const R̃e = ũm * H̃ / ν̃ # normalized Reynolds number matches the physical one
@assert ũm < 1.0 "Population moves MORE than one cell in a timestep"

using StaticArrays
using Einsum
# constants
const D        = 2 # spatial dimensions
const Q        = 9 # discrete velocities 
const ex       = @SVector[0, 1, 1, 0,-1,-1,-1, 0, 1]
const ey       = @SVector[0, 0, 1, 1, 1, 0,-1,-1,-1]
const west     = @SVector[2, 3, 9]
const east     = @SVector[5, 6, 7]
const north    = @SVector[3, 4, 5]
const south    = @SVector[7, 8, 9]
const weights  = @SVector[4/9,1/9,1/36,1/9,1/36,1/9,1/36,1/9,1/36]
const opposite = @SVector[1, 6, 7, 8, 9, 2, 3, 4, 5]

# populations on lattice (PDFs)
const f  = ones(H̃, W̃, Q) / Q # normalized distribution
const f′ = similar(f)        # post-collision distribution

# derived quantities on lattice
const rho = zeros(H̃, W̃)
const ux  = zeros(H̃, W̃)
const uy  = zeros(H̃, W̃)

const fi  = zeros(H̃, W̃)
function stream()
    for i in 1:Q
        copyto!(fi,view(f, :, :, i))
        circshift!(view(f, :, :, i), fi, (ex[i], ey[i]))
    end

    return nothing
end

function collide()
    for n in 1:H̃, m in 1:W̃
        vy = τ * g̃ * rho[n, m]
        uu = ux[n, m]^2 + uy[n, m]^2 + vy^2
        
        for i in 1:Q
            eu = (ex[i] * ux[n,m] + ey[i] * uy[n,m] + ey[i] * vy)
            fi = f[n, m, i]
            feq = rho[n, m] * weights[i] *
                (1. + 3. * eu + 4.5 * eu^2 - 1.5 * uu)
            f[n, m, i] -= (1.0 / τ) * (fi - feq)
        end
    end

    return nothing
end

function boundary()
    for n=1, m in 1:W̃, i in east
        f′[n, m, i] = f[n, m, opposite[i]]
    end # bounce-back from east wall

    for n=H̃, m in 1:W̃, i in west
        f′[n, m, i] = f[n, m, opposite[i]]
    end # bounce-back from west wall

    for n=1, m in 1:W̃, i in east
        f[n, m, i] = f′[n, m, i]
    end

    for n=H̃, m in 1:W̃, i in west
        f[n, m, i] = f′[n, m, i]
    end

    return nothing
end

function step()
    collide()
    stream()
    boundary()

    @einsum rho[n, m] = f[n, m, i]
    @einsum ux[n, m]  = f[n, m, i] * ex[i] / rho[n, m]
    @einsum uy[n, m]  = f[n, m, i] * ey[i] / rho[n, m]
    
    return nothing
end

# simulation and post-processing
using BenchmarkTools
using ProgressMeter
using UnicodePlots

@showprogress for t=1:500 step() end
show(lineplot(Cu * view(uy, :, 5),
              title="Physical y-velocity profile"))
show(extrema(Cu * uy))