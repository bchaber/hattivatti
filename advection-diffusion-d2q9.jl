# physical input parameters
#const ν = 0.1 # kinematic viscosity [m^2/s]
const H  = 1.0 # channel height [m]
const ρ  = 1.0 # density [kg/m^3]
const D  = 1.5 # diffiusion [m^2/s]
# derived parameters
const τg = 3.0D + 1/2

# simulation parameters
const NT = 200 # number of timesteps [1]
const H̃  = 512 # height of the channel [l. u.]
const W̃  = 512 # width of the channel [l. u.]
const ρ̃  = 1.0 # average density [1]
const σ0 = 10. # initial width [1]
const Δx̃ = 1.0 # assumed (not used)
const Δt̃ = 1.0 # assumed (not used)
const τ̃  = τg/Δt̃ # relaxation time [1]
const λ̃  = τg

# conversion factors
CH = H / H̃                      # [m]
Cρ = ρ / ρ̃                      # [kg / m^3] 

# normalized parameters
const ũ = 0.0
@assert ũ < 1.0 "Population moves MORE than one cell in a timestep"

# constants
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
const ux = [ũ for n=1:H̃, m=1:W̃]
const uy = [ũ for n=1:H̃, m=1:W̃]
# derived, physical quantities on lattice
const rho = ones(H̃, W̃) * ρ

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
        uu = ux[n,m]^2 + uy[n,m]^2
        @inbounds for i in 1:Q
            eu  = ex[i] * ux[n, m] + ey[i] * uy[n, m]
            feq = weights[i] * rho[n, m] * (1. + 3. * eu + 9/2 * eu^2 - 3/2 * uu)
            fi  = f[n, m, i]
            Si  = 0.0
            f′[n, m, i] = fi + ω̃ * (feq - fi) + Δt̃ * Si
        end
    end

    return nothing
end

using Einsum
function moments()
    @einsum rho[n, m] = f[n, m, i]
    
    return nothing
end

function gauss(dx)
    return exp(-0.5*dx^2/σ0^2)
end

function init()
    @inbounds for n in 1:H̃, m in 1:W̃
        rho[n, m] = gauss(n - 200) * gauss(m - 200)
        uu = ũ^2 + ũ^2
        @inbounds for i in 1:Q
            eu  = ex[i] * ũ + ey[i] * ũ
            f[n, m, i] = weights[i] * ρ̃ * rho[n, m] * (1. + 3. * eu + 9/2 * eu^2 - 3/2 * uu)
        end
    end

    return nothing
end

function step()
    moments()
    collide()
    stream()
    
    return nothing
end

function evaluate()
    @show findmax(rho)
    @show C0 = σ0^2 / (σ0^2 + 2.0D * 200)
end