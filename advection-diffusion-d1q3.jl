# based on a C code (GaussianHill) from Alexandr Kuzmin
# https://github.com/shurikkuzmin/LatticeBoltzmannMethod

const NX = 512
const NT = 200
const σ0 = 10.0 
const Q = 3
const f  = zeros(NX, Q) # population distribution
const f′ = zeros(NX, Q) # post-collision distributio
const rho = zeros(NX)   # 0-th moment

const ũx = 0.1 # 0.1 # 0.0

const weights = [2/3 1/6 1/6];
const compliment = [1 3 2]
const cx = [0 1 -1]

const τg = 0.513 # 0.513 # 5.0
const Λ = 1/6 # 1/12  # 1/6
const D = (1/3) * (τg - 1/2)

rhos = zeros(NX, NT)
function saverho(t)
    rhos[:, t] .= rho
end

function writerho(fname)
    open(fname, "w") do f
        for n=1:NX
            print(f, rho[n], " ")
        end
        println(f)
    end
    
    return nothing
end

function gauss(dx)
    return exp(-0.5*dx^2/σ0^2)
end

function init()
    for n=1:NX
        rho[n] = gauss(n - 200)
    end
    
    for n = 1:NX, i=1:Q
        f[n, i] = weights[i]*rho[n]*(1.0 + 3.0 * cx[i]*ũx + 4.5 * (cx[i]^2 - 1/3)*ũx^2)
    end
    
    return nothing
end

const g⁺   = zeros(Q)
const g⁻   = zeros(Q)
const geq  = zeros(Q)
const geq⁺ = zeros(Q)
const geq⁻ = zeros(Q)  
function collide()
    for n=1:NX
        
        # for i=1:Q
        #     g⁺[i] = 0.5*(f[n, i] + f[n, compliment[i]])
        #     g⁻[i] = 0.5*(f[n, i] - f[n, compliment[i]])
        # end
        
        for i=1:Q
            geq[i] = weights[i]*rho[n]*(1.0+3.0*cx[i]*ũx+4.5*(cx[i]^2 - 1/3)*ũx^2);			
        end
        
        # for i=1:Q
        #     geq⁺[i] = 0.5*(geq[i] + geq[compliment[i]]);
        #     geq⁻[i] = 0.5*(geq[i] - geq[compliment[i]]);
        # end
        
        # ω_plus_rho  = 1.0/τg
        ω = 1.0 / τg
        for i=1:Q
            f′[n, i] = f[n, i] - ω*(f[n, i]-geq[i])
        end
    end
    
    return nothing
end

fi = zeros(NX)
function stream()
    @inbounds for i in 1:Q
        copyto!(fi,view(f′, :, i))
        circshift!(view(f,  :, i), fi, (cx[i]))
    end
    return nothing
end

using Einsum
using BenchmarkTools
function moments()
    @einsum rho[n] = f[n, i]

    return nothing
end

function step()
    moments()
    collide()
    stream()
    
    return nothing
end

function main()
    init()
    for t=1:NT
        step()
        saverho(t)
    end
    println("done.")
    return 0;
end

using UnicodePlots
plt = lineplot([gauss(n - 200) for n in 1:NX], name="analytical")
init()
plt = lineplot!(plt, rho, name="initial")
main()
plt = lineplot!(plt, rho, name="final")
show(plt)

# benchmark
#t = @belapsed step()
#MLUps = NX / t * 1e-6
#@show MLUps
