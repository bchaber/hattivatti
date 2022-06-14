# based on a C code (GaussianHill) from Alexandr Kuzmin
# https://github.com/shurikkuzmin/LatticeBoltzmannMethod

const NX = 80
const NT = 20
const xInit = 40 + 1
const σ = 8


const g = zeros(NX, 3)
const g2= zeros(NX, 3)
const phase = zeros(NX)

const ux = 0.5

const weights = [2/3 1/6 1/6];
const compliment =[1 3 2]
const cx = [0 1 -1]

const Λ = 1.0/6.0
const ω =0.3


function writephase(fname)
    open(fname, "w") do f
        for iX=1:NX
            print(f, phase[iX], " ")
        end
        println(f)
    end
    
    return nothing
end

function gauss(dx)
    return exp(-0.5*dx^2/σ^2)
end

function init()
    for iX=1:NX
        if abs(iX - xInit) <= 3σ
            phase[iX] = gauss(iX-xInit)
        else
            phase[iX] = 0.0
        end
    end
    
    for iX = 1:NX, k=1:3
        g[iX, k] = weights[k]*phase[iX]*(1.0 + 3.0 * cx[k]*ux + 4.5 * (cx[k]*cx[k]-1.0/3.0)*ux*ux)
    end
    
    return nothing
end

function collide()
    for iX=1:NX
        phase[iX]=0.0;
        for iPop=1:3
            phase[iX]+=g[iX, iPop]
        end
    end
    
    for iX=1:NX
        g_plus    = zeros(3)
        g_minus   = zeros(3)
        geq       = zeros(3)
        geq_plus  = zeros(3)
        geq_minus = zeros(3)  
        
        for k=1:3
            g_plus[k]  = 0.5*(g[iX, k]+g[iX, compliment[k]])
            g_minus[k] = 0.5*(g[iX, k]-g[iX, compliment[k]])
        end
        
        for k=1:3
            geq[k]=weights[k]*phase[iX]*(1.0+3.0*cx[k]*ux+4.5*(cx[k]*cx[k]-1.0/3.0)*ux*ux);			
        end
        
        
        for k=1:3
            geq_plus[k]  = 0.5*(geq[k]+geq[compliment[k]]);
            geq_minus[k] = 0.5*(geq[k]-geq[compliment[k]]);
        end
        
        ω_minus_phase = 1.0/(Λ/(1.0/ω-0.5)+0.5);
        ω_plus_phase  = 1.0/(Λ/(1.0/ω-0.5)+0.5);
        
        for k=1:3
            g2[iX, k] = g[iX, k]-ω_plus_phase*(g_plus[k]-geq_plus[k])-ω_minus_phase*(g_minus[k]-geq_minus[k]);
        end
    end
    
    return nothing
end

function stream()
    for iX=1:NX
        for iPop=1:3
            iX2 = (iX - cx[iPop] + NX) % NX + 1;
            g[iX, iPop] = g2[iX2, iPop];
        end
    end
end

function main()
    for t=0:NT
        collide()
        stream()
        if t % 10 == 0
            writephase("phase$t.dat")
        end
    end
    println("done.")
    return 0;
end

using UnicodePlots
plt = lineplot([gauss(dx) for dx = -30:+50], name="analytical")
init()
plt = lineplot!(plt, phase, name="initial")
main()
plt = lineplot!(plt, phase, name="final")
show(plt)