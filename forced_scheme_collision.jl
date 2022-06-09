
function collide2019() # [Khazaeli et al. 2019]
    @inbounds for n in 1:H̃, m in 1:W̃
        Δu = Δt̃ * F̃ / rho[n, m]
        uu = ux[n, m]^2 + uy[n, m]^2
        vv = uu + 2.0 * uy[n, m] * 0.5Δu + 0.25Δu^2
        @inbounds for i in 1:Q
            ev   = ex[i] * ux[n,m] + ey[i] * uy[n,m] + ey[i] * 0.5Δu
            ueq  = uy[n,m] + Δu
            feq  = weights[i] * rho[n, m] * (1. + 3. * ev + 9/2 * ev^2 - 3/2 * vv)
            fi   = f[n, m, i]
            Si   = weights[i] * (1.0 - 1.0 / 2.0τ̃) *
                    (3. * (ey[i] - ueq) + 9/2 * ey[i]^2 * ueq) * F̃
            f[n, m, i] +=  ω̃ * (feq - fi)
            f[n, m, i] += Δt̃ * Si
        end
    end

    return nothing
end
function collide1993() # [Shan and Chen 1993]
    @inbounds for n in 1:H̃, m in 1:W̃
        Δu = τ̃ * Δt̃ * F̃ / rho[n, m]
        uu = ux[n, m]^2 + uy[n, m]^2
        vv = uu + 2.0 * uy[n, m] * Δu + Δu^2
        @inbounds for i in 1:Q
            ev   = ex[i] * ux[n,m] + ey[i] * uy[n,m] + ey[i] * Δu
            feq  = weights[i] * rho[n, m] * (1. + 3. * ev + 9/2 * ev^2 - 3/2 * vv)
            fi   = f[n, m, i]
            Si  = 0.0
            f[n, m, i] +=  ω̃ * (feq  - fi)
            f[n, m, i] += Δt̃ * Si
        end
    end

    return nothing
end