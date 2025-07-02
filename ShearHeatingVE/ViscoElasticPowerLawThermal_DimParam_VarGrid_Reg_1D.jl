using Plots, LinearAlgebra, Statistics, StaticArrays, BackwardDifferenceFormula

@views function ViscoElastic_CrankNicolson(;θM=1.0, θT=1.0, nresy=1, nrest=1, adapt_dy=false, adapt_dt=false, viz=false, noisy=false)

    CFL_M    = 0.5
    CFL_T    = 0.5
    c_fact   = 1.8 # Damping
    rel      = 0.5
    bdf      = 4   # if 0 then NO BDF is used (CN or BE depending on θM, θT) 
    debug    = false
    nitermax = 3e5
    nresmax  = 128
    spY      = 3600*24*365.25

    # Scales
    sc = (L=1e3, t=1e7, σ=1e6, T=1000)
    m  = sc.σ * sc.L * sc.t^2.0
    J  = m * sc.L^2.0 / sc.t^2.0
    W  = J/sc.t

    # For this setting BDF schemes predict early runaway (at 5<t<10) 
    n      = 3.5
    L      = 20e3     / sc.L
    G      = 120e9    / sc.σ
    Δt_ref = 1.28*spY / sc.t
    t_end  = 600*spY  / sc.t
    ε̇      = 4e-13    * sc.t
    H_R    = 63748.0  / sc.T
    η_reg  = 1e16     / sc.σ / sc.t
    τ0     = 1e9      / sc.σ
    Tini   = 973.15   / sc.T
    C0     = 1.1e-16  / (sc.σ^(-n)/sc.t) #Pa^-ndis s^-1
    k      = 3.0      / (W/sc.L/sc.T) 
    ρ      = 3469.0   / (m/sc.L^3)
    cp     = 1000.0   / (J/m/sc.T)
    Nt     = 6000*nresmax
    ε̇0     = abs(ε̇)
   
    σ      = 0.05 * L
    ωmax   = 10.0
    hL     = 0.02
    η0     = 0.5 * C0^(-1.0/n) * ε̇0^(1.0/n - 1.0) * exp(H_R/(n*Tini))
    y      = (min = -L/2, max   = L/2)
    ncy    = 100*nresy

    # time stepping
    if adapt_dt # adaptive
        Δt_max = Δt_ref / 2.0
        Δt_min = Δt_ref / 128.0
        dT_ref = 10.0   / sc.T
        dτ_ref = 50e6   / sc.σ
        Δt     = copy(Δt_ref)
    else
        Δt     = Δt_ref/nrest
        Nt     = Int64(ceil(t_end/Δt))
    end
    println("Timestep: $(Δt/spY) yr")

    # 1D
    if adapt_dy
        y0     = (y.max + y.min)/2
        σy     = 0.5*(nresy)^(1/3)                 # the larger, the finer the grid
        yv0    = LinRange(y.min, y.max, ncy+1) 
        ymin1  = (sinh.( σy.*(yv0[1]  .-y0) ))
        ymax1  = (sinh.( σy.*(yv0[end].-y0) ))
        sy     = (yv0[end]-yv0[1])/(ymax1-ymin1)
        yv     = (sinh.( σy.*(yv0.-y0) )) .* sy  .+ y0  
        Δy_c   = diff(yv)
        yc     = 0.5*(yv[1:end-1] .+ yv[2:end])
        yce    = zeros(ncy+2); yce[2:end-1] .= yc; yce[1] = yv[1]-Δy_c[1]/2; yce[end] = yv[end]+Δy_c[end]/2
        Δy_v   = diff(yce)
        ncy_eq = Int64(ceil((y.max - y.min) / minimum(Δy_c)))
        # display(plot(yc0, log10.(Δy_c./Δy_c[1])))
        # display(scatter(yv[2:end-1], max.(Δy_c[1:end-1], Δy_c[2:end])./min.(Δy_c[1:end-1], Δy_c[2:end]), title="ratio"))
        # error()
    else
        Δy      = (y.max - y.min) / ncy
        yv      = LinRange(y.min, y.max, ncy+1) 
        yc      = LinRange(y.min+Δy/2, y.max-Δy/2, ncy) 
        yce     = LinRange(y.min-Δy/2, y.max+Δy/2, ncy+2) 
        Δy_c    = Δy*ones(ncy)
        Δy_v    = Δy*ones(ncy+1)
        ncy_eq  = ncy
    end

    Vx      = zeros(ncy+2)
    Rx      = zeros(ncy+2)
    Rx0     = zeros(ncy+2)
    ∂Vx∂τ   = zeros(ncy+2)
    τxy     = zeros(ncy+1)
    τxy0    = zeros(ncy+1)
    τxy00   = zeros(ncy+1)
    τxy000  = zeros(ncy+1)
    τxy0000 = zeros(ncy+1)
    τxy_reg = zeros(ncy+1)
    τII_dis = zeros(ncy+1)
    f0_M    = zeros(ncy+1) 
    ε̇xy     = zeros(ncy+1)
    ε̇xy_ve  = zeros(ncy+1)
    ε̇xy_v   = zeros(ncy+1)
    η_ve    = zeros(ncy+1)
    a_ve    = zeros(ncy+1)
    T       = zeros(ncy+2)
    T0      = zeros(ncy+2)
    T00     = zeros(ncy+2)
    T000    = zeros(ncy+2)
    T0000   = zeros(ncy+2)
    Hs      = zeros(ncy+0)
    τII     =1*ones(ncy+1) 
    ηv      =  ones(ncy+1) 
    ηv_true =  ones(ncy+1) 
    RT      = zeros(ncy+2)
    RT0     = zeros(ncy+2)
    ∂T∂τ    = zeros(ncy+2)
    qT      = zeros(ncy+1)
    Tv      = zeros(ncy+1)
    f0_T    = zeros(ncy+0)
    dum     = zeros(ncy+2)
    dumv    = zeros(ncy+1)
 
    T      .= Tini #.+ 50*exp.(-(yce.-mean(yce)).^2/σ.^2)

    # set rheology anomaly
    C       = ones(ncy+1) .* C0
    σ       = hL*L/(2*sqrt(2*log(2)))
    ω       = 1.0 .+ (ωmax-1.0) * exp.(-0.5*((yv)./σ).^2)
    C     .*= ω

    # set initial velocity
    Vx     .= ε̇  .* yce

    # set initial stress
    τxy    .= τ0
    τxy0   .= τ0
    τII    .= τ0

    # Analytics
    Δtvec    = zeros(Nt)
    tvec     = zeros(Nt)
    τvec_ana = zeros(Nt)
    τvec     = zeros(Nt)
    err      = (T = zeros(Nt), V = zeros(Nt))
    probes   = (vmax=zeros(Nt),)

    a, b, c, d, e = 0.0, 0.0, 0.0, 0.0, 0.0
    Δt0, Δt00 = Δt, Δt
    it = 0
    t  = 0.

    # Max number of time steps
    nitmax = debug ? 2 : 1e4

    # Time integration
    @time twall = @elapsed while t < t_end && it<nitmax
        it += 1

        # adapt time step
        if adapt_dt && it > 1
            dum .= abs.(T .- T0)
            max_ΔT = maximum(dum)
            dumv .= abs.(τxy .- τxy0)
            max_Δτ = maximum(dumv)
            Δt_T   = Δt * dT_ref / max_ΔT
            Δt_τ   = Δt * dτ_ref / max_Δτ
            Δt     = max(Δt_min, min(Δt_T, Δt_τ, Δt * 1.25, Δt_max))
        end
        #println("Timestep: $(Δt/spY) yr")
        t += Δt

        τxy0000 .= τxy000
        τxy000  .= τxy00
        τxy00   .= τxy0
        τxy0    .= τxy

        T0000   .= T000
        T000    .= T00
        T00     .= T0
        T0      .= T

        Δt000 = Δt00
        Δt00  = Δt0
        Δt0   = Δt

        if bdf==1 || it<1+1 || bdf==0
            a = 1/Δt
            b =-1/Δt
            c = 0.
            d = 0.
            e = 0.
        end

        if bdf==2 && it>1
            coeffs  = bdf_coefficients(@SVector([ -Δt-Δt0, -Δt, 0]))
            a, b, c = coeffs[1], coeffs[2], coeffs[3]
            d = 0.
            e = 0.
        end

        if bdf==3 && it>2
            coeffs  = bdf_coefficients(@SVector([-Δt-Δt0-Δt00, -Δt-Δt0, -Δt, 0]))
            a, b, c, d = coeffs[1], coeffs[2], coeffs[3], coeffs[4]
            e = 0.
        end

        if bdf==4 && it>3
            coeffs  = bdf_coefficients(@SVector([-Δt-Δt0-Δt00-Δt000, -Δt-Δt0-Δt00, -Δt-Δt0, -Δt, 0]))
            a, b, c, d, e = coeffs[1], coeffs[2], coeffs[3], coeffs[4], coeffs[5]
        end

        f0_M .= (2*G.*ε̇xy .- G./ηv.*τxy)
        f0_T .= (qT[2:end]  .- qT[1:end-1])./Δy_c .- Hs
        
        if bdf==0
            ηe    = θM*G*Δt
        else
            ηe    = G/a
        end
        Tv   .= 0.5.*(T[2:end] .+ T[1:end-1])
        ηv   .= 0.5.*C.^(-1).*τII.^(1-n).*exp.(H_R./Tv) .+ η_reg
        η_ve .= (1 ./ ηv .+ 1 ./ ηe).^(-1)

        # PT params 1
        λmax  = 4*maximum(η_ve)/minimum(Δy_c)^2
        λmin  = 0.
        Δτx   = 2/sqrt(λmax) * CFL_M
        cx    = sqrt(λmin) * 0.5
        αx    = (2 - cx*Δτx) / (2 + cx*Δτx)
        βx    = 2Δτx^2 / (2 + cx*Δτx)

        λmax  = 4*k/minimum(Δy_c)^2 + ρ*cp/Δt
        λmin  = 0.
        ΔτT   = 2/sqrt(λmax) *CFL_T
        cT    = sqrt(λmin) * 0.5
        αT    = (2 - cT*ΔτT) / (2 + cT*ΔτT)
        βT    = 2ΔτT^2 / (2 + cT*ΔτT)

        nRx0, nRT0 = 1.0, 1.0
        tol        = 1e-5

        for iter=1:nitermax

            Rx0         .= Rx
            RT0         .= RT

            Vx[1]        = 2*ε̇0*yv[1]   - Vx[2]
            Vx[end]      = 2*ε̇0*yv[end] - Vx[end-1]
            ε̇xy         .= (Vx[2:end] .- Vx[1:end-1])./Δy_v
            ε̇xy_v       .= ε̇xy .- (τxy .- τxy0) ./ (2 * ηe)
            τxy_reg     .= 2 .* η_reg .* ε̇xy_v
            τII_dis     .= abs.(τxy .- τxy_reg)
            Tv          .= 0.5.*(T[2:end] .+ T[1:end-1])
            ηv_true     .= 0.5.*C.^(-1).*τII_dis.^(1-n).*exp.(H_R./Tv) .+ η_reg
            ηv          .= exp.( rel.*log.(ηv_true) .+ (1-rel).*log.(ηv))
            # ηv_dis      .= exp.( rel.*log.(ηv_true) .+ (1-rel).*log.(ηv_dis))
            # ηv          .= ηv_dis .+ η_reg
            η_ve        .= (1 ./ ηv .+ 1 ./ ηe).^(-1)
            # # # Vx[1]        = 2*ε̇0*yv[1]   - Vx[2]
            # # # Vx[end]      = 2*ε̇0*yv[end] - Vx[end-1]
            # # # ε̇xy         .= (Vx[2:end] .- Vx[1:end-1])./Δy_v
            # # # Tv          .= 0.5.*(T[2:end] .+ T[1:end-1])
            # # # ηv_true     .= 0.5.*C.^(-1).*τII.^(1-n).*exp.(H_R./Tv)
            # # # ηv          .= exp.( rel.*log.(ηv_true) .+ (1-rel).*log.(ηv))
            # # # η_ve        .= (1 ./ ηv .+ 1 ./ ηe).^(-1)
            if bdf==0
                a_ve        .= (1-θM) .* ηv.*Δt./(ηv .+ ηe)
                ε̇xy_ve      .= ε̇xy .+ τxy0./(2 .* ηe)
                τxy         .= 2 .*η_ve.* ε̇xy_ve
            else
                τxy         .= 2 .* η_ve .* (   ε̇xy  .- (b.*τxy0 .+ c.*τxy00 .+ d.*τxy000 .+ e.*τxy0000)./2 ./G) 
            end
            τII         .= abs.(τxy)
            Rx[2:end-1] .= (τxy[2:end] .- τxy[1:end-1])./Δy_c
            ∂Vx∂τ       .= Rx .+ αx.*∂Vx∂τ
            Vx         .+= βx.*∂Vx∂τ
            
            T[1]         = T[2]
            T[end]       = T[end-1]
            qT          .= -k.*(T[2:end].-T[1:end-1])./Δy_v
            Hs          .= 0.5.*(τII[1:end-1].^2 ./ ηv[1:end-1] .+ τII[2:end].^2 ./ ηv[2:end])
            if bdf==0
                RT[2:end-1] .= .-(ρ.*cp.*(T[2:end-1].-T0[2:end-1])./Δt .+ θT.*((qT[2:end] .- qT[1:end-1])./Δy_c .- Hs) .+ (1-θT).*f0_T)
            else
                RT[2:end-1] .= .-(ρ.*cp.*(a.*T[2:end-1].+b.*T0[2:end-1].+c.*T00[2:end-1].+d.*T000[2:end-1].+e.*T0000[2:end-1]) .+ ((qT[2:end] .- qT[1:end-1])./Δy_c .- Hs))
            end
            ∂T∂τ        .= RT .+ αT.*∂T∂τ
            T          .+= βT.*∂T∂τ

            if iter==1 || mod(iter, 100)==0
                if iter==1 nRx0, nRT0 = norm(Rx)/ncy, norm(RT)/ncy end
                errV = max(norm(Rx)/ncy_eq, norm(Rx)/ncy_eq/nRx0)
                errT = max(norm(RT)/ncy_eq, norm(RT)/ncy_eq/nRT0)
                # @show it, iter, errV, errT
                isnan(errV) && error("NaN @ it = $(iter)")
                if (errT<tol)  && (errV<tol)
                    err.T[it] = errV
                    err.V[it] = errT
                    noisy ? println("$(iter), $(iter/ncy), $(Δt)") : nothing
                    break
                end
                # PT params 2
                λmax  = 4*maximum(η_ve)/minimum(Δy_c)^2
                Δτx   = 2/sqrt(λmax) * CFL_M
                dum  .= βx.*∂Vx∂τ.*(Rx.-Rx0)
                λmin  = abs.(sum(dum))
                dum  .= βx.^2 .*∂Vx∂τ.^2
                λmin /= abs.(sum(dum))
                cx    = sqrt(λmin) * c_fact
                αx    = (2 - cx*Δτx) / (2 + cx*Δτx)
                βx    = 2Δτx^2 / (2 + cx*Δτx)
                dum  .= βT.*∂T∂τ.*(RT.-RT0)
                λmin  = abs.(sum(dum))
                dum  .= βT.^2 .*∂T∂τ.^2
                λmin /= abs.(sum(dum))
                cT    = sqrt(λmin) * c_fact
                αT    = (2 - cT*ΔτT) / (2 + cT*ΔτT)
                βT    = 2ΔτT^2 / (2 + cT*ΔτT) 
            end
        end
        Δtvec[it]       = Δt
        tvec[it]        = t
        τvec_ana[it]    = 2η0*ε̇*(1 .- exp.(-G/η0.*t))
        τvec[it]        = mean(τII)
        probes.vmax[it] = maximum(Vx)/abs(ε̇0*yv[end])
    
        if (mod(it, 200)==0 || it==Nt || it==nitmax) && viz
            # Visualise
            p1 = plot(legend=:bottomright)
            p1 = plot!(tvec[1:it] .* sc.t ./ spY, τvec[1:it] .* sc.σ ./ 1e6)
            #p1 = plot!(tvec[1:it] ./ spY, τvec_ana[1:it] ./ 1e6, label="ref", ls=:dash )
            p2 = plot(T[2:end-1] .* sc.T, yc .* sc.L ./ 1e3, label="T")
            # p3 = plot(log10.(ηv), yv, label="η")
            p3 = scatter(tvec[1:it] .* sc.t ./ spY, probes.vmax[1:it], label="Vmax")
            p4 = plot(Vx[2:end-1] .* sc.L/sc.t, yc .* sc.L ./ 1e3, label="Vx")
            p5 = plot(legend=:bottomright)
            p5 = scatter!(tvec[1:it] .* sc.t ./ spY, log10.(err.T[1:it]), label="err T")
            p5 = scatter!(tvec[1:it] .* sc.t ./ spY, log10.(err.V[1:it]), label="err V")
            p6 = plot(yc, log10.(Δy_c./Δy_c[1]), label = "Δy")
            display(plot(p1,p2,p3,p4,p5,p6, layout=(3, 2)))
        end
    end 
    vmax,imax = findmax(probes.vmax) 
    @info vmax
    @info tvec[imax]
    @info twall
    dt_vals = (minimum((Δtvec[1:it])), maximum((Δtvec[1:it])))
    dx_vals = (minimum((Δy_c)), maximum((Δy_c)))
    return vmax, tvec[imax], twall, dt_vals, dx_vals 
end
 
# let
#     # n = [2 4 8 16 32]
#     n = [64]
#     v     = zeros(length(n))
#     t     = zeros(length(n))
#     dtmin = zeros(length(n))
#     dtmax = zeros(length(n))
   
#     for ires in eachindex(n)
#         @show n[ires]
#         v[ires], t[ires], dtmin[ires], dtmax[ires] = ViscoElastic_CrankNicolson(;θM=1/2, θT=1/2, nres=n[ires], adapt_dt=false, viz=false)
#     end
#     @show v
#     @show t
#     @show dtmin
#     @show dtmax
# end


ViscoElastic_CrankNicolson(;θM=1.0, θT=1.0, nresy=1, nrest=1, adapt_dy=true, adapt_dt=true, viz=true, noisy=false)

# v = [119.61046691382055, 76.78011622622239, 76.14098484053791, 76.15711043461168, 76.17197880423882]
# t = [14.944999999999169, 14.95875000000304, 14.959374999995456, 14.959374999990926, 14.959218749988661]
# dtmin = [0.0025, 0.00125, 0.000625, 0.0003125, 0.00015625]
# dtmax = [0.0025, 0.00125, 0.000625, 0.0003125, 0.00015625]


# # tougher (rho=5)
# v = [204.52516689855722, 155.96885327473433, 140.9638203888602, 136.36902800372212, 134.11591522309152]
# t = [20.600000000001312, 20.616249999999372, 20.623124999990306, 20.626249999985774, 20.627968749983506]
# dtmin = [0.0025, 0.00125, 0.000625, 0.0003125, 0.00015625]
# dtmax = [0.0025, 0.00125, 0.000625, 0.0003125, 0.00015625]
# Adapt
# [ Info: 134.00921452749344
# [ Info: 20.623782224164753
# (134.00921452749344, 20.623782224164753, 0.00015625, 0.005)


# Tougher 2: ρ      = 90.0
# n[ires] = 2
# [ Info: 353.26600619557814
# [ Info: 36.745000000002825
# [ Info: 76.563806334
# n[ires] = 4
# [ Info: 335.32053921912154
# [ Info: 36.7549999999847
# [ Info: 133.376973375
# n[ires] = 8
# [ Info: 328.8749764527673
# [ Info: 36.759374999975634
# [ Info: 238.71558325
# n[ires] = 16
# [ Info: 325.8475035241317
# [ Info: 36.7618749999711
# [ Info: 371.055890625
# n[ires] = 32
# [ Info: 324.34702621873
# [ Info: 36.76312500007713
# [ Info: 278.030052541
# n[ires] = 64
# 306.746844 seconds
# [ Info: 323.593683036392
# [ Info: 36.763750000439806
# [ Info: 306.746843459
# Adapt 32
# [ Info: 324.1710873602131
# [ Info: 36.75833708002466
# [ Info: 339.209569459
# Adapt 64
# [ Info: 323.4110038456593
# [ Info: 36.758670603675
# [ Info: 82.90927225
# (324.1710873602131, 36.75833708002466, 0.00015625, 0.005, 339.209569459)