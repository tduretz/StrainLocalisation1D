using Plots, LinearAlgebra, Statistics, StaticArrays, BackwardDifferenceFormula

@views function ViscoElastic_CrankNicolson(;θM=1.0, θT=1.0, nres=1, adapt_dt=false, viz=false)

    CFL_M = 0.1
    CFL_T = 0.5
    rel   = 0.2
    bdf   = 1   # if 0 then NO BDF is used (CN or BE depending on θM, θT) 

    nresmax = 32
    spY     = 3600*24*365.25

    # For this setting BDF schemes predict early runaway (at 5<t<10) 
    G      = 120e9
    Δt_ref = 1.28*spY 
    t_end  = 800*spY
    Nt     = 6000*nresmax
    ε̇      = 4e-13
    H_R    = 63748.0
    ε̇0     = abs(ε̇)
    τ0     = 1e9
    Tini   = 973.15
    n      = 3.5
    C0     = 1.1e-16   #Pa^-ndis s^-1
    k      = 3.0
    ρ      = 3469.0
    cp     = 1000.0
    L      = 20e3
    σ      = 0.05 * L
    ωmax   = 10.0
    hL     = 0.02
    η0     = 0.5 * C0^(-1.0/n) * ε̇0^(1.0/n - 1.0) * exp(H_R/(n*Tini))

    # time stepping
    if adapt_dt # adaptive
        Δt_max = Δt_ref / 2.0
        Δt_min = Δt_ref / 64.0
        dT_ref = 10.0
        dτ_ref = 50e6
        Δt     = copy(Δt_ref)
    else
        Δt     = Δt_ref/nresmax
        Nt     = Int64(ceil(t_end/Δt))
    end
    println("Timestep: $(Δt/spY) yr")

    # 1D
    y       = (min = 0, max   = L)
    ncy     = 100
    Δy      = (y.max - y.min) / ncy
    yv      = LinRange(y.min, y.max, ncy+1) 
    yc      = LinRange(y.min+Δy/2, y.max-Δy/2, ncy) 
    yce     = LinRange(y.min-Δy/2, y.max+Δy/2, ncy+2) 
    Vx      = zeros(ncy+2)
    Rx      = zeros(ncy+2)
    Rx0     = zeros(ncy+2)
    ∂Vx∂τ   = zeros(ncy+2)
    τxy     = zeros(ncy+1)
    τxy0    = zeros(ncy+1)
    τxy00   = zeros(ncy+1)
    τxy000  = zeros(ncy+1)
    τxy0000 = zeros(ncy+1)
    f0_M    = zeros(ncy+1) 
    ε̇xy     = zeros(ncy+1)
    ε̇xy_ve  = zeros(ncy+1)
    η_ve    = zeros(ncy+1)
    a_ve    = zeros(ncy+1)
    T       = zeros(ncy+2)
    T0      = zeros(ncy+2)
    T00     = zeros(ncy+2)
    T000    = zeros(ncy+2)
    T0000   = zeros(ncy+2)
    Hs      = zeros(ncy+0)
    ε̇II_v   = zeros(ncy+1) 
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
    ω       = 1.0 .+ (ωmax-1.0) * exp.(-0.5*((yv .- L/2.0)./σ).^2)
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
    probes   = (vmax=zeros(Nt),)

    a, b, c, d, e = 0, 0, 0, 0, 0
    Δt0, Δt00 = Δt, Δt
    it = 0
    t  = 0.

    # scaling for residual
    Tscal = Tini
    Vscal = maximum(Vx) - minimum(Vx)
    
    # Time integration
    @time twall = @elapsed while t < t_end  #&& it< 2
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
        # @show Δt_min, Δt_T, Δt_τ, Δt * 1.25, Δt_max
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
        f0_T .= (qT[2:end]  .- qT[1:end-1])./Δy .- Hs
        
        if bdf==0
            ηe    = θM*G*Δt
        else
            ηe    = G/a
        end
        Tv   .= 0.5.*(T[2:end].+T[1:end-1])
        ηv   .= 0.5.*C.^(-1).*τII.^(1-n).*exp.(H_R./Tv)
        η_ve .= (1 ./ ηv .+ 1 ./ ηe).^(-1)

        # PT params 1
        λmax  = 4*maximum(η_ve)/Δy^2
        λmin  = 0.
        Δτx   = 2/sqrt(λmax) * CFL_M
        cx    = sqrt(λmin) * 0.5
        αx    = (2 - cx*Δτx) / (2 + cx*Δτx)
        βx    = 2Δτx^2 / (2 + cx*Δτx)

        λmax  = 4*k/Δy^2 + ρ*cp/Δt
        λmin  = 0.
        ΔτT   = 2/sqrt(λmax) *CFL_T
        cT    = sqrt(λmin) * 0.5
        αT    = (2 - cT*ΔτT) / (2 + cT*ΔτT)
        βT    = 2ΔτT^2 / (2 + cT*ΔτT)

        nRx0, nRT0 = 1.0, 1.0
        tol        = 1e-7

        for iter=1:100000

            Rx0         .= Rx
            RT0         .= RT

            Vx[1]        = 2*ε̇0*yv[1]   - Vx[2]
            Vx[end]      = 2*ε̇0*yv[end] - Vx[end-1]
            ε̇xy         .= (Vx[2:end] .- Vx[1:end-1])./Δy
            Tv          .= 0.5.*(T[2:end].+T[1:end-1])
            ηv_true     .= 0.5.*C.^(-1).*τII.^(1-n).*exp.(H_R./Tv)
            ηv          .= exp.( rel.*log.(ηv_true) .+ (1-rel).*log.(ηv))
            η_ve        .= (1 ./ ηv .+ 1 ./ ηe).^(-1)
            if bdf==0
                a_ve        .= (1-θM) .* ηv.*Δt./(ηv .+ ηe)
                ε̇xy_ve      .= ε̇xy .+ τxy0./(2 .* ηe)
                τxy         .= 2 .*η_ve.* (   ε̇xy  .- (b.*τxy0 .+ c.*τxy00 .+ d.*τxy000 .+ e.*τxy0000)./2 ./G)
            else
                τxy         .= 2 .*η_ve.* (   ε̇xy  .- (b.*τxy0 .+ c.*τxy00 .+ d.*τxy000 .+ e.*τxy0000)./2 ./G)
            end
            τII         .= abs.(τxy)
            Rx[2:end-1] .= (τxy[2:end] .- τxy[1:end-1])./Δy
            ∂Vx∂τ       .= Rx .+ αx.*∂Vx∂τ
            Vx         .+= βx.*∂Vx∂τ
            
            T[1]         = T[2]
            T[end]       = T[end-1]
            qT          .= -k.*(T[2:end].-T[1:end-1])./Δy
            Hs          .= 0.5.*(τII[1:end-1].^2 ./ ηv[1:end-1] .+ τII[2:end].^2 ./ ηv[2:end])
            if bdf==0
                RT[2:end-1] .= .-(ρ.*cp.*(T[2:end-1].-T0[2:end-1])./Δt .+ θT.*((qT[2:end] .- qT[1:end-1])./Δy .- Hs) .+ (1-θT).*f0_T)
            else
                RT[2:end-1] .= .-(ρ.*cp.*(a.*T[2:end-1].+b.*T0[2:end-1].+c.*T00[2:end-1].+d.*T000[2:end-1].+e.*T0000[2:end-1]) .+ ((qT[2:end] .- qT[1:end-1])./Δy .- Hs))
            end
            ∂T∂τ        .= RT .+ αT.*∂T∂τ
            T          .+= βT.*∂T∂τ

            if iter==1 || mod(iter, 100)==0
                errV = Rx ./ Vscal
                errT = RT ./ Tscal
                #@show it, iter, norm(errV), norm(errT)
                if iter==1 nRx0, nRT0 = norm(errV)/ncy, norm(errT)/ncy end
                isnan(norm(errV)) && error("NaN")
                if (norm(errT)<tol || norm(errT)/nRT0<tol)  && (norm(errV)<1e-6 || norm(errV)/nRx0<1e-6)
                #if (norm(errT)<tol)  && (norm(errV)<1e-6)
                    println("$(iter), $(iter/ncy), $(Δt)")
                    break
                end
                # PT params 2
                λmax  = 4*maximum(η_ve)/Δy^2
                Δτx   = 2/sqrt(λmax) * CFL_M
                dum  .= βx.*∂Vx∂τ.*(Rx.-Rx0)
                λmin  = abs.(sum(dum))
                dum  .= βx.^2 .*∂Vx∂τ.^2
                λmin /= abs.(sum(dum))
                cx    = sqrt(λmin) * 0.5
                αx    = (2 - cx*Δτx) / (2 + cx*Δτx)
                βx    = 2Δτx^2 / (2 + cx*Δτx)
                dum  .= βT.*∂T∂τ.*(RT.-RT0)
                λmin  = abs.(sum(dum))
                dum  .= βT.^2 .*∂T∂τ.^2
                λmin /= abs.(sum(dum))
                cT    = sqrt(λmin) *2
                αT    = (2 - cT*ΔτT) / (2 + cT*ΔτT)
                βT    = 2ΔτT^2 / (2 + cT*ΔτT) 
            end
        end
        Δtvec[it]       = Δt
        tvec[it]        = t
        τvec_ana[it]    = 2η0*ε̇*(1 .- exp.(-G/η0.*t))
        τvec[it]        = mean(τII)
        probes.vmax[it] = maximum(Vx)/abs(ε̇0*yv[end])
    
        if mod(it, 10)==0 && viz
            # Visualise
            p1 = plot(legend=:bottomright)
            p1 = plot!(tvec[1:it] ./ spY, τvec[1:it] ./ 1e6)
            #p1 = plot!(tvec[1:it] ./ spY, τvec_ana[1:it] ./ 1e6, label="ref", ls=:dash )
            p2 = plot(T[2:end-1], yc ./ 1e3, label="T")
            # p3 = plot(log10.(ηv), yv, label="η")
            p3 = scatter(tvec[1:it] ./ spY, probes.vmax[1:it], label="Vmax")
            p4 = plot(Vx[2:end-1], yc ./ 1e3, label="Vx")
            display(plot(p1,p2,p3,p4))
        end
    end 
    vmax,imax = findmax(probes.vmax) 
    @info vmax
    @info tvec[imax]
    @info twall
    return vmax, tvec[imax], minimum(Δtvec), maximum(Δtvec), twall 
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


ViscoElastic_CrankNicolson(;θM=1.0, θT=1.0, nres=1, adapt_dt=true, viz=true)

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