using Plots, JLD2, Printf, DelimitedFiles

nc      = 100
dtref   = 0.01
P       = readdlm("BFD4/Param.csv")
numRuns = size(P, 1)


runList = []
nanList = []

for i = 1 : numRuns; push!(runList, @sprintf("BFD4/BFD4_%02d.jld2", i)); end

spY   = 3600*24*365.25
sc    = (L=1e3, t=1e7, σ=1e6, T=1000)
L     = 20e3
Tmax  = zeros(numRuns)
vmax  = zeros(numRuns)
SZW   = zeros(numRuns)
Cells = zeros(numRuns)
Δt    = zeros(numRuns)

for i = eachindex(runList)
    name = runList[i]
    if P[i,1] > 0
        Cells[i] = Int64(P[i,1] * nc)
    else
        Cells[i] = 6400
    end
    if P[i,2] > 0
        Δt[i]    = P[i,2] * dtref
    else
        Δt[i]    = dtref * 1.5
    end

    if isfile(name)
        d = load(name)
        tend = maximum(d["t"]) * sc.t / spY
        if tend < 600
            println("$(name) did not finish. ($(round(tend)) yr)")
            Tmax[i], vmax[i], SZW[i], Cells[i], Δt[i] = NaN, NaN, NaN, NaN, NaN
            push!(nanList, i)
        else
            ind     = findall(d["t"] .> 0)
            Tmax[i] = maximum(d["evo"].Tmax[ind]) * sc.T
            vmax[i] = maximum(d["evo"].vmax[ind])
            if i < 43
                SZW[i]  = min(minimum(d["evo"].szw[ind]) * L / (P[i,1] * nc), L)
            else
                nc_sz   = minimum(d["evo"].szw[ind])
                Δy      = d["Δy"]
                ind     = Int64(floor((ncy+1-nc_sz)/2))
                SZW[i]  = sum(Δy[ind+1:ind+nc_sz])
            end
        end
    else
        println("$(name) does not exist.")
        Tmax[i], vmax[i], SZW[i], Cells[i], Δt[i] = NaN, NaN, NaN, NaN, NaN
        push!(nanList, i)
    end
end

ti_Cells   = [100, 200, 400, 800, 1600, 3200, 6400]
ti_Δt      = [1.5, 3.0, 6.0, 12.0, 24.0, 48.0, 96.0] .* dtref
la_Cells   = ["100", "200", "400", "800", "1600", "3200", "Refined"]
la_Δt      = ["Adaptive", "0.03", "0.06", "0.12", "0.24", "0.48", "0.96"]

p1 = scatter(Cells, Δt, marker_z = Tmax, colorbar_title="Tmax [K]",                             label=:none, ms=20, title="BDF 4", xlabel="# of cells", ylabel="Δt [yr]", xaxis=:log, yaxis=:log, xticks=(ti_Cells, la_Cells), yticks=(ti_Δt, la_Δt), frame=:box)
p2 = scatter(Cells, Δt, marker_z = vmax, colorbar_title="vmax []",                              label=:none, ms=20, title="BDF 4", xlabel="# of cells", ylabel="Δt [yr]", xaxis=:log, yaxis=:log, xticks=(ti_Cells, la_Cells), yticks=(ti_Δt, la_Δt), frame=:box)
p3 = scatter(Cells, Δt, marker_z = log10.(SZW), colorbar_title="log10(SZW [m])", clim=(1.7, 3), label=:none, ms=20, title="BDF 4", xlabel="# of cells", ylabel="Δt [yr]", xaxis=:log, yaxis=:log, xticks=(ti_Cells, la_Cells), yticks=(ti_Δt, la_Δt), frame=:box)

display(p1)
display(p2)
display(p3)

savefig(p1, "BDF4_T.png")
savefig(p2, "BDF4_V.png")
savefig(p3, "BDF4_SZW.png")