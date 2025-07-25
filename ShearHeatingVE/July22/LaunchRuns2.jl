using DelimitedFiles, Printf
include("TargetTestRefineGrid.jl")

P = readdlm("Param.csv", '\t')

istart = parse(Int, ARGS[1])
istop  = parse(Int, ARGS[2])

for i = istart : istop
    name = @sprintf("BFD4_%02d", i)
    nresy = Int64(P[i,1])
    tfac  = P[i,2]
    if nresy > 0 && tfac > 0
        ViscoElastic_CrankNicolson(name; θM=1.0, θT=1.0, nresy=nresy, nrest=128.0/tfac, adapt_dy=false, adapt_dt=false, viz=false, noisy=true)
    elseif nresy > 0
        ViscoElastic_CrankNicolson(name; θM=1.0, θT=1.0, nresy=nresy, nrest=1, adapt_dy=false, adapt_dt=true, viz=false, noisy=true)
    elseif tfac > 0
        ViscoElastic_CrankNicolson(name; θM=1.0, θT=1.0, nresy=1, nrest=128.0/tfac, adapt_dy=true, adapt_dt=false, viz=false, noisy=true)
    else
        ViscoElastic_CrankNicolson(name; θM=1.0, θT=1.0, nresy=1, nrest=1, adapt_dy=true, adapt_dt=true, viz=false, noisy=true)
    end
end