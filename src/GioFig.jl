module GioFig

greet() = print("Hello World!")

import Ripserer
import Plots
import Graphs
import SimpleWeightedGraphs as SWG
import Distances

using Printf
using LinearAlgebra


sample_data() = [Float64.((cos(θ),sin(θ))) for θ=map(k->2π*k/8, 0:7)]

function sample_data2()
    data = sample_data()
    foo = [(0.6*cos(θ), 0.6*sin(θ)) for θ=map(k->2π*k/8, 0:7)]
    data2 = [data; map(x->x.+(1,1),foo)]
    return data2
end
    
function plot_rips_subcomplexes(data::AbstractVector)
    rips = Ripserer.Rips(data; metric=Distances.Euclidean(1e-6))
    # thresholds = sort(collect(Set(Float16.((Ripserer.birth).(Ripserer.edges(rips))))))
    # T = (typeof(rips).parameters)[2]    
    # _plot_subcomplexes(data, rips, T.(thresholds))

    T = (typeof(rips).parameters)[2]
    cuts = Vector{T}()
    rawthresholds = sort(collect(Set((Ripserer.birth).(Ripserer.edges(rips)))))
    print(rawthresholds, "\n")    
    push!(cuts, pop!(rawthresholds))
    
    while length(rawthresholds) > 0
        cur_t = pop!(rawthresholds)
        if !isapprox(cur_t, cuts[end])
            push!(cuts, cur_t)
        end
        if length(rawthresholds) == 0
            break
        end
    end
    print(cuts,"\n")
    _plot_subcomplexes(data, rips, cuts)


    result = Ripserer.ripserer(rips)
    Plots.plot(result[2])
    Plots.png("pd1")
    return nothing
end


function _plot_subcomplexes(
    data::AbstractVector,
    rips::Ripserer.AbstractRipsFiltration{I,T},
    thresholds::Vector{T}
    ) where {I,T}
    sorted_thresholds = sort(thresholds)
    
    # figure out 2-cliques (triangles)
    relevant_edges = [x for x in Ripserer.edges(rips) if x.birth <= sorted_thresholds[end]]
    st = map(Ripserer.vertices, relevant_edges)
    
    g = SWG.SimpleWeightedGraph(
        getfield.(st, 1),
        getfield.(st, 2),
        map(Ripserer.birth, relevant_edges)
    )

    triangles = Vector{Set{Int64}}()
    for e in Graphs.edges(g)
        common = intersect(Graphs.neighbors(g, Graphs.src(e)), Graphs.neighbors(g, Graphs.dst(e)))
        for vert in common
            push!(triangles, Set([Graphs.src(e),Graphs.dst(e), vert]))
        end
    end
    triangles = Set(triangles)
    relevant_triangles = map(x->Ripserer.simplex(rips, Val(2), Tuple(x)), collect(triangles))

    simplices = sort([relevant_triangles; relevant_edges]; by=Ripserer.birth, rev=true)
    Plots.scatter(data; markersize=8, legend=false, aspect_ratio=1, axis=([], false))
    Plots.png(@sprintf("%.12f", 0.0))
    for t in sorted_thresholds
        while length(simplices) > 0 && Ripserer.birth(simplices[end]) <= t
            spx = pop!(simplices)
            print(spx, Ripserer.dim(spx), "\n")
            if Ripserer.dim(spx) == 1
                Plots.plot!(spx, data; linewidth=4, linealpha=1, linecolor=:black)
            else
                Plots.plot!(spx, data; linewidth=4, linealpha=1, linecolor=:black, fill=(0,0.5,:blue))
            end            
        end
        Plots.png(@sprintf("%.12f",t))
        # Plots.png(string(t))
    end
    return nothing
end

end # module GioFig
