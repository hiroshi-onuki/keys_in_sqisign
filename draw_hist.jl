using Plots

function read_file(filename::String)
    fp = open(filename, "r")
    lines = readlines(fp)
    data = Int[]
    for l in lines
        push!(data, parse(Int, l))
    end
    close(fp)
    return data
end

N = 251
for N in [211, 223, 227, 239, 251, 1019, 2003]
    ps = []
    for ex in [0, 5, 11]
        data = read_file("data/secret_$(N)_$(ex)_full.txt")
        push!(ps, histogram(data, bin=range(1,N+1, step=1), ylims=(0,250),
            label="", title="ε=$(ex)", linewidth=0, xtickfontsize=11, ytickfontsize=11))
    end
    plot(ps[1],ps[2],ps[3], layout=(1,3), size=(1000,250))
    savefig("fig/hist_$(N).png")
end
