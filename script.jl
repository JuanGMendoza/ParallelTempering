using Plots
using JLD2
include("tools.jl")

timesteps = 250
fileNames = Vector{String}()
println(timesteps)
q_average::Vector{Float64} = zeros(timesteps)
q_list = Vector{Vector{Float64}}(undef, 10)
i = 0

for index in (3:10)
	push!(fileNames, "./History_Files/all_historyh1_" * string(index) * ".jld2")
end


for ID in (1:10)
	q_list[ID] = autocorrelation(fileNames, UInt8(ID))

end

for t in (1:timesteps)

	for ID in (1:10)
		q_average[t] += q_list[ID][t]
	end
	q_average[t] = q_average[t]/10 + 1
end

plot(q_average, seriestype = :scatter, yaxis=:log, xaxis=:log)
savefig("autocorrelationLogx.png")