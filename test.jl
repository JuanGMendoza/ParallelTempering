include("tools.jl")
include("ParallelTempering.jl")

using Plots

k_B = 1#1.380649 * 10^-23

function magnetization_expec_value(h, T)

	stateMax = (2^size(h.J,1)) - 1
	z = 0
	expectation = 0
	unnormal_probabilities = Vector{Float64}(undef, stateMax+1)
	i = 1
	magnetizations = Vector{Float64}(undef, stateMax+1)
	for state in (0:stateMax)
		stateVector = split(bitstring(state)[64 - size(h.J,1) + 1:end],"")
		stateVectorInt = Vector{UInt8}(undef, length(stateVector))
		

		for j in (1:length(stateVector))
			stateVectorInt[j] = parse(UInt8, stateVector[j])
		end
		
		unnormal_probabilities[i] = exp(-evaluate_energy(stateVectorInt, h)/(k_B*T))
		magnetizations[i] = magnetization(stateVectorInt)

		i += 1
	end

	for k in (1:stateMax+1)

		z += unnormal_probabilities[k]

	end
	
	for l in (1:stateMax+1)
		expectation += unnormal_probabilities[l]*magnetizations[l]
	end
	return expectation/z

end

h1 = Hamiltonian(1, zeros(4,4))

main(h1, "unit-test1.jld2")

J = ones(Float64, (4,4))
h2 = Hamiltonian(2, J)

println("50% Complete")
main(h2, "unit-test2.jld2")

println("H1 Test:")

for i in (1:10)

	println("T= ", i)
	PT = load_and_calc_expectation("unit-test1.jld2", Float64(i))[2]
	@assert PT != 0
	difference = magnetization_expec_value(h1, i) - PT
	
	println(difference)

end


println("H2 Test:")
for i in (1:10)


	println("T= ", i)
	PT = load_and_calc_expectation("unit-test2.jld2", Float64(i))[2]
	@assert PT != 0
	difference = magnetization_expec_value(h2, i) - PT
	
	println(difference)

end
