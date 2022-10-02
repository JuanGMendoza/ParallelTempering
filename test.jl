include("tools.jl")

k_B = 1#1.380649 * 10^-23

function energy_expec_value(h, T)

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

h = Hamiltonian(1, zeros(4,4))

for i in (1:10)
replicas = load_T_history("test_history_2.jld2", UInt8(i))
println("T= ",i)
println("Ex: ",energy_expec_value(h, i))
println("PT: ",calculate_expectation(magnetization, replicas), '\n')

end

