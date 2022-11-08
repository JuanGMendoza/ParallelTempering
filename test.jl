include("tools.jl")
include("ParallelTempering.jl")
using Statistics

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


function run_tests()

	h1 = Hamiltonian(1, zeros(4,4))
	J = ones(Float64, (4,4))
	h2 = Hamiltonian(2, J)
	println("\n")
	for i in (1:10)

		#main(h1, "unit-test1_" * string(i) * ".jld2")
		#main(h2, "unit-test2_" * string(i) * ".jld2")

		
		print("Running simulations...[" *string(i) *"0%]\r")
		
	end
	



	numTests = 0
	passedTests = 0

	println("\nH1 Test:\n")

	for i in (1:10)

		samples::Vector{Float64} = Vector{Float64}(undef, 10)
		for sample in (1:10)
			samples[sample] = load_and_calc_expectation("unit-test1_" * string(sample) * ".jld2", Float64(i))[2]
		end
		PT = mean(samples)
		sigma = std(samples)
		Exact = magnetization_expec_value(h1, i)
		difference = Exact - PT
		

		if difference > sigma
			print("FAILED ")		
		else
			print("PASSED ")
			passedTests += 1
		end
		numTests += 1

		print("| T= ", i)
		print(" | Ex = ", Exact, " | PT = ",PT, " ± ", sigma, "\n\n")


	end


	println("\nH2 Test:\n")
	for i in (1:10)

		samples::Vector{Float64} = Vector{Float64}(undef, 10)
		for sample in (1:10)
			samples[sample] = load_and_calc_expectation("unit-test2_" * string(sample) * ".jld2", Float64(i))[2]
		end
		
		PT = mean(samples)
		sigma = std(samples)
		Exact = magnetization_expec_value(h2, i)
		difference = abs(Exact - PT)

		
		

		if difference > sigma
			print("FAILED ")		
		else
			print("PASSED ")
			passedTests += 1
		end
		numTests += 1

		print("| T= ", i)
		print(" | Ex = ", Exact, " | PT = ",PT, " ± ", sigma, "\n\n")

	
	end
	println(passedTests, "/", numTests, " tests passed")
end

run_tests()