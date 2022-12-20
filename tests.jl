include("tools.jl")
include("ParallelTempering.jl")
using Statistics

k_B = 1#1.380649 * 10^-23


function test1(temperature=1, correct=-3.0463766238230594, hamiltonian="1")

	outcome::Bool = false
	feedback::String = ""
	samples::Vector{Float64} = Vector{Float64}(undef, 10)
	PT::Float64 = 0
	sigma::Float64 = 0
	difference::Float64 = 0

	if temperature == 1 && hamiltonian == "1"
		println("\nComplete simulation tests:\n")
	end

	for sample in (1:10)
		samples[sample] = load_and_calc_expectation("unit-test" * hamiltonian * "_" * string(sample) * ".jld2", Float64(temperature))[2]
	end

	PT = mean(samples)
	sigma = std(samples)
	difference = abs(correct - PT)

	
	if difference < sigma
		outcome = true		
	end

	if hamiltonian == "1"
		feedback = " | B field 1.0, no inter. | <M> | T= " * string(temperature) * " | Ex = " * string(correct) * " | PT = " * string(PT) * " ± " * string(sigma) 
	else
		feedback = " | B field 2.0, 2body inter.(all to all) 1.0 | <M> | T= " * string(temperature) * " | Ex = " * string(correct) * " | PT = " * string(PT) * " ± " * string(sigma) 
	end
	return outcome, feedback

end

function test2()
	return test1(2, -1.8484686290400385)
end

function test3()
	return test1(5, -0.789501280899616)
end

function test4()
	return test1(9, -0.4426244420989519)
end

function test5()
	return test1(10, -0.3986719784998232)
end

function test6()
	return test1(1, -1.7254383196404988, "2")
end

function test7()
	return test1(2, -1.4024877844260495, "2")
end

function test8()
	return test1(5, -0.9401356167065131, "2")
end

function test9()
	return test1(9, -0.6475840976857671, "2")
end

function test10()
	return test1(10, -0.600244381967262, "2")
end

#Individual Function Tests tools.jl

#function energy_difference(ID::UInt8, different_spin::UInt8, state_matrix::Vector{UInt128}, hamiltonian::Hamiltonian)


#2 body interaction test
function test11()


	passed = false
	hami = Hamiltonian(0, [ [[UInt64(2)]], [[UInt64(1)]], [[]]], [[2.0],[2.0], []])
	replicas = 2
	spins = 3
	state_matrix::Vector{UInt128} = Vector{UInt128}(undef, 3)

	state_matrix[1] = UInt128(2)
	state_matrix[2] = UInt128(3)
	state_matrix[3] = UInt128(2)

	different_spin = UInt64(1)

	output = energy_difference(UInt8(2), different_spin, state_matrix, hami)
	correct = -4.0

	if output == correct

		passed = true
	end


	return passed, " | input 000 | bond 1,2 (2.0) | flip 1 | output " * string(output) * " | should be -4.0" 
end

#3 body interaction test
function test12()


	passed = false
	hami = Hamiltonian(0, [ [[2,3]], [[1,3]], [[1,2]]], [[2.0],[2.0], [2.0]])
	replicas = 2
	spins = 3
	state_matrix::Vector{UInt128} = Vector{UInt128}(undef, 3)

	state_matrix[1] = UInt128(2)
	state_matrix[2] = UInt128(3)
	state_matrix[3] = UInt128(2)

	different_spin = UInt64(1)

	output = [energy_difference(UInt8(2), different_spin, state_matrix, hami), energy_difference(UInt8(1), different_spin, state_matrix, hami)]
	correct = [4.0, 4.0]

	if output == correct

		passed = true
	end


	return passed, " | input -1-1-1 and 1-11 | bond 1,2,3 (2.0) | flip 1 | output " * string(output) * " | should be " * string(correct) 
end



#External B test
function test12()


	passed = false
	hami = Hamiltonian(15.5, [ [[]], [[]], [[]]], [[],[], []])
	replicas = 2
	spins = 3
	state_matrix::Vector{UInt128} = Vector{UInt128}(undef, 3)

	state_matrix[1] = UInt128(2)
	state_matrix[2] = UInt128(3)
	state_matrix[3] = UInt128(2)

	different_spin = UInt64(1)

	output = [energy_difference(UInt8(1), different_spin, state_matrix, hami), energy_difference(UInt8(1), different_spin+1, state_matrix, hami)]
	correct = [-2*15.5, 2*15.5]

	if output == correct

		passed = true
	end


	return passed, " | input -11-1 | ext. B 15.5 | flip 1 then 2 | output " * string(output) * " | should be " * string(correct) 
end

#3 body interaction and external B test
function test13()


	passed = false
	hami = Hamiltonian(3, [ [[2,3]], [[1,3]], [[1,2]]], [[2.0],[2.0], [2.0]])
	replicas = 2
	spins = 3
	state_matrix::Vector{UInt128} = Vector{UInt128}(undef, 3)

	state_matrix[1] = UInt128(2)
	state_matrix[2] = UInt128(2)
	state_matrix[3] = UInt128(2)

	different_spin = UInt64(1)

	output = [energy_difference(UInt8(2), different_spin, state_matrix, hami), energy_difference(UInt8(1), different_spin, state_matrix, hami)]
	correct = [10.0, -10.0]

	if output == correct

		passed = true
	end


	return passed, " | input -1-1-1 and 111 | bond 1,2,3 (2.0) + ext. B 3.0 | flip 1 | output " * string(output) * " | should be " * string(correct) 
end
