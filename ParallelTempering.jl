using Random
using Bits
include("tools.jl")


"""
    parallel_tempering(h::Hamiltonian, temperatures::Vector{Float64}, fileName::String, measurement::Function, name::String, steps::UInt64)

Run parallel tempering given all the parameters, measure a specific operator for every replica in each step, and store the value on disk.

# Arguments

- `n::Integer`: the number of elements to compute.

- `dim::Integer=1`: the dimensions along which to perform the computation.

- `h::Hamiltonian`: the hamiltonian of the system to be simulated. More details on this structure in tools.jl.

- `temperatures::Vector{Float64}`: a sorted list (increasing) of the temperatures to assign to each corresponding replica. This determines the number of replicas.

- `fileName::String`: the name to be assigned to the file on disk containing measurement outcomes. Omit any file extensions.

- `measurement::Function`: the measurement to be calculated at every iteration. A function that takes in a state vector and returns a scalar. See tools.jl/magnetization.

- `name::Function`: the name of the measurement function. Ex: magnetization

- `steps::UInt64`: the number of steps to run the simulation for.


"""
function parallel_tempering(h::Hamiltonian, temperatures::Vector{Float64}, fileName::String, measurement::Function, name::String, steps::UInt64)


	#Size of the states
	size::UInt8 = length(h.bonds)

	#Number of Replicas
	N = length(temperatures)

	replica_list::Vector{Replica} = Vector{Replica}(undef, N)
	state_matrix::Vector{UInt128} = Vector{UInt128}(undef, size)
	row::UInt128 = 0
	
	energy::Float64 = 0

	jldsave(fileName * "_meas.jld2", measurement=name ,temps=temperatures)
	indices::StepRange{Int64, Int64} = (0:0)

	#Generate Replicas
	for i = (1:size)

		row = rand(0:(UInt128(2)^N)-1)
		state_matrix[i] = row
	end

	for k in (1:N)
		energy = evaluate_energy(UInt8(k), h, state_matrix)
		replica_list[k] = Replica(temperatures[k], temperatures[k].^-1, k, Queue{UInt64}(), energy)
		refill_random_bits!(replica_list[k], size)
	end

	#This variable decides which pairs are considered for exchange
	#it is functionally a boolean
	toggle::UInt8 = 0

	for j in (1:steps)

		toggle = toggle ⊻ 1

		#This appropriately produces the indices of the left-most element
		#of the pairs to be proposed for exchange

		if N ⊻ 1 == 1

			indices = (1 + toggle : 2 : N - 1 + toggle)
		else
			indices = (1 + toggle : 2 : N - toggle)
		end

		save_measurements(fileName, replica_list, j, measurement, state_matrix)
		#save_all_history(fileName, replica_list, j)

		for replica in replica_list

			evolve!(replica, h, state_matrix)

		end


		for i in indices

			if propose_exchange(replica_list[i], replica_list[i+1], h, state_matrix)

				exchange!(replica_list, UInt8(i))

			end
		end
		
	end


end



