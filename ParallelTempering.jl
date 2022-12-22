using Random
using Bits
include("tools.jl")


function parallel_tempering(h::Hamiltonian, temperatures::Vector{Float64}, fileName::String)


	#Size of the states
	size::UInt8 = length(h.bonds)

	#Number of Replicas
	N = length(temperatures)
	replica_list::Vector{Replica} = Vector{Replica}(undef, N)
	state_matrix::Vector{UInt128} = Vector{UInt128}(undef, size)
	row::UInt128 = 0
	
	energy::Float64 = 0
	jldsave(fileName, measurement="magnetization",temps=temperatures)
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

	max_steps::UInt16 = 1000
	for j in (1:max_steps)

		toggle = toggle ⊻ 1

		#This appropriately produces the indices of the left-most element
		#of the pairs to be proposed for exchange

		if N ⊻ 1 == 1

			indices = (1 + toggle : 2 : N - 1 + toggle)
		else
			indices = (1 + toggle : 2 : N - toggle)
		end

		save_measurements(fileName, replica_list, j, magnetization, state_matrix)
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



