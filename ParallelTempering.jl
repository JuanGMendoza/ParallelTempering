using Random
using Bits
include("tools.jl")


function parallel_tempering(h::Hamiltonian, fileName::String)


	
	#Number of replicas
	N::UInt8 = 10
	#Size of the states
	size::UInt8 = length(h.bonds)

	replica_list::Vector{Replica} = Vector{Replica}(undef, N)

	state_matrix::Vector{UInt128} = Vector{UInt128}(undef, size)
	row::UInt128 = 0
	#Defining these is what Tameem suggested we research
	temperatures::Vector{Float64} = Vector{Float64}(1:N)
	jldsave(fileName, measurement="magnetization",temps=temperatures)

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

	j = 1
	while j <= 1000

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
			printState(replica, state_matrix)
			evolve!(replica, h, state_matrix)
			printState()
		end

#=
		for i in indices

			if propose_exchange(replica_list[i], replica_list[i+1], h, state_matrix)

				exchange!(replica_list, UInt8(i))
				true

			end
		end
		=#
		j += 1
		
	end


end



