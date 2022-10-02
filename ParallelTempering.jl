using Random

include("tools.jl")


function main(h::Hamiltonian, size::UInt8)


	history = "test_history_2.jld2"

	#Number of replicas
	N::UInt8 = 10

	replica_list::Vector{Replica} = Vector{Replica}(undef, N)

	state::Vector{UInt8} = Vector{UInt8}(undef, size)

	#Defining these is what Tameem suggested we research
	temperatures::Vector{UInt8} = Vector{UInt8}(1:N)

	#Generate Replicas
	for i = (1:N)

		state = rand([0,1], size)
		replica_list[i] = Replica(temperatures[i], temperatures[i].^-1, state, i)

	end

	#This variable decides which pairs are considered for exchange
	#it is functionally a boolean
	toggle::UInt8 = 0

	j = 1
	while j <= 10000

		toggle = toggle ⊻ 1

		#This appropriately produces the indices of the left-most element
		#of the pairs to be proposed for exchange

		if N ⊻ 1 == 1

			indices = (1 + toggle : 2 : N - 1 + toggle)
		else
			indices = (1 + toggle : 2 : N - toggle)
		end

		save_history(history, replica_list, j)

		for replica in replica_list
			evolve!(replica, h)
		end

		for i in indices

			if propose_exchange(replica_list[i], replica_list[i+1], h)

				exchange!(replica_list, UInt8(i))

			end
		end
		
		j += 1
		
	end

	for replica in replica_list

		println(replica.state)
	end


end

h = Hamiltonian(1, zeros(4,4))
main(h, UInt8(2))





#println(states)