using Random

include("tools.jl")


function main(h::Hamiltonian, size::UInt8)

	#Temperature Path
	B_path = zeros(40)


	#Number of replicas
	N::UInt8 = 10

	replica_list::Vector{Replica} = Vector{Replica}(undef, N)


	state::Vector{UInt8} = Vector{UInt8}(undef, size)

	#Defining these is what Tameem suggested we research
	temperatures = Vector{UInt8}(1:N)

	#Generate Replicas
	for i = (1:N)

		state = rand([0,1], size)
		replica_list[i] = Replica(temperatures[i].^-1, state)

	end

	replica_track = replica_list[1]

	#This variable decides which pairs are considered for exchange
	#it is functionally a boolean
	toggle::UInt8 = 0

	j = 1
	while j <= 40

		toggle = toggle ⊻ 1

		#This appropriately produces the indices of the left-most element
		#of the pairs to be proposed for exchange

		if N % 2 == 1

			indices = (1 + toggle : 2 : N - 1 + toggle)
		else
			indices = (1 + toggle : 2 : N - toggle)
		end

		
		for replica in replica_list
			evolve!(replica, h)
		end

		for i in indices

			if propose_exchange(replica_list[i], replica_list[i+1], h)

				exchange!(replica_list, UInt8(i))

			end
		end
		

		
		B_path[j] = replica_track.B
		j += 1
		
	end

	for replica in replica_list

		println(replica.state)
	end

	#myplot = plot((1:j),B_path, seriestype = :scatter) 
	#display(myplot)
	
	#readline()

end

h = Hamiltonian(1, zeros(4,4))
main(h, UInt8(4))





#println(states)