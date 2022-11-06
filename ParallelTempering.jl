using Random

include("tools.jl")


function main(h::Hamiltonian, fileName::String)

	
	#Number of replicas
	N::UInt8 = 10
	#Size of the states
	size::UInt8 = length(h.J[1,:])

	replica_list::Vector{Replica} = Vector{Replica}(undef, N)

	state::Vector{UInt8} = Vector{UInt8}(undef, size)

	#Defining these is what Tameem suggested we research
	temperatures::Vector{Float64} = Vector{Float64}(1:N)
	#jldsave(fileName, measurement="magnetization",temps=temperatures)

	#Generate Replicas
	for i = (1:N)

		state = rand([0,1], size)
		replica_list[i] = Replica(temperatures[i], temperatures[i].^-1, state, i, Queue{UInt64}())
		refill_random_bits!(replica_list[i], size)

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

		#save_measurements(fileName, replica_list, j, magnetization)
		save_all_history(fileName, replica_list, j)

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


end


J = ones(Float64, (4,4))
h2 = Hamiltonian(2, J)

for i in (3:10)
	println(i)
	main(h2, "all_historyh2_"* string(i) * ".jld2")
end