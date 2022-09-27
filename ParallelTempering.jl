using Random
using Plots

#How many bits define a state
#If you change this variable you must change the data type for all other state variables
STATE_LENGTH = 8

mutable struct Hamiltonian

	#Local Field
	h::Float64

	#Interaction Terms
	J::Array{Float64}

end

mutable struct Replica

	#Inverse Temperature
	B::Float64

	location::UInt8

end


#Markov chain evolution
function evolve!(state::Vector{UInt8}, temperature::Float64 ,h::Hamiltonian )

	new_state = flip_random_bit(state)

	#VERIFY this algorithm is correct
	criterion = exp(replica.B*(evaluate_energy(replica.state, h) - evaluate_energy(new_state, h)))
	
	if criterion < 1

	 	if rand() < criterion

	 		replica.state = new_state
	 	end
	
	else
		replica.state = new_state
	end

end

function flip_random_bit!(state::Vector{UInt8})

	flip::UInt8 = rand(1:length(state))

	state[flip] = state[flip] ⊻ 1


end

function evaluate_energy(state::UInt8, h::Hamiltonian)

	E::Float64 = 0

	string = bitstring(state)

	for spin in string

			E = E - (-1)^Int8(spin)*h.h

	end

	for i in (1:STATE_LENGTH-1)

		for j in (i+1:STATE_LENGTH)

			E = E + h.J[i,j]*(-1)^((UInt8(string[i]) ⊻ UInt8(string[j])))

			#println(string[i],string[j], ' ',h.J[i,j]*(-1)^((UInt8(string[i]) + UInt8(string[j])) % 2) )
		end
	end
	return E
end


function propose_exchange(replica1::Replica, replica2::Replica, h::Hamiltonian)

	delta = (replica1.B - replica2.B)*(evaluate_energy(replica1.state, h) - evaluate_energy(replica2.state, h))

	if delta < 0

		return true

	elseif delta >= 0

		W = exp(-delta)

	end

	if rand() < W 
		return true 
	else 
		return false
	end

end

#exchanges the replicas (but not their temperature), replica_list[index] <-> replica_list[index + 1]
function exchange!(replica_list::Vector{Replica}, index::UInt8)

	replica1 = replica_list[index]

	temp_B = replica_list[index].B

	replica_list[index].B = replica_list[index + 1].B 

	replica_list[index + 1].B = temp_B

	replica_list[index] = replica_list[index + 1]

	replica_list[index + 1] = replica1

end

#Discuss equation 4.4 on Hukushima with Tameem
function autocorrelation()

	p = rand()
	if p < .05
		return true
	end

	return false
end

function brute_force_ground_state(h::Hamiltonian)

	temp::Float64 = 0
	min::Float64 = typemax(Float64)
	minState::UInt8 = 0

	for i in (0:2^7)
		temp = evaluate_energy(UInt8(i), h)
		if temp < min
			min = temp
			minState = i
		end
	end
	return min, minState
end

function main(h::Hamiltonian)

	#Temperature Path
	B_path = zeros(40)


	#Number of replicas
	N::UInt8 = 10

	replica_list = Array{Replica}(undef, N)

	#Defining these is what Tameem suggested we research
	temperatures = Array{UInt8}(1:N)

	#2-D array where each row represents the bitstring state of a replica
	#states::Array{Int8} = zeros(num_replicas, STATE_LENGTH)

	#Generate Replicas
	for i = (1:N)

		replica = Replica(temperatures[i].^-1, rand(0:2^(STATE_LENGTH)-1))

		replica_list[i] = replica
	end

	replica_track = replica_list[1]

	#This variable decides which pairs are considered for exchange
	#it is functionally a boolean
	toggle::UInt8 = 0

	j = 1
	while j < 40

		toggle = (toggle + 1) % 2

		#This appropriately produces the indices of the left-most element
		#of the pairs to be proposed for exchange

		if N % 2 == 1

			indices = (1 + toggle : 2 : N - 1 + toggle)
		else
			indices = (1 + toggle : 2 : N - toggle)
		end

		
		for i in indices

			if propose_exchange(replica_list[i], replica_list[i+1], h)

				exchange!(replica_list, UInt8(i))

			end
		end
		

		for replica in replica_list
			evolve!(replica, h)
		end
		
		B_path[j] = replica_track.B
		j += 1

		stop = autocorrelation()
		
	end

	for replica in replica_list

		println(bitstring(replica.state))
	end

	myplot = plot((1:j),B_path, seriestype = :scatter) 
	#savefig(myplot, "plot.png")
	display(myplot)
	
	readline()

end

#h::Hamiltonian = Hamiltonian(3,ones(8,8))

#println(bitstring(brute_force_ground_state(h)[2]))

#main(h)

num_replicas = 3

states::Array{UInt8} = zeros(num_replicas, STATE_LENGTH)


state = states[1,:]

state1::Vector{UInt8} = [1,1,1,1]

println(state1)
state[1] = 1

flip_random_bit!(state1)
println(state1)


#println(states)