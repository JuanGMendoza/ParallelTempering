#spin -1 -> 0
#spin 1 -> 1


#H = h * spins + Ji,j * spin_i * spin_j
mutable struct Hamiltonian

	#Local Field
	h::Float64

	#Interaction Terms
	J::Array{Float64}

end

mutable struct Replica

	#Inverse Temperature
	B::Float64

	state::Vector{UInt8}

end


#Markov chain evolution
function evolve!(replica::Replica ,h::Hamiltonian)

	random_bit::UInt8 = rand(1:length(replica.state))
	new_state = copy(replica.state)
	new_state[random_bit] = new_state[random_bit] ⊻ 1
	
	#VERIFY this algorithm is correct
	criterion = exp(replica.B*(evaluate_energy(replica.state, h) - evaluate_energy(new_state, h)))
	
	println("criterion: ", criterion)
	if criterion < 1

	 	if rand() < criterion

	 		replica.state[random_bit] = replica.state[random_bit] ⊻ 1
	 	end
	
	else
		replica.state[random_bit] = replica.state[random_bit] ⊻ 1
	end

	
end

function flip_random_bit(state::Vector{UInt8})

	state_copy::Vector{UInt8} = copy(state)

	flip::UInt8 = rand(1:length(state))

	state_copy[flip] = state_copy[flip] ⊻ 1

	return state_copy

end

function evaluate_energy(state::Vector{UInt8}, h::Hamiltonian)

	E::Float64 = 0

	for spin in state

			E = E - (-1)^Int8(spin)*h.h

	end

	for i in (1:length(state)-1)

		for j in (i+1:length(state))

			E = E + h.J[i,j]*(-1)^(state[i] ⊻ state[j])

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