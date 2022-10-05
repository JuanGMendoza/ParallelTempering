using JLD2

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

	#Temperature
	T::UInt8

	#Inverse Temperature
	B::Float64

	state::Vector{UInt8}

	#Number to identify the replica as it moves through B space
	ID::UInt8

end


#Markov chain evolution
function evolve!(replica::Replica ,h::Hamiltonian)

	random_bit::UInt8 = rand(1:length(replica.state))
	new_state = copy(replica.state)
	new_state[random_bit] = new_state[random_bit] ⊻ 1
	
	#VERIFY this algorithm is correct
	criterion = exp(replica.B*(evaluate_energy(replica.state, h) - evaluate_energy(new_state, h)))

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

#Returns the energy corresponding to the input state
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

	delta::Float64 = -(replica1.B - replica2.B)*(evaluate_energy(replica1.state, h) - evaluate_energy(replica2.state, h))

	#println(delta)
	if delta < 0

		return true

	elseif delta >= 0

		W::Float64 = exp(-delta)

	end

	if rand() < W 
		return true 
	else 
		return false
	end

end

#exchanges the replicas (but not their temperature), replica_list[index] <-> replica_list[index + 1]
function exchange!(replicaList::Vector{Replica}, index::UInt8)

	#println(replicaList)
	#=
	replica1::Replica = replica_list[index]

	temp_T::Float64 = replica_list[index].T

	replica_list[index].T = replica_list[index + 1].T

	replica_list[index + 1].T = temp_T

	replica_list[index] = replica_list[index + 1]

	replica_list[index + 1] = replica1

	replica_list[index].B = 1/replica_list[index].T
	replica_list[index + 1].B = 1/replica_list[index + 1].T 

	=#

	#Discuss efficiency of this
	temp_state::Vector{UInt8} = replicaList[index].state
	temp_ID::UInt8 = replicaList[index].ID

	replicaList[index].ID = replicaList[index + 1].ID
	replicaList[index].state = replicaList[index + 1].state

	replicaList[index + 1].ID = temp_ID
	replicaList[index + 1].state = temp_state


	#println(replicaList)

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


function save_history(filename::String, replicas::Vector{Replica}, t::Int64)

	jldopen(filename, "a+") do file
		
		i = 1
		groupString = "t" * string(t) * "/" * "replica" 
		#println(groupString)

		for replica in replicas

			file[groupString * string(i)] = replica

			i += 1

		end
		
	end
end

#Returns all replicas at temperature T from file fileName, ordered by timestep
function load_T_history(fileName::String, T::UInt8)

	
	#characters in the word replica + 1
	replicaLength = 8
	replicaList::Vector{Replica} = []
	jldopen(fileName, "r") do file

		timesteps = length(keys(file))
		replicasPerTimestep = parse(UInt8, last(keys(file["t1"]))[replicaLength:end])
		replicaList = Vector{Replica}(undef, timesteps)
		indexOfDesiredT = 0

		for j in (1:replicasPerTimestep)

			 if file["t1/replica" * string(j)].T == T
			 	indexOfDesiredT = j
			 	break
			 end
		end 
		

		for k in (1:timesteps)
			
			replicaList[k] = file["t" * string(k) * "/replica" * string(indexOfDesiredT)]
			
		end
	end
	
	return replicaList
end

function load_ID_history(filename::String, ID::UInt8)
	
	#characters in the word replica + 1
	replicaLength = 8	

	replicaList::Vector{Replica} = []

	jldopen(filename) do file
		timesteps = length(keys(file))
		replicaList = Vector{Replica}(undef, timesteps)
		replicasPerTimestep = parse(UInt8, last(keys(file["t1"]))[replicaLength:end])
		for i in (1:timesteps)
			for j in (1:replicasPerTimestep)
				replica = file["t"*string(i)]["replica"*string(j)]
				if replica.ID == ID
					replicaList[i] = replica 
				end
			end
		end

	end
	return replicaList
end
#load_T_history("test_history.jld2", UInt8(10))
function calculate_expectation(operator::Function, replicas::Vector{Replica})

	average::Float64 = 0

	for replica in replicas
		average += operator(replica.state)
	end
	return average/length(replicas)

end

function magnetization(state::Vector{UInt8})

	m = 0 
	for spin in state

		m += (-1)^(spin ⊻ 1)
	end
	return m
end

