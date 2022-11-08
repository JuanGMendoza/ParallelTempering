using JLD2
using DataStructures
using Random


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

	bitsToFlip::Queue{UInt64}


end


#Markov chain evolution
function evolve!(replica::Replica ,h::Hamiltonian)

	
	if isempty(replica.bitsToFlip)
		refill_random_bits!(replica, length(h.J[1,:]))
	end
	random_bit::UInt8 = dequeue!(replica.bitsToFlip)
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

	#Discuss efficiency of this
	temp_state::Vector{UInt8} = replicaList[index].state
	temp_ID::UInt8 = replicaList[index].ID

	replicaList[index].ID = replicaList[index + 1].ID
	replicaList[index].state = replicaList[index + 1].state

	replicaList[index + 1].ID = temp_ID
	replicaList[index + 1].state = temp_state


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


#Creates a file containing measurements for input operator for every temperature
#Receives a replica list whose order matches the temperature list
# t represents the timestep of the simulation
function save_measurements(fileName::String, replicas::Vector{Replica}, t::Int64, operator::Function)

	groupString = "t" * string(t) * "/" 

	jldopen(fileName, "a+") do file

		i::UInt8 = 1
		for replica in replicas

			
			file[groupString * string(i) * "/measurement"] = operator(replica.state)
			i = i + 1
		end
	end
end

#returns the expectation value of all measurements stored in fileName at temperature T
function load_and_calc_expectation(fileName::String, T::Float64)

	expectation::Float64 = 0
	timesteps::UInt64 = 0
	name::String = ""
	jldopen(fileName, "r") do file
		name = file["measurement"]

		timesteps = length(keys(file)) - 2


		indexOfDesiredT = 0
		replicasPerTimestep = length(file["temps"])

		for j in (1:replicasPerTimestep)

			if file["temps"][j] == T
			 	indexOfDesiredT = j
			 	break
			end
		end 

		if indexOfDesiredT == 0

			println("The requested T, was not found amongst the temperatures available in " * filename * string(file["temps"]))

			expectation = 0

		

		else
			for j in (1:timesteps)
				
				expectation += file["t" * string(j) * "/" * string(indexOfDesiredT) * "/measurement"]
			end
		end
	end

	return name, expectation/timesteps
end


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


#This function takes in a replica and the state size and it
#replaces the replica.bitsToFlip queue to a new random one
function refill_random_bits!(replica::Replica, size)

	randomList::Vector{UInt64} = collect(1:size)

	shuffle!(randomList)

	bitsToFlip::Queue{UInt64} = Queue{UInt64}()

	for bit in randomList

		enqueue!(bitsToFlip, bit)

	end

	replica.bitsToFlip = bitsToFlip

end

#The functions below work with different file formatting#
#########################################################


function save_all_history(fileName::String, replicas::Vector{Replica}, t::Int64)

	jldopen(fileName, "a+") do file
		
		i = 1
		groupString = "t" * string(t) * "/" * "replica" 
		#println(groupString)

		for replica in replicas

			file[groupString * string(i)] = replica

			i += 1

		end
		
	end
end

function autocorrelation(fileNames::Vector{String}, ID::UInt8)

	numFiles = length(fileNames)
	replicas::Vector{Vector{Replica}} = Vector{Vector{Replica}}(undef,length(fileNames))
	timesteps::UInt64 = 0
	stateLength::UInt64 = 0
	average::Float64 = 0
	sum::Float64 = 0
	j::UInt8 = 1

	for fileName in fileNames
		#println(fileName)
		replicas[j] = load_ID_history(fileName, ID)
		j += 1
	end

	timesteps = length(replicas[1])
	println("t",timesteps)
	stateLength = length(replicas[1][1].state)
	
	q::Vector{Float64} = Vector{Float64}(undef, timesteps)


	for t in (1:timesteps)
		sum = 0
		for i in (1:stateLength)
			average = 0
			for fileIndex in (1:numFiles)
				average +=  (-1)^(replicas[fileIndex][1].state[i] ⊻ replicas[fileIndex][t].state[i])
			end
			average = average/numFiles
			sum += average
		end
		q[t] = sum/ stateLength

	end

	return q	
end

#This function tracks a replica through its evolution in temperature space
#It returns the replica structure at every timestep
#fileName must have been produced using save_all_history
function load_ID_history(fileName::String, ID::UInt8)
	
	#characters in the word replica + 1
	replicaLength = 8	

	replicaList::Vector{Replica} = []


	jldopen(fileName) do file
		#println("h ", length(keys(file)))
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


#Returns all replicas at temperature T from file fileName, ordered by timestep
#fileName must contain all replica history, created with save_all_history()
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


