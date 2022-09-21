using Random

#How many bits define a state
#If you change this variable you must change the data type for all other state variables
STATE_LENGTH = 8

mutable struct Hamiltonian

	placeHolder::UInt8
end

mutable struct Replica

	#Size
	L::UInt8

	#Inverse Temperature
	B::Float64

	state::UInt8

end

#Begin definition of all the Replica class functions
function change_temperature(replica::Replica, new_B::Float64)

	replica.B = new_B

end

#Markov chain evolution
function evolve(replica::Replica, h::Hamiltonian)

	new_state = flip_random_bit(replica.state)

	#VERIFY this algorithm is correct
	criterion = exp(-replica.B*(evaluate_energy(replica.state, h) - evaluate_energy(new_state, h)))
	println(criterion)
	
	if criterion < 1

	 	if rand() < criterion

	 		replica.state = new_state
	 	end
	
	else
		replica.state = new_state
	end

end

function flip_random_bit(bits::UInt8)

	flip::UInt8 = 2^rand(0:STATE_LENGTH-1)

	#XOR
	return bits ⊻ flip

end

#To finish later
function evaluate_energy(state::UInt8, h::Hamiltonian)

return rand(1:5)

end


function propose_echange(replica1, replica2, h)

	delta = (replica1.B - replica2.B)*(evaluate_energy(replica1.state, h) - evaluate_energy(replica2.state, h))

	if delta > 0

		W = exp(-delta)

	elseif delta <= 0

		W = 1

	end

	if rand() < W 
		return true 
	else 
		return false
	end

end

function main()

	#Number of replicas
	N::UInt8 = 10


	replica_list = Array{Replica}(undef, N)

	#Come back and define a way to fill this
	system_sizes = Array{UInt8}(1:N)

	#Defining these is what Tameem suggested we research
	temperatures = Array{UInt8}(1:N)

	#Generate Replicas
	for i = (1:N)

		replica = Replica(system_sizes[i], temperatures[i], rand(0:2^(STATE_LENGTH)-1))

		replica_list[i] = replica
	end

	stop::Bool = false

	#This variable decides which pairs are considered for exchange
	#it is functionally a boolean
	toggle::UInt8 = 0
	N = 6
	while stop == false

		toggle = (toggle + 1) % 2

		#This appropriately produces the indeces of the left-most element
		#of the pairs to be proposed for exchange

		if N % 2 == 1

			indeces = (1 + toggle : 2 : N - 1 + toggle)
		else
			indeces = (1 + toggle : 2 : N - toggle)
		end

		for i in indeces
			println(i)
			readline()
		end
	end

	
end

main()