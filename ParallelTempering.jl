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


function main()

	#Number of replicas
	N::UInt8 = 10


	replica_list = Array{Replica}(undef, N)

	#Come back and define a way to fill this
	system_sizes = Array{UInt8}(1:N)

	temperatures = Array{UInt8}(1:N)

	

	
end

main()