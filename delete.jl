
include("tools.jl")

hami = Hamiltonian(2, [ [[]], [[]], [[]], [[]] ], [[],[], []])

state_matrix = Vector{UInt128}(undef, 4)

state_matrix[1] = UInt8(0)
state_matrix[2] = UInt8(0)
state_matrix[3] = UInt8(0)
state_matrix[4] = UInt8(0)

rep = Replica(1,1,1, Queue{UInt64}(), 8)

println(state_matrix)
println(rep.curEnergy)
evolve!(rep, hami, state_matrix)
println(state_matrix)
println(rep.curEnergy)