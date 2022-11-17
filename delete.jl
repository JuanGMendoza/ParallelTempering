using Bits

size = 4

state_matrix::Vector{UInt128} = Vector{UInt128}(undef, size)

state_matrix[1] = 2

println(length(bits(state_matrix[1])))