#using Plots
using Random


x = 1:10; y = rand(10); # These are the plotting data
#display(plot(x, y))

#readline()

system_sizes = Array{UInt8}(1:10)

println(system_sizes.^-1)
