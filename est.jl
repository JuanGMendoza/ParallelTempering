module test

include("test.jl")

end

println(names(test, all=true, imported=true))