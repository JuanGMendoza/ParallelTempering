module tests
	include("tests.jl")
end

include("tools.jl")
include("ParallelTempering.jl")


function custom_cmp(x::String)
   number_idx = findfirst(isdigit, x)
   str, num = SubString(x, 1, number_idx-1), SubString(x, number_idx, length(x))
   return str, parse(Int, num)
end

#Needs to be updated for new data structures
function generate_required_files()

		h1 = Hamiltonian(1, [[[]],[[]],[[]],[[]]], [[],[], [], []])
		bonds = [[[2], [3], [4]] ,[[1], [3], [4] ] ,[[1], [2], [4] ] ,[[1], [2], [3]]]
		h2 = Hamiltonian(2, bonds, [[1, 1, 1], [1,1,1], [1,1,1], [1,1,1]])
		temperatures::Vector{Float64} = Vector{Float64}(1:10)
		println("\n")
		#print("Running simulations...[00%]\r")
		print("Running simulations[0%]  " * repeat(".", 50 ) *"\r")
		t = @elapsed begin
			for i in (1:10)

				
				parallel_tempering(h1, temperatures, "./unit_tests/unit-test1_" * string(i), magnetization, "magnetization", UInt64(1000))
				parallel_tempering(h2, temperatures, "./unit_tests/unit-test2_" * string(i), magnetization, "magnetization", UInt64(1000))
				print("Running simulations[" *string(i) *"0%]  "* repeat("#", i*5) * repeat(".", 50 - i*5) *"\r")
				
				
			end
		end

		println("\ntime = ", t, " s")

end


function run_tests()
	totalTests = 0
	passedTests = 0

	functions = filter( x->getproperty(tests,x) isa Function && occursin("test", String(x)), names( tests, all = true ))
	functions = sort(string.(functions), by=custom_cmp)
	generate_required_files()

	for test in functions


		outcome, feedback = tests.eval(Meta.parse(String(test)))()

		if outcome
			passedTests += 1
			print("PASSED | ")
		else
			print("FAILED | ")
		end
		print(test, feedback * "\n\n")
		totalTests += 1


	end

	println(passedTests, "/", totalTests, " tests passed")
end

run_tests()