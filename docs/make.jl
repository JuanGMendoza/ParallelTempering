nside make.jl
push!(LOAD_PATH,"../src/")
using ParallelTempering
using Documenter
makedocs(
	          sitename = "ParallelTempering.jl",
		           modules  = [VegaGraphs],
			            pages=[
					                   "Home" => "index.md"
							                  ])
deploydocs(;
	       repo="github.com/JuanGMendoza/ParallelTempering.jl",
	       )
