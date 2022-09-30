using JLD2

function save()

	name::String = "test.jld2"
	jldopen(name, "a+") do file


		file["t0/replica1/state"] = [0,0,0,0]
		file["t0/replica2/state"] = [1,1,0,0]
		file["t0/replica1/T"] = 10
		file["t0/replica2/T"] = 20

		
	end
end


save()