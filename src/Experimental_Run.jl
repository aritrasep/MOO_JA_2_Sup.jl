###############################################################################
#                                                                             #
#  This file is part of the julia module for Multi Objective Optimization     #
#  (c) Copyright 2017 by Aritra Pal                                           #
#                                                                             #
# Permission is hereby granted, free of charge, to any person obtaining a     # 
# copy of this software and associated documentation files (the "Software"),  #
# to deal in the Software without restriction, including without limitation   #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,    #
# and/or sell copies of the Software, and to permit persons to whom the       #
# Software is furnished to do so, subject to the following conditions:        #
#                                                                             #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
#                                                                             #
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                   #
#                                                                             #
# Every publication and presentation for which work based on the Program or   #
# its output has been used must contain an appropriate citation and           #
# acknowledgment of the author(s) of the Program.                             #
###############################################################################

@everywhere include("Libs_and_Funcs.jl")

@inbounds function generate_approximate_true_frontier()
	data = readtable("../results/Experimental_Results.csv")
	data = @where(data, :Problem_Type .== "MOMBP")
	unique_instances = unique(collect(data[:Instance_ID]))
	for i in 1:length(unique_instances)
		tmp = @where(data, :Instance_ID .== unique_instances[i])
		if nrow(tmp) >= 1
			true_frontier = readdlm("../results/non_dom_sols/True_Frontier/MOMBP_Aritra_$(unique_instances[i])_True_Frontier.txt")
			for j in 1:nrow(tmp)
				algorithm = tmp[:Algorithm][j]!="V5"?tmp[:Algorithm][j]:"$(tmp[:Algorithm][j]) T$(tmp[:Total_Threads][j])"
				tmp2 = readdlm("../results/non_dom_sols/$(algorithm)/MOMBP_Aritra_$(unique_instances[i])_$(tmp[:LP_Solver][j]).txt")
				true_frontier = select_non_dom_sols(vcat(true_frontier, tmp2))
			end
			writedlm("../results/non_dom_sols/True_Frontier/MOMBP_Aritra_$(unique_instances[i])_True_Frontier.txt", sort_non_dom_sols(true_frontier))
		end
	end
end

@inbounds function update_experimental_results()
	generate_approximate_true_frontier()
	data = readtable("../results/Experimental_Results.csv")
	for i in 1:nrow(data)
		if data[:Problem_Type][i] .== "MOMBP"
			println("Earlier")
			println(data[i, :])
			true_frontier = readdlm("../results/non_dom_sols/True_Frontier/MOMBP_Aritra_$(data[:Instance_ID][i])_True_Frontier.txt")
			algorithm = data[:Algorithm][i]!="V5"?data[:Algorithm][i]:"$(data[:Algorithm][i]) T$(data[:Total_Threads][i])"
			non_dom_sols = readdlm("../results/non_dom_sols/$(algorithm)/MOMBP_Aritra_$(data[:Instance_ID][i])_$(data[:LP_Solver][i]).txt")
			hypervolume_gap, cardinality, max_coverage, avg_coverage, uniformity = compute_quality_of_norm_apprx_frontier(non_dom_sols, true_frontier)
			data[:Hypervolume_Gap][i] = hypervolume_gap
			data[:Cardinality][i] = cardinality
			data[:Coverage][i] = avg_coverage
			data[:Uniformity][i] = uniformity
			println("Current")
			println(data[i, :])
		end
	end
	writetable("../results/Experimental_Results.csv", data)
end

@inbounds function run_experiment()
	for i in 1:10
		if length(procs()) <= 4
			addprocs(5 - length(procs()))
			@everywhere include("Libs_and_Funcs.jl")
		end
		queue = generate_queue()
		ind = Int64[]
		for i in 1:length(queue)
			if queue[i][:total_threads] == 1
				push!(ind, i)
			end
		end
		inds = generate_parallel_queues(ind)
		procs_::Vector{Int64} = procs()
		@sync begin
        	for p in procs_
        		if p != myid() || length(procs_) == 1
        		    @async begin
        		    	remotecall_fetch(solve_instance, p, queue, sort(inds[p]))
        		    end
        		end	
    		end
    	end
    end
    println("Finished Parallel Queue")
    queue = generate_queue()
	ind = Int64[]
	for i in 1:length(queue)
		if queue[i][:total_threads] > 1
			push!(ind, i)
		end
	end
    for i in ind
		try
			solve_instance(queue[i])
		catch
		end
	end
	println("Finished Serial Queue")
	update_experimental_results()
end

run_experiment()
