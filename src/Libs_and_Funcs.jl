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

using Modof, Modolib, JuMP, FPBH, FPBHCPLEX, GLPKMathProgInterface, Clp, SCIP, Gurobi, CPLEX, Match, DataFrames, DataFramesMeta

try
	model = ModoModel()
	@variable(model, x[1:4], Bin)
	objective!(model, 1, :Min, x[1]+x[2]+x[3]+x[4])
	objective!(model, 2, :Min, -x[1]-x[2]-x[3]-x[4])
	@constraint(model, x[1]+x[2]+x[3]+x[4] <= 3)

	@time solution = fpbhcplex(model)
	@time solution = fpbh(model, lp_solver=GLPKSolverLP())
	@time solution = fpbh(model, lp_solver=SCIPSolver())
	@time solution = fpbh(model, lp_solver=CplexSolver())

	model = ModoModel()
	@variable(model, x[1:4], Bin)
	objective!(model, 1, :Min, x[1]+x[2]+x[3]+x[4])
	objective!(model, 2, :Min, -x[1]-x[2]-x[3]-x[4])
	objective!(model, 3, :Min, x[1]+2x[2]+3x[3]+4x[4])
	@constraint(model, x[1]+x[2]+x[3]+x[4] <= 3)

	@time solution = fpbhcplex(model)
	@time solution = fpbh(model, lp_solver=GLPKSolverLP())
	@time solution = fpbh(model, lp_solver=ClpSolver())
	@time solution = fpbh(model, lp_solver=GurobiSolver())

	model = ModoModel()
	@variable(model, x[1:2], Bin)
	@variable(model, y[1:2] >= 0.0)
	objective!(model, 1, :Min, x[1] + x[2] + y[1] + y[2])
	objective!(model, 2, :Min, -x[1] - x[2] - y[1] - y[2])
	@constraint(model, x[1] + x[2] <= 1) 
	@constraint(model, 2.72y[1] + 7.39y[2] >= 1) 

	@time solution = fpbhcplex(model)
	@time solution = fpbh(model, lp_solver=GLPKSolverLP())

	model = ModoModel()
	@variable(model, x[1:2], Bin)
	@variable(model, y[1:2] >= 0.0)
	objective!(model, 1, :Min, x[1] + x[2] + y[1] + y[2])
	objective!(model, 2, :Min, -x[1] - x[2] - y[1] - y[2])
	objective!(model, 3, :Min, x[1] + 2x[2] + y[1] + 2y[2])
	@constraint(model, x[1] + x[2] <= 1) 
	@constraint(model, 2.72y[1] + 7.39y[2] >= 1)
 
	@time solution = fpbhcplex(model)
	@time solution = fpbh(model, lp_solver=GLPKSolverLP())
catch
end

function solve_instance(params)
	t0::Float64, t1::Float64 = 0.0, 0.0
	if params[:problem_type] in ["BOMBLP", "BOUFLP"]
		instance, true_non_dom_sols = eval(params[:read_instance])(params[:read_instance_input])
		println("Started Solving Instance = $(params[:read_instance_input])")	
	else
		instance = @match params[:problem_type] begin
			"MOAP" || "MOKP" => MOBPInstance()
			"MOMBP" => MOMBLPInstance()
		end
		true_non_dom_sols::Array{Float64, 2} = zeros(1, 1)
		try
			instance, true_non_dom_sols = eval(params[:read_instance])(params[:read_instance_input]...)
			println("Started Solving Instance = $(params[:read_instance_input]...)")
		catch
			return
		end
	end	
	if params[:algorithm] == "MDLS" && params[:problem_type] == "MOKP"
		t0 = time()
		non_dom_sols = mdls_kp(instance)
		t1 = time()
	else
		if params[:lpsolver] == "CPLEXEfficient"
			if params[:algorithm] != "V5"
				t0 = time()
				non_dom_sols = @match params[:algorithm] begin
					"V1" => fpbhcplex(instance, obj_fph=false, local_search=false, decomposition=false, solution_polishing=false)
					"V2" => fpbhcplex(instance, local_search=false, decomposition=false, solution_polishing=false)
					"V3" => fpbhcplex(instance, decomposition=false, solution_polishing=false)
					"V4" => fpbhcplex(instance, solution_polishing=false)
				end 
				t1 = time()
			else
				t0 = time()
				non_dom_sols = @match params[:total_threads] begin
					1 => fpbhcplex(instance)
					_ => fpbhcplex(instance, threads=params[:total_threads], parallelism=true)
				end 
				t1 = time()
			end
			println("Finished")
		else
			t0 = time()
			non_dom_sols = @match params[:lpsolver] begin
				"GLPK" => fpbh(instance, lp_solver=GLPKSolverLP(tol_bnd=1.0e-9))
				"Clp" => fpbh(instance, lp_solver=ClpSolver(PrimalTolerance=1.0e-9))
				"SCIP" => fpbh(instance, lp_solver=SCIPSolver("lp/threads", 1, "display/lpinfo", false, "numerics/lpfeastol", 1e-9))
				"Gurobi" => fpbh(instance, lp_solver=GurobiSolver(Threads=1, OutputFlag=0, FeasibilityTol=1.0e-9))
				"CPLEX" => fpbh(instance, lp_solver=CplexSolver(CPX_PARAM_THREADS=1, CPX_PARAM_SCRIND=0, CPX_PARAM_EPRHS=1.0e-9))
			end
			t1 = time()
		end
	end
	
	non_dom_sols = check_feasibility(non_dom_sols, instance)
	if params[:algorithm] != "V5"
		if params[:algorithm] == "MDLS"
			write_nondominated_frontier(non_dom_sols, "../results/non_dom_sols/$(params[:algorithm])/$(params[:problem_type])_$(params[:instance_type])_$(params[:instance_id]).txt")
		else
			write_nondominated_frontier(non_dom_sols, "../results/non_dom_sols/$(params[:algorithm])/$(params[:problem_type])_$(params[:instance_type])_$(params[:instance_id])_$(params[:lpsolver]).txt")
		end
	else
		write_nondominated_frontier(non_dom_sols, "../results/non_dom_sols/$(params[:algorithm]) T$(params[:total_threads])/$(params[:problem_type])_$(params[:instance_type])_$(params[:instance_id])_$(params[:lpsolver]).txt")
	end
	writedlm("../results/non_dom_sols/True_Frontier/$(params[:problem_type])_$(params[:instance_type])_$(params[:instance_id])_True_Frontier.txt", true_non_dom_sols)
	
	if params[:problem_type] in ["BOMBLP", "BOUFLP"]
		hypervolume_gap, cardinality, max_coverage, avg_coverage, uniformity = compute_quality_of_norm_apprx_frontier(wrap_sols_into_array(non_dom_sols), true_non_dom_sols, true)
	else
		hypervolume_gap, cardinality, max_coverage, avg_coverage, uniformity = compute_quality_of_norm_apprx_frontier(wrap_sols_into_array(non_dom_sols), true_non_dom_sols)
	end	
	println("Algorithm = $(params[:algorithm]), LP Solver = $(params[:lpsolver]), Hypervolume Gap = $hypervolume_gap Cardinality = $cardinality Average Coverage = $avg_coverage Uniformity = $uniformity")
	
	if params[:total_threads] == 1
		filename = "../results/Experimental_Results_$(myid()).csv"
	else
		filename = "../results/Experimental_Results_$(params[:total_threads]).csv"
	end
	obj::Int64 = 2
	if typeof(instance) == MOBPInstance || typeof(instance) == MOMBLPInstance
		obj = size(instance.c)[1]
	end
	if isfile(filename) == false
		results::DataFrame = DataFrame(Problem_Type=params[:problem_type], Instance_Type=params[:instance_type], Objectives=obj, Variables=size(instance.A)[2], Constraints=size(instance.A)[1], Instance_ID=string(params[:instance_id]), Algorithm=params[:algorithm], LP_Solver=params[:lpsolver], Total_Threads=params[:total_threads], Num_Non_Dom_Pts_Found=length(non_dom_sols), Hypervolume_Gap=hypervolume_gap, Cardinality=cardinality, Coverage=avg_coverage, Uniformity=uniformity, Time=round(t1-t0, 7))
	else
		results = readtable(filename)
		results[:Instance_ID] = string.(results[:Instance_ID])
		push!(results, [params[:problem_type], params[:instance_type], obj, size(instance.A)[2], size(instance.A)[1], string(params[:instance_id]), params[:algorithm], params[:lpsolver], params[:total_threads], length(non_dom_sols), hypervolume_gap, cardinality, avg_coverage, uniformity, round(t1-t0, 7)])
	end
	writetable(filename, results)
end

function solve_instance(params, inds::Vector{Int64})
	for i in inds
		try
			solve_instance(params[i])
		catch
		end
	end
end

#####################################################################
# Generating the Queue for the Experimental Run                     #
#####################################################################

function generate_queue()
	queue = Dict{Any, Any}[]
	
	algorithms = ["V1", "V2", "V3", "V4", "V5"]
	lpsolvers = ["GLPK", "Clp", "SCIP", "Gurobi", "CPLEX", "CPLEXEfficient"]
	threads = [1:4...]
	
	for i in 1:25, algorithm in algorithms, lpsolver_ in lpsolvers, thread in threads
		if algorithm != "V5" && lpsolver_ in ["GLPK", "Clp", "SCIP", "Gurobi", "CPLEX"]
			continue
		end
		if algorithm != "V5" && thread >= 2
			continue
		end
		if lpsolver_ in ["GLPK", "Clp", "SCIP", "Gurobi", "CPLEX"] && thread >= 2
			continue
		end
		params = Dict()
    	params[:read_instance] = :read_bomip_hadi
		params[:read_instance_input] = i
		params[:problem_type] = "BOMBLP"
		params[:instance_type] = "Hadi"
		params[:instance_id] = i
		params[:algorithm] = algorithm
		params[:lpsolver] = lpsolver_
		params[:total_threads] = thread
		push!(queue, params)
	end
	
	algorithms = ["V1", "V2", "V3", "V4", "V5"]
	lpsolvers = ["Clp", "Gurobi", "CPLEX", "CPLEXEfficient"]
	threads = [1:4...]
	
	for i in 1:12, algorithm in algorithms, lpsolver_ in lpsolvers, thread in threads
		if algorithm != "V5" && lpsolver_ in ["GLPK", "Clp", "SCIP", "Gurobi", "CPLEX"]
			continue
		end
		if algorithm != "V5" && thread >= 2
			continue
		end
		if lpsolver_ in ["GLPK", "Clp", "SCIP", "Gurobi", "CPLEX"] && thread >= 2
			continue
		end
    	params = Dict()
    	params[:read_instance] = :read_bouflp_hadi
		params[:read_instance_input] = i
		params[:problem_type] = "BOUFLP"
		params[:instance_type] = "Hadi"
		params[:instance_id] = i
		params[:algorithm] = algorithm
		params[:lpsolver] = lpsolver_
		params[:total_threads] = thread
		push!(queue, params)
	end
	
	algorithms = ["MDLS", "V1", "V2", "V3", "V4", "V5"]
	lpsolvers = ["None", "GLPK", "Clp", "SCIP", "Gurobi", "CPLEX", "CPLEXEfficient"]
	threads = [1:4...]
	
	for i in 3:5, j in [10:10:100...], k in 1:10, algorithm in algorithms, lpsolver_ in lpsolvers, thread in threads
		if (algorithm == "MDLS" && lpsolver_ != "None") || (algorithm != "MDLS" && lpsolver_ == "None")
			continue
		end
		if algorithm != "V5" && lpsolver_ in ["GLPK", "Clp", "SCIP", "Gurobi", "CPLEX"]
			continue
		end
		if algorithm != "V5" && thread >= 2
			continue
		end
		if lpsolver_ in ["GLPK", "Clp", "SCIP", "Gurobi", "CPLEX"] && thread >= 2
			continue
		end
		if i == 4 && j >= 50
			continue
		end
		if i == 5 && j >= 30
			continue
		end
    	params = Dict()
    	params[:read_instance] = :read_mokp_kirlik
		params[:read_instance_input] = [i, j, k]
		params[:problem_type] = "MOKP"
		params[:instance_type] = "Kirlik"
		params[:instance_id] = "$(i)_$(j)_$(k)"
		params[:algorithm] = algorithm
		params[:lpsolver] = lpsolver_
		params[:total_threads] = thread
		push!(queue, params)
	end
	
	algorithms = ["V1", "V2", "V3", "V4", "V5"]
	lpsolvers = ["GLPK", "Clp", "SCIP", "Gurobi", "CPLEX", "CPLEXEfficient"]
	threads = [1:4...]
	
    for i in 3, j in [5:5:50...], k in 1:10, algorithm in algorithms, lpsolver_ in lpsolvers, thread in threads
		if algorithm != "V5" && lpsolver_ in ["GLPK", "Clp", "SCIP", "Gurobi", "CPLEX"]
			continue
		end
		if algorithm != "V5" && thread >= 2
			continue
		end
		if lpsolver_ in ["GLPK", "Clp", "SCIP", "Gurobi", "CPLEX"] && thread >= 2
			continue
		end
    	params = Dict()
    	params[:read_instance] = :read_moap_kirlik
		params[:read_instance_input] = [i, j, k]
		params[:problem_type] = "MOAP"
		params[:instance_type] = "Kirlik"
		params[:instance_id] = "$(i)_$(j)_$(k)"
		params[:algorithm] = algorithm
		params[:lpsolver] = lpsolver_
		params[:total_threads] = thread
		push!(queue, params)
	end
	
	algorithms = ["V1", "V2", "V3", "V4", "V5"]
	lpsolvers = ["GLPK", "Clp", "SCIP", "Gurobi", "CPLEX", "CPLEXEfficient"]
	threads = [1:4...]
	
    for p in 3:5, c in [80, 160, 320, 640, 1280], i in 1:5, algorithm in algorithms, lpsolver_ in lpsolvers, thread in threads
		if algorithm != "V5" && lpsolver_ in ["GLPK", "Clp", "SCIP", "Gurobi", "CPLEX"]
			continue
		end
		if algorithm != "V5" && thread >= 2
			continue
		end
		if lpsolver_ in ["GLPK", "Clp", "SCIP", "Gurobi", "CPLEX"] && thread >= 2
			continue
		end
    	params = Dict()
    	params[:read_instance] = :read_mombp_aritra
		params[:read_instance_input] = [p, c, i]
		params[:problem_type] = "MOMBP"
		params[:instance_type] = "Aritra"
		params[:instance_id] = "$(p)_$(c)_$(i)"
		params[:algorithm] = algorithm
		params[:lpsolver] = lpsolver_
		params[:total_threads] = thread
		push!(queue, params)
	end
	
	try
		files = readdir("../results/")
		files = files[[contains(files[i], "Experimental_Results") for i in 1:length(files)]]
		data = readtable("../results/$(files[1])")
		data[:Instance_ID] = string.(data[:Instance_ID])
		for i in 2:length(files)
    		tmp = readtable("../results/$(files[i])")
    		tmp[:Instance_ID] = string.(tmp[:Instance_ID])
    		data = vcat(data, tmp)
		end
		data = data[completecases(data), :]
		data = data[names(data)[1:15]]
		unique!(data)
		sort!(data)
		#[rm("../results/$(files[i])") for i in 1:length(files)]
		writetable("../results/Experimental_Results.csv", data)
		ind1 = Int64[]
		for i in 1:length(queue)
    		tmp = data[data[:Problem_Type] .== queue[i][:problem_type],:]
    		tmp = tmp[tmp[:Instance_Type] .== queue[i][:instance_type],:]
    		tmp = tmp[tmp[:Instance_ID] .== string(queue[i][:instance_id]),:]
    		tmp = tmp[tmp[:Algorithm] .== queue[i][:algorithm],:]
    		tmp = tmp[tmp[:LP_Solver] .== queue[i][:lpsolver],:]
    		tmp = tmp[tmp[:Total_Threads] .== queue[i][:total_threads],:]
    		if nrow(tmp) == 0
        		push!(ind1, i)
        		continue
    		end
		end
		return queue[ind1]
	catch
		return queue
	end
end

function generate_parallel_queues(inds::Vector{Int64})
	procs_::Vector{Int64} = setdiff(procs(), myid())
	inds_ = Dict()
	shuffle!(inds)
	num = ceil(Int64, length(inds)/length(procs_))
	count::Int64 = 1
    for p in procs_
    	try
    		inds_[p] = inds[num*(count-1)+1:num*(count)]
    	catch
    		inds_[p] = inds[num*(count-1)+1:end]
    	end
    	count += 1
    end
    inds_
end


