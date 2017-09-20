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

using DataFrames, DataFramesMeta, Match, PyCall, PyPlot, Modoplots

@pyimport seaborn as sns
sns.set()
color_palette="Set2"
sns.set(context="paper", style="darkgrid", palette=color_palette, font_scale=4)

#####################################################################
# Collecting the Experimental Results in one place                  #
#####################################################################

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
[rm("../results/$(files[i])") for i in 1:length(files)]
writetable("../results/Experimental_Results.csv", data)

#####################################################################
# Removing MOMBP Instances with less than 80 variables              #
#####################################################################

data = readtable("../results/Experimental_Results.csv")
inds = Int64[]
for i in 1:nrow(data)
	if data[:Problem_Type][i] == "MOMBP" && data[:Variables][i] in [20, 40, 1280] 
	else
		push!(inds, i)
	end
end
data = data[inds, :]
sort!(data)
writetable("../results/Experimental_Results.csv", data)

#####################################################################
# Generating Tables                                                 #
#####################################################################

with_instance_id = unique(data[[:Problem_Type, :Objectives, :Variables, :Constraints, :Instance_ID]])
without_instance_id = unique(data[[:Problem_Type, :Objectives, :Variables, :Constraints]])
num = Int64[]
for i in 1:nrow(without_instance_id)
    tmp = with_instance_id[with_instance_id[:Problem_Type] .== without_instance_id[:Problem_Type][i], :]
    tmp = tmp[tmp[:Objectives] .== without_instance_id[:Objectives][i], :]
    tmp = tmp[tmp[:Variables] .== without_instance_id[:Variables][i], :]
    tmp = tmp[tmp[:Constraints] .== without_instance_id[:Constraints][i], :]
    push!(num, nrow(tmp))
end
without_instance_id[:Instances] = num
cont_vars = zeros(Int64,nrow(without_instance_id))
bin_vars = zeros(Int64, nrow(without_instance_id))
for i in 1:nrow(without_instance_id)
    if without_instance_id[:Problem_Type][i] in ["BOMBLP", "BOUFLP", "MOMBP"]
        if without_instance_id[:Problem_Type][i] in ["BOMBLP", "MOMBP"]
            cont_vars[i] = bin_vars[i] = without_instance_id[:Variables][i]/2
        else
            cont_vars[i], bin_vars[i] = @match without_instance_id[:Variables][i] begin
                816 => 800, 16
                1275 => 1250, 25
                2550 => 2500, 50
            end
        end
    else
        bin_vars[i] = without_instance_id[:Variables][i]
    end
end
without_instance_id[:Continuous_Variables] = cont_vars
without_instance_id[:Binary_Variables] = bin_vars
without_instance_id = without_instance_id[[:Problem_Type, :Objectives, :Continuous_Variables, :Binary_Variables, :Constraints, :Instances]]
writetable("../results/Instances_Table.csv", without_instance_id)

data_1 = readtable("../results/Experimental_Results.csv")
data_1 = @where(data_1, :Algorithm .== "V5", :Total_Threads .== 1)
data_2 = readtable("../results/Instances_Table.csv")
lp_solvers = ["Clp", "GLPK", "SCIP", "Gurobi", "CPLEX", "CPLEXEfficient"]
count = zeros(Int64, nrow(data_2), length(lp_solvers))
for i in 1:nrow(data_2)
	tmp = data_1[data_1[:Problem_Type] .== data_2[:Problem_Type][i], :]
	tmp = tmp[tmp[:Objectives] .== data_2[:Objectives][i], :]
    tmp = tmp[tmp[:Variables] .== data_2[:Continuous_Variables][i] + data_2[:Binary_Variables][i], :]
    tmp = tmp[tmp[:Constraints] .== data_2[:Constraints][i], :]
	for j in 1:length(lp_solvers)
		tmp2 = tmp[tmp[:LP_Solver] .== lp_solvers[j], :]
		count[i, j] = nrow(tmp2)
	end
end
for i in 1:length(lp_solvers)
    data_2[Symbol(lp_solvers[i])] = count[:, i]
end
writetable("../results/Instances_Table_2.csv", data_2)

#####################################################################
# Summarizing Results for Biobjective Mixed Integer &               #
# Uncapacitated Facility Location Problems                          #
# Single Thread                                                     #
#####################################################################

data = readtable("../results/Experimental_Results.csv")
data = @where(data, :Instance_Type .== "Hadi", :LP_Solver .== "CPLEXEfficient")
sort!(data)
for i in 1:nrow(data)
	if data[:Algorithm][i] == "V5"
		data[:Algorithm][i] = "V5 T$(data[:Total_Threads][i])"
	end
end

f = plt[:figure](figsize=(16, 12))
yscale("log")
sns.boxplot(x=data[:Algorithm], y=data[:Num_Non_Dom_Pts_Found], palette=color_palette)
legend()
xlabel("Algorithms")
ylabel("Local nondominated points")
plt[:tight_layout]()
f[:savefig]("../plots/BOMIP_Local_Non_Dom_Pts.eps")

f = plt[:figure](figsize=(16, 12))
yscale("log")
sns.boxplot(x=data[:Algorithm], y=data[:Hypervolume_Gap], palette=color_palette)
legend()
xlabel("Algorithms")
ylabel("Hypervolume gap")
plt[:tight_layout]()
f[:savefig]("../plots/BOMIP_Hypervolume_Gap.eps")

f = plt[:figure](figsize=(16, 12))
yscale("log")
sns.boxplot(x=data[:Algorithm], y=data[:Coverage], palette=color_palette)
legend()
xlabel("Algorithms")
ylabel("Coverage")
plt[:tight_layout]()
f[:savefig]("../plots/BOMIP_Coverage.eps")

f = plt[:figure](figsize=(16, 12))
yscale("log")
sns.boxplot(x=data[:Algorithm], y=data[:Uniformity], palette=color_palette)
legend()
xlabel("Algorithms")
ylabel("Uniformity")
plt[:tight_layout]()
f[:savefig]("../plots/BOMIP_Uniformity.eps")

f = plt[:figure](figsize=(16, 12))
yscale("log")
sns.boxplot(x=data[:Algorithm], y=data[:Time], palette=color_palette)
legend()
xlabel("Algorithms")
ylabel("Computing time (Sec.)")
plt[:tight_layout]()
f[:savefig]("../plots/BOMIP_Computing_Time.eps")

#####################################################################
# Summarizing Results for Sensitivity of using Different LP Solvers #
#####################################################################

data = readtable("../results/Experimental_Results.csv")
data = @where(data, :Instance_Type .== "Hadi", :Algorithm .== "V5", :Total_Threads .== 1)
sort!(data)
solvers = ["Clp", "GLPK", "SCIP", "Gurobi", "CPLEX", "CPLEXEfficient"]
solvers_ = ["Clp", L"GLPK^*", L"SCIP^*", "Gurobi", "CPLEX", "CPLEXEfficient"]

f = plt[:figure](figsize=(16, 12))
yscale("log")
sns.boxplot(x=data[:LP_Solver], y=data[:Num_Non_Dom_Pts_Found], palette=color_palette, order=solvers)
legend()
xticks([0:5...], solvers_, rotation="vertical")
xlabel("LP Solver")
ylabel("Local nondominated points")
plt[:tight_layout]()
f[:savefig]("../plots/BOMIP_LP_Solver_Local_Non_Dom_Pts.eps")

f = plt[:figure](figsize=(16, 12))
yscale("log")
sns.boxplot(x=data[:LP_Solver], y=data[:Hypervolume_Gap], palette=color_palette, order=solvers)
legend()
xticks([0:5...], solvers_, rotation="vertical")
xlabel("LP Solver")
ylabel("Hypervolume gap")
plt[:tight_layout]()
f[:savefig]("../plots/BOMIP_LP_Solver_Hypervolume_Gap.eps")

f = plt[:figure](figsize=(16, 12))
yscale("log")
sns.boxplot(x=data[:LP_Solver], y=data[:Coverage], palette=color_palette, order=solvers)
legend()
xticks([0:5...], solvers_, rotation="vertical")
xlabel("LP Solver")
ylabel("Coverage")
plt[:tight_layout]()
f[:savefig]("../plots/BOMIP_LP_Solver_Coverage.eps")

f = plt[:figure](figsize=(16, 12))
yscale("log")
sns.boxplot(x=data[:LP_Solver], y=data[:Uniformity], palette=color_palette, order=solvers)
legend()
xticks([0:5...], solvers_, rotation="vertical")
xlabel("LP Solver")
ylabel("Uniformity")
plt[:tight_layout]()
f[:savefig]("../plots/BOMIP_LP_Solver_Uniformity.eps")

f = plt[:figure](figsize=(16, 12))
yscale("log")
sns.boxplot(x=data[:LP_Solver], y=data[:Time], palette=color_palette, order=solvers)
legend()
xticks([0:5...], solvers_, rotation="vertical")
xlabel("LP Solver")
ylabel("Computing time (Sec.)")
plt[:tight_layout]()
f[:savefig]("../plots/BOMIP_LP_Solver_Time.eps")

#####################################################################
# Approximate Frontier of a Biobjective Mixed Integer Program       #
#####################################################################

num = 8
instance_id = "BOMBLP_Hadi_$(num)"
non_dom_sols = [readdlm("../results/non_dom_sols/True_Frontier/$(instance_id)_True_Frontier.txt")]
push!(non_dom_sols, readdlm("../results/non_dom_sols/V1/$(instance_id)_CPLEXEfficient.txt"))
push!(non_dom_sols, readdlm("../results/non_dom_sols/V2/$(instance_id)_CPLEXEfficient.txt"))
push!(non_dom_sols, readdlm("../results/non_dom_sols/V3/$(instance_id)_CPLEXEfficient.txt"))
push!(non_dom_sols, readdlm("../results/non_dom_sols/V4/$(instance_id)_CPLEXEfficient.txt"))
push!(non_dom_sols, readdlm("../results/non_dom_sols/V5 T1/$(instance_id)_CPLEXEfficient.txt"))
algorithms = ["True Frontier", "V1", "V2", "V3", "V4", "V5"]
plt_non_dom_frntr_bomip(non_dom_sols[[2, 3, 4]], algorithms[[2, 3, 4]], "../plots/$(instance_id)_Frontier_V1_V2_V3.eps")
plt_non_dom_frntr_bomip(non_dom_sols[[4, 5]], algorithms[[4, 5]], "../plots/$(instance_id)_Frontier_V3_V4.eps")
plt_non_dom_frntr_bomip(non_dom_sols[[5, 6]], algorithms[[5, 6]], "../plots/$(instance_id)_Frontier_V4_V5.eps")
plt_non_dom_frntr_bomip(non_dom_sols[[1, 6]], algorithms[[1, 6]], "../plots/$(instance_id)_Frontier_True_Frontier_V5.eps")

#####################################################################
# Summarizing Results for Multiobjective Pure Integer Programs      #
#####################################################################

data = readtable("../results/Experimental_Results.csv")
data = @where(data, :Instance_Type .== "Kirlik", :LP_Solver .== "CPLEXEfficient")
sort!(data)
for i in 1:nrow(data)
	if data[:Algorithm][i] == "V5"
		data[:Algorithm][i] = "V5 T$(data[:Total_Threads][i])"
	end
end

f = plt[:figure](figsize=(16, 12))
yscale("log")
sns.boxplot(x=data[:Algorithm], y=data[:Num_Non_Dom_Pts_Found], palette=color_palette)
legend()
xlabel("Algorithms")
ylabel("Local nondominated points")
plt[:tight_layout]()
f[:savefig]("../plots/MOBP_Local_Non_Dom_Pts.eps")

f = plt[:figure](figsize=(16, 12))
yscale("log")
sns.boxplot(x=data[:Algorithm], y=data[:Hypervolume_Gap], palette=color_palette)
legend()
xlabel("Algorithms")
ylabel("Hypervolume gap")
plt[:tight_layout]()
f[:savefig]("../plots/MOBP_Hypervolume_Gap.eps")

f = plt[:figure](figsize=(16, 12))
yscale("log")
sns.boxplot(x=data[:Algorithm], y=data[:Cardinality], palette=color_palette)
legend()
xlabel("Algorithms")
ylabel("Cardinality")
plt[:tight_layout]()
f[:savefig]("../plots/MOBP_Cardinality.eps")

f = plt[:figure](figsize=(16, 12))
yscale("log")
sns.boxplot(x=data[:Algorithm], y=data[:Coverage], palette=color_palette)
legend()
xlabel("Algorithms")
ylabel("Coverage")
plt[:tight_layout]()
f[:savefig]("../plots/MOBP_Coverage.eps")

f = plt[:figure](figsize=(16, 12))
yscale("log")
sns.boxplot(x=data[:Algorithm], y=data[:Uniformity], palette=color_palette)
legend()
xlabel("Algorithms")
ylabel("Uniformity")
plt[:tight_layout]()
f[:savefig]("../plots/MOBP_Uniformity.eps")

f = plt[:figure](figsize=(16, 12))
yscale("log")
sns.boxplot(x=data[:Algorithm], y=data[:Time], palette=color_palette)
legend()
xlabel("Algorithms")
ylabel("Computing time (Sec.)")
plt[:tight_layout]()
f[:savefig]("../plots/MOBP_Computing_Time.eps")

#####################################################################
# Summarizing Results for MDLS on Multiobjective Knapsack Problems  #
#####################################################################

data = readtable("../results/Experimental_Results.csv")
data = @where(data, :Problem_Type .== "MOKP", :LP_Solver .!= "Clp", :Total_Threads .== 1)
data = @where(data, :LP_Solver .!= "GLPK")
data = @where(data, :LP_Solver .!= "SCIP")
data = @where(data, :LP_Solver .!= "Gurobi")
data = @where(data, :LP_Solver .!= "CPLEX")
sort!(data)
for i in 1:nrow(data)
	if data[:Algorithm][i] == "V5"
		data[:Algorithm][i] = "V5 T$(data[:Total_Threads][i])"
	end
end

f = plt[:figure](figsize=(16, 12))
yscale("log")
sns.boxplot(x=data[:Algorithm], y=data[:Num_Non_Dom_Pts_Found], palette=color_palette)
legend()
xlabel("Algorithms")
ylabel("Local nondominated points")
plt[:tight_layout]()
f[:savefig]("../plots/MOKP_Local_Non_Dom_Pts.eps")

f = plt[:figure](figsize=(16, 12))
yscale("log")
sns.boxplot(x=data[:Algorithm], y=data[:Hypervolume_Gap], palette=color_palette)
legend()
xlabel("Algorithms")
ylabel("Hypervolume gap")
plt[:tight_layout]()
f[:savefig]("../plots/MOKP_Hypervolume_Gap.eps")

f = plt[:figure](figsize=(16, 12))
yscale("log")
sns.boxplot(x=data[:Algorithm], y=data[:Cardinality], palette=color_palette)
legend()
xlabel("Algorithms")
ylabel("Cardinality")
plt[:tight_layout]()
f[:savefig]("../plots/MOKP_Cardinality.eps")

f = plt[:figure](figsize=(16, 12))
yscale("log")
sns.boxplot(x=data[:Algorithm], y=data[:Coverage], palette=color_palette)
legend()
xlabel("Algorithms")
ylabel("Coverage")
plt[:tight_layout]()
f[:savefig]("../plots/MOKP_Coverage.eps")

f = plt[:figure](figsize=(16, 12))
yscale("log")
sns.boxplot(x=data[:Algorithm], y=data[:Uniformity], palette=color_palette)
legend()
xlabel("Algorithms")
ylabel("Uniformity")
plt[:tight_layout]()
f[:savefig]("../plots/MOKP_Uniformity.eps")

f = plt[:figure](figsize=(16, 12))
yscale("log")
sns.boxplot(x=data[:Algorithm], y=data[:Time], palette=color_palette)
legend()
xlabel("Algorithms")
ylabel("Computing time (Sec.)")
plt[:tight_layout]()
f[:savefig]("../plots/MOKP_Computing_Time.eps")

#####################################################################
# Summarizing Results for Sensitivity of using Different LP Solvers #
#####################################################################

data = readtable("../results/Experimental_Results.csv")
data = @where(data, :Instance_Type .== "Kirlik", :Algorithm .== "V5", :Total_Threads .== 1)
sort!(data)
solvers = ["Clp", "GLPK", "SCIP", "Gurobi", "CPLEX", "CPLEXEfficient"]

f = plt[:figure](figsize=(16, 12))
yscale("log")
sns.boxplot(x=data[:LP_Solver], y=data[:Num_Non_Dom_Pts_Found], palette=color_palette, order=solvers)
legend()
xticks(rotation="vertical")
xlabel("LP Solver")
ylabel("Local nondominated points")
plt[:tight_layout]()
f[:savefig]("../plots/MOBP_LP_Solver_Local_Non_Dom_Pts.eps")

f = plt[:figure](figsize=(16, 12))
yscale("log")
sns.boxplot(x=data[:LP_Solver], y=data[:Hypervolume_Gap], palette=color_palette, order=solvers)
legend()
xticks(rotation="vertical")
xlabel("LP Solver")
ylabel("Hypervolume gap")
plt[:tight_layout]()
f[:savefig]("../plots/MOBP_LP_Solver_Hypervolume_Gap.eps")

f = plt[:figure](figsize=(16, 12))
yscale("log")
sns.boxplot(x=data[:LP_Solver], y=data[:Cardinality], palette=color_palette, order=solvers)
legend()
xticks(rotation="vertical")
xlabel("LP Solver")
ylabel("Cardinality")
plt[:tight_layout]()
f[:savefig]("../plots/MOBP_LP_Solver_Cardinality.eps")

f = plt[:figure](figsize=(16, 12))
yscale("log")
sns.boxplot(x=data[:LP_Solver], y=data[:Coverage], palette=color_palette, order=solvers)
legend()
xticks(rotation="vertical")
xlabel("LP Solver")
ylabel("Coverage")
plt[:tight_layout]()
f[:savefig]("../plots/MOBP_LP_Solver_Coverage.eps")

f = plt[:figure](figsize=(16, 12))
yscale("log")
sns.boxplot(x=data[:LP_Solver], y=data[:Uniformity], palette=color_palette, order=solvers)
legend()
xticks(rotation="vertical")
xlabel("LP Solver")
ylabel("Uniformity")
plt[:tight_layout]()
f[:savefig]("../plots/MOBP_LP_Solver_Uniformity.eps")

f = plt[:figure](figsize=(16, 12))
yscale("log")
sns.boxplot(x=data[:LP_Solver], y=data[:Time], palette=color_palette, order=solvers)
legend()
xticks(rotation="vertical")
xlabel("LP Solver")
ylabel("Computing time (Sec.)")
plt[:tight_layout]()
f[:savefig]("../plots/MOBP_LP_Solver_Time.eps")

#####################################################################
# Approximate Frontier of a Multiobjective Knapsack Problem         #
#####################################################################

p, n, instance = 3, 10, 5
instance_id = "MOKP_Kirlik_$(p)_$(n)_$(instance)"
non_dom_sols = [readdlm("../results/non_dom_sols/True_Frontier/$(instance_id)_True_Frontier.txt")]
push!(non_dom_sols, readdlm("../results/non_dom_sols/V1/$(instance_id)_CPLEXEfficient.txt"))
push!(non_dom_sols, readdlm("../results/non_dom_sols/V2/$(instance_id)_CPLEXEfficient.txt"))
push!(non_dom_sols, readdlm("../results/non_dom_sols/V3/$(instance_id)_CPLEXEfficient.txt"))
push!(non_dom_sols, readdlm("../results/non_dom_sols/V4/$(instance_id)_CPLEXEfficient.txt"))
push!(non_dom_sols, readdlm("../results/non_dom_sols/V5 T1/$(instance_id)_CPLEXEfficient.txt"))
algorithms = ["True Frontier", "V1", "V2", "V3", "V4", "V5"]
plt_discrete_non_dom_frntr(non_dom_sols, algorithms, true, "../plots/$(instance_id)_Frontier")

#####################################################################
# Summarizing Results for Multiobjective Mixed Integer Programs     #
#####################################################################

data = readtable("../results/Experimental_Results.csv")
data = @where(data, :Problem_Type .== "MOMBP", :LP_Solver .== "CPLEXEfficient")
sort!(data)
for i in 1:nrow(data)
	if data[:Algorithm][i] == "V5"
		data[:Algorithm][i] = "V5 T$(data[:Total_Threads][i])"
	end
end

f = plt[:figure](figsize=(16, 12))
yscale("log")
sns.boxplot(x=data[:Algorithm], y=data[:Num_Non_Dom_Pts_Found], palette=color_palette)
legend()
xlabel("Algorithms")
ylabel("Local nondominated points")
plt[:tight_layout]()
f[:savefig]("../plots/MOMBP_Local_Non_Dom_Pts.eps")

f = plt[:figure](figsize=(16, 12))
yscale("log")
sns.boxplot(x=data[:Algorithm], y=data[:Hypervolume_Gap], palette=color_palette)
legend()
xlabel("Algorithms")
ylabel("Hypervolume gap")
plt[:tight_layout]()
f[:savefig]("../plots/MOMBP_Hypervolume_Gap.eps")

f = plt[:figure](figsize=(16, 12))
yscale("log")
sns.boxplot(x=data[:Algorithm], y=data[:Cardinality], palette=color_palette)
legend()
xlabel("Algorithms")
ylabel("Cardinality")
plt[:tight_layout]()
f[:savefig]("../plots/MOMBP_Cardinality.eps")

f = plt[:figure](figsize=(16, 12))
yscale("log")
sns.boxplot(x=data[:Algorithm], y=data[:Coverage], palette=color_palette)
legend()
xlabel("Algorithms")
ylabel("Coverage")
plt[:tight_layout]()
f[:savefig]("../plots/MOMBP_Coverage.eps")

f = plt[:figure](figsize=(16, 12))
yscale("log")
sns.boxplot(x=data[:Algorithm], y=data[:Uniformity], palette=color_palette)
legend()
xlabel("Algorithms")
ylabel("Uniformity")
plt[:tight_layout]()
f[:savefig]("../plots/MOMBP_Uniformity.eps")

f = plt[:figure](figsize=(16, 12))
yscale("log")
sns.boxplot(x=data[:Algorithm], y=data[:Time], palette=color_palette)
legend()
xlabel("Algorithms")
ylabel("Computing time (Sec.)")
plt[:tight_layout]()
f[:savefig]("../plots/MOMBP_Computing_Time.eps")

#####################################################################
# Summarizing Results for Sensitivity of using Different LP Solvers #
#####################################################################

data = readtable("../results/Experimental_Results.csv")
data = @where(data, :Problem_Type .!= "MOMBP", :Algorithm .== "V5", :Total_Threads .== 1)
sort!(data)
solvers = ["Clp", "GLPK", "SCIP", "Gurobi", "CPLEX", "CPLEXEfficient"]

f = plt[:figure](figsize=(16, 12))
yscale("log")
sns.boxplot(x=data[:LP_Solver], y=data[:Num_Non_Dom_Pts_Found], palette=color_palette, order=solvers)
legend()
xticks(rotation="vertical")
xlabel("LP Solver")
ylabel("Local nondominated points")
plt[:tight_layout]()
f[:savefig]("../plots/MOMBP_LP_Solver_Local_Non_Dom_Pts.eps")

f = plt[:figure](figsize=(16, 12))
yscale("log")
sns.boxplot(x=data[:LP_Solver], y=data[:Hypervolume_Gap], palette=color_palette, order=solvers)
legend()
xticks(rotation="vertical")
xlabel("LP Solver")
ylabel("Hypervolume gap")
plt[:tight_layout]()
f[:savefig]("../plots/MOMBP_LP_Solver_Hypervolume_Gap.eps")

f = plt[:figure](figsize=(16, 12))
yscale("log")
sns.boxplot(x=data[:LP_Solver], y=data[:Cardinality], palette=color_palette, order=solvers)
legend()
xticks(rotation="vertical")
xlabel("LP Solver")
ylabel("Cardinality")
plt[:tight_layout]()
f[:savefig]("../plots/MOMBP_LP_Solver_Cardinality.eps")

f = plt[:figure](figsize=(16, 12))
yscale("log")
sns.boxplot(x=data[:LP_Solver], y=data[:Coverage], palette=color_palette, order=solvers)
legend()
xticks(rotation="vertical")
xlabel("LP Solver")
ylabel("Coverage")
plt[:tight_layout]()
f[:savefig]("../plots/MOMBP_LP_Solver_Coverage.eps")

f = plt[:figure](figsize=(16, 12))
yscale("log")
sns.boxplot(x=data[:LP_Solver], y=data[:Uniformity], palette=color_palette, order=solvers)
legend()
xticks(rotation="vertical")
xlabel("LP Solver")
ylabel("Uniformity")
plt[:tight_layout]()
f[:savefig]("../plots/MOMBP_LP_Solver_Uniformity.eps")

f = plt[:figure](figsize=(16, 12))
yscale("log")
sns.boxplot(x=data[:LP_Solver], y=data[:Time], palette=color_palette, order=solvers)
legend()
xticks(rotation="vertical")
xlabel("LP Solver")
ylabel("Computing time (Sec.)")
plt[:tight_layout]()
f[:savefig]("../plots/MOMBP_LP_Solver_Time.eps")
