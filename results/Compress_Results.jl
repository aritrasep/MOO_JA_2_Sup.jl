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

using DataFrames

files = readdir()
files = files[[contains(files[i], "Experimental_Results") for i in 1:length(files)]]
data = readtable("$(files[1])")
data[:Instance_ID] = string.(data[:Instance_ID])
for i in 2:length(files)
	tmp = readtable("$(files[i])")
	tmp[:Instance_ID] = string.(tmp[:Instance_ID])
	data = vcat(data, tmp)
end
data = data[completecases(data), :]
data = data[names(data)[1:15]]
unique!(data)
sort!(data)
writetable("Experimental_Results.csv", data)
