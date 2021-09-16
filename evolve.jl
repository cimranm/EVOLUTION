#! /usr/bin/env julia

using DelimitedFiles

using CSV

# TODO: 
# Use like MyModule.jl  
# have a script run.jl   that runs MyModule.main(args ...)  
# this can serve as the 'gray box' CLI tool 

# TODO: 
# perhaps use classes so you can do 
# e = EvolutionHandler(substitution_matrix, other, parameters ...)
# foreach sequence: evolve sequence 
#
# for AA in seq 
#   AA.evolve() 
# end 


# TODO have a way of using AA as index i.e. 'AA_F' or 'AA_C' but it's actually just an int (enum?)



# QUESTIONS:
# - can we have multiple files / modules?
# - what is last line of matrix for? 
# - can we have command line arg parsing? 
# - what libraries can I import for julia?

const GENERATIONS = 500


function getmatrix()

    filename = "matrix.txt"
    matrix = readdlm(filename, ',')

    # Get the first 'header' row of AA letters (will be used to assign index values)
    header = matrix[1, 2:end]

    # Strip matrix of the first col / first row (we only need the numbers)
    return matrix[2:end, 2:end]
end


# SUBSTITUTION MATRIX

matrix = """ 
,A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V
A,9867,2,9,10,3,8,17,21,2,6,4,2,6,2,22,35,32,0,2,18
R,1,9914,1,0,1,10,0,0,10,3,1,19,4,1,4,6,1,8,0,1
N,4,1,9822,36,0,4,6,6,21,3,1,13,0,1,2,20,9,1,4,1
D,6,0,42,9859,0,6,53,6,4,1,0,3,0,0,1,5,3,0,0,1
C,1,1,0,0,9973,0,0,0,1,1,0,0,0,0,1,5,1,0,3,2
Q,3,9,4,5,0,9876,27,1,23,1,3,6,4,0,6,2,2,0,0,1
E,10,0,7,56,0,35,9865,4,2,3,1,4,1,0,3,4,2,0,1,2
G,21,1,12,11,1,3,7,9935,1,0,1,2,1,1,3,21,3,0,0,5
H,1,8,18,3,1,20,1,0,9913,0,1,1,0,2,3,1,1,1,4,1
I,2,2,3,1,2,1,2,0,0,9871,9,2,12,7,0,1,7,0,1,33
L,3,1,3,0,0,6,1,1,4,22,9947,2,45,13,3,1,3,4,2,15
K,2,37,25,6,0,12,7,2,2,4,1,9924,20,0,3,8,11,0,1,1
M,1,1,0,0,0,2,0,0,0,5,8,4,9875,1,0,1,2,0,0,4
F,1,1,1,0,0,0,0,1,2,8,6,0,4,9944,0,2,1,3,28,0
P,13,5,2,1,1,8,3,2,5,1,2,2,1,1,9924,12,4,0,0,2
S,28,11,34,7,11,4,6,16,2,2,1,7,4,3,17,9840,38,5,2,2
T,22,2,13,4,1,3,2,2,1,11,2,8,6,1,5,32,9869,0,2,9
W,0,2,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,9976,1,0
Y,1,0,3,0,3,0,1,0,4,1,1,0,0,21,0,1,1,2,9947,1
V,13,2,1,1,3,2,2,3,3,57,11,1,17,1,3,2,10,0,2,9901
"""

# TODO enum? i.e. number 1 - 20 (exactly) that corresponds to a particular AA

# Map indexes to AAs (have a way to not hardcode this; but get it from the matrix CSV)
function getaminomap()

    d = Dict(
        "A" => 1
    )


    return d
end

# Get index from a given AA character 
function aminoindex(aa)

    #TODO
    return 1



end



# Randomly sample discrete variable
function sample(items, weights)

    # STEPS
    # get all weights.  'spread' weights over some distance (using cumsum)
    # iterate through all array indexes; find array index SUCH THAT:
    #     rand() <= weights[index]
    #
    # This will be the index that we use (as it means that our random variable 
    # has 'fallen' on this region of the distribution)

    # TODO this may not necessarily add up to 1.0 (so make sure we scale it so it does)
    distribution = cumsum(weights)
    random = rand()

    for index in 1:length(items)
        if random <= distribution[index]
            return items[index]
        end
    end

end 


# Input from stdin 
# check if FASTA format 
# TODO


# Benchmark

# TODO: able to give weights as relative numbers 
# i.e. function will automatically add them to total, and divide by total to give prob (so don't 
# the user doesn't have to make sure that they sum to 1.0 exactly ) 

# Make dict?      Dict(x->f(x) for ....)

function main()
    items = ["A", "B", "C", "D"]

    weights = [0.1, 0.2, 0.1, 0.6]



    # Initialise dict 
    count = Dict()
    for item in items
        count[item] = 0
    end 

    # Sample randomly

    total = 10^6

    for i in 1:total
        x = sample(items, weights)
        count[x] += 1

    end

    # Display 
    for item in items
        c = count[item] * 100 / total 
        println("$item:  $c%")
    end

    println()
    
    matrix = getmatrix()

    

end


main()
