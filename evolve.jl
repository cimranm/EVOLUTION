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


# TODO: maybe use dict that maps item=> weight ; instead of 2 arrays? 
function sample(items, weights)

    # STEPS
    # get all weights.  'spread' weights over some distance (using cumsum)
    # iterate through all array indexes; find array index SUCH THAT:
    #     rand() <= weights[index]
    #
    # This will be the index that we use (as it means that our random variable 
    # has 'fallen' on this region of the distribution)

    # TODO this may not necessarily add up to 1.0 (so make sure we scale it so it does)


    # Scale the weights so that they add to 1.0 
    f(x, y) = x / y
    weights = f.(weights, sum(weights))         # TODO modify in place with broadcast? 
    


    # Get the index that our random number 'falls' on
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

    weights = [2, 4, 2, 12]



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
