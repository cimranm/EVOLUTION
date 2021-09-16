#! /usr/bin/env julia

const GENERATIONS = 500


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




# Benchmark

# TODO: able to give weights as relative numbers 
# i.e. function will automatically add them to total, and divide by total to give prob (so don't 
# the user doesn't have to make sure that they sum to 1.0 exactly ) 

# Make dict?      Dict(x->f(x) for ....)

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


