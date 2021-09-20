"""
View statistics for various test cases
"""
function probability_analysis()

    # Initialise dict 
    count = Dict()
    for item in items
        count[item] = 0
    end 

    # Sample randomly

    total = 10^3

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
end