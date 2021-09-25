#! /usr/bin/env julia
"""
# AMINO ACID SEQUENCE -- EVOLUTIONARY MUTATION MODELLER 
# -----------------------------------------------------

# usage:  evolve < INFILE > OUTFILE 

julia version 1.6.2

Written by Cam McMenamie
"""

using FASTX     # BioJulia / FASTX.jl       https://biojulia.net 

"""
TODO 
Detect AAs that aren't real 
"""

"""
Global variables
"""
const GENERATIONS = 500
const SUBSTITUTIONS = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
const SUBMATRIX = [
                    # A R N D C Q E G H I L K M F P S T W Y V
                    9867 2 9 10 3 8 17 21 2 6 4 2 6 2 22 35 32 0 2 18
                    1 9914 1 0 1 10 0 0 10 3 1 19 4 1 4 6 1 8 0 1
                    4 1 9822 36 0 4 6 6 21 3 1 13 0 1 2 20 9 1 4 1
                    6 0 42 9859 0 6 53 6 4 1 0 3 0 0 1 5 3 0 0 1
                    1 1 0 0 9973 0 0 0 1 1 0 0 0 0 1 5 1 0 3 2
                    3 9 4 5 0 9876 27 1 23 1 3 6 4 0 6 2 2 0 0 1
                    10 0 7 56 0 35 9865 4 2 3 1 4 1 0 3 4 2 0 1 2
                    21 1 12 11 1 3 7 9935 1 0 1 2 1 1 3 21 3 0 0 5
                    1 8 18 3 1 20 1 0 9913 0 1 1 0 2 3 1 1 1 4 1
                    2 2 3 1 2 1 2 0 0 9871 9 2 12 7 0 1 7 0 1 33
                    3 1 3 0 0 6 1 1 4 22 9947 2 45 13 3 1 3 4 2 15
                    2 37 25 6 0 12 7 2 2 4 1 9924 20 0 3 8 11 0 1 1
                    1 1 0 0 0 2 0 0 0 5 8 4 9875 1 0 1 2 0 0 4
                    1 1 1 0 0 0 0 1 2 8 6 0 4 9944 0 2 1 3 28 0
                    13 5 2 1 1 8 3 2 5 1 2 2 1 1 9924 12 4 0 0 2
                    28 11 34 7 11 4 6 16 2 2 1 7 4 3 17 9840 38 5 2 2
                    22 2 13 4 1 3 2 2 1 11 2 8 6 1 5 32 9869 0 2 9
                    0 2 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 9976 1 0
                    1 0 3 0 3 0 1 0 4 1 1 0 0 21 0 1 1 2 9947 1
                    13 2 1 1 3 2 2 3 3 57 11 1 17 1 3 2 10 0 2 9901 ]


"""
Get index where a given AA occurs in SUBSTITUTIONS 

TODO: maybe use map instead of this?
"""
function aminoindex(AA)

    for i in 1:length(SUBSTITUTIONS)
        if AA == SUBSTITUTIONS[i]
            return i 
        end
    end 
end


"""
Take a random sample of a probability distribution 
(Discrete outcomes, but with individual weights)

STEPS
    # Get all weights (i.e. relative probabilities)
    # Scale weights by their sum (so they all add up to 1.0)
    # 'Spread' weights over some distance (using cumsum)
    # Iterate through all array indexes; find array index SUCH THAT:
    #     rand() <= weights[index]
    #
    # This will be the index that we use (as it means that our random variable 
    # has 'fallen' on this region of the distribution)

"""
function sample(items, weights)

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

"""
For a given Amino Acid, get the weights associated for all possible transitions 
(As per the substitution matrix)
"""
function getweights(AA)

    index = aminoindex(AA)              # Get the column number for this AA 
    return SUBMATRIX[1:end, index]      # Return the column from the substitution matrix
end


"""
Evolve a sequence by 1 generation.  Mutable i.e. changes the given array.
"""
function evolve!(seq::Vector{Char})

    # For each Amino Acid in the sequence, get all the transition states 
    # i.e. all the possible outcomes (subsequent amino acids) and their
    # respective weights.  Then place these weights on a distribution
    # from which a sample is taken randomly.  The amino acid at this 
    # point of the sequence is set to be the result of this 'sampling'. 
    for AA in 1:length(seq)
        weights = getweights(seq[AA])               # Get weights for this AA transitioning (i.e. column in submatrix)
        seq[AA] = sample(SUBSTITUTIONS, weights)    # 'Mutate' this particular AA randomly 
    end 
end

"""
Check that AA sequence is valid
"""
function validate(seq::Vector{Char})

    for c in seq
        if c ∉ SUBSTITUTIONS
            error("Not a valid amino acid: $c")
        end
    end
end

"""
Main function.
"""
function main()

    # Read FASTA formatted sequence from STDIN 
    r = FASTA.Reader(stdin)

    # Store all FASTA records 
    records = []
    for record in r 
        push!(records, record)
    end 

    # Handle invalid input
    if length(records) < 1
        error("No sequence given.")
    end 
    if length(records) > 1 
        error("More than one sequence given.")
    end

    # Save attributes from the single FASTA record 
    rec = records[1]
    id = identifier(rec)
    seq = sequence(rec)
    seq = convert(String, seq)  # Represent our sequence as a String 

    # Writer for sending FASTA formatted sequences to STDOUT
    w = FASTA.Writer(stdout, width=60)

    # Write initial sequence numbered `0`
    rec = FASTA.Record("$id  0", seq)
    write(w, rec)

    # Convert sequence from String** to an array of characters (so each character can be mutated individually)
    #   ** Strings are immutable 
    seq = Vector{Char}(seq)

    validate(seq)



    # For each generation, evolve the sequence and output it as a new FASTA record
    for i in 1:GENERATIONS
        evolve!(seq)    # Evolve the sequence by 1 generation 

        # Make a FASTA record containing the description, the generation number, and the new sequence.
        # Write to STDOUT. 
        rec = FASTA.Record("$id  $i", join(seq))
        write(w, rec)
    end

    close(w)
    close(r)
end

# Run main function
main()





