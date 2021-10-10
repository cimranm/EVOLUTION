#! /usr/bin/env julia

"""TODO

- don't show stacktrace when throwing error
- convert to uppercase?

"""



"""
# ALIGNMENT
# implementation of Global Alignment 
# Needleman-Wunsch algorithm 
# -----------------------------------------------------

# usage:  align < INFILE > OUTFILE 
# 
# Treats first sequence as 'reference' sequence; compares all sequences to this first sequence
# Outputs alignm

julia version 1.6.2

Written by Cam McMenamie
"""

using FASTX     # BioJulia / FASTX.jl       https://biojulia.net 

"""
Global variables
"""
const SUBSTITUTIONS = ['A', 'B', 'C', 'D', 'E',	'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'X', 'Y', 'Z']
const BLOSUM62 = [

    #A 	B 	C 	D 	E 	F 	G 	H 	I 	K 	L 	M 	N 	P 	Q 	R 	S 	T 	V 	W 	X 	Y 	Z
    4 	-2 	0 	-2 	-1 	-2 	0 	-2 	-1 	-1 	-1 	-1 	-2 	-1 	-1 	-1 	1 	0 	0 	-3 	-1 	-2 	-1
    -2 	6 	-3 	6 	2 	-3 	-1 	-1 	-3 	-1 	-4 	-3 	1 	-1 	0 	-2 	0 	-1 	-3 	-4 	-1 	-3 	2
    0 	-3 	9 	-3 	-4 	-2 	-3 	-3 	-1 	-3 	-1 	-1 	-3 	-3 	-3 	-3 	-1 	-1 	-1 	-2 	-1 	-2 	-4
    -2 	6 	-3 	6 	2 	-3 	-1 	-1 	-3 	-1 	-4 	-3 	1 	-1 	0 	-2 	0 	-1 	-3 	-4 	-1 	-3 	2
    -1 	2 	-4 	2 	5 	-3 	-2 	0 	-3 	1 	-3 	-2 	0 	-1 	2 	0 	0 	-1 	-2 	-3 	-1 	-2 	5
    -2 	-3 	-2 	-3 	-3 	6 	-3 	-1 	0 	-3 	0 	0 	-3 	-4 	-3 	-3 	-2 	-2 	-1 	1 	-1 	3 	-3
    0 	-1 	-3 	-1 	-2 	-3 	6 	-2 	-4 	-2 	-4 	-3 	0 	-2 	-2 	-2 	0 	-2 	-3 	-2 	-1 	-3 	-2
    -2 	-1 	-3 	-1 	0 	-1 	-2 	8 	-3 	-1 	-3 	-2 	1 	-2 	0 	0 	-1 	-2 	-3 	-2 	-1 	2 	0
    -1 	-3 	-1 	-3 	-3 	0 	-4 	-3 	4 	-3 	2 	1 	-3 	-3 	-3 	-3 	-2 	-1 	3 	-3 	-1 	-1 	-3
    -1 	-1 	-3 	-1 	1 	-3 	-2 	-1 	-3 	5 	-2 	-1 	0 	-1 	1 	2 	0 	-1 	-2 	-3 	-1 	-2 	1
    -1 	-4 	-1 	-4 	-3 	0 	-4 	-3 	2 	-2 	4 	2 	-3 	-3 	-2 	-2 	-2 	-1 	1 	-2 	-1 	-1 	-3
    -1 	-3 	-1 	-3 	-2 	0 	-3 	-2 	1 	-1 	2 	5 	-2 	-2 	0 	-1 	-1 	-1 	1 	-1 	-1 	-1 	-2
    -2 	1 	-3 	1 	0 	-3 	0 	1 	-3 	0 	-3 	-2 	6 	-2 	0 	0 	1 	0 	-3 	-4 	-1 	-2 	0
    -1 	-1 	-3 	-1 	-1 	-4 	-2 	-2 	-3 	-1 	-3 	-2 	-2 	7 	-1 	-2 	-1 	-1 	-2 	-4 	-1 	-3 	-1
    -1 	0 	-3 	0 	2 	-3 	-2 	0 	-3 	1 	-2 	0 	0 	-1 	5 	1 	0 	-1 	-2 	-2 	-1 	-1 	2
    -1 	-2 	-3 	-2 	0 	-3 	-2 	0 	-3 	2 	-2 	-1 	0 	-2 	1 	5 	-1 	-1 	-3 	-3 	-1 	-2 	0
    1 	0 	-1 	0 	0 	-2 	0 	-1 	-2 	0 	-2 	-1 	1 	-1 	0 	-1 	4 	1 	-2 	-3 	-1 	-2 	0
    0 	-1 	-1 	-1 	-1 	-2 	-2 	-2 	-1 	-1 	-1 	-1 	0 	-1 	-1 	-1 	1 	5 	0 	-2 	-1 	-2 	-1
    0 	-3 	-1 	-3 	-2 	-1 	-3 	-3 	3 	-2 	1 	1 	-3 	-2 	-2 	-3 	-2 	0 	4 	-3 	-1 	-1 	-2
    -3 	-4 	-2 	-4 	-3 	1 	-2 	-2 	-3 	-3 	-2 	-1 	-4 	-4 	-2 	-3 	-3 	-2 	-3 	11 	-1 	2 	-3
    -1 	-1 	-1 	-1 	-1 	-1 	-1 	-1 	-1 	-1 	-1 	-1 	-1 	-1 	-1 	-1 	-1 	-1 	-1 	-1 	-1 	-1 	-1
    -2 	-3 	-2 	-3 	-2 	3 	-3 	2 	-1 	-2 	-1 	-1 	-2 	-3 	-1 	-2 	-2 	-2 	-1 	2 	-1 	7 	-2
    -1 	2 	-4 	2 	5 	-3 	-2 	0 	-3 	1 	-3 	-2 	0 	-1 	2 	0 	0 	-1 	-2 	-3 	-1 	-2 	5]

const NUCLEOTIDES = ['A', 'C', 'G', 'T']
const MATRIX = [
    2 0 0 0
    0 2 0 0
    0 0 2 0
    0 0 0 2]

const GAP_OPEN = 2
const GAP_EXTEND = 1

"""
Get index for matrix from letter.  TODO Change back to AA
"""
function getindex(a::Char)
    for i in 1:length(NUCLEOTIDES)
        if NUCLEOTIDES[i] == a 
            return i
        end 
    end
end

"""
Get score for match / mismatch.
"""
function matchscore(a::Char, b::Char)

    i = getindex(a)
    j = getindex(b)
    return MATRIX[i, j]
end


"""
Align two sequences. 
"""
function align(a::Vector{Char}, b::Vector{Char})

    m = length(a) + 1
    n = length(b) + 1

    println("M: $m")


    # Create matrices 
    S =     Matrix(undef, m, n)
    Ix =    Matrix(undef, m, n)
    Iy =    Matrix(undef, m, n)

    # Initialise to minus infinity to ensure 'max' is chosen (certain cells in matrices might be undefined)
    for i in 1:m 
        for j in 1:n 
            S[i, j] = -Inf
            Ix[i, j] = -Inf
            Iy[i, j] = -Inf 
        end
    end

    # INITIALISE ORIGIN
    S[1, 1] = 0 
    Ix[1, 1] = 0
    Iy[1, 1] = 0
    
    # Initialise left-hand column
    S[2, 1] = -GAP_OPEN 
    Iy[2, 1] = -GAP_OPEN
    for i in 3:m 
        val = S[i-1, 1] - GAP_EXTEND 
        S[i, 1] = val 
        Iy[i, 1] = val
    end

    # Initialise top row
    S[1, 2] = -GAP_OPEN
    Ix[1, 2] = -GAP_OPEN
    for j in 3:n 
        val = S[1, j-1] - GAP_EXTEND 
        S[1, j] = val
        Ix[1, j] = val
    end



    
    for i in 2:m 
        for j in 2:n

            # Set S to be max(continue match, end insertion in x, end insertion in y)

             # index a and b with -1 to account for index offset
            match = matchscore(a[i-1], b[j-1])
            S[i, j] =   max(S[i-1, j-1] + match,
                        Ix[i-1, j-1] + match,
                        Iy[i-1, j-1] + match)  
                        
            # Set Ix 
            Ix[i, j] =  max(S[i-1, j] - GAP_OPEN,
                            Ix[i-1, j] - GAP_EXTEND)

            Iy[i, j] =  max(S[i, j-1] - GAP_OPEN, 
                            Iy[i, j-1] - GAP_EXTEND)

        end
    end

    display(S)


    

    #=
    for i in 1:m 
        for j in 1:n 
            S[i, j] = -1
        end
    end
    =#


    """

    For i 
        for j 

            down =      # work out affine; i.e. if previous had gap opening or what
            across = 

            S[i, j] = min(diagonal + score(i, j), down, up)
            
            assign x, y


    """


    # 
    """
    Traceback 
    
    i = m+1
    j = n+1
    k = 0   # alignment length 
    id = 0  # count of identical residues matches

    while (i > 1 or j > 1) 
        iprev = x[i, j]
        jprev = y[i, j]

        if iprev == i       
            j = jprev

        elseif jprev == j   # Gap 
            i = iprev

        else                # match

            # check if residue exact match
            i = iprev
            j = jprev
        end

        k += 1

    end

    percentid = 100 * (id / k)
    """
end

"""
Check that AA sequence is valid. 
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

    # Read FASTA formatted sequences from STDIN 
    r = FASTA.Reader(stdin)

    # Store all FASTA records 
    records = []
    sequences = []
    for record in r 
        push!(records, record)
        seq = sequence(record)
        seq = convert(String, seq)
        seq = Vector{Char}(seq)
        push!(sequences, seq)
    end 

    # Handle invalid input
    if length(records) < 1
        error("No sequence given.")
    end 

    # Save the 'reference' sequence (first one input)
    rec = records[1]
    id = identifier(rec)
    seq = sequence(rec)
    seq = convert(String, seq)  # Represent our sequence as a String 

    REFSEQ = Vector{Char}(seq)
    validate(REFSEQ)           # Check for invalid letter codes

    #=
    for seq in sequences

        alignment = align(seq, REFSEQ)

        # Output
        n = 1
        score = 2
        id = 0.80
        println("$n\t$score\t$id")
        
    end
    =#


    a = "GAACGT"
    a = Vector{Char}(a)
    b = "GAT"
    b = Vector{Char}(b)
    align(a, b)

        # TODO make 'alignment' object? so you can do 
        # alignment.getscore() or something?
end

# Run main function
main()