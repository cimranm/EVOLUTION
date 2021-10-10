#! /usr/bin/env julia

"""
# ALIGNMENT
# implementation of Global Alignment 
# Needleman-Wunsch algorithm 
# -----------------------------------------------------

usage:  align < INFILE > OUTFILE

evolve < INFILE | align > OUTFILE 



# Treats first sequence as 'reference' sequence; compares all sequences to this first sequence
# Outputs alignment number, score, and percentage identity (tab-delimited)

julia version 1.6.2

Written by Cam McMenamie
"""

using FASTX     # BioJulia / FASTX.jl       https://biojulia.net 

"""
Global variables
"""
const SUBSTITUTIONS = ['A', 'B', 'C', 'D', 'E',	'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'X', 'Y', 'Z']
const MATRIX = [

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

const GAP_OPEN = 7
const GAP_EXTEND = 1

# 'Enum' for pointer matrix 
const POINTER_DIAGONAL = 1
const POINTER_INSERTX  = 2
const POINTER_INSERTY   = 3

"""
Get index for matrix from letter.  TODO Change back to AA
"""
function getindex(a::Char)
    for i in 1:length(SUBSTITUTIONS)
        if SUBSTITUTIONS[i] == a 
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

    # Create matrices 
    S =     Matrix(undef, m, n)
    Ix =    Matrix(undef, m, n)
    Iy =    Matrix(undef, m, n)
    Px  =   Matrix(undef, m, n)
    Py  =   Matrix(undef, m, n)
    P  =   Matrix(undef, m, n)

    # Initialise to minus infinity to ensure 'max' is chosen (certain cells in matrices might be undefined)
    for i in 1:m 
        for j in 1:n 
            S[i, j]     = -Inf
            Ix[i, j]    = -Inf
            Iy[i, j]    = -Inf 
            Px[i, j]    = 0
            Py[i, j]    = 0
            P[i, j]     = 0
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

            paths = [S[i-1, j-1] + match,   # DIAGONAL
                    Ix[i-1, j-1] + match,       # INSERTION IN X
                    Iy[i-1, j-1] + match]       # INSERTION IN Y

            # Store the co-ordinates from which we 'arrived' at this cell (i.e. save the best path we took that gave us the 
            # max score)

            pathindex = Int(argmax(paths))
            S[i, j] = paths[pathindex] # store max value in S 

            P[i, j] = pathindex
            # 'Diagonal' case
            if pathindex == 1

                

                Px[i, j] = i-1 
                Py[i, j] = j-1
            
            # Insertion in X (gap in Y)
            # We decrement back along X (sequence A) while leaving a gap in Y
            elseif pathindex == 2
                Px[i, j] = i-1 
                Py[i, j] = j
            
            # Insertion in Y (gap in X)
            # We decrement back along Y (sequence B) while leaving a gap in X
            elseif pathindex == 3
                Px[i, j] = i 
                Py[i, j] = j-1

            end
                
                        
            # Set Ix 
            Ix[i, j] =  max(S[i-1, j] - GAP_OPEN,
                            Ix[i-1, j] - GAP_EXTEND)

            Iy[i, j] =  max(S[i, j-1] - GAP_OPEN, 
                            Iy[i, j-1] - GAP_EXTEND)

        end
    end

    #display(S)

    # Backtracing
    aout = []
    bout = []

    
    i = m 
    j = n 

    score = S[i, j]

    k = 0   # Alignment length 
    id = 0  # Number of exact matches 
    
    while i > 1 && j > 1 

        iprev = i  
        jprev = j 

        if P[i, j] == POINTER_DIAGONAL 

            i -= 1
            j -= 1

            push!(aout, a[iprev-1])
            push!(bout, b[jprev-1])

            # Check if match 
            if a[iprev-1] == b[jprev-1]
                id += 1
            end 

        elseif P[i, j] == POINTER_INSERTX 

            i -= 1 
            push!(aout, a[iprev-1])
            push!(bout, '-')

        elseif P[i, j] == POINTER_INSERTY 

            j -= 1
            push!(bout, b[jprev-1])
            push!(aout, '-')

        end

        k += 1  # increment counter 

    end 

    percentid = 100 * (id / k) 
    return score, percentid

    


    #=
    i = m
    j = n 
    while i > 1 && j > 1

        iprev = i 
        jprev = j

        if i == iprev
            print("-")
        elseif j == jprev 
            print(a[i-1])
        else 
            print(a[i-1])
        end
        
        if i == 0 || j == 0
            break
        end
        i = Px[i, j]
        j = Py[i, j]
    end
    =#


        





    

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

    n = 0 
    for seq in sequences

        score, id = align(seq, REFSEQ)

        # Output
        println("$n\t$score\t$id")
        
        n += 1
        
    end    
end

# Run main function
main()