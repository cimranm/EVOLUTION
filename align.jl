#! /usr/bin/env julia
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


"""
Align two sequences. 
"""
function align(a::Vector{Char}, b::Vector{Char})

    m = length(a)
    n = length(b)
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
        push!(sequences, sequence(record))
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


    for seq in sequences

        #alignment = align(seq, REFSEQ)

        # Output
        n = 1
        score = 2
        id = 0.80
        println("$n\t$score\t$id")
        
    end

        # TODO make 'alignment' object? so you can do 
        # alignment.getscore() or something?
end

# Run main function
main()