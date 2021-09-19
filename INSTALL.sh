#! /usr/bin/env bash
filename="evolve.jl"

# install julia TODO
# Using julia version 1.6.2 



chmod +x "$filename"

# install FASTX package
julia -e 'using Pkg; Pkg.add("FASTX")'
