# EVOLUTION

Simulate evolving an amino acid sequence

[Spec](http://www.cse.unsw.edu.au/~bi3020/21T3/spec11.html)


```bash 
wget http://www.cse.unsw.edu.au/~bi3020/21T3/sample_seqs.zip
unzip sample_seqs.zip 
```

## Hardcoding matrix.txt to Julia matrix 

```bash
# Getting from text to matrix in julia
cat matrix.txt | cut -d ',' -f2-22 | tr , ' ' >> evolve.jl 
cat BLOSUM62.txt | cut -f2-24 >> align.jl

# Display words per line; should be 20 (quick check to make sure we haven't left out any AAs) 
awk '{print (NF, $0)}' 
```


## DEPENDENCIES 

Uses FASTX package from BioJulia.  Can be installed from Julia using REPL 
```bash
$ julia 
``` 

```julia 
using Pkg 
Pkg.add("FASTX")
```


## Julia syntax to keep in mind 

### Broadcasting 

functions can 'broadcast' to each element of an array e.g. 

```julia
f(x) = x + 1
array = [1, 2, 3, 4]

f(array)    # Perform function on whole array
f.(array)   # Perform function on each element in array 
```

### Indexing 

Indexes start at `1`, not `0`. 

### Mutable functions 

Functions that 'mutate' or change their inputs are named with a `!` by convention.  For example, 
```
pop!(list)  # Remove the last element
```

