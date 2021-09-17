# EVOLUTION

Simulate evolving a sequence

[Spec](http://www.cse.unsw.edu.au/~bi3020/21T3/spec11.html)


```
$ wget http://www.cse.unsw.edu.au/~bi3020/21T3/sample_seqs.zip
$ unzip sample_seqs.zip 
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

