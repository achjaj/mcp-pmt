# How to run the script

The best way how to run the script is by this command:
```
julia -L 3d.jl
```

This will load the script into Julia shell. When you want to run the simulation again just type `run()` into the shell. This way you do not have to wait for the slow initialization after the first run.

# Plotting the result
Type this command into the running Julia shell, after the simulation finished

```julia
drawSteps(steps)[1]
```

# Changing parameters "on the fly"
If you want to change some parameters after the first run and run the simulation again, then you need to call the `run` function with the appropriate keyword argument. For example, if you want to change the bias voltage to 2 MV, type

```
run(V_in=2e6)
``` 
into the Julia shell. Or if you want to change the channel length to 1 cm and diameter to 3 mm type

```
run(L_in = 1e-2, R_in = 3e-3)
``` 

See line 68 in file [3d.jl](3d.jl) for list of available parameters.

# Second simulation
To run the second script, run this command
```
julia -L 3d2.jl
```

The rest is the same as for the first simulation.

# Plotting the result of the second simulation
Type 
```julia
drawPositions(positions)[1]
```
into the Julia shell

# Changing parameters
Same as for the first simulation