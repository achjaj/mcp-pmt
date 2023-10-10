# TLM model

## Files description

* `tlm.jl`: the source code of the TLM model written in
Julia language
* `sim-param.json`: configuration file of the simulation
* `tlm-ploting.jl`: create plots from the generated data
* `tlm-explore.jl`: small utility to explore the generated
data

## Configuration of the simulation

The configuration is done using the `sim-param.json` file.
It allows to set the parameters of the MCP-PMT detector and
other parameters of the calculations.

### List of the parameters:

These are the parameters of the simulated MCP-PMT (see the
[Gudicottis paper](dx.doi.org/10.1016/j.nima.2011.07.017) for more information):

* `Vs`: MCP bias voltage
* `Vm`: MCP crossover voltage
* `C`: anode capacitance
* `R`: anode resistance
* `Lmcp`: length of the MCP/microchannels
* `k`: gain coefficient

Other parameters are:

* `xcoords`: The coordinates of the mesh. The gain is 
calculated in these points. The value of this parameter can 
a list or a range. The range is defined by
    - `min`: the minimum value
    - `max`: the maximum value
    - `len`: the length of the range
    - `step`: alternatively, the difference between two values
    can be defined instead of the `len` parameter

* `times`: The times in which the gain is calculated. This can be defined by a list or range like `xcoords`

* `i0`: The input current. It can be defined as a constant or a function of `times`

* `psi0`: The initial value of the voltage equation. It can also be defined as a constant or as a function

* `outputDir`: Directory to which the output is written

* `output`: What should be written to the output
    - `"all"`: everything
    - `[name&[step]]`: array of variable names