# mcp-pmt

## guest. simulation output

The output file from the simulation is in JSON format. It conatins JSON array which contains arrays that represents electrons.
The electron array has this format: `[T₀, Tc, ϕ, θ, l₀, l, idx]`.
  - `T₀` - array of initial energies; the firts element of this array is the energy of the electron at emission, the rest are eneries of the electron after collisions
  - `Tc` - array of collision energies
  - `ϕ` - array of initial angles measured from normal
  - `θ` - array of collision angles measured from normal
  - `l₀`- distance from the beggining of channel where the electron was emitted
  - `l` - array of distances where the electron collided with the channel; first elemnt is l₀
  - `idx` - unique ID of the electron for tracking

### Units
 - Energies: electronvolts
 - Angles: radians
 - Distances: meters

### Example output
``` json
[
  [[100, 105], [-1, 150], [0.7, 0.75], [-1.0, 1.2], [0], [0, 1e-5], [1]],
  [[21, 18], [-1, 27], [0.65, 0.6], [-1, 1.5], [1e-5], [1e-5, 5e-5], [2]]
]

```
The negative ones should be ignored as they are used to initialize the arrays. I shall remove them in the future.
