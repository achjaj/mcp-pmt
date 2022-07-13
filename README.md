# mcp-pmt

## guest. simulation output

The output from the simulation is a JSON array of arrays. These inner arrays represent individual
electrons inside the channel. The electron array is also an array of arrays
with this structure: `[T₀, Tc, ϕ, θ, l₀, l, idx]`.
  - `T₀` - array of initial energies; the first element of this array is the energy of the electron at emission, the rest are eneries of the electron after collisions
  - `Tc` - array of collision energies; it's the energy the electron has when it collides with the wall
  - `ϕ` - array of initial angles measured from normal
  - `θ` - array of collision angles measured from normal
  - `l₀`- distance between channel opening and point of emission
  - `l` - array of distances where the electron collided with the channel; first element is l₀
  - `idx` - unique ID of the electron for tracking

### Units
 - Energies: electronvolts
 - Angles: radians
 - Distances: meters

### Example output
``` json
[
  [[100, 105], [-1, 150], [0.7, 0.75], [-1.0, 1.2], [0], [0, 1e-5], [1]] <-- this is one electron,
  [[21, 18], [-1, 27], [0.65, 0.6], [-1, 1.5], [1e-5], [1e-5, 5e-5], [2]] <-- this is another electron
]

```
The negative ones should be ignored as they are used to initialize the arrays. I shall remove them in the future.
