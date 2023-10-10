# COMSOL

This folder contains all the model files I created. It also
contains the source code of the C library. Unfortunately I could
not include the raw data, because the files were too big. Feel
free to contact me if you want the data. I am not sure anymore what
are the differences between the model files.

However, I used the model (and its copies) named
`particle.mph` for testing new concepts, for example
different geometries.

The numbers in the model names refer to the maximum number
of secondary particles the simulation can produce.

The models with the text `pic<some number>` are models in
which one simulated particle represents the specified numbers
of real electrons.

All the files may contain absolute file paths valid for my system.
You will have to change them.

## Bash scripts

Here is the description of the three bash scripts in the
folder:

 - `getframes.sh`: This script looks for all exported
 animations (in MP4 format) and exports the frames of the
 animation in PNG format into a folder named after the
 animation file name.

 - `launch_comsol`: This script sets up the environment for
 the COMSOL, and then it launches COMSOL. It also runs the
 `oom-kill-stop` script.

 - `oom-kill-stop`: This script looks for all the running 
 processes that have the string `comsol` in their name.
 Then, it tries to prevent the [OOM killer](https://rakeshjain-devops.medium.com/linux-out-of-memory-killer-31e477a45759) from
 killing the processes. I am not sure how well the script
 works. I think that if the system is almost completely out
 of memory, the OOM killer will kill the processes anyway.
 From what I know, this script should be safe to use, but it
 on your own risk.

 ## Julia scripts

 The folder also contains some Julia scripts. Here is a brief
 description of what they do:

 - `parseExport.jl`: This script defines `struct Particle`
 and a function to parse the exported data. The data can
 then be easily accessed through the `struct`. This script 
 is used in other Julia scripts.

 - `particleHistogram.jl`: This script generates histograms
 from the exported data.

 - `qyAnalysis.jl`: Old script that I used to analyse the
 exported positions and velocities of the particles.

 - `outputVelocities.jl`: This script generates the `pivi/animations/outVel.mp4` animation.

 ## The C library

 The folder `c_library` contains the source code of the C
 library. This library implements the Furman-Pivi model.
 The library is loaded using the `External` node. See the
 comments in the source code for more information. Also
 read the [documentation](https://doc.comsol.com/6.1/docserver/#!/com.comsol.help.comsol/application_programming_guide.14.54.html)
 and the page 109 from [here](https://doc.comsol.com/6.0/doc/com.comsol.help.comsol/COMSOL_ProgrammingReferenceManual.pdf).