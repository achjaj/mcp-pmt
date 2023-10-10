# Microchannel plate photomutiplier simulations

This repository contains source codes and model files for
simulations I tried to develop as part of my masters' thesis.
The goal was to develop a full simulation of MCP-PMT with an
ALD coating. I managed to write very simple simulations of
a single MCP channel without ALD. None of the simulations
provides realistic results. Nevertheless, I hope that my work
will be useful as a starting point for someone's future work.

I worked on four types of simulations:

* Macroscopic TLM model based on a model developed by
[Guidicotti](dx.doi.org/10.1016/j.nima.2011.07.017)
* "Quasyanalitical" microscopic model that exploited
homogeneity of the electric field inside a channel to simplify
the simulation
* Particle-In-Cell simulation that calculated positions and
velocities of electrons in discrete time steps
* Model developed using COMSOL Multiphysics software

### TLM model

The TLM model uses Transmission Line Modeling methods to
derive an implicit equation for MCP-PMT gain saturated due to
wall charge effects. The equation can be solved numerically
in an iterative manner. The advantage of this simulation is
its performance. It works with macroscopic quantities such as
electric current and does not have to keep information about
every single electron. This means that it has low memory and
CPU requirements. A single simulation can be finished within a
few minutes. For more information, see the original [Guidicottis](dx.doi.org/10.1016/j.nima.2011.07.017)
paper, my thesis or the source codes.

### Quasianalitical simulation

This model tracks the positions and energies of all electrons
that are in a microchannel. This requires a lot of memory, so
I tried to simplify and optimise the simulation. The
electric field created by the bias voltage should be
homogeneous and time-invariant inside a microchannel. This
means that the equation of motion of an electron in the channel
can be derived and solved beforehand. From the EOM, I derived
equations for collision energies and collision angles; these
quantities are important for calculating the secondary
emission yield. The simulation then works by simply
calculating the relevant quantities from derived equations.
See the thesis and the source codes for more information.

### PIC model

This model also tracks the positions and energies of every
single electron. However, this simulation is not optimised
and instead calculates the energies and positions in discrete time steps from velocity and position vectors
calculated in a previous step. It uses a lot of memory and
CPU time. It is a slow simulation that can be improved by
utilising multithreading and multiprocessing. I tried to take
this approach, but I could not get it fully working.
See the thesis and the source codes for more information.

### COMSOL Multiphysics

After three failed attempts to develop a working simulation
I tried to simplify my work by using commercial simulation
software. COMSOL allows one to create simulations without
the need to write one's own code. Most of the things can be done 
using a graphical interface. However, I had specific 
requirements for secondary emission, and the simple model 
that COMSOL uses was not enough for me. I had to resort to
writing my own C library that can be used inside COMSOL.
This simulation is very similar to the PIC model. However, it
should perform much better as the underlying COMSOL code was
written by professional programmers who know how to optimise
such simulations. See the thesis and the source codes for
more information.

### Secondary emission

I used two models of secondary emission. At the beginning, I
used a simple model used by [Guest](https://psec.uchicago.edu/Papers/Guest.pdf). The model works by
calculating a mean yield after a collision. This mean is
then put into Poisson distribution, and a random number of
secondary electrons is generated.

The second model I used was developed by Furman and Pivi.
This model is much more complicated but more realistic.
This model was also harder to use in COMSOL, and I had to
resort to writing my own C library.

### Failures

The TLM model looks very promising however it is a
macroscopic model that does not allow to calculate the
influence of microscopic effects as fringe fields or
interactions between the electrons. The model also uses
assumptions that I find unphysical.

The quasianalitical and the PIC models are poorly written
and are very slow. They are unable to work with a high number
of electrons. The mistake I made was my choice of
programming language. The code is written in Julia. The
language uses Just-In-Time compilation, which means that the
code should perform as well as a C code. The language is also
easy to code and work with and offers a lot of packages
similar to Python. However, the language uses a garbage
collector, which significantly reduces the performance of the
code I wrote. I tried to get around it, but I was not able
to. But this does not mean that the language is the problem.
It just means that I am not skilled enough to be able to
write good code.

The COMSOL models showed the best results. They were much
faster than my codes, and at the beginning, they also produced
reasonable data. The problem is the secondary emission.
When the C library with the Furman-Pivi model is used, the
electron cloud that is created in the simulated microchannel
have strange behaviour. It does not propagate. I think that
the problem is in the library, and my implementation of the
Furman-Pivi model.

## Directory structure

* `literature`: contains a copy of the thesis and a list of
used literature 

* `simulations/comsol`: all the COMSOL models I developed.
It includes the source code of the C library and some highly
experimental models that did not end up in the thesis

* `simulations/monte-carlo`: source codes of the PIC Montr
Carlo simulation

* `simulations/tlm`: source codes of the TLM model and some
results

## Contact
If you have any questions, feel free to contact me via email:
bucko.jakub@proton.me