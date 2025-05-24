# 2D Ising Model

This repo contains a small code to simulate an Ising model on a rectangular 2D domain. You can find some explanations in the '[note\_ising](https://github.com/roomate/2D-Ising-Model/blob/master/note_ising.md)' markdown.

The command line allows you to tune out some parameters, such as the dimensions of the lattice, or the initial conditions. More specifically, the parser admits the following flags:

- --B. External magnetic field applied in Tesla. Default Value: 0 T.
- --T. Temperature of the thermal bath in kelvin. Default Value: 300K.
- --J. Exchange Interaction in Joule. Default Value: $4e-21$ J.
- --bound\_cond. Boundary conditions of the domain. Default Value: "periodic".
- --Initialization. How to initialize the domain. Default Value: "random".
- --cmap. Color map to display the lattice. Default Value: "viridis".

To run the simulation, run in the folder of your choice:
`mkdir Ising-Simulation; git clone git@github.com:roomate/2D-Ising-Model.git`

Then, to run the simulation, type the prompt:
`python3 ising_model`
followed by the potential flags desired, as described above.
