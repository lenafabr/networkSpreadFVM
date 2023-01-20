# networkSpreadFVM
Finite volume simulation of concentration profiles spreading over network structures

Public code is currently under construction and more examples and implementation details will be provided at a future time.

# Basic starter guide:

1) Acquire the gfortran compiler
2) Go into source/ and time "make".
An executable file netmeshdynamicsFVM.exe will be generated
3) Go into examples/ and run one of the example parameter files.

Example 1:
Run with

../netmeshdynamicsFVM.exe example1_ca

This will simulate local calcium release from an example honeycomb network.

Two relevant output files:

example1_ca.out contains the total flux of calcium out of the network at each time point
example1_ca.snap.txt contains the concentration of free calcium and total buffer protein at all mesh cells, at each snapshot time

To parse this data and generate plots, an example matlab script is provided in scripts/example1_ca.m

Example 2:
Run with
../netmeshdynamicsFVM.exe example_PA

This will simulate spreading of a 'photoactivated' bolus of particles from an initial starter region.

To visualize the spread, use the matlab script in scripts/example_PA.m
