Programmed population control
-----------------------------

Design an extension of the bacterial population control system
developed by You et al that shows better performances. Indeed, the only demonstrated way to
tune the cell density at steady state of the induced system is by changing the pH of the growth
media. Firstly, this is not a very practical solution. Secondly, the dynamic range –defined here as
the ratio of the highest to the lowest cell density that one can get at steady state– is rather limited
this way.
To tentatively address these limitations, we investigate whether the additional use of an inducible
promoter to control the intracellular levels of I and/or R proteins could lead to a system with a
larger dynamic range for cell densities. 
We will do this in three steps. 
Firstly, we will model and simulate the existing system, and compute its predicted dynamic range. 
Secondly, we will propose various extensions of the existing system, and model them. 
Thirdly, we will simulate their behaviors and compute their predicted dynamic ranges. 

********************************************************

Folders:
--------

Data:
-----
Contains the data generated randomly and during optimization using CMAES

Plots:
------
Contains the plots generated during analysis

Scripts:
---------
Contains MATLAB scripts that are used for:
a) generate data
b) differential equations to define the model
c) cmaes for parameter optimization
d) computing cost while optimizing parameters
e) computing steady state parameters



