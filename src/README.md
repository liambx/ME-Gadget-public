ME-Gadget

Jiajun Zhang
contact me by email: liamzhang@sjtu.edu.cn

ME-Gadget is a modified version of the publicly avaiable N-body simulation code Gadget2.
ME is short for dark Matter and dark Energy.
It also means you may run N-body simulations with arbitrary cosmological models.
I have used this code to run Interacting Dark Energy (IDE) model, effective field theory model and modified gravity model.
Please feel free to use this code to run cosmological N-body simulations for your non-standard cosmological model.
just remember to cite the following papers:

https://ui.adsabs.harvard.edu/abs/2019MNRAS.tmp.1961A/abstract
https://ui.adsabs.harvard.edu/abs/2019ApJ...875L..11Z/abstract
https://ui.adsabs.harvard.edu/abs/2018PhRvD..98j3530Z/abstract

The details and tests about the code is given in:

https://ui.adsabs.harvard.edu/abs/2018PhRvD..98j3530Z/abstract

Dependency:

The code need the following library to be installed:

a linux system in your computer
OpenMPI v2.0 or higher
GSL v2.4 or higher
FFTW v2.1.5

You may install the packages using my one-click script:
https://github.com/liambx/libinstall


What to provide:
1. a table about how the universe expand, containing one column of scale factor a and another column of H/H_0, calculated by Friedmann equation. FILE: hubble_table.txt
2. a table about how the dark matter particle mass change, containing one column of scale factor a and another column of the simulation particle mass in the unit of 10^10 M_sun/h. FILE: dmmass_table.txt
3. a table about how the dark matter particle feels drag force which is only determined by the velocity of the particle, containing one column of scale factor a and another column of -3(xi_1+xi_2/r)Ha for DMDEI model, which is obtained from the Euler equation. the unit should be in Gyr^-1 FILE drag_table.txt
4. a table about how the gravitational constant changes comparing to newtonian gravitational constant, containing one column of scale factor a and another column of G_eff/G_newtonian. FILE: varg_table.txt 

Warning:
The code has only passed cosmological N-body simulation with periodic boundary condition functional test.
I don't garantee any hydrodynamic use or non-periodic boundary condition or zoom-in use.

How to use:
First of all, you need to know how to use the original Gadget2 code to do a cosmological N-body simulation with periodic boundary condition.

In the Makefile, there are five additional options, they are:
just uncomment the needed OPT, you can use the related functions, as long as you prepared the needed tables in the right name in the code running folder.

#OPT   +=  -DHUBBLE_TABLE 
#OPT   +=  -DDMMASS_TABLE 
#OPT   +=  -DDRAG_TABLE
#OPT   +=  -DDE_TABLE 
#OPT   +=  -DVARG_TABLE 

the DE_TABLE function is not available in the public version, if you need it, as decribed in my paper, please contact me.
Any more questions, please contact me, I am mostly availabe and happy to listen to comments and suggestions.




