# CalculatePolarization

The code calculates the spin polarization of hyperons, using the spin Cooper-Frye formula that associates
the momentum-space distribution of the hyperon polarization with the position-space vorticity of the fluid.

The repository contains a modified version of MUSIC code(MUSIC_POLZ) which can generate required format of hypersurface file
which will go inside the "CalculatePolarization" code for polarization calculation.

Currently, the "init.cpp" in "MUSIC_POLZ" is comaptible to take SARJ initial profile. One can modify the function 
"init.cpp/initial_IPGlasma_XY_TRIBHUBAN()"  in order to accomodate any arbitrary initial profile.
