Julia Code:
neogrowth_script_deterministic.jl
- Main code that runs the model and plots policy functions
- This is the code file that you should run to run the Julia code. Just make sure neogrowth_model_deterministic.jl is in the same folder.
neogrowth_model_deterministic.jl
- Code that contains auxiliary functions used to solve value function iteration

Fortran Code:
neogrowth_deterministic.f90
- I've included instructions in the file for how to run the fortran code.
- The most important thing is you'll need a compiler (I use gfortran)
- If you don't want to install a fortran compiler on your local device, linstat has gfortan installed.
- If you're struggling to figure out how to read/understand/compile/run the fortran code, don't worry. We'll go over how to do so in our first meeting. 
