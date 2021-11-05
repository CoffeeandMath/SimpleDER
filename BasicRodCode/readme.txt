#############################################################################
### README
#############################################################################

This readme is intended as a guide to aid in setting up and running a basic version of the code implemented in the rest of the paper. Note that this code was developed and run on Ubuntu 20.04.3. Two libraries necessary for running the code are the Eigen3 (vector and tensor manipulation) library and GSL (GNU Scientific Library). CMake is used to control software compilation. Results are outputted as .vtk files and visualized using Paraview. I will write commands to be written in a linux terminal to manipulate folders and files.

#############################################################################
###### COMPILATION
#############################################################################

The first step is to create the build folder. In the source directory (BasicBeamCode), create a folder called "build". This is where the code will be compiled. Then navigate into that folder.

	mkdir build
	cd build

We then need to run CMake to set up the compiler. 

	cmake ..

If CMake is successful, we then want to make to code to start, this is done by running

	make

Next, we would like to create the folder for the data to be output in. We also want to import a data file that will have some basic simulation parameters. This file is included to be able to change basic parameters without having to recompile the code. This is done by running

	mkdir data
	cp ../inputFile.in .

All the vtk data files will be placed in this folder and be numbered. We will discuss later on how the distination folder and file names can be changed.

At this point, the code can be running the command

	./run

It will then ask for an input file. By default this is given by
	
	inputFile.in

At this point, the code should run and output files should appear in the data folder.


#############################################################################
###### SIMULATION SETUP
#############################################################################
In this basic simulation, we are modeling an anisotropic Kirchhoff rod of length 2 using 80 nodes. The first boundary is fixed in place and not allowed to rotate. To make the simulation interesting, we apply the following to the beam in increasing magnitude.

	1) Internal curvature
	2) Constant global gravitational load
	3) Force on the end point

This problem is non-trivial due to the 3D nature of the deformation. The rod will attempt to curl due to the internal curvature; however, the gravitational load and force on the end point will cause twist. The bending/twisting moduli are defined as in the Overcurvature example, and we give a width and height of the rod. We start with an initially straight rod and incrementally increase the magnitudes of the three items listed previously. This is done over 1000 steps; however, much smaller steps can be done as well.



#############################################################################
###### INPUT FILE
#############################################################################
The input file (inputFile.in) feeds in a few parameters in order to not have to recompile the code if we change a simple quantity such as number of nodes or number of iterations. The parameters given in this file are

Nn - Number of nodes
Nsteps - Number of solve increments
w - Width of the rod
h - Height of the rod
L - Length of the rod
YoungsMod - Youngs Modulus of the rod
rho - Unscaled Graviational Force
fend - Force on the end of the rod
SolverRead - Level of detail displayed in the terminal readout
foldername - Name of folder data is being stored to (data by default)
filename - Name of files being stored (dataout by default)

By changing any of these parameters, we can rerun the code without having to recompile.

#############################################################################
###### MAIN.CPP
#############################################################################

The main file is structured as follows:

	Importing data from inputFile.in
	Constructing the bending modulus vector
	Constructing vector of nodes
	Fixing the first two nodes in place
	Constructing vector of frames
	Fixing the rotation of the first frame
	Constructing vector of elements
	Constructing vector of material objects
	Constructing vector of local constructors
	Constructing rod object and associating previous vectors with it
	Constructing vector of constraints, particularly inextensibility
	Within the for loop:
		Iteratively increasing curvature, gravity, and end point load
		Pre-calculating Hessian (r1.Update(2))
		Solving the nonlinear equilibrium problem through Newton-Raphson (s1.Solve())
		Update the D frame using the newly solved d frames (r1.setdtoD());
		Saving output file


One thing to note in gradient and hessian calulations is the command r1.Update(k) where k is some integer. k denotes the order of update to do. k = 0 updates just the energy, k = 1 updates the energy and gradient, and k = 2 updates the energy, gradient, and hessian. This is done to eliminate some unneccesary calculations in the Newton-Raphson steps.


#############################################################################
###### OUTPUT/VISUALIZATION
#############################################################################
The output data is stored as iterative .vtk files in the data folder. It stores the nodal positions as well as the frames, projected onto the nodes. This allows for visualization of the deformation as well as twists induced. We use Paraview to visualize the results. 


#############################################################################
###### CHANGING THE STRAIN ENERGY FUNCTION
#############################################################################
The local strain energy function is defined in the /src/Material.cpp file. We strain energy functions can be written in their continuous forms. We can change the strain energy function and its derivatives by changing the Efunction, dEfunction, and ddEfunction, referring to the objective, gradient, and hessian, respectively. Rescaling for the discrete elastic rod formulation is done automatically in the implementation.