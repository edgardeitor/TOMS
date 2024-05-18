# TOMS
Villar-Sepúlveda, E., & Champneys, A. (2023). Computation of Turing Bifurcation Normal Form for n-Component Reaction-Diffusion Systems. ACM Transactions on Mathematical Software, 49(4), 1-24.

Scripts to compute the criticality of the Turing bifurcation.

There are two folders. In the folder 'Python + Mathematica', you need to run two scripts to get the bifurcation curves. Though this seems to be more complicated, the combination of scripts can be more efficient than the script in the other folder that does all the calculation in Mathematica. You are free to run both or make use of the one you feel most comfortable with.

INSTRUCTIONS IF YOU WANT TO RUN THE CODE IN FOLDER 'Python + Mathematica'

The files 'Criticality.py', 'functions.py' and 'Plotter.nb' must be placed into the same folder to run the script.

Create different directories for different models and place a file called 'name_of_the_directory.py' with all the required variables indicated in it.

The only thing that remains to be done is running the script 'Criticality.py' and write the name of the directory as the first input when you are asked to enter it. After that, get into the folder and run the Mathematica file 'Plotter.nb' that will be already copied (provided that all the variables were well-provided) into the directory of your system for you.

INSTRUCTIONS IF YOU WANT TO RUN THE CODE IN FOLDER 'Mathematica only'

Create different directories for different models and place the file 'Plotter.nb' inside of them. Input the variables explained in the file and run the script. The script will output some files. For that reason, it is recommended to create a different folder for each model that is to be analyzed.
