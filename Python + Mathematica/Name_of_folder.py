parameters={} '''Dictionary with the form 'name_of_parameter': value. It is important to put the name
              between quotation marks and not in LaTeX form. They will be transformed automatically
              by the script.'''

var=[] '''List of variables written as 'name_of_variable'. It is important to write the name of the variables
       between quotation marks.'''

diffmatrix=[] '''List of lists with the entries of the diffusion matrix of your reaction-diffusion equation.
              All symbolic terms must be written between quotation marks while numerical values do not
              require them. This matrix must be either a constant or it must depend on the parameters provided.'''

kinetics=[] '''List of kinetics of the reaction.diffusion equation. Again, each symbolic function must be
            written using quotation marks. This list must depend on the parameters and variables provided.'''

tol=1e-7 '''Tolerance to check whether there is a Turing bifurcation at a given point.'''

parameter_functions={} '''Dictionary to set parameters as functions of other parameters of the system
                       if necesary.'''

phiunit='n' '''This variable must be equal to either 'y' or 'n', depending on whether you would like
            \phi_1^1 to be a unit vector or not, respectively. By default, it is set as 'n'.'''

crosscoef='n' '''This variable must be equal to either 'y' or 'n', depending on whether you would like
              to get the cross-order coefficient, respectively. If you want to get it, you must provide
              a parameter with respect to which you would like to compute the coefficient and an
              expression of the equilibrium below, At least depending on the parameter you are to
              compute the cross-order coefficient with respect to.'''
             
numberofzerotemporalderivatives=0 '''This variable represents the number of time derivatives that are
                                    equal to zero in the equation. These equations must be placed at the end
                                    of the system.'''

crosspar='' '''This variable is a string that must have the name of the cross-order parameter if you want
               to compute the cross-order coeffiient.'''

fifthcoef='n' '''This variable must be equal to either 'y' or 'n', depending on whether you would like
              to get the expression of the fifth-order coefficient, respecively. If your system is large,
              the proess of obtaining this in terms of dummy variables could be quick but the saving
              and evaluating processes ould take ages. We only recommend you to set it as 'y' if the system
              is small or it does not have big expressions for the low-order vectors.'''

orthogonal='n' '''This variable must be equal to either 'y' or 'n', depending on whether you would like
               Q_1^3 to be orthogonal to \phi_1^1 or not, respectively. This variable is 'n' by default.'''

cod2='y' '''This variable must be equal to either 'y' or 'n', depending on whether you would like Mathematica
         to find codimension-two bifuration points on each of the bifurcation curves. If the expression of
         the third-order oefficient is big, this operation can take a long time. By default thisvariable
         is set as 'y'.'''
         
considerC3='n' '''This variable must be equal to either 'y' or'n', depending on whether you want to put the third
                order coefficient in all the calculations from third to fifth order.'''
         
alphaval='n' '''This variable must be equal to either 'y' or 'n', depending on whether you want to compute
                the coefficients \alpha_i. If this is set as 'y', then the calculation will be carried out if
                and only if fifthcoeff=='y' '''

plot2d='y' '''This variable must be equal to either 'y' or 'n', depending on whether you would like to get
           the output required by Mathematica to plot the bifurcation diagram in terms of two parameters.
           or not, respectively. If you set it as 'y', you must provide the names of the parameters that
           will be on the axes of the diagram, the intervals and lines to search, explained below.'''

parameters_on_axes=[] '''This is a list with two elements that are parameter names written between quotation
                      marks. These are the parameters that will use the [x,y]-axes in the diagram, respectively.'''

names_of_parameters=[] '''Python has some limitations with names. For instance, lambda cannot be a name of a
                       parameter beause it is a variable used for anoher thing. The vector names_of_parameters
                       lets you give your parameters the name you wish in the bifurcation diagram.
                       Be mindful that these must be strings written in LaTeX directly. For instance,
                       if you want the symbol lambda to be the name of one of your axes, then you can
                       input ['\\lambda','another_parameter'].'''

intervalx=[] '''List with two numerical values, a lower and an upper bound for the parameter you
             will have in the x-axis.'''

intervaly=[] '''List with two numerical values, a lower and an upper bound for the parameter you
             will have in the y-axis.'''

lines_to_search={} '''Dictionary with the format {'name_of_parameter':[list of values]}. Again, it is
                   important that the name of the parameter is written between quotation marks.
                   If the list of values has only one element, it an be written as a simple
                   numerical value instead of a list.'''

plot3d='n' '''This variable must be equal to either 'y' or 'n', depending on whether you would like to get
           the output required by Mathematica to plot the codimension-two curves in terms of three parameters.
           or not, respectively. This will be done only if plot2d is equal to 'y'.
           Furthermore, If you set it as 'y', you must provide the names of the parameters that
           will be on the axes of the diagram, the extra interval, and parameter values of
           the extra parameter for the Turing bifuration curves you would like to add to the diagram.'''

parameters_on_axes3=[] '''This is a list with three elements that are parameter names written between quotation
                       marks. These are the parameters that will use the [x,y,z]-axes in the diagram, respectively.'''

name_of_extra_parameter='' '''The analog thing to names_of_parameters but only for the extra one that
                           will appear in the 3D diagram.'''

extrainterval=[] '''Same thing as for intervalx and intervaly for the extra parameter.'''

extraturing=[] '''By default, after you run the script in Mathematica, you will get some codimension-two
               curves if they exist. Nevertheless, you may want to add a few Turing nbifurcation curves
               for some values of the extra parameter, that is, some curves restricted to some planes.
               Extraturing is a list with particular numerical values for which you would like to get the
               extra Turing bifuration curves. These extra curves will start from the codimension-two
               points so they will be associated with the codimension-two curve.'''