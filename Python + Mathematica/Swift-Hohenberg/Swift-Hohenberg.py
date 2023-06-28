parameters={'r': 0,
        'q_c': 1,
        'b3': 1,
        'b5': 1
        }

var=['u',
     'v']

diffmatrix=[[0, -1],
            [-1, 0]]

kinetics=['r*u - q_c**4*u - 2*q_c**2*v + b3*u**3 - b5*u**5',
          'v']

tol=1e-7

parameter_functions={}

phiunit='n'

crosscoef='n'

numberofzerotemporalderivatives=1

crosspar=''

equilibrium=[]

fifthcoef='y'

considerC3='n'

orthogonal='n'

cod2='n'

alphaval='y'

plot2d='n'

parameters_on_axes=[]

names_of_parameters=[]

intervalx=[]

intervaly=[]

lines_to_search={}

plot3d='n'

parameters_on_axes3=[]

name_of_extra_parameter=''

extrainterval=[]

extraturing=[]