parameters={'a': 0.7591863073854013,
            'b': 0.5,
            'c': 1,
            'd': 0,
            'h': 0,
            'delta':0.07,
            'h_1':0.5,
            'h_2':-0.001}

var=['u','v']

diffmatrix=[['delta^2', 'h_1'],
            ['h_2', 1]]

kinetics=['a - c*u + d*v + u**2*v',
          'b - d*v + h*u - u**2*v']

tol=1e-7

parameter_functions={}

phiunit='n'

crosscoef='n'

numberofzerotemporalderivatives=0

crosspar=''

fifthcoef='y'

considerC3='n'

orthogonal='n'

cod2='y'

alphaval='n'

plot2d='y'

parameters_on_axes=['b','a']

names_of_parameters=[]

intervalx=[0,10]

intervaly=[0,5]

lines_to_search={'b':0.5}

plot3d='n'

parameters_on_axes3=[]

name_of_extra_parameter=''

extrainterval=[]

extraturing=[]