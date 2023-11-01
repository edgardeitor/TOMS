parameters={"a": 0.0310302,
    "b": 1.00311,
    "c": 0.01,
    "eta": 0.05,
    "epsilon": 0.025,
    "lambda_0": 0,
    "delta": 0.2,
    "delta_2": 0.05,
    "mu": 0.19,
    "mu_2": 1}

var=['u','v','w','z']

diffmatrix=[['delta*lambda_0 + delta_2**2*(1 - lambda_0)', 0, 0, 0],
            [0, 1, 0, 0],
            [0, 0, 'delta*lambda_0*mu + delta_2**2*mu_2*(1 - lambda_0)', 0],
            [0, 0, 0, 'lambda_0*mu + mu_2*(1 - lambda_0)']]

kinetics=[
    "lambda_0*(epsilon*(-u + w) + eta*(a + u**2*v - u*(b + 1))) + (1 - lambda_0)*(c*(-u + w) + eta*(a + u**2*v - u*(b + 1)))",
    "lambda_0*(epsilon*(-v + z) + eta*(b*u - u**2*v)) + (1 - lambda_0)*(c*(-v + z) + eta*(b*u - u**2*v))",
    "lambda_0*(epsilon*(u - w) + eta*(a + w**2*z - w*(b + 1))) + (1 - lambda_0)*(c*(u - w) + eta*(a + w**2*z - w*(b + 1)))",
    "lambda_0*(epsilon*(v - z) + eta*(b*w - w**2*z)) + (1 - lambda_0)*(c*(v - z) + eta*(b*w - w**2*z))"
]

tol=1e-7

parameter_functions={}

phiunit='n'

crosscoef='n'

numberofzerotemporalderivatives=0

crosspar=''

fifthcoef='n'

considerC3='n'

orthogonal='n'

cod2='y'

alphaval='n'

plot2d='y'

parameters_on_axes=['b','a']

names_of_parameters=[]

intervalx=[0,2]

intervaly=[0,8]

lines_to_search={'b':1.00311}

plot3d='y'

parameters_on_axes3=['b','lambda_0','a']

name_of_extra_parameter=''

extrainterval=[0,1]

extraturing=[0,0.5,1]