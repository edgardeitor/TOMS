from IPython import get_ipython
get_ipython().run_line_magic('reset', '-sf')

import os
import shutil
import math
from sys import exit
from sympy import *
init_printing()
from sympy.solvers import solve
from mpmath import findroot

'''The script runs the file 'name_of_model.py'.'''

modelname = input('Enter the name of the system you would like to analyze: ')

if not os.path.isdir(modelname):
    print('The directory of that system was not found. Create it first and place ' + 
          'the data file inside named in the same way.')
    exit()
else:
    os.chdir(os.getcwd() + '\\' + modelname)

try:
    exec(open(modelname + '.py').read())
except:
    print('The file ' + modelname + '.py could not be run')
    exit()
    
nvar = len(var)

'''The script runs the file 'functions.py' '''

try:
    exec(open(os.path.dirname(os.path.realpath(__file__)) + '\\functions.py').read())
except:
    print('File functions.py is not in the same folder as the script you are running')
    exit()

file = open('List of Variables.txt', 'w')

'''The script defines the vector of variables as symbols in sympy.'''

try:
    var
except:
    print('Variables were not provided')
    exit()


for varnum in range(nvar):
    try:
        exec(var[varnum] + ' = symbols(var[varnum], real=True)')
        var[varnum]=eval(var[varnum])
    except:
        print('The script could not define ' + (var[varnum]) + ' as a variable')
        exit()
    file.write(latex(var[varnum]) + '\n')
    
'''The script defines the parameters as symbols in sympy.'''

try:
    parameters
except:
    print('Parameters were not provided. The script will assume that there are no parameters.')
    parameters = []

npar=len(parameters)

if npar>0:
    newparameters = dict()
    for key in parameters.keys():
        try:
            exec(key + ' = symbols(key, real=True)')
            newparameters[eval(key)] = parameters[key]
        except:
            print('The script could not define your variable ' + key + ' as a variable')
            exit()
    parameters = newparameters
    
'''It defines the diffusion matrix as a symbolic matrix in terms of the parameters of the system.'''
    
try:
    diffmatrix
except:
    print('Diffusion matrix was not provided.')
    exit()
    
try:
    for row in range(nvar):
        for col in range(nvar):
            diffmatrix[row][col] = eval(str(diffmatrix[row][col]).replace('^','**'))
except:
    print('The diffusion matrix is not a function of the parameters of the system')
    exit()
        
diffmatrix = Matrix(diffmatrix)

# Kinetics

'''The script defines the list of kinetics as a symbolic matrix in terms of the parameters of the system.'''

try:
    kinetics
except:
    print('Kinetics were not provided.')
    exit()
    
for functionnumber in range(nvar):
    try:
        kinetics[functionnumber] = eval(str(kinetics[functionnumber]).replace('^','**'))
    except:
        print('The expression ' + kinetics[functionnumber] + 'is not a function of the parameters of your system')
        exit()

kinetics = Matrix(kinetics)

'''The script writes the list of kinetics functions in LaTeX in a text file.'''

file = open('Kinetics.txt','w')
for functionnumber in range(nvar):
    file.write(latex(kinetics[functionnumber]) + '\n')
file.close()

jacobianmat = kinetics.jacobian(var)

# Equilibria

'''The script evaluates the matrix of kinetics and searches for the equilibria at the parameters provided.'''

kineticsevaluated = kinetics
kineticsevaluated = kineticsevaluated.subs(parameters)
eq=solve(kineticsevaluated,var)

'''If the script does not find any equilibrium point, it asks you whether you would like to 
change a parameter to find it and, if yes, then you must provide a new parameter value.'''

while eq==[]:
    print('This script could not find a real equilibrium point for the given parameter values shown below: ')
    print(parameters)
    parameterchange=input('Would you like to change a parameter value? [y/n] ')
    while parameterchange!='y' and parameterchange!='n':
        parameterchange=input('You must enter your answer in the format [y/n]. ' +
                              'Would you like to change a parameter value? ')
    if parameterchange=='y':
        whichparam=input('Which parameter would you like to change? ')
        while True:
            try:
                whichparam=eval(whichparam)
                while whichparam not in parameters.keys():
                    whichparam=input('That is not a parameter of the system. ' +
                                     'Which parameter would you like to change? ')
                break
            except:
                whichparam=input('That is not a parameter of the system. Which parameter would you like to change? ')
        parameters[whichparam]=input('Enter a value of ' + str(whichparam) + ': ')
        while not isfloat(parameters[whichparam]):
            parameters[whichparam]=input('What you entered before is not a number. ' +
                                         'Enter a value of ' + whichparam + ': ')
        parameters[whichparam]=eval(parameters[whichparam])
        kineticsevaluated=kinetics
        kineticsevaluated = kineticsevaluated.subs(parameters)
        eq=solve(kineticsevaluated,var)
    else:
        break
    
'''If there is more than one real equilibrium, the script asks you to choose one to analyze.
If there is only one equilibrium, then the script picks it automatically.'''

if isinstance(eq,list) and len(eq)>1:
    print('Your system has ' + str(len(eq)) + ' real equilibria for the given parameter values given by:')
    for eqnum in range(len(eq)):
        print(str(eqnum+1) + '.- ' + str(eq[eqnum]) + '\n')
    eqnumber = input('Enter the number of the equilibrium point you want to consider: ')
    while not eqnumber.isnumeric() or int(eqnumber)==0:
        eqnumber=input('What you entered before is not a positive integer. ' +
                       'Enter the number of the equilibrium point you want to consider: ')
    eqnumber = int(eqnumber)
    eq = Matrix(list(eq[eqnumber - 1]))
elif isinstance(eq,list) and len(eq)==1:
    print('Your system has only one real equilibrium point given by:')
    eqnumber=0
    print(eq[0])
    eq=Matrix(eq[0])

# Turing conditions

'''The script defines symbols and standard matrices that will be used to find C3.'''

muNF = symbols('mu_NF')

eigval = symbols('eigval', real=True)

for counter in range(4):
    if counter==1:
        coefmat1 = matrix('coef1mat', nvar, Add(jacobianmat, Mul(-1, muNF, diffmatrix), Mul(-1, eigval, eye(nvar))))
    else:
        exec(f'coefmat{counter} = matrix("coef{counter}mat", nvar, Add(jacobianmat, Mul(- Pow(counter, 2), muNF, diffmatrix)))')

'''The Jacobian matrix can be complicated algebraically so the script computes the determinant of a
dummy matrix to simplify the process.'''

jacobianmatdet = coefmat1.dummy.det()

'''The script then replaces the actual entries of the dummy matrix in the expression of the determinant.'''

for row in range(nvar):
    for col in range(nvar):
        jacobianmatdet = jacobianmatdet.subs(coefmat1.dummy[row, col], coefmat1.actualcoord[row, col])

nonzero = diff(jacobianmatdet, eigval).subs(eigval, 0)
coefmat1.actualcoord = coefmat1.actualcoord.subs(eigval, 0)
jacobianmatdet = jacobianmatdet.subs(eigval, 0)

determinantderivative = diff(jacobianmatdet, muNF)

determinanteval = jacobianmatdet
derivativeeval = determinantderivative

'''If the script finds at least one equilibrium, it evaluates the determinant and the Jacobian
matrix at the parameters provided.'''

if eq!=[]:
    for varnum in range(nvar):
        try:
            determinanteval = determinanteval.subs(var[varnum], eq[varnum])
            derivativeeval = derivativeeval.subs(var[varnum], eq[varnum])
        except:
            determinanteval = determinanteval.subs(var[varnum], eq[var[varnum]])
            derivativeeval = derivativeeval.subs(var[varnum], eq[var[varnum]])
    determinanteval = determinanteval.subs(parameters)
    derivativeeval = derivativeeval.subs(parameters)
    getout = 0
    ksquared = nan
    
    '''It checks whether there is a value of \mu that solves both equations
    that define the Turing bifurcation at the same time. If not, it sets a
    value of k=1 to continue the calculation'''
    
    try:
        mucritical = solve(derivativeeval, muNF)
        for muvalue in mucritical:
            muvalue = complex(muvalue).real
            if abs(N(determinanteval.subs(muNF, muvalue)))<5*tol:
                ksquared = muvalue
                getout = 1
                break
            if muvalue==mucritical[-1]:
                print('There is no Turing bifurcation for the parameter values you provided. ' +
                      'The calculation will be carried out symbolically.')
    except:
        pass
    
'''The script saves the functions to find a bifurcation into text files.'''
    
file = open('Determinant of the Jacobian matrix.txt', 'w')
file.write(latex(jacobianmatdet))
file.close()
file = open('Derivative of the Determinant.txt', 'w')
file.write(latex(determinantderivative))
file.close()
file = open("Non-zero variable.txt", 'w')
file.write(latex(nonzero))
file.close()

# Normal form

'''The script defines a vector that will be used to find the solutions to all the linear systems
of equations and finds the derivatives of the kinetics up to order five.'''

negativeRHS = Vector('negativeRHS')

firstorderderivatives = list()
secondorderderivatives = list()
thirdorderderivatives = list()
fourthorderderivatives = list()
fifthorderderivatives = list()
for counter1 in range(nvar):
    firstorderderivatives.append(diff(kinetics,var[counter1]))
    secondorderderivatives.append(list())
    thirdorderderivatives.append(list())
    fourthorderderivatives.append(list())
    fifthorderderivatives.append(list())
    for counter2 in range(nvar):
        secondorderderivatives[counter1].append(diff(firstorderderivatives[counter1],var[counter2]))
        thirdorderderivatives[counter1].append(list())
        fourthorderderivatives[counter1].append(list())
        fifthorderderivatives[counter1].append(list())
        for counter3 in range(nvar):
            thirdorderderivatives[counter1][counter2].append(diff(secondorderderivatives[counter1][counter2],
                                                                  var[counter3]))
            fourthorderderivatives[counter1][counter2].append(list())
            fifthorderderivatives[counter1][counter2].append(list())
            for counter4 in range(nvar):
                fourthorderderivatives[counter1][counter2][counter3].append(
                    diff(thirdorderderivatives[counter1][counter2][counter3], var[counter4]))
                fifthorderderivatives[counter1][counter2][counter3].append(list())
                for counter5 in range(nvar):
                    fifthorderderivatives[counter1][counter2][counter3][counter4].append(
                        diff(fourthorderderivatives[counter1][counter2][counter3][counter4], var[counter5]))

phiNF = Vector('phiNF')
W02NF = Vector('W02NF')
W22NF = Vector('W22NF')
psiNF = Vector('psiNF')
W13NF = Vector('W13NF')
W33NF = Vector('W33NF')
W04NF = Vector('W04NF')
W24NF = Vector('W24NF')

getout = 0

'''The script looks for an invertible (n-1) x (n-1) submatrix of the Jacobian matrix of the system.'''

if eq!=[] and ksquared!=nan:
    for row in range(nvar):
        for col in range(nvar):
            submatrixrows = list(range(nvar))
            submatrixcols = list(range(nvar))
            submatrixrows.remove(row)
            submatrixcols.remove(col)
            invertiblesubmatrix = coefmat1.actualcoord.extract(submatrixrows, submatrixcols)
            submatrixeval = invertiblesubmatrix
            submatrixeval = submatrixeval.subs(parameters)
            if eq!=[]:
                for varnum in range(nvar):
                    submatrixeval = submatrixeval.subs(var[varnum], eq[varnum])
                submatrixeval = submatrixeval.subs(muNF, ksquared)
                if abs(N(submatrixeval.det()))>tol:
                    phiNF.actualcoord[col] = 1
                    criticalrow = row
                    criticalcol = col
                    getout = 1
                    break
            else:
                phiNF.actualcoord[0] = 1
                criticalrow = 0
                criticalcol = 0
                break
        if getout==1:
            break
else:
    submatrixrows = list(range(nvar))
    submatrixcols = list(range(nvar))
    submatrixrows.remove(0)
    submatrixcols.remove(0)
    invertiblesubmatrix=coefmat1.actualcoord.extract(submatrixrows, submatrixcols)
    phiNF.actualcoord[0] = 1
    criticalrow = 0
    criticalcol = 0
    
coefsubmatrix = matrix('dummysubmatrix',nvar-1,invertiblesubmatrix)

'''The script defines the vector \phi_1^1 that spans the kernel of the Jacobian matrix with dummy matrices.'''

auxiliaryterm, = linsolve(Add(Mul(coefsubmatrix.dummy, Matrix(phiNF.actualcoord).extract(submatrixcols, [0])),
                              coefmat1.dummy.extract(submatrixrows, [criticalcol])),
                          list(Matrix(phiNF.actualcoord).extract(submatrixcols, [0])))
    
phiNF.actualcoord[0:criticalcol] = auxiliaryterm[0:criticalcol]
phiNF.actualcoord[criticalcol+1:nvar] = auxiliaryterm[criticalcol:nvar-1]

phiNF.actualcoord = Matrix(phiNF.actualcoord)
    
'''The script replaces the dummy variables by the actual values.'''

for row in range(nvar):
    for col in range(nvar):
        phiNF.actualcoord = phiNF.actualcoord.subs(coefmat1.dummy[row,col],coefmat1.actualcoord[row,col])
        if row<nvar-1 and col<nvar-1:
            phiNF.actualcoord = phiNF.actualcoord.subs(coefsubmatrix.dummy[row,col],coefsubmatrix.actualcoord[row,col])

'''The script normalizes the \phi_1^1 if requested.'''

if phiunit=='y':
    phiNF.actualcoord = Mul(Pow(sqrt(phiNF.actualcoord.dot(phiNF.actualcoord)), -1), phiNF.actualcoord)
    
hatphiNF = Vector('hatphiNF')
hatphiNF.dummy = Matrix(phiNF.dummy[0:nvar-numberofzerotemporalderivatives]
                        + [0]*numberofzerotemporalderivatives)

phiNF_eval = evaluation_dict(phiNF)

print('First-order ready')

'''The script solves the linear equations to find the second-order vectors.'''

DS_phiphi = secondorderapplied(phiNF, phiNF)

negativeRHS.actualcoord = DS_phiphi

W02NF = linearsolver(W02NF, negativeRHS, coefmat0)

W22NF = linearsolver(W22NF, negativeRHS, coefmat2)

W02NF_eval = evaluation_dict(W02NF)
W22NF_eval = evaluation_dict(W22NF)
        
print('Second-order ready')

'''The script defines a few terms to find the third-order coefficient.'''

DS_phiW02 = secondorderapplied(phiNF, W02NF)
DS_phiW22 = secondorderapplied(phiNF, W22NF)

TS_phiphiphi = thirdorderapplied(phiNF, phiNF, phiNF)
    
Tdummycoefsubmatrix = transpose(coefsubmatrix.dummy)
Tdummycoefmat1 = transpose(coefmat1.dummy)

psiNF.actualcoord[criticalrow] = 1

'''The script finds the vector \psi_0^1 that spans the Kernel of the adjoint of the Jacobian matrix.'''

auxiliaryterm, = linsolve(Add(Mul(Tdummycoefsubmatrix,Matrix(psiNF.actualcoord).extract(submatrixrows, [0])),
                              Tdummycoefmat1.extract(submatrixcols, [criticalrow])),
                          list(Matrix(psiNF.actualcoord).extract(submatrixrows, [0])))
    
psiNF.actualcoord[0:criticalrow] = auxiliaryterm[0:criticalrow]
psiNF.actualcoord[criticalrow+1:nvar] = auxiliaryterm[criticalrow:nvar-1]

psiNF.actualcoord = Matrix(psiNF.actualcoord)

for row in range(nvar):
    for col in range(nvar):
        psiNF.actualcoord = psiNF.actualcoord.subs(coefmat1.dummy[row,col],coefmat1.actualcoord[row,col])
        if row<nvar-1 and col<nvar-1:
            psiNF.actualcoord = psiNF.actualcoord.subs(coefsubmatrix.dummy[row,col],coefsubmatrix.actualcoord[row,col])

psiNF_eval = evaluation_dict(psiNF)

'''The script gets the third-order oefficient.'''

denominator = hatphiNF.dummy.dot(psiNF.dummy)

C3 = Mul(Pow(denominator, -1), psiNF.dummy.dot(Add(Mul(4, DS_phiW02), Mul(2, DS_phiW22),
                                                   Mul(3, TS_phiphiphi))))

'''The script gets the cross-order coefficient if requested.'''

if crosscoef=='y':
    try:            
        crosspar = eval(crosspar)
        
        pre_C11 = Mul(Pow(hatphiNF.dummy.dot(psiNF.dummy), -1),
                  psiNF.dummy.dot(Mul(Add(jacobianmat, Mul(-1, muNF, diffmatrix)), phiNF.dummy)))
        
        # for varnum in range(nvar):
        #     if varnum<nvar - 1:
        #         file.write(latex(phiNF.actualcoord[varnum]) + ',')
        #     else:
        #         file.write(latex(phiNF.actualcoord[varnum]) + '\n')
        # for varnum in range(nvar):
        #     if varnum<nvar - 1:
        #         file.write(latex(psiNF.actualcoord[varnum]) + ',')
        #     else:
        #         file.write(latex(psiNF.actualcoord[varnum]))
        # file.close()
        
        file = open('To get cross-order coefficient.txt', 'w')
        file.write(latex(crosspar) + '\n')
        file.write(latex(pre_C11))
        file.close()
    except:
        print('The parameter crosspar is not a parameter of the system.')

'''The script continues with the following orders to find C5 if requested.'''

if fifthcoef=='y':
    '''The script finds W_1^3 using an analogous approach to the one used to find \phi_1^1.'''
    
    negativeRHS.actualcoord = Add(Mul(4, DS_phiW02), Mul(2, DS_phiW22), Mul(3, TS_phiphiphi))
    
    if considerC3=='y':
        negativeRHS.actualcoord = Add(negativeRHS.actualcoord, Mul(-1, C3, hatphiNF.dummy))
    
    W13NF = critical_linearsolver(W13NF, negativeRHS, criticalcol, coefsubmatrix, submatrixrows, submatrixcols)
    
    '''The script orthogonalizes W_1^3 with respect to \phi_1^1 if requested.'''
    
    if orthogonal=='y':
        if phiunit=='y':
            W13NF.actualcoord = Add(W13NF.actualcoord, Mul(-1, W13NF.actualcoord.dot(phiNF.dummy), phiNF.dummy))
        else:
            W13NF.actualcoord = Add(W13NF.actualcoord, Mul(-1, W13NF.actualcoord.dot(phiNF.dummy),
                                                           Pow(phiNF.dummy.dot(phiNF.dummy), -1), phiNF.dummy))
        
    negativeRHS.actualcoord = Add(Mul(2, DS_phiW22), TS_phiphiphi)
    
    W33NF = linearsolver(W33NF, negativeRHS, coefmat3) # The script finds W_3^3
    
    W13NF_eval = evaluation_dict(W13NF)
    W33NF_eval = evaluation_dict(W33NF)
    
    print('Third-order ready')
    
    '''The code finds the fourth-order vectors.'''
        
    DS_W02W02 = secondorderapplied(W02NF, W02NF)
    DS_W22W22 = secondorderapplied(W22NF, W22NF)
    DS_phiW13 = secondorderapplied(phiNF, W13NF)
    
    DS_W02W22 = secondorderapplied(W02NF, W22NF)
    DS_phiW33 = secondorderapplied(phiNF, W33NF)
    
    TS_phiphiW02 = thirdorderapplied(phiNF, phiNF, W02NF)
    TS_phiphiW22 = thirdorderapplied(phiNF, phiNF, W22NF)
    
    Q4S_phiphiphiphi = fourthorderapplied(phiNF, phiNF, phiNF, phiNF)
    
    hatW02NF = Vector('hatW02NF')
    hatW02NF.dummy = Matrix(W02NF.dummy[0:nvar-numberofzerotemporalderivatives]
                            + [0]*numberofzerotemporalderivatives)
    
    negativeRHS.actualcoord = Add(Mul(2, DS_W02W02), DS_W22W22, Mul(2, DS_phiW13), Mul(6, TS_phiphiW02),
                                  Mul(3, TS_phiphiW22), Mul(3, Q4S_phiphiphiphi))
    
    if considerC3=='y':
        negativeRHS.actualcoord = Add(negativeRHS.actualcoord, Mul(-2, C3, hatW02NF.dummy))
    
    W04NF = linearsolver(W04NF, negativeRHS, coefmat0)
    
    hatW22NF = Vector('hatW22NF')
    hatW22NF.dummy = Matrix(W22NF.dummy[0:nvar-numberofzerotemporalderivatives]
                            + [0]*numberofzerotemporalderivatives)
    
    negativeRHS.actualcoord = Add(Mul(4, DS_W02W22), Mul(2, DS_phiW13), Mul(2, DS_phiW33),
                                  Mul(6, TS_phiphiW02), Mul(6, TS_phiphiW22), Mul(4, Q4S_phiphiphiphi))
    
    if considerC3=='y':
        negativeRHS.actualcoord = Add(negativeRHS.actualcoord, Mul(-2, C3, hatW22NF.dummy))
    
    W24NF = linearsolver(W24NF, negativeRHS, coefmat2)
    
    W04NF_eval = evaluation_dict(W04NF)
    W24NF_eval = evaluation_dict(W24NF)
        
    print('Fourth-order ready')
    
    '''The script computes C5.'''
    
    DS_phiW04 = secondorderapplied(phiNF, W04NF)
    DS_phiW24 = secondorderapplied(phiNF, W24NF)
    DS_W02W13 = secondorderapplied(W02NF, W13NF)
    DS_W22W13 = secondorderapplied(W22NF, W13NF)
    DS_W22W33 = secondorderapplied(W22NF, W33NF)
    
    TS_phiphiW13 = thirdorderapplied(phiNF, phiNF, W13NF)
    TS_phiphiW33 = thirdorderapplied(phiNF, phiNF, W33NF)
    TS_phiW02W02 = thirdorderapplied(phiNF, W02NF, W02NF)
    TS_phiW02W22 = thirdorderapplied(phiNF, W02NF, W22NF) 
    TS_phiW22W22 = thirdorderapplied(phiNF, W22NF, W22NF)
    
    Q4S_phiphiphiW02 = fourthorderapplied(phiNF, phiNF, phiNF, W02NF)
    Q4S_phiphiphiW22 = fourthorderapplied(phiNF, phiNF, phiNF, W22NF)
    
    Q5S_phiphiphiphiphi = fifthorderapplied(phiNF, phiNF, phiNF, phiNF, phiNF)
    
    hatW13NF = Vector('hatW13NF')
    hatW13NF.dummy = Matrix(W13NF.dummy[0:nvar-numberofzerotemporalderivatives]
                            + [0]*numberofzerotemporalderivatives)

    C5 = Mul(Pow(denominator, -1),
             psiNF.dummy.dot(Add(Mul(4, DS_phiW04), Mul(2, DS_phiW24), Mul(4, DS_W02W13),
                                 Mul(2, DS_W22W13), Mul(2, DS_W22W33), Mul(9, TS_phiphiW13),
                                 Mul(3, TS_phiphiW33), Mul(12, TS_phiW02W02), Mul(12, TS_phiW02W22),
                                 Mul(6, TS_phiW22W22), Mul(24, Q4S_phiphiphiW02),
                                 Mul(16, Q4S_phiphiphiW22), Mul(10, Q5S_phiphiphiphiphi))))
    
    if considerC3=='y':
        C5 = Add(C5, Mul(Pow(denominator, -1), psiNF.dummy.dot(Mul(-3, C3, hatW13NF.dummy))))
    
    print('The calculation of the fifth-order coefficient has been carried out successfully.')
    
# C3 = C3.subs(W02NF_eval)
# C3 = C3.subs(W22NF_eval)

# C3 = C3.subs(phiNF_eval)
# C3 = C3.subs(psiNF_eval)

C3 = simplify(C3)

file = open('Third-order coefficient.txt', 'w')
file.write(latex(C3))
file.close()

file = open('Vectors.txt', 'w')
write_vector(phiNF, file)
write_vector(psiNF, file)
write_vector(W02NF, file)
write_vector(W22NF, file)
file.close()

print('The third-order coefficient was computed and saved into a text file.')
    
if fifthcoef=='n':
    if crosscoef=='y':
        print('Everything but the fifth-order coefficient was computed and saved.')
    else:
        print('Everything but the fifth-order and cross-order coefficients was computed and saved.')
if fifthcoef=='y':
    file = open('Vectors.txt', 'a')
    write_vector(W13NF, file)
    write_vector(W33NF, file)
    write_vector(W04NF, file)
    write_vector(W24NF, file)
    file.close()
    
    file = open('Fifth-order coefficient.txt','w')
    file.write(latex(C5))
    file.close()
    
    print('The fifth-order coefficient was computed and saved into a text file.')
    
if alphaval=='y' and fifthcoef=='y':
    W123NF = Vector('W123^NF')
    W024NF = Vector('W024^NF')
    W224NF = Vector('W224^NF')
    
    negativeRHS.actualcoord = Mul(2, sqrt(muNF), diffmatrix, phiNF.dummy)
    
    '''The script finds W_{1,2}^3 using an analogous approach to the one used to find \phi_1^1.'''
    
    W123NF = critical_linearsolver(W123NF, negativeRHS, criticalcol, coefsubmatrix, submatrixrows, submatrixcols)
    
    DS_phiW123 = secondorderapplied(phiNF, W123NF)
    
    negativeRHS.actualcoord = Mul(2, DS_phiW123)
    
    W024NF = linearsolver(W024NF, negativeRHS, coefmat0)
    
    negativeRHS.actualcoord = Add(Mul(2, DS_phiW123), Mul(8, sqrt(muNF), diffmatrix, W22NF.dummy))
    
    W224 = linearsolver(W224NF, negativeRHS, coefmat2)
    
    W123NF_eval = evaluation_dict(W123NF)
    W024NF_eval = evaluation_dict(W024NF)
    W224NF_eval = evaluation_dict(W224NF)
    
    DS_phiW024 = secondorderapplied(phiNF, W024NF)
    DS_phiW224 = secondorderapplied(phiNF, W224NF)
    DS_W02W123 = secondorderapplied(W02NF, W123NF)
    
    DS_W22W123 = secondorderapplied(W22NF, W123NF)
    
    TS_phiphiW123 = thirdorderapplied(phiNF, phiNF, W123NF)
    
    alpha1 = Mul(- 2, sqrt(muNF), Pow(denominator, -1), psiNF.dummy.dot(Mul(diffmatrix, W123NF.dummy)))
    
    alpha2 = Mul(2, Pow(denominator, -1), psiNF.dummy.dot(Add(DS_phiW024, DS_phiW224,
                                                              Mul(2, DS_W02W123), Mul(3, TS_phiphiW123),
                                                              Mul(2, sqrt(muNF), diffmatrix, W13NF.dummy))))
    
    alpha3 = Mul(Pow(denominator, -1), psiNF.dummy.dot(Add(Mul(- 2, DS_phiW024),
                                                          Mul(- 2, DS_W22W123), Mul(- 3, TS_phiphiW123),
                                                          Mul(2, sqrt(muNF), diffmatrix, W13NF.dummy))))
    
    # alpha1 = alpha1.subs(W123NF_eval)
    # alpha1 = alpha1.subs(phiNF_eval)
    # alpha1 = alpha1.subs(psiNF_eval)
    
    # alpha2 = alpha2.subs(W024NF_eval)
    # alpha2 = alpha2.subs(W224NF_eval)
    # alpha2 = alpha2.subs(W13NF_eval)
    # alpha2 = alpha2.subs(W123NF_eval)
    # alpha2 = alpha2.subs(W02NF_eval)
    # alpha2 = alpha2.subs(W22NF_eval)
    # alpha2 = alpha2.subs(phiNF_eval)
    # alpha2 = alpha2.subs(psiNF_eval)
    
    # alpha3 = alpha3.subs(W024NF_eval)
    # alpha3 = alpha3.subs(W13NF_eval)
    # alpha3 = alpha3.subs(W123NF_eval)
    # alpha3 = alpha3.subs(W02NF_eval)
    # alpha3 = alpha3.subs(W22NF_eval)
    # alpha3 = alpha3.subs(phiNF_eval)
    # alpha3 = alpha3.subs(psiNF_eval)
    
    file = open('alpha1.txt', 'w')
    file.write(latex(alpha1))
    file.close()
    
    file = open('alpha2.txt', 'w')
    file.write(latex(alpha2))
    file.close()
    
    file = open('alpha3.txt', 'w')
    file.write(latex(alpha3))
    file.close()
    
file = open('Find codimension-two bifurcation points.txt','w')
try:
    cod2
    file.write(str(cod2))
except:
    file.write('n')
file.close()

''' The sript saves text files with all the information required by Mathematica to find and plot
the bifurcation diagram'''
    
if plot2d=='y':
    try:
        for parnum in range(2):
            parameters_on_axes[parnum]=eval(parameters_on_axes[parnum])
    except:
        print('The variable párameters_on_axes is not well defined.')
        exit()

    try:
        file = open('Initial conditions for Turing bifurcation curves.txt','w')
        
        for key in lines_to_search.keys():
            if isinstance(lines_to_search[key],list):
                for initialsolnum in range(len(lines_to_search[key])):
                    if isfloat(str(lines_to_search[key][initialsolnum])):
                        file.write(latex(eval(key)) + ',' + 
                                   latex(eval(str(lines_to_search[key][initialsolnum]))) + '\n')
            else:
                if isfloat(str(lines_to_search[key])):
                    file.write(latex(eval(key)) + ',' + latex(eval(str(lines_to_search[key]))) + '\n')
                
        file.close()    
    except:
        print('The variable lines_to_search is not well defined.')
        exit()
    
    try:
        auxpar=dict()
        
        for key in parameter_functions.keys():
            auxpar[eval(key)] = eval(parameter_functions[key])
        parameter_functions=auxpar
    except:
        parameter_functions=dict()
    
    file = open('Fixed parameter values.txt', 'w')
    for key in parameters.keys():
        if key not in parameters_on_axes:
            if key not in parameter_functions.keys():
                file.write(latex(key) + ',' + latex(parameters[key]) + '\n')
            else:
                file.write(latex(key) + ',' + latex(parameter_functions[key]) + '\n')
            
    file.close()
        
    file = open('Parameters on axes.txt', 'w')
    file.write(latex(parameters_on_axes[0]) + ',' + latex(intervalx[0]) + ',' + latex(intervalx[1]) + '\n')
    file.write(latex(parameters_on_axes[1]) + ',' + latex(intervaly[0]) + ',' + latex(intervaly[1]) + '\n')
    file.close()

    try:
        if len(names_of_parameters)==0:
            for parnum in range(2):
                names_of_parameters[parnum] = latex(parameters_on_axes[parnum])
    except:
        names_of_parameters = parameters_on_axes
        for parnum in range(2):
            names_of_parameters[parnum] = latex(names_of_parameters[parnum])
        
    file = open('Actual names of parameters.txt','w')
    file.write(names_of_parameters[0] + ',' + names_of_parameters[1])
    file.close()
    
    if not os.path.isfile('Plotter.nb'):
        shutil.copyfile(os.path.dirname(os.path.realpath(__file__)) + '\\Plotter.nb', 'Plotter.nb')
    
    print('The variables to plot the bifurcation diagram in Mathematica were correctly saved')

if plot2d=='y' and plot3d=='y':
    try:
        for parnum in range(3):
            parameters_on_axes3[parnum] = eval(parameters_on_axes3[parnum])
    except:
        print('The variable párameters_on_axes3 is not well defined.')
        exit()
    
    interval3 = []
    counter = 0
    auxpar = parameters_on_axes.copy()
    for parnum in range(2):
        auxpar[parnum] = eval(auxpar[parnum])
    for parnum in range(3):
        if counter>1:
            print('The variable parameter_on_axes3 is not well defined.')
            exit()
        elif parameters_on_axes3[parnum] in auxpar:
            ind = auxpar.index(parameters_on_axes3[parnum])
            if ind==0:
                interval3.append(intervalx)
            elif ind==1:
                interval3.append(intervaly)
        else:
            interval3.append(extrainterval)
            extraparnum = parnum
            counter+= 1
    
    try:
        file = open('Extra Turing curves in codimension-two bifurcation diagram.txt','w')
        for parnum in range(len(extraturing)):
            file.write(latex(extraturing[parnum]) + '\n')
        file.close()
    except:
        pass
            
    try:
        if len(name_of_extra_parameter)==0:
            name_of_extra_parameter=latex(parameters_on_axes3[extraparnum])
        counter = 0
        names_of_parameters3=[]
        for parnum in range(3):
            if parnum!=extraparnum:
                names_of_parameters3.append(names_of_parameters[counter])
                counter+=1
            else:
                names_of_parameters3.append(name_of_extra_parameter)
    except:
        name_of_extra_parameter = latex(parameters_on_axes3[extraparnum])
        counter = 0
        names_of_parameters3 = []
        for parnum in range(3):
            if parnum!=extraparnum:
                names_of_parameters3.append(names_of_parameters[counter])
                counter+=1
            else:
                names_of_parameters3.append(name_of_extra_parameter)
    
    file = open('Codimension-two actual names of parameters.txt','w')
    file.write(names_of_parameters3[0] + ',' + names_of_parameters3[1] + ',' + names_of_parameters3[2])
    file.close()
    
    file = open('Codimension-two parameters on axes.txt','w')
    file.write(latex(parameters_on_axes3[0]) + ',' + latex(interval3[0][0]) + ',' + latex(interval3[0][1]) + '\n')
    file.write(latex(parameters_on_axes3[1]) + ',' + latex(interval3[1][0]) + ',' + latex(interval3[1][1]) + '\n')
    file.write(latex(parameters_on_axes3[2]) + ',' + latex(interval3[2][0]) + ',' + latex(interval3[2][1]) + '\n')
    file.close()
    
    print('The variables to plot the codimension-two bifurcation diagram in Mathematica were correctly saved')