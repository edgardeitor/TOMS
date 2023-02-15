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

'''The script runs the file 'data.py'.'''

modelname=input('Enter the name of the system you would like to analyze: ')

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

file=open('List of Variables.txt','w')

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
    parameters=[]

npar=len(parameters)

if npar>0:
    newparameters=dict()
    for key in parameters.keys():
        try:
            exec(key + ' = symbols(key, real=True)')
            newparameters[eval(key)]=parameters[key]
        except:
            print('The script could not define your variable ' + key + ' as a variable')
            exit()
    parameters=newparameters
    
'''It defines the diffusion matrix as a symbolic matrix in terms of the parameters of the system.'''
    
try:
    diffmatrix
except:
    print('Diffusion matrix was not provided.')
    exit()
    
try:
    for row in range(nvar):
        for col in range(nvar):
            diffmatrix[row][col]=eval(str(diffmatrix[row][col]).replace('^','**'))
except:
    print('The diffusion matrix is not a function of the parameters of the system')
    exit()
        
diffmatrix=Matrix(diffmatrix)

# Kinetics

'''The script defines the list of kinetics as a symbolic matrix in terms of the parameters of the system.'''

try:
    kinetics
except:
    print('Kinetics were not provided.')
    exit()
    
for functionnumber in range(nvar):
    try:
        kinetics[functionnumber]=eval(str(kinetics[functionnumber]).replace('^','**'))
    except:
        print('The expression ' + kinetics[functionnumber] + 'is not a function of the parameters of your system')
        exit()

kinetics=Matrix(kinetics)

'''The script writes the list of kinetics functions in LaTeX in a text file.'''

file=open('Kinetics.txt','w')
for functionnumber in range(nvar):
    file.write(latex(kinetics[functionnumber]) + '\n')
file.close()

jacobianmat=kinetics.jacobian(var)

# Equilibria

'''The script evaluates the matrix of kinetics and searches for the equilibria at the parameters provided.'''

kineticsevaluated=kinetics
for key in parameters:
    kineticsevaluated=kineticsevaluated.subs(key,parameters[key])
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
        for key in parameters:
            kineticsevaluated=kineticsevaluated.subs(key,parameters[key])
        eq=solve(kineticsevaluated,var)
    else:
        break
    
'''If there is more than one real equilibrium, the script asks you to choose one to analyze.
If there is only one equilibrium, then the script picks it automatically.'''

if isinstance(eq,list) and len(eq)>1:
    print('Your system has ' + str(len(eq)) + ' real equilibria for the given parameter values given by:')
    for eqnum in range(len(eq)):
        print(str(eqnum+1) + '.- ' + str(eq[eqnum]) + '\n')
    eqnumber=input('Enter the number of the equilibrium point you want to consider: ')
    while not eqnumber.isnumeric() or int(eqnumber)==0:
        eqnumber=input('What you entered before is not a positive integer. ' +
                       'Enter the number of the equilibrium point you want to consider: ')
    eqnumber=int(eqnumber)
    eq=Matrix(list(eq[eqnumber-1]))
elif isinstance(eq,list) and len(eq)==1:
    print('Your system has only one real equilibrium point given by:')
    eqnumber=0
    print(eq[0])
    eq=Matrix(eq[0])

# Turing conditions

'''The script defines symbols and standard matrices that will be used to find C3.'''

muNF = symbols('mu_NF')

eigval = symbols('eigval', real=True)

coefmat0 = matrix('coef0mat', nvar, jacobianmat)
coefmat1 = matrix('coef1mat', nvar, Add(jacobianmat, Mul(-1, muNF, diffmatrix), Mul(-1, eigval, eye(nvar))))
coefmat2 = matrix('coef2mat', nvar, Add(jacobianmat, Mul(-1, 4, muNF, diffmatrix)))
coefmat3 = matrix('coef3mat', nvar, Add(jacobianmat, Mul(-1, 9, muNF, diffmatrix)))

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
            determinanteval=determinanteval.subs(var[varnum],eq[varnum])
            derivativeeval=derivativeeval.subs(var[varnum],eq[varnum])
        except:
            determinanteval=determinanteval.subs(var[varnum],eq[var[varnum]])
            derivativeeval=derivativeeval.subs(var[varnum],eq[var[varnum]])
    for key in parameters:
        determinanteval=determinanteval.subs(key,parameters[key])
        derivativeeval=derivativeeval.subs(key,parameters[key])
    getout=0
    
    '''It tries to check whether there is a value of \mu that solves both equations
    that define the Turing bifurcation at the same time.'''
    
    try:
        mucritical=solve(derivativeeval,muNF)
        for muvalue in mucritical:
            muvalue=complex(muvalue).real
            if abs(N(determinanteval.subs(muNF,muvalue)))<5*tol:
                ksquared=muvalue
                getout=1
                break
    except:
        pass
    
    '''If there is no Turing bifurcation, the script asks you to provide a parameter that can be
    changed in order to find the bifurcation.'''
    
    if getout==0:
        determinanteval=jacobianmatdet
        derivativeeval=determinantderivative
        print('There is no Turing bifurcation for the parameters provided.')
        print('The parameters of the system are:')
        for key in parameters.keys():
            print(key)
        parametertochange=input('Enter a parameter that can be changed in order to find the wavenumber ' +
                                '(Enter 0 if you do not want Python to find the bifurcation condition): ')
        counter=0
        while True:
            
            '''The script checks whether the equations change after changing the parameter provided.'''
            
            try:
                if parametertochange=='0':
                    ksquared=0
                    getout=1
                    break
                if counter==1:
                    parametertochange=input('Enter a parameter that can be changed in order to find the wavenumber ' +
                                            '(Enter 0 if you do not want Python to find the bifurcation condition): ')
                parametertochange=eval(parametertochange)
                if parametertochange not in parameters.keys():
                    print('You did not enter a parameter of the system.')
                    continue
                else:
                    if (determinanteval==determinanteval.subs(parametertochange,parameters[parametertochange])
                        and derivativeeval==derivativeeval.subs(parametertochange,parameters[parametertochange])):
                            print('Your parameter does not produce any changes in the ' +
                                  'equations for the Turing bifurcation.')
                            continue
                break
            except:
                counter=1
                
        if getout==0:
            kineticsevaluated=kinetics
            
            '''The script evaluates the kinetics to find the new equilibrium.'''
            
            for key in parameters:
                if key!=parametertochange:
                    kineticsevaluated=kineticsevaluated.subs(key,parameters[key])
                    determinanteval=determinanteval.subs(key,parameters[key])
                    derivativeeval=derivativeeval.subs(key,parameters[key])
                else:
                    initialpartochange=parameters[key]
            initialeq=eq
            eq=solve(kineticsevaluated,var)[eqnumber]
            for varnum in range(nvar):
                determinanteval=determinanteval.subs(var[varnum],eq[varnum])
                derivativeeval=derivativeeval.subs(var[varnum],eq[varnum])
                
            '''The script tries to find a solution to both equations numerically starting from \mu=1
            and the parameter value provided previously.'''
                
            try:
                zerofunction=[lambda mutofind, parametertofind:
                              determinanteval.subs([(muNF,mutofind),(parametertochange,parametertofind)]),
                                  lambda mutofind, parametertofind:
                                      derivativeeval.subs([(muNF,mutofind),(parametertochange,parametertofind)])]
                [muvalue,newparval]=findroot(zerofunction,(1,initialpartochange))
                muvalue=complex(muvalue).real
                newparval=complex(newparval).real
                if (abs(simplify(determinanteval.subs([(muNF,muvalue),(parametertochange,newparval)])))<tol
                    and abs(simplify(derivativeeval.subs([(muNF,muvalue),(parametertochange,newparval)])))<tol):
                    ksquared=muvalue
                    parameters[parametertochange]=newparval
                    print('The Turing bifurcation was found for ' + str(parametertochange) + '=' + str(newparval))
                    getout=1
            except:
                pass
        if getout==0:
            
            '''The script tries to find the solution by first solving one equation to get \mu
            in terms of the parameter and then solves the other equation to find it.'''
            
            try:
                mucritical=solve(derivativeeval,muNF)
                for muvalue in mucritical:
                    try:
                        newparval=float(findroot(lambda parametertofind:
                                                 determinanteval.subs([(muNF,muvalue),
                                                                       (parametertochange,parametertofind)]),
                                                     initialpartochange))
                        muvalue=complex(simplify(muvalue.subs(parametertochange,newparval))).real
                        newparval=complex(newparval).real
                        if simplify(determinanteval.subs([(muNF,muvalue),(parametertochange,newparval)]))<tol:
                            ksquared=muvalue
                            parameters[parametertochange]=newparval
                            print('The Turing bifurcation was found for ' + str(parametertochange) + '='
                                  + str(newparval))
                            getout=1
                            break
                    except:
                        continue
            except:
                pass
            
        if getout==0:
            
            '''The script finds the resultant between the functions that define the equations and then it finds the
            value of the parameter by solving resultant=0 with the value provided preiously as an
            initial ondition.'''
            
            try:
                intersection=resultant(determinanteval, derivativeeval, muNF)
                newparval=float(findroot(lambda parametertofind:
                                         intersection.subs(parametertochange,parametertofind),initialpartochange))
                mucritical=solve(derivativeeval.subs(parametertochange,newparval),muNF)
                newparval=complex(newparval).real
                for muvalue in mucritical:
                    muvalue=complex(muvalue).real
                    if abs(simplify(determinanteval.subs([(muNF,muvalue),(parametertochange,newparval)])))<tol:
                        ksquared=muvalue
                        parameters[parametertochange]=newparval
                        print('The Turing bifurcation was found for ' + str(parametertochange) + '=' + str(newparval))
                        getout=1
                        break
            except:
                pass
            
            '''If none of the previous methods worked, then the code will ask you to provide a wavenumber.'''
            
        if getout==0:
            eq=initialeq
            kval=input('This script could not find the value of k. ' +
                       'Make sure that you have a Turing bifurcation for the parameter values provided. ' +
                           'If you are at a Turing bifurcation. Enter a value of k you want to consider: ')
            while not isfloat(kval):
                kval=input('What you entered before is not a number. ' +
                           'Enter a value of k you want to consider: ')
            kval=eval(kval)
            ksquared=Pow(kval,2)
else:
    ksquared=0
    
'''The script saves the functions to find a bifurcation into text files.'''
    
file = open('Determinant of the Jacobian matrix.txt', 'w')
file.write(latex(jacobianmatdet))
file.close()
file = open('Derivative of the Determinant.txt', 'w')
file.write(latex(determinantderivative))
file.close()
file = open("Non-zero variable.txt",'w')
file.write(latex(nonzero))
file.close()

# Normal form

'''The script defines a vector that will be used to find the solutions to all the linear systems
of equations and finds the derivatives of the kinetics up to order five.'''

negativeRHS=Vector('negativeRHS')

firstorderderivatives=list()
secondorderderivatives=list()
thirdorderderivatives=list()
fourthorderderivatives=list()
fifthorderderivatives=list()
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

phiNF = Vector('phi^NF')
Q02NF = Vector('Q02^NF')
Q22NF = Vector('Q22^NF')
psiNF = Vector('psi^NF')
Q13NF = Vector('Q13^NF')
Q33NF = Vector('Q33^NF')
Q04NF = Vector('Q04^NF')
Q24NF = Vector('Q24^NF')

getout=0

'''The script looks for an invertible (n-1) x (n-1) submatrix of the Jacobian matrix of the system.'''

if eq!=[]:
    for row in range(nvar):
        for col in range(nvar):
            submatrixrows=list(range(nvar))
            submatrixcols=list(range(nvar))
            submatrixrows.remove(row)
            submatrixcols.remove(col)
            invertiblesubmatrix=coefmat1.actualcoord.extract(submatrixrows,submatrixcols)
            submatrixeval=invertiblesubmatrix
            for key in parameters:
                submatrixeval=submatrixeval.subs(key,parameters[key])
            if eq!=[]:
                for varnum in range(nvar):
                    submatrixeval=submatrixeval.subs(var[varnum],eq[varnum])
            submatrixeval=submatrixeval.subs(muNF,ksquared)
            if abs(N(submatrixeval.det()))>tol:
                phiNF.actualcoord[col]=1
                criticalrow=row
                criticalcol=col
                getout=1
                break
        if getout==1:
            break
else:
    submatrixrows = list(range(nvar))
    submatrixcols = list(range(nvar))
    submatrixrows.remove(0)
    submatrixcols.remove(0)
    invertiblesubmatrix=coefmat1.actualcoord.extract(submatrixrows,submatrixcols)
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
    phiNF.actualcoord = Mul(Pow(sqrt(phiNF.actualcoord.dot(phiNF.actualcoord)),-1), phiNF.actualcoord)
    
hatphiNF = Vector('hatphiNF')
hatphiNF.dummy = Matrix(phiNF.dummy[0:nvar-numberofzerotemporalderivatives]
                        + [0]*numberofzerotemporalderivatives)

phiNF_eval = evaluation_dict(phiNF)

print('First-order ready')

'''The script solves the linear equations to find the second-order vectors.'''

DS_phiphi = secondorderapplied(phiNF, phiNF)

negativeRHS.actualcoord = DS_phiphi

Q02NF = linearsolver(Q02NF, negativeRHS, coefmat0)

Q22NF = linearsolver(Q22NF, negativeRHS, coefmat2)

Q02NF_eval = evaluation_dict(Q02NF)
Q22NF_eval = evaluation_dict(Q22NF)
        
print('Second-order ready')

'''The script defines a few terms to find the third-order coefficient.'''

DS_phiQ02 = secondorderapplied(phiNF, Q02NF)
DS_phiQ22 = secondorderapplied(phiNF, Q22NF)

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

C3 = Mul(Pow(denominator, -1), psiNF.dummy.dot(Add(Mul(4, DS_phiQ02), Mul(2, DS_phiQ22),
                                                   Mul(3, TS_phiphiphi))))

'''The script gets the cross-order coefficient if requested.'''

if crosscoef=='y':
    try:
        for varnum in range(nvar):
            equilibrium[varnum] = eval(equilibrium[varnum])
            
        crosspar = eval(crosspar)
            
        SS_phi = crossorderapplied(phiNF, crosspar, equilibrium)
        
        C11 = Mul(Pow(hatphiNF.dummy.dot(psiNF.dummy),-1),
                  psiNF.dummy.dot(Add(SS_phi, Mul(-1, muNF, diff(diffmatrix, crosspar), phiNF.dummy))))
        
        for varnum in range(nvar):
            C11 = C11.subs(phiNF.dummy[varnum], phiNF.actualcoord[varnum])
            C11 = C11.subs(psiNF.dummy[varnum], psiNF.actualcoord[varnum])
        
        file = open('Cross-order coefficient.txt','w')        
        file.write(latex(C11))
        file.close()
    except:
        print('The equilibrium you introduced is not defined in terms of the parameters of the system' +
              'or the parameter crosspar is not a parameter of the system.')

'''The script continues with the following orders to find C5 if requested.'''

if fifthcoef=='y':
    Q13NF.actualcoord[criticalcol] = 0
    
    '''The script finds Q_1^3 using an analogous approach to the one used to find \phi_1^1.'''
    
    negativeRHS.actualcoord = Add(Mul(-1, C3, hatphiNF.dummy), Mul(4, DS_phiQ02), Mul(2, DS_phiQ22),
                                  Mul(3, TS_phiphiphi))
    
    auxiliaryterm, = linsolve(Add(Mul(coefsubmatrix.dummy,Matrix(Q13NF.actualcoord).extract(submatrixcols,[0])),
                                  negativeRHS.dummy.extract(submatrixrows,[0])),
                              list(Matrix(Q13NF.actualcoord).extract(submatrixcols,[0])))
        
    Q13NF.actualcoord[0:criticalcol] = auxiliaryterm[0:criticalcol]
    Q13NF.actualcoord[criticalcol+1:nvar] = auxiliaryterm[criticalcol:nvar-1]
    
    Q13NF.actualcoord = Matrix(Q13NF.actualcoord)
    
    for row in range(nvar):
        Q13NF.actualcoord = Q13NF.actualcoord.subs(negativeRHS.dummy[row], negativeRHS.actualcoord[row])
        for col in range(nvar-1):
            if row<nvar-1:
                Q13NF.actualcoord = Q13NF.actualcoord.subs(coefsubmatrix.dummy[row, col],
                                                           coefsubmatrix.actualcoord[row, col])
    
    '''The script orthogonalizes Q_1^3 with respect to \phi_1^1 if requested.'''
    
    if orthogonal=='y':
        if phiunit=='y':
            Q13NF.actualcoord = Add(Q13NF.actualcoord, Mul(-1, Q13NF.actualcoord.dot(phiNF.dummy), phiNF.dummy))
        else:
            Q13NF.actualcoord = Add(Q13NF.actualcoord, Mul(-1, Q13NF.actualcoord.dot(phiNF.dummy),
                                                           Pow(phiNF.dummy.dot(phiNF.dummy), -1), phiNF.dummy))
        
    negativeRHS.actualcoord = Add(Mul(2, DS_phiQ22), TS_phiphiphi)
    
    Q33NF = linearsolver(Q33NF, negativeRHS, coefmat3) # The script finds Q_3^3
    
    Q13NF_eval = evaluation_dict(Q13NF)
    Q33NF_eval = evaluation_dict(Q33NF)
    
    print('Third-order ready')
    
    '''The code finds the fourth-order vectors.'''
        
    DS_Q02Q02 = secondorderapplied(Q02NF, Q02NF)
    DS_Q22Q22 = secondorderapplied(Q22NF, Q22NF)
    DS_phiQ13 = secondorderapplied(phiNF, Q13NF)
    
    DS_Q02Q22 = secondorderapplied(Q02NF, Q22NF)
    DS_phiQ33 = secondorderapplied(phiNF, Q33NF)
    
    TS_phiphiQ02 = thirdorderapplied(phiNF, phiNF, Q02NF)
    TS_phiphiQ22 = thirdorderapplied(phiNF, phiNF, Q22NF)
    
    Q4S_phiphiphiphi = fourthorderapplied(phiNF, phiNF, phiNF, phiNF)
    
    hatQ02NF = Vector('hatQ02NF')
    hatQ02NF.dummy = Matrix(Q02NF.dummy[0:nvar-numberofzerotemporalderivatives]
                            + [0]*numberofzerotemporalderivatives)
    
    negativeRHS.actualcoord = Add(Mul(-2, C3, hatQ02NF.dummy), Mul(2, DS_Q02Q02),
                                  DS_Q22Q22, Mul(2, DS_phiQ13), Mul(6, TS_phiphiQ02),
                                  Mul(3, TS_phiphiQ22), Mul(3, Q4S_phiphiphiphi))
    
    Q04NF = linearsolver(Q04NF, negativeRHS, coefmat0)
    
    hatQ22NF = Vector('hatQ22NF')
    hatQ22NF.dummy = Matrix(Q22NF.dummy[0:nvar-numberofzerotemporalderivatives]
                            + [0]*numberofzerotemporalderivatives)
    
    negativeRHS.actualcoord = Add(Mul(-2, C3, hatQ22NF.dummy), Mul(4, DS_Q02Q22), Mul(2, DS_phiQ13),
                                  Mul(2, DS_phiQ33), Mul(6, TS_phiphiQ02), Mul(6, TS_phiphiQ22),
                                  Mul(4, Q4S_phiphiphiphi))
    
    Q24NF = linearsolver(Q24NF, negativeRHS, coefmat2)
    
    Q04NF_eval = evaluation_dict(Q04NF)
    Q24NF_eval = evaluation_dict(Q24NF)
        
    print('Fourth-order ready')
    
    '''The script computes C5.'''
    
    DS_phiQ04 = secondorderapplied(phiNF, Q04NF)
    DS_phiQ24 = secondorderapplied(phiNF, Q24NF)
    DS_Q02Q13 = secondorderapplied(Q02NF, Q13NF)
    DS_Q22Q13 = secondorderapplied(Q22NF, Q13NF)
    DS_Q22Q33 = secondorderapplied(Q22NF, Q33NF)
    
    TS_phiphiQ13 = thirdorderapplied(phiNF, phiNF, Q13NF)
    TS_phiphiQ33 = thirdorderapplied(phiNF, phiNF, Q33NF)
    TS_phiQ02Q02 = thirdorderapplied(phiNF, Q02NF, Q02NF)
    TS_phiQ02Q22 = thirdorderapplied(phiNF, Q02NF, Q22NF) 
    TS_phiQ22Q22 = thirdorderapplied(phiNF, Q22NF, Q22NF)
    
    Q4S_phiphiphiQ02 = fourthorderapplied(phiNF, phiNF, phiNF, Q02NF)
    Q4S_phiphiphiQ22 = fourthorderapplied(phiNF, phiNF, phiNF, Q22NF)
    
    Q5S_phiphiphiphiphi = fifthorderapplied(phiNF, phiNF, phiNF, phiNF, phiNF)
    
    hatQ13NF = Vector('hatQ13NF')
    hatQ13NF.dummy = Matrix(Q13NF.dummy[0:nvar-numberofzerotemporalderivatives]
                            + [0]*numberofzerotemporalderivatives)

    C5 = Mul(Pow(denominator, -1),
             psiNF.dummy.dot(Add(Mul(-3, C3, hatQ13NF.dummy), Mul(4, DS_phiQ04), Mul(2, DS_phiQ24),
                                 Mul(4, DS_Q02Q13), Mul(2, DS_Q22Q13), Mul(2, DS_Q22Q33),
                                 Mul(9, TS_phiphiQ13), Mul(3, TS_phiphiQ33), Mul(12, TS_phiQ02Q02),
                                 Mul(12, TS_phiQ02Q22), Mul(6, TS_phiQ22Q22), Mul(24, Q4S_phiphiphiQ02),
                                 Mul(16, Q4S_phiphiphiQ22), Mul(10, Q5S_phiphiphiphiphi))))
    
    print('The calculation of the fifth-order coefficient has been carried out successfully. ' +
          'The saving process could take longer.')
    
C3 = C3.subs(Q02NF_eval)
C3 = C3.subs(Q22NF_eval)

C3 = C3.subs(phiNF_eval)
C3 = C3.subs(psiNF_eval)
    
file=open('Third-order coefficient.txt','w')
file.write(latex(C3))
file.close()

print('The third-order coefficient was computed and saved into a text file.')
    
if fifthcoef=='n':
    if crosscoef=='y':
        print('Everything but the fifth-order coefficient was computed and saved.')
    else:
        print('Everything but the fifth-order and cross-order coefficients was computed and saved.')
if fifthcoef=='y':
    C5 = C5.subs(Q04NF_eval)
    C5 = C5.subs(Q24NF_eval)
    
    C5 = C5.subs(Q13NF_eval)
    C5 = C5.subs(Q33NF_eval)
    
    C5 = C5.subs(Q02NF_eval)
    C5 = C5.subs(Q22NF_eval)
    
    C5 = C5.subs(phiNF_eval)
    C5 = C5.subs(psiNF_eval)
    
    file = open('Fifth-order coefficient.txt','w')
    file.write(latex(C5))
    file.close()
    
    print('The fifth-order coefficient was computed and saved into a text file.')
    
file=open('Find codimension-two bifurcation points.txt','w')
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
        file=open('Initial conditions for Turing bifurcation curves.txt','w')
        
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
    
    file=open('Fixed parameter values.txt','w')
    for key in parameters.keys():
        if key not in parameters_on_axes:
            if key not in parameter_functions.keys():
                file.write(latex(key) + ',' + latex(parameters[key]) + '\n')
            else:
                file.write(latex(key) + ',' + latex(parameter_functions[key]) + '\n')
            
    file.close()
        
    file=open('Parameters on axes.txt','w')
    file.write(latex(parameters_on_axes[0]) + ',' + latex(intervalx[0]) + ',' + latex(intervalx[1]) + '\n')
    file.write(latex(parameters_on_axes[1]) + ',' + latex(intervaly[0]) + ',' + latex(intervaly[1]) + '\n')
    file.close()

    try:
        if len(names_of_parameters)==0:
            for parnum in range(2):
                names_of_parameters[parnum] = latex(parameters_on_axes[parnum])
    except:
        names_of_parameters=parameters_on_axes
        for parnum in range(2):
            names_of_parameters[parnum] = latex(names_of_parameters[parnum])
        
    file=open('Actual names of parameters.txt','w')
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
        file=open('Extra Turing curves in codimension-two bifurcation diagram.txt','w')
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