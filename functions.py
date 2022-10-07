class Vector:
    '''This class defines a vector with two components. A dummy version of it and
    its actual coordinates'''
    def __init__(self,name):
        self.dummy=[]
        self.actualcoord=[]
        for varnum in range(nvar):
            self.dummy.append(symbols(name + '_' + str(varnum)))
            self.actualcoord.append(symbols(name + '_' + str(varnum)))
        self.dummy=Matrix(self.dummy)
            
class matrix:
    '''This class defines a sympy Matrix with two components. A dummy version of it and
    its actual coordinates. Here you can provide the actual matrix to be saves into the
    latter one'''
    def __init__(self,name,ncol=nvar,refmatrix=ones(nvar)):
        self.dummy=zeros(ncol)
        if refmatrix!=ones(nvar):
            self.actualcoord=refmatrix
        for row in range(ncol):
            for col in range(ncol):
                if refmatrix[row,col]!=0:
                    self.dummy[row,col]=symbols(name + '_' + str(row+1) + '^' + str(col+1))
                    
def isfloat(value):
    '''This function tells you whether a string is a float number or not'''
    try:
      float(eval(value))
      return True
    except:
      return False
                    
def crossorderapplied(u,parameter,equilibrium):
    '''This function is a coded version of F_{1,1}'''
    SS=zeros(nvar,1)
    firstordereval=firstorderderivatives
    for counter1 in range(nvar):
        for counter2 in range(nvar):
            firstordereval[counter1]=Matrix(firstordereval[counter1]).subs(var[counter2], equilibrium[counter2])
    for counter1 in range(nvar):
        SS=Add(SS,Mul(u.dummy[counter1], \
                      diff(firstordereval[counter1],parameter)))
    return SS
    
def secondorderapplied(u,v):
    '''This function is a coded version of F_2'''
    DS=zeros(nvar,1)
    for counter1 in range(nvar):
        for counter2 in range(nvar):
            DS=Add(DS,Mul(u.dummy[counter1],v.dummy[counter2], \
                          secondorderderivatives[counter1][counter2]))
    return DS
    
def thirdorderapplied(u,v,w):
    '''This function is a coded version of F_3'''
    TS=zeros(nvar,1)
    for counter1 in range(nvar):
        for counter2 in range(nvar):
            for counter3 in range(nvar):
                TS=Add(TS,Mul(u.dummy[counter1],v.dummy[counter2],w.dummy[counter3], \
                              thirdorderderivatives[counter1][counter2][counter3]))
    return TS

def fourthorderapplied(u,v,w,r):
    '''This function is a coded version of F_4'''
    Q4S=zeros(nvar,1)
    for counter1 in range(nvar):
        for counter2 in range(nvar):
            for counter3 in range(nvar):
                for counter4 in range(nvar):
                    Q4S=Add(Q4S,Mul(u.dummy[counter1],v.dummy[counter2],w.dummy[counter3],r.dummy[counter4], \
                                    fourthorderderivatives[counter1][counter2][counter3][counter4]))
    return Q4S

def fifthorderapplied(u,v,w,r,s):
    '''This function is a coded version of F_5'''
    Q5S=zeros(nvar,1)
    for counter1 in range(nvar):
        for counter2 in range(nvar):
            for counter3 in range(nvar):
                for counter4 in range(nvar):
                    for counter5 in range(nvar):
                        Q5S=Add(Q5S,Mul(u.dummy[counter1],v.dummy[counter2], \
                                        w.dummy[counter3],r.dummy[counter4],s.dummy[counter5], \
                                        fifthorderderivatives[counter1][counter2][counter3][counter4][counter5]))
    return Q5S

def dummyvareval(vector,negativeRHS,coefmat):
    '''This function evaluates all the dummy variables used to solve a linear
    system in a simple way'''
    for row in range(nvar):
        vector=vector.subs(negativeRHS.dummy[row],negativeRHS.actualcoord[row])
        for col in range(nvar):
            vector=vector.subs(coefmat.dummy[row,col],coefmat.actualcoord[row,col])
    return vector

def linearsolver(vector,negativeRHS,coefmat):
    '''This function solves a linear system with dummy variables and then uses
    the previous function to evaluate the dummy variables'''
    vector.actualcoord=linsolve(Add(Mul(coefmat.dummy,vector.dummy),negativeRHS.dummy),list(vector.dummy))
    vector.actualcoord=transpose(Matrix(list(vector.actualcoord)))
    vector.actualcoord=dummyvareval(vector.actualcoord,negativeRHS,coefmat)
    return vector.actualcoord