import numpy as np
import matplotlib.pyplot as plt
import itertools


def symplectic_stabilizer(stabList):
    """
    Input: (str) list of stabilizers in the form "XXI;IXX;XXX"
    Output: ([np.array, np.array]) symplectic representation matrix
    """
    stab = stabList.split(";")
    stab = [i.replace(" ","") for i in stab]
    
    #test
    lenStabs = [len(i) for i in stab]
    eqLens = all(i==lenStabs[0] for i in lenStabs)
    if eqLens==False:
        print("Stab list no equal length")
        return 0
    
    matrixX = np.zeros((len(stab),lenStabs[0]))
    matrixZ = np.zeros((len(stab),lenStabs[0]))

    for indxgen, gen in enumerate(stab):
        for indxop,op in enumerate(gen):
            if op=='X':
                matrixX[indxgen, indxop] = 1
            elif op=='Z':
                matrixZ[indxgen, indxop] = 1
            elif op=='Y':
                matrixX[indxgen, indxop] = 1
                matrixZ[indxgen, indxop] = 1
    
    return [matrixX, matrixZ]

def symplectic_innerProd(op1, op2):
    """
    Input: opi = [binary vector, binary vector]
    Output: (0 or 1) symplectic innerProduct between op1 and op2
    """
    op1 = np.asarray(op1)
    op2 = np.asarray(op2)
    return np.sum(op1[0]*op2[1] + op1[1]*op2[0])%2

def error_syndrome(stabList, error):
    """
    Input:  symplectic matrix of stabilizer
    Output: np.array with error syndrome:
    """
    error_synd = np.zeros(stabList[0].shape[0])
    for i in range(stabList[0].shape[0]):
        error_synd[i] = symplectic_innerProd([stabList[0][i], stabList[1][i]], error)
    
    return error_synd.astype(int)

def classify_errors(stabList):
    """
    Input: stabilizerList
    Output: (dict) with key syndrome and value list of all posible errors in symplectic form
    """
    class_errors = {}
    symplectic_matrix = symplectic_stabilizer(stabList)
    
    numQubits = symplectic_matrix[0].shape[1]
    numStabilizers = symplectic_matrix[0].shape[0]
    
    for error in itertools.product([0,1], repeat=numStabilizers):
        class_errors[tuple(error)] = []
    
    oneError = [np.array(i) for i in itertools.product([0, 1], repeat=numQubits)]
    errors = [[x,y] for x,y in itertools.product(oneError, oneError)]
    
    for error in errors:
        syndrome = error_syndrome(symplectic_matrix, error)
        class_errors[tuple(syndrome)].append(error)

    
    return class_errors
    

def classify_errors_weight(stabList, weight):
    """
    Input: stabilizerList
    Output: (dict) with key syndrome and value list of all posible errors in symplectic form
    """
    class_errors = {}
    symplectic_matrix = symplectic_stabilizer(stabList)
    
    numQubits = symplectic_matrix[0].shape[1]
    numStabilizers = symplectic_matrix[0].shape[0]
    
    for error in itertools.product([0,1], repeat=numStabilizers):
        class_errors[tuple(error)] = []
    
    oneError = [np.array(i) for i in itertools.product([0, 1], repeat=numQubits)]
    errors = [[x,y] for x,y in itertools.product(oneError, oneError) if np.sum(np.maximum(x,y))==weight]
    
    for error in errors:
        syndrome = error_syndrome(symplectic_matrix, error)
        class_errors[tuple(syndrome)].append(error)

    
    return class_errors
    
    
def minimum_weight_error(errors):
    """
    Input: dictionary with errors in syndromes
    Output: dictionary with lowest weight error in syndromes
    """
    minWeight = {}
    for syndr in errors:
        minW = errors[syndr][0][0].shape[0]
        detected = 0
        for error in errors[syndr]:
            x,y = error
            #weight of error
            w = np.sum(np.maximum(x,y))
            if w<minW:
                minW = w
                detected = [error]
            elif w==minW:
                detected.append(error)
        minWeight[syndr] = detected
    return minWeight 

def symplectic2txt(vector):
    """
    Input: a symplectic vector
    Output: txt operators
    """
    x,z = vector
    operator = ""
    for i in range(x.shape[0]):
        if x[i]==1 and z[i]==0:
            operator+="X"
        elif x[i]==0 and z[i]==1:
            operator+="Z"
        elif x[i]==1 and z[i]==1:
            operator+="Y"
        elif x[i]==0 and z[i]==0:
            operator+="I"
            
    return operator

def syndromes_symplectic2txt(syndrom_dict):
    """
    Input: dictionary of syndromes in symplectic
    Output: dictionary of syndromes in txt
    """
    syndromes = syndrom_dict.copy()
    for i in syndromes:
        if len(syndromes[i])==1:
            syndromes[i] = symplectic2txt(*syndromes[i])
        else:
            newtxt = []
            for j in syndromes[i]:
                newtxt.append(symplectic2txt(j))
            syndromes[i] = newtxt
    return syndromes

def distance_of_code(syndromes):
    """
    Input: dict of syndromes in symplectic form
    Output: distance of code
    """
    commutingOps = syndromes[list(syndromes.keys())[0]]
    
    dist = commutingOps[0][0].shape[0]
    for op in commutingOps:
        x,z = op
        weight = np.sum(np.maximum(x,z))
        
        if weight!=0 and weight<dist:
            dist = weight
    return dist

