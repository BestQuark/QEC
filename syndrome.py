import numpy as np
import matplotlib.pyplot as plt
import itertools
import pickle

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

def classify_errors_small(stabList):
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

def dump_chunk_errors(errorList, iterat):
    with open("{iterat}_error.p".format(iterat=str(iterat)), 'wb') as fp:
        pickle.dump(errorList, fp, protocol=pickle.HIGHEST_PROTOCOL)
        
def load_chunk_errors(iterat):
    with open("{iterat}_error.p".format(iterat=str(iterat)), 'rb') as fp:
        chunkErrors = pickle.load(fp)
    return chunkErrors

def classify_errors_chunks(stabList,  chunkSize = 50000):
    """
    Input: stabilizerList
    Output: (dict) with key syndrome and value list of all posible errors in symplectic form
    """
    class_errors = {}
    symplectic_matrix = symplectic_stabilizer(stabList)

    numQubits = symplectic_matrix[0].shape[1]
    numStabilizers = symplectic_matrix[0].shape[0]
    
    perm1 = itertools.product([0,1], repeat=numStabilizers)

    for error in perm1:
        class_errors[tuple(error)] = []
    
    perm2 = itertools.product([0, 1], repeat=numQubits)
    oneError = []
    for i in perm2:
        oneError.append(np.array(i))
    
    errors = []
    perm3 = itertools.product(oneError, oneError)
    iteration = 0
    for x,y in perm3:
        errors.append([x,y])
        if len(errors)>=chunkSize:
            iteration+=1
            dump_chunk_errors(errors, iteration)
            errors = []
                        
    while iteration>=0:
        
        for error in errors:
            syndrome = error_syndrome(symplectic_matrix, error)

            class_errors[tuple(syndrome)].append(error)

        if iteration!=0:
            errors = load_chunk_errors(iteration)
        iteration-=1

    return class_errors

def classify_errors_weight(stabList, weights):
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
    
    for weight in weights:
        oneError = [np.array(i) for i in itertools.product([0, 1], repeat=numQubits)]
        errors = [[x,y] for x,y in itertools.product(oneError, oneError) if np.sum(np.maximum(x,y))==weight]

        for error in errors:
            syndrome = error_syndrome(symplectic_matrix, error)
            class_errors[tuple(syndrome)].append(error)

    return class_errors
    
    
def classify_errors(stabList, chunkSize = 0, weights = []):
    """
    Input: Symplectic stabilizer list
    Output:(dict) Classification of errors requested
    """
    
    if chunkSize == 0 and (not weights):
        return classify_errors_small(stabList)
    
    elif chunkSize>0 and (not weights):
        return classify_errors_chunks(stabList,  chunkSize = chunkSize)
    
    elif chunkSize==0 and (not not weights):
        return classify_errors_weight(stabList, weights)
    else:
        print("Classification not possible")
        
def minimum_weight_error(errors):
    """
    Input: dictionary with errors in syndromes
    Output: dictionary with lowest weight error in syndromes
    """
    minWeight = {}
    for syndr in errors:
        if not errors[syndr]:
            continue
        minW = errors[syndr][0][0].shape[0]
        detected = []
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
    try:
        dist = commutingOps[0][0].shape[0]
    except IndexError:
        print("Cannot calculate the distance of the code because the normalizer of the stabilizer is empty")
        return -999
    for op in commutingOps:
        x,z = op
        weight = np.sum(np.maximum(x,z))
        
        if weight!=0 and weight<dist:
            dist = weight
    return dist


def fidelity(minWeightErrors, probX = 1/3, probY = 1/3, probZ = 1/3):
    """
    Input: dictionary of minWeightErrors (only works with one encoded qubit)
    Output: calculates fidelity of state after error correction on channel (1- probX - probY - probZ) I + probX (X rho X) + probY (Y rho Y)     probZ (Z rho Z)
    """
    probI = 1 - (probX + probY + probZ) 
    fid = 0
    for error in minWeightErrors:
        if len(minWeightErrors[error])==1:
            probFix = 1
            x,z = minWeightErrors[error][0]        
            for qubit in range(x.shape[0]):
                if x[qubit]==0 and z[qubit]==0:
                    probFix = probFix*probI
                elif x[qubit]==1 and z[qubit]==0:
                    probFix = probFix*probX
                elif x[qubit]==0 and z[qubit]==1:
                    probFix = probFix*probZ
                elif x[qubit]==1 and z[qubit]==1:
                    probFix = probFix*probY
            fid += probFix
    return fid

def save_syndromes(nombreArchivo, errorSyndrome) -> None:
    """
    Input: (str, errorSyndrome) saves syndrome in three files, min weights, allErrors and properties
    """
    minErrors = minimum_weight_error(errorSyndrome)
    minErrorsTxt = syndromes_symplectic2txt(minErrors)
    distance = distance_of_code(errorSyndrome)
    
    f = open("{fnombre}_Properties.txt".format(fnombre=nombreArchivo), "w")
    f.write("Distance of code: "+ str(distance))
    f.write("\nCorrectable Errors: ")
    for i in minErrorsTxt:
        f.write("\nSyndrome: "+str(i)+ "; Error: " + str(minErrorsTxt[i]))
    f.close()
    
    with open("{fnombre}_AllErrors.p".format(fnombre=str(nombreArchivo)), 'wb') as fp:
        pickle.dump(errorSyndrome, fp, protocol=pickle.HIGHEST_PROTOCOL)
        
    with open("{fnombre}_MinWeightErrors.p".format(fnombre=str(nombreArchivo)), 'wb') as fp:
        pickle.dump(minErrors, fp, protocol=pickle.HIGHEST_PROTOCOL)
        
def open_syndromes(nombreArchivo):
    """
    Input: (str) name of code 
    Output: error syndrome and minimum weight errors
    """
    with open("{fnombre}_AllErrors.p".format(fnombre=str(nombreArchivo)), 'rb') as fp:
        all_errors = pickle.load(fp)
        
    with open("{fnombre}_MinWeightErrors.p".format(fnombre=str(nombreArchivo)), 'rb') as fp2:
        minWeight_errors = pickle.load(fp2)
    return all_errors, minWeight_errors









