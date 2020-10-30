import os
import numpy as np
import numpy as np
from sklearn.multiclass import OneVsRestClassifier
from sklearn.svm import SVC
from sklearn.utils import shuffle


def map2binary(): # Map 2 and -2 to 1 and 1,-1,0 to 0 for copy states
    
    xstate = np.load('dataVectors/xState.npy') #Call state vector
    
    #Order is important! If firstly 2 and -2 will map, then that 1's also will be turn to 0.
    
    xstate[(xstate == 1) | (xstate == -1)] = 0 #Map 1 and 1,-1,0 to 0
    xstate[(xstate == 2) | (xstate == -2)] = 1 #Map 2 and -2 to 1
    
    
    return xstate
    

def map2binary_(): # Map 2, -2, 1, -1 to 1 -> It means all values except 0 will be 1
    
    xstate = np.load('dataVectors/xState.npy') #Call state vector
    xstate = np.where(xstate != 0,1,xstate) #Map all values except 0 to 1
    return xstate

def map2binary3(): # {-1, 0,1} where CNV -2, -1 are mapped to -1 and CNV 1, 2  are mapped to 1        
    
    xstate = np.load('dataVectors/xState.npy') #Call state vector
    
    
    xstate[(xstate == -2)] = -1 #Map -2 to -1
    xstate[(xstate == 2)] = 1 #Map 2 to 1
    
    
    return xstate

def map2binary4(): # {-1, 0,1}  where CNV -2 mapped to -1 and CNV 2 mapped to 1. All others -1,0,1 are mapped to 0.
    
    xstate = np.load('dataVectors/xState.npy') #Call state vector
    #Order is important
    xstate[(xstate == -1)] = 0 #Map -1 to 0
    xstate[(xstate == 1)] = 0 #Map 1 to 0
    
    xstate[(xstate == -2)] = -1 #Map 1 and 1,-1,0 to 0
    xstate[(xstate == 2)] = 1 #Map 2 and -2 to 1
    
    
    return xstate
    


def prepareInput(stateVec): #Concatenate mutation and copy state vectors
    mutVec = np.load('dataVectors/xMutation.npy') #Call mutation vector
    
    return np.concatenate((mutVec,stateVec),axis=1)
 

def getInputVectors(mapType):
    
    yTrain = np.load('dataVectors/yTrain.npy') #Call label vector
    
    if (mapType == "map1"):
        
        stateVec = map2binary()
        xTrain = prepareInput(stateVec)
       
        return xTrain, yTrain
        
    elif (mapType == "map2"):
        
        stateVec = map2binary_()
        xTrain = prepareInput(stateVec)
        
        return xTrain, yTrain
    
    elif (mapType == "map3"):
        
        stateVec = map2binary3()
        xTrain = prepareInput(stateVec)
        
        return xTrain, yTrain
    
    elif (mapType == "map4"):
        
        stateVec = map2binary4()
        xTrain = prepareInput(stateVec)
        
        return xTrain, yTrain
        
        
    else:
        raise Exception("Invalid map type! map1 or map2 can be used.")
    

print("We will use map: where CNV features -2, -1 are mapped to -1 and CNV features 1, 2  are mapped to 1     ")


#mapType = input("Map type: ")
mapType = "map3"
X,y = getInputVectors(mapType)

print()
print("Data is ready to use!")
X, y = shuffle(X, y, random_state=5)



def extractOVRModel(ovr):
    if not os.path.exists('modelCoef'):
        os.makedirs('modelCoef')
    np.savetxt("modelCoef/coef.csv", ovr.coef_, delimiter=",",fmt="%-.8f",)
    np.savetxt("modelCoef/rho.csv", ovr.intercept_, delimiter=",",fmt="%-.8f",)
    


path = os.getcwd()
dataPath = path
featuresFile2 = dataPath+"/tfeatures/Features.txt"

with open(featuresFile2, "r") as f:

            for line in f:
                    line = line[1:len(line)-1]
                    line = line.split(',')
                    featureArray = np.asarray(line,dtype='int')
                    
featureArray = list(np.sort(np.asarray(featureArray)))
X_train = X[:,featureArray]
print("We have that number of features: ",X_train.shape[1])
y_train = y
svmModel = OneVsRestClassifier(SVC(kernel='linear',probability=True))
result = svmModel.fit(X_train, y_train)
extractOVRModel(svmModel)

print("Model parameters were saved in modelCoef/coef.csv and modelCoef/rho.csv")

