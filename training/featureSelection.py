import os
import numpy as np
import shap
from xgboost import XGBClassifier

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


print("We have that number of features : ", X.shape[1])
print("Feature selection is starting...")




def shapped_data(shap_values,X_,threshold):
    feature_list = []
    for class_ in range(11):
                shap_values_class = shap_values[class_]
                featureImportance = {}

                for features in range(X_.shape[1]):

                        feature_name = features
                        feature_value = abs(shap_values_class[:,features]).mean()
                        if (feature_value >threshold):
                            if feature_name not in feature_list:
                                feature_list.append(feature_name)
    return feature_list

if not os.path.exists('tfeatures'):
    os.makedirs('tfeatures')

featureFile = open("tfeatures/Features.txt","w") 

print("Our first classifier is XGBOOST")
print("===========Start XGBOOST Classifier Section===========")

featureWriter = []

tree = 100
depth = 2

from sklearn.utils import shuffle
X, y = shuffle(X, y, random_state=5)
 
xgbModel = XGBClassifier(n_estimators= tree,max_depth= depth, objective='multi:softprob')
 
result = xgbModel.fit(X, y)
             
shap_values_train = shap.TreeExplainer(result).shap_values(X)
 
 
featureList = shapped_data(shap_values_train,X,0.1)
 
featureWriter.append(str(featureList))
 
featureFile.writelines(featureWriter)
featureFile.close() 

print("Important features were saved in tfeatures/Features.txt")