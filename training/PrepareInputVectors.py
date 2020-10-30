import numpy as np
import os


#path = os.getcwd()
dataPath = input("Enter your data path: ") 

#dataPath = path+"/Challenge/"
#sampleName = "/challenge_sample_sites.txt"


blaaderVartxt = "/Bladder_challenge_variants.txt"
breastVartxt = "/Breast_challenge_variants.txt"
bronchVartxt = "/Bronchusandlung_challenge_variants.txt"
cervixVartxt = "/Cervixuteri_challenge_variants.txt"
colonVartxt = "/Colon_challenge_variants.txt"
corpusuterVartxt = "/Corpusuteri_challenge_variants.txt"
kidneyVartxt = "/Kidney_challenge_variants.txt"
liverVartxt = "/Liverandintrahepaticbileducts_challenge_variants.txt"
ovaryVartxt = "/Ovary_challenge_variants.txt"
skinVartxt = "/Skin_challenge_variants.txt"
stomachVartxt = "/Stomach_challenge_variants.txt"



blaaderVarpath = dataPath + blaaderVartxt
breastVarpath = dataPath + breastVartxt
bronchVarpath = dataPath + bronchVartxt
cervixVarpath = dataPath + cervixVartxt
colonVarpath = dataPath + colonVartxt
corpusuterVarpath = dataPath + corpusuterVartxt
kidneyVarpath = dataPath + kidneyVartxt
liverVarpath = dataPath + liverVartxt
ovaryVarpath = dataPath + ovaryVartxt
skinVarpath = dataPath + skinVartxt
stomachVarpath = dataPath + stomachVartxt

def patient2Gene(index,dataPath,patientDict,totalGeneDict,locationName):
    patientGene = {}
    with open(dataPath, "r") as f:
        
        for line in f:
            
            patientId = line.split()[0]
            gene = line.split()[1]
            if gene not in totalGeneDict:
                totalGeneDict[gene] = index # Create gene dictionary
                index += 1
                #print(index)
            if patientId not in patientGene:
                patientGene[patientId] = [gene]
            else:
                if gene not in patientGene[patientId]:
                    patientGene[patientId].append(gene)
        
        patientDict[locationName] = dict(patientGene)
        
    return index




index = 0



patientDict = {}

totalGeneDict = {}
index = patient2Gene(index,blaaderVarpath,patientDict,totalGeneDict,'Bladder')

index = patient2Gene(index,breastVarpath,patientDict,totalGeneDict,'Breast')


index = patient2Gene(index,bronchVarpath,patientDict,totalGeneDict,'Bronchusandlung')


index = patient2Gene(index,cervixVarpath,patientDict,totalGeneDict,'Cervixuteri')


index = patient2Gene(index,colonVarpath,patientDict,totalGeneDict,'Colon')


index = patient2Gene(index,corpusuterVarpath,patientDict,totalGeneDict,'Corpusuteri')


index = patient2Gene(index,kidneyVarpath,patientDict,totalGeneDict,'Kidney')


index = patient2Gene(index,liverVarpath,patientDict,totalGeneDict,'Liverandintrahepaticbileducts')


index = patient2Gene(index,ovaryVarpath,patientDict,totalGeneDict,'Ovary')


index = patient2Gene(index,skinVarpath,patientDict,totalGeneDict,'Skin')


index = patient2Gene(index,stomachVarpath,patientDict,totalGeneDict,'Stomach')




blaaderCNtxt = "/Bladder_challenge_CNs.txt"
breastCNtxt = "/Breast_challenge_CNs.txt"
bronchCNtxt = "/Bronchusandlung_challenge_CNs.txt"
cervixCNtxt = "/Cervixuteri_challenge_CNs.txt"
colonCNtxt = "/Colon_challenge_CNs.txt"
corpusuterCNtxt = "/Corpusuteri_challenge_CNs.txt"
kidneyCNtxt = "/Kidney_challenge_CNs.txt"
liverCNtxt = "/Liverandintrahepaticbileducts_challenge_CNs.txt"
ovaryCNtxt = "/Ovary_challenge_CNs.txt"
skinCNtxt = "/Skin_challenge_CNs.txt"
stomachCNtxt = "/Stomach_challenge_CNs.txt"




blaaderCNpath = dataPath + blaaderCNtxt
breastCNpath = dataPath + breastCNtxt
bronchCNpath = dataPath + bronchCNtxt
cervixCNpath = dataPath + cervixCNtxt
colonCNpath = dataPath + colonCNtxt
corpusuterCNpath = dataPath + corpusuterCNtxt
kidneyCNpath = dataPath + kidneyCNtxt
liverCNpath = dataPath + liverCNtxt
ovaryCNpath = dataPath + ovaryCNtxt
skinCNpath = dataPath + skinCNtxt
stomachCNpath = dataPath + stomachCNtxt




def copyNumberStates2dict(variantsPath,loc2Cn,cnsGeneDict,location):
    patientList = []
    patient2Dict = {}
    geneDict = []
    with open(variantsPath, "r") as f:
            
            firstLine = True
            
            for line in f:
                if firstLine == True:
                    #geneName = line.split()[:2]
                    #geneName =  ' '.join(geneName)
                    patientIds = line.split()[2:]
                    #patientIds.insert(0,geneName)
                    for item in patientIds:
                        patientList.append(item)

                    firstLine = False
                    continue
                
                geneName = line.split()[0]
                copyLevel = line.split()[1:]
                geneDict.append(geneName)
                
                index = 0
                for item in copyLevel:
                    patientId = patientList[index]
                    index += 1
                    
                    if patientId not in patient2Dict:
                        patient2Dict[patientId] = [item]
                    else:
                        patient2Dict[patientId].append(item)
                    
    cnsGeneDict[location] = geneDict
    loc2Cn[location] = dict(patient2Dict)
            
           
        


cnsDict = {}
cnsGeneDict = {}
copyNumberStates2dict(blaaderCNpath,cnsDict,cnsGeneDict,'Bladder')

copyNumberStates2dict(breastCNpath,cnsDict,cnsGeneDict,'Breast')


copyNumberStates2dict(bronchCNpath,cnsDict,cnsGeneDict,'Bronchusandlung')


copyNumberStates2dict(cervixCNpath,cnsDict,cnsGeneDict,'Cervixuteri')


copyNumberStates2dict(colonCNpath,cnsDict,cnsGeneDict,'Colon')


copyNumberStates2dict(corpusuterCNpath,cnsDict,cnsGeneDict,'Corpusuteri')


copyNumberStates2dict(kidneyCNpath,cnsDict,cnsGeneDict,'Kidney')


copyNumberStates2dict(liverCNpath,cnsDict,cnsGeneDict,'Liverandintrahepaticbileducts')


copyNumberStates2dict(ovaryCNpath,cnsDict,cnsGeneDict,'Ovary')


copyNumberStates2dict(skinCNpath,cnsDict,cnsGeneDict,'Skin')


copyNumberStates2dict(stomachCNpath,cnsDict,cnsGeneDict,'Stomach')





def dict2vec(statesDict):
    
    for location in statesDict:
        for patientId in statesDict[location]:
            
            statesDict[location][patientId] = np.array(statesDict[location][patientId],dtype='int8')
    return statesDict




cnsVec = dict2vec(cnsDict)




##################Start Of Select Feature Genes###############




def mutationGene2Index(totalGeneDict): #This is for creating features txt file
    myDictGene2Index = {}
    for gene in totalGeneDict:
        if totalGeneDict[gene] in myDictGene2Index:
            print(gene)
        myDictGene2Index[totalGeneDict[gene]] = gene
    return myDictGene2Index
    




def variationGene2Index(cnsGeneDict):
    myDictGene2Index = {}
    index = 0
    for location in cnsGeneDict:
        geneList = cnsGeneDict[location]
        for i in range(len(geneList)):
            myDictGene2Index[index + len(totalGeneDict)] = geneList[i]
        
            index += 1
        break
    return myDictGene2Index
        





variationGene2Index = variationGene2Index(cnsGeneDict)





mutationGene2Index = mutationGene2Index(totalGeneDict)





##################End Of Select Feature Genes###############





def determineLabel(labelName):
    if labelName == 'Bronchusandlung':
        return 0
    elif labelName == 'Bladder':
        return 1
    elif labelName == 'Colon':
        return 2
    elif labelName == 'Skin':
        return 3
    elif labelName == 'Stomach':
        return 4
    elif labelName == 'Corpusuteri':
        return 5
    elif labelName == 'Liverandintrahepaticbileducts':
        return 6
    elif labelName == 'Ovary':
        return 7
    elif labelName == 'Kidney':
        return 8
    elif labelName == 'Cervixuteri':
        return 9
    elif labelName == 'Breast':
        return 10





def mutationVector(patientDict,totalGeneDict,patientId,location):
    
    totalGeneMutate = len(totalGeneDict)
    
    myVec = np.zeros((1,totalGeneMutate))
    
    for gene in patientDict[location][patientId]:
            index = totalGeneDict[gene]
            myVec[0,index] = 1
            
    return myVec





def stateVector(cnsDict,location,patientId):
    
    locationName = list(cnsDict.keys())[0]
    totalGeneStates = len(cnsDict[locationName][list(cnsDict[locationName].keys())[0]])
    
    if patientId not in cnsDict[location]:
        
            return np.zeros((1,totalGeneStates))
        
    else:
            cnsVec = dict2vec(cnsDict)
            myVec = cnsVec[location][patientId]
            return np.reshape(myVec,(1,myVec.shape[0]))
        
        





def xVector(patientDict,totalGeneDict,cnsDict,patientId,location,state):

    if state == "Mutation":
        return mutationVector(patientDict,totalGeneDict,patientId,location)
    
    elif state == "State":
        return stateVector(cnsDict,location,patientId)
    
   





def data2vector(patientDict,totalGeneDict,cnsDict):
    firstPass = True

    for location in patientDict:
        for patientId in patientDict[location]:
            #PatientOrder.append(patientId)
            xMutation = xVector(patientDict,totalGeneDict,cnsDict,patientId,location,"Mutation")
            xStates = xVector(patientDict,totalGeneDict,cnsDict,patientId,location,"State")

            label = determineLabel(location)
            y = np.asarray([label],dtype='int8')

            if firstPass == True:
                yVec = y
                xmutVec = xMutation
                xstateVec = xStates
                firstPass = False
            else:
                xmutVec = np.concatenate((xmutVec,xMutation),axis=0)
                xstateVec = np.concatenate((xstateVec,xStates),axis=0)
                yVec = np.concatenate((yVec,y),axis=0)
        
                
    return xmutVec,xstateVec,yVec





xmutVec, xstateVec, yVec = data2vector(patientDict,totalGeneDict,cnsDict)



if not os.path.exists('dataVectors'):
    os.makedirs('dataVectors')

np.save('dataVectors/xMutation', xmutVec) # save the file as "outfile_name.npy" 
np.save('dataVectors/xState', xstateVec) # save the file as "outfile_name.npy" 
np.save('dataVectors/yTrain', yVec) # save the file as "outfile_name.npy"


print("Input vectors were created!")
