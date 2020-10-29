# In[]:
import numpy as np
from sklearn.metrics import roc_curve, auc

def micro_auc(y_test,y_score):
    # Compute ROC curve and ROC area for each class
    fpr = dict()
    tpr = dict()
    roc_auc = dict()

    # Compute micro-average ROC curve and ROC area
    fpr["micro"], tpr["micro"], _ = roc_curve(y_test.ravel(), y_score.ravel())
    roc_auc["micro"] = auc(fpr["micro"], tpr["micro"])

    return roc_auc["micro"]


def flat(Y):
    Y = np.asarray(Y,dtype='int8')
    newArray = np.zeros((Y.size, Y.max()+1))

    newArray[np.arange(Y.size),Y] = 1

    return newArray

# In[]:
yy_true = np.load('Y.npy')
yy_pred = np.loadtxt('results.csv', delimiter=',')

# In[]:
yy_true.shape
# In[]:
myauc = micro_auc(flat(yy_true)[0:1000,], yy_pred)

# %%
print(myauc)

# %%
