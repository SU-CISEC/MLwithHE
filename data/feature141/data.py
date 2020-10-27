# In[]:
import numpy as np

X = np.load('X.npy') # loads your saved array into variable a.
Y = np.load('Y.npy') # loads your saved array into variable a.

# In[]:
X.shape
# %%
f = open("FeaturesLast.txt")
line = f.readline()
# %%
features = [int(x) for x in line.strip('[]').split(',')]
# %%
features
# %%
