# In[]:
import numpy as np

X = np.load('X.npy') # loads your saved array into variable a.
Y = np.load('Y.npy') # loads your saved array into variable a.

# In[]:
X.shape
# %%
f = open("../../model/FeaturesLastSorted.txt")
line = f.readline()
# %%
features = [int(x) for x in line.strip('[]').split(',')]
# %%
features.sort()
print(features)
# %%
X_sorted = np.take(X, features, axis=1)
# %%
X_sorted.shape
# %%
np.savetxt("X_test.csv", X_sorted, delimiter=",",fmt="%-.1f")
# %%
