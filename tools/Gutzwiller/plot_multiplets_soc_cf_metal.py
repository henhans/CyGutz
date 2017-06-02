'''
Standard plot of the histogram.
'''

import numpy as np
import h5py

prepath = "/Metal_soccf"

f = h5py.File("multiplets.h5", 'r')
weights = f["weights"][...]
deg_label = f["deg_label"][...]
j_label = f["j_label"][...]
n_label = f["n_label"][...]

chi_label = []; i = 0
while True:
  try:
    chi_label.append(f["/chi_label"+str(i)][...])
  except:
    break
  i += 1
f.close()
weights_1 = weights/deg_label
f_weights_1 = -np.log(weights_1)
f_weights_1 = f_weights_1 - f_weights_1[0]

# order
idx = f_weights_1.argsort()
f_weights_1 = f_weights_1[idx]
j_label = np.array(j_label)[idx]
weights = weights[idx]
n_label = np.array(n_label)[idx]
deg_label = deg_label[idx]
chi_label = np.array(chi_label)[idx]

np.set_printoptions(precision=2, suppress=True, linewidth=150)
print 'N: \n', n_label
print 'J: \n', j_label
print 'Degeneracy: \n', deg_label
print 'f_weight: \n', f_weights_1
print 'chi_label:'
for i,label in enumerate(chi_label):
  print i, np.array(label).real

#f = h5py.File("../../../UO2_prl_data.h5", 'r+')
#
#try:
#  del f[prepath]
#except:
#  pass
#
#n_min = int(np.min(n_label)+0.1)
#n_max = int(np.max(n_label)+0.1)
#
#for i in range(n_min, n_max+1):
#  idx = np.where(abs(n_label-i) < 0.1)
#  f[prepath+"/fn_gsoccf_"+str(i)] = f_weights_1[idx]
#  f[prepath+"/wn_gsoccf_"+str(i)] = weights[idx]
#
#f.close()

import matplotlib.pyplot as plt
plt.figure(0)
plt.xlabel('$f_{n}$')
plt.ylabel('$e^{-f_{n}}d_{n}$')
plt.plot(f_weights_1, weights,  'o')
plt.show()


