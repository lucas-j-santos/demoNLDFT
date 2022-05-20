import numpy as np 
from matplotlib import pyplot as plt

dft = np.loadtxt("saida.dat", dtype=float) 
#mc = np.loadtxt("MC.dat", dtype=float) 

def cm_to_inch(value):
    return value/2.54

plt.figure(figsize=(cm_to_inch(13.5), cm_to_inch(10.5)))

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.size": 14
})

plt.plot(dft[:,0], dft[:,1])
#plt.plot(mc[:,0]/3.575, mc[:,1]*3.575**3, 'ko', mfc = 'none', label='GCMC')
#plt.xlim([-0.1,5])
#plt.ylim([-0.2,8.5])
plt.xlabel(r'$ z/ \sigma $', fontsize=15)
plt.ylabel(r'$\rho \sigma^3 $', fontsize=15)
#plt.legend(fontsize=12, frameon=False)
plt.savefig('figure.pdf', bbox_inches='tight')
plt.show()