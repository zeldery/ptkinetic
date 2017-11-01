'''
The code to simulate kinetic
'''

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm

class Kinetic:
    def __init__(self,delta = 0.0001,criterion = 0.01):
        self.current = 0
        self.chemical = []
        self.concentrate = []
        self.constant = []
        self.k = []
        self.inp = []
        self.outp = []
        self.delta = delta
        self.criterion = criterion
        
    def add_chemical(self,name,concentrate,constant = False):
        self.chemical.append(name)
        self.concentrate.append([concentrate])
        self.constant.append(constant)
    
    def add_reaction(self,inp,outp,k):
        self.k.append(k)
        self.inp.append([self.chemical.index(name) for name in inp])
        self.outp.append([self.chemical.index(name) for name in outp])
        
    def run(self,n):
        for i in range(n):
            for c in self.concentrate:
                c += [c[self.current]]
            v = []
            for j in range(len(self.k)):
                temp = self.k[j]
                for m in self.inp[j]:
                    temp *= self.concentrate[m][self.current]
                v += [temp]
            self.current += 1
            for j in range(len(v)):
                v[j] *= self.delta
            for j in range(len(self.k)):
                for inp in self.inp[j]:
                    if not self.constant[inp]:
                        if self.concentrate[inp][self.current] >= 0.000001 and v[j] / self.concentrate[inp][self.current] >= self.criterion:
                            print(self.current)
                            print(inp)
                            print(j)
                            raise ValueError('The change in one step is too large, decrease the delta')
                        self.concentrate[inp][self.current] -= v[j]
                for outp in self.outp[j]:
                    if not self.constant[outp]:
                        self.concentrate[outp][self.current] += v[j]
    
    def reset(self):
        self.current = 0
        for i in range(len(self.concentrate)):
            self.concentrate[i] = self.concentrate[i][:1]
    
    def plot(self, chemical = None):
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        if chemical == None:
            ind = [i for i in range(len(self.concentrate))]
        else:
            ind = [self.chemical.index(name) for name in chemical]
        color_map = cm.rainbow(np.linspace(0,1,len(ind)))
        for i,j in enumerate(ind):
            ax.plot(np.arange(self.current+1)*self.delta,self.concentrate[j],c=color_map[i])
        fig.show()
        
    
    
