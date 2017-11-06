'''
The code to simulate kinetic

'''

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm

class Kinetic:
    def __init__(self,delta = 0.0001,criterion = 0.01):
        self.current = 0
        self.concentrate = pd.DataFrame(index = [0])
        self.constant = None
        self.k = None
        self.inp = []
        self.outp = []
        self.delta = delta
        self.criterion = criterion
        self.fig = []
        
    def add_chemical(self,name,concentrate,constant = False):
        self.concentrate[name] = concentrate
        if isinstance(self.constant,pd.Series):
            self.constant[name] = constant
        else:
            self.constant = pd.Series([constant],index=[name])
    
    def add_reaction(self,inp,outp,k):
        if isinstance(self.k,pd.Series):
            self.k[len(self.k)] = k
        else:
            self.k = pd.Series([k])
        self.inp.append(inp)
        self.outp.append(outp)
        
    def run(self,n):
        if not isinstance(n,int):
            n = int(n)
        temp = pd.DataFrame(columns = self.concentrate.columns, index = np.arange(self.current + n + 1))
        temp.iloc[:(self.current+1),:] = self.concentrate.iloc[:(self.current+1),:]
        self.concentrate = temp
        for i in range(n):
            self.current += 1
            self.concentrate.iloc[self.current,:] = self.concentrate.iloc[self.current-1,:]
            v = self.k * self.delta
            for j in range(len(self.k)):
                for inp in self.inp[j]:
                    v[j] *= self.concentrate.loc[self.current,inp]
            for j in range(len(self.k)):
                for inp in self.inp[j]:
                    if not self.constant[inp]:
                        if self.concentrate.loc[self.current,inp] >= 0.000001 and v[j] / self.concentrate.loc[self.current,inp] >= self.criterion:
                            print(self.current)
                            print(inp)
                            print(j)
                            raise ValueError('The change in one step is too large, decrease the delta')
                        self.concentrate.loc[self.current,inp] -= v[j]
                for outp in self.outp[j]:
                    if not self.constant[outp]:
                        self.concentrate.loc[self.current,outp] += v[j]
    
    def reset(self):
        self.current = 0
        self.concentrate = self.concentrate.iloc[-1,:]
        for fig in self.fig:
            plt.close(fig)
            
    def save(self,file_name):
        self.concentrate.to_csv(file_name)
        
    def plot(self, chemical = None):
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        if chemical == None:
            ind = self.concentrate.columns
        else:
            ind = chemical
        color_map = cm.rainbow(np.linspace(0,1,len(ind)))
        for i,j in enumerate(ind):
            ax.plot(np.arange(self.current+1)*self.delta,self.concentrate[j],c=color_map[i])
        fig.show()
        self.fig += [fig]
        
    
    
