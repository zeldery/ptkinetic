'''
The code to simulate kinetic

'''

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm

class Kinetic:
    def __init__(self,delta = None, check = True, criterion = None):
        self.delta = delta
        self.criterion = criterion
        self.check = True
        self.status = 0 # 1 for add chemical, 2 for add reaction, 3 for initialization, 4 for already run
        # Chemical properties
        self.chemicals = []
        self.concentrations = []
        self.stables = []
        # Reaction properties
        self.reaction_inputs = []
        self.reaction_outputs = []
        self.reaction_constants = []
        # Other properties
        self.fig = []
        
        
    def add_chemical(self,name,concentration,stable = False):
        if self.status == 0:
            self.status = 1
        if self.status == 1:
            if name in self.chemicals:
                raise ValueError('The chemical name is already in the Kinetic object')
            self.chemicals += [name]
            self.concentrations += [concentration]
            self.stables += [stable]
        else:
            raise RuntimeError('Cannot add chemical after adding reaction')
    
    
    def add_reaction(self, inp, outp, constant):
        if self.status == 1:
            self.status = 2
        if self.status == 2:
            for name in inp:
                if name not in self.chemicals:
                    raise ValueError('The chemical name has not been added')
            for name in outp:
                if name not in self.chemicals:
                    raise ValueError('The chemical name has not been added')
            self.reaction_inputs += [inp]
            self.reaction_outputs += [outp]
            self.reaction_constants += [constant]
        elif self.status == 0:
            raise RuntimeError('Cannot add reaction before adding chemical')
        else:
            raise RuntimeError('Cannot add reaction after running simulation')
    
    def init(self, delta = None, check = None, criterion = None):
        if self.status == 2:
            if delta != None:
                self.delta = delta
            if criterion != None:
                self.criterion = criterion
            if check != None:
                self.check = check
            if self.delta == None or self.criterion == None:
                raise ValueError('Missing step size or criterion')
            self.status = 3
            n_chem = len(self.chemicals)
            n_react = len(self.reaction_constants)
            self.k = np.array(self.reaction_constants) * self.delta
            self.inp_marker = np.zeros((n_react,n_chem))
            for i in range(n_react):
                for name in self.reaction_inputs[i]:
                    self.inp_marker[i,self.chemicals.index(name)] += 1
            self.outp_marker = np.zeros((n_react,n_chem))
            for i in range(n_react):
                for name in self.reaction_outputs[i]:
                    self.outp_marker[i,self.chemicals.index(name)] += 1
            self.data = np.array(self.concentrations).reshape((1,n_chem))
            self.stables = np.array(self.stables)
        else:
            if self.status < 2:
                raise RuntimeError('Cannot initialize before adding reaction')
            else:
                raise RuntimeError('Already initialized before')
    
    def run(self,times):
        if self.status == 3:
            n_react, n_chem = self.inp_marker.shape
            self.status = 4
            self.data = np.zeros((times+1, n_chem))
            self.data[0,:] = np.array(self.concentrations)
            for i in range(times):
                rate = ( self.data[i,:] ** self.inp_marker ).prod(axis=1) * self.k
                change_outp = (self.outp_marker.T * rate).sum(axis = 1)
                change_inp = (self.inp_marker.T * rate).sum(axis = 1)
                change = np.where(self.stables, 0.0, change_outp - change_inp)
                if self.check and np.where(self.data[i,:] != 0, np.abs(change_inp/self.data[i,:]) , 0.0).max() > self.criterion:
                    raise ValueError('The change in one step is too large, decrease the step size')
                self.data[i+1,:] = self.data[i,:] + change
        elif self.status == 4:
            n_react, n_chem = self.inp_marker.shape
            data = np.zeros((times+1,n_chem))
            data[0,:] = self.data[-1,:]
            for i in range(times):
                rate = ( data[i,:] ** self.inp_marker ).prod(axis=1) * self.k
                change_outp = (self.outp_marker.T * rate).sum(axis = 1)
                change_inp = (self.inp_marker.T * rate).sum(axis = 1)
                change = np.where(self.stables, 0.0, change_outp - change_inp)
                if self.check and np.where(data[i,:] !=0 , np.abs(change_inp/data[i,:]) , 0.0).max() > self.criterion:
                    raise ValueError('The change in one step is too large, decrease the step size')
                data[i+1,:] = data[i,:] + change
            self.data = np.concatenate((self.data,data[1:,:]),axis = 0)
        if self.status < 3:
            raise RuntimeError('Cannot run before initialization')
    
    def reset(self):
        self.data = self.data[-1,:]
        self.status = 3
        for fig in self.fig:
            plt.close(fig)
            
    def save(self,file_name):
        data = pd.DataFrame(self.data,columns = self.chemicals)
        data.to_csv(file_name)
        
    def plot(self, chemical = None):
        if self.status < 4:
            raise RuntimeError('Cannot plot before run')
        n_run, n_chem = self.data.shape
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        if chemical == None:
            ind = range(n_chem)
        else:
            ind = [self.chemicals.index(x) for x in chemical]
        color_map = cm.rainbow(np.linspace(0,1,len(ind)))
        for i,j in enumerate(ind):
            temp, = ax.plot(self.delta*np.arange(n_run), self.data[:,j], c=color_map[i],label = self.chemicals[j])
        ax.legend()
        fig.show()
        self.fig += [fig]
        
    
    
