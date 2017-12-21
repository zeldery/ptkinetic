'''
ptkinetic package
The code for kinetic simulation
Written by: Thien-Phuc Tu-Nguyen
Last modified: December 2017
Source: https://github.com/zeldery/ptkinetic
'''

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm

class Kinetic:
    '''The main class for Kinetic
'''
    def __init__(self):
        '''Variable declaration'''
        self.delta = 0.0 # the time step, must be small enough
        self.check = False # Whether check the criterion or not
        self.criterion = 1.0 # the criterion of the ratio between the change and current concentration of any chemicals
        self.num_of_steps = 1
        
        # Flag variables
        self.flag_init = False
        self.flag_run = False
        # Chemical properties
        self.chemicals = [] # Contain the name of reaction
        self.concentrations = [] # Contain the initial concentration
        self.chemical_stables = [] # True if the concentration does not change during the simulation
        # Reaction properties
        self.reaction_inputs = [] # List of reactants
        self.reaction_outputs = [] # List of products
        self.reaction_constants = [] # List of rate constants
        
        # Data variables
        self.times = None
        self.data = None
        self.k = None
        # Other properties
        self.fig = [] # List of figure object, for clean-up only
        
        
    def add_chemical(self,name,concentration,stable = False):
        '''Add chemical nomenclature to the model'''
        if self.flag_init:
            raise RuntimeError('Cannot add chemical after initialized')
        if name in self.chemicals:
            raise ValueError('The chemical name is already in the Kinetic object')
        self.chemicals += [name]
        self.concentrations += [concentration]
        self.chemical_stables += [stable]
    
    
    def add_reaction(self, inp, outp, constant):
        '''Add reactions to model, including reactants, products and rate constant'''
        if self.flag_init:
            raise RuntimeError('Cannot add reaction after initialized')
        for name in inp:
            if name not in self.chemicals:
                raise ValueError('The chemical name has not been added')
        for name in outp:
            if name not in self.chemicals:
                raise ValueError('The chemical name has not been added')
        self.reaction_inputs += [inp]
        self.reaction_outputs += [outp]
        self.reaction_constants += [constant]
    
    def init(self, num_of_steps = 1, check = False, criterion = None):
        '''Initialize the variable necessary in the simulation'''
        if self.flag_init:
            raise RuntimeError('Cannot reinitialize')
        if check and criterion == None:
            raise ValueError('Criterion is not given to check')
        self.check = check
        self.criterion = criterion
        self.flag_init = True
        self.num_of_steps = num_of_steps
        n_chem = len(self.chemicals)
        n_react = len(self.reaction_constants)
        self.k = np.array(self.reaction_constants) # The constant with the delta multiplied in advance
        self.inp_marker = np.zeros((n_react,n_chem)) # A matrix, each row is an reaction, each column is an chemical
        # Mark the reactant 1 for each time they appear
        for i in range(n_react):
            for name in self.reaction_inputs[i]:
                self.inp_marker[i,self.chemicals.index(name)] += 1
        # The same with products
        self.outp_marker = np.zeros((n_react,n_chem))
        for i in range(n_react):
            for name in self.reaction_outputs[i]:
                self.outp_marker[i,self.chemicals.index(name)] += 1
        self.times = np.array([0.0])
        self.data = np.array(self.concentrations).reshape((1,n_chem)) # Matrix with each row is run, each column is a chemical
        self.stables = np.array(self.chemical_stables) # Convert to numpy object
        

    
    def run(self, cycle, delta):
        ''' Run function '''
        if not self.flag_init:
            raise RuntimeError('You have to initialize before run')
        self.flag_run = True
        k = self.k * delta / self.num_of_steps
        current, n_chem = self.data.shape
        temp = self.times[-1] + np.arange(1,cycle+1)*delta
        self.times = np.concatenate((self.times, temp), axis = 0)
        temp = np.zeros((cycle, n_chem))
        self.data = np.concatenate((self.data, temp), axis = 0)
        temp = self.data[current-1,:]
        for i in range(cycle):
            for j in range(self.num_of_steps):
                rate = (temp ** self.inp_marker).prod(axis=1) * k
                change_outp = (self.outp_marker.T * rate).sum(axis = 1)
                change_inp = (self.inp_marker.T * rate).sum(axis = 1)
                change = np.where(self.stables, 0.0, change_outp - change_inp) # Do not change if the concentration is stable
                # Check if the criterion is violated
                # Only check for consumed amount, not forming amount
                if self.check and np.where(temp != 0, np.abs(change_inp/temp) , 0.0).max() > self.criterion:
                    raise ValueError('The change in one step is too large, decrease the step size')
                temp = temp + change
            self.data[current,:] = temp
            current += 1
    
    def reset(self):
        '''Free the data and figure, keep only the last concentration'''
        self.data = self.data[-1,:]
        for fig in self.fig:
            plt.close(fig)
            
    def save(self,file_name):
        '''Save the data to a csv file'''
        data = pd.DataFrame(self.data,columns = self.chemicals)
        data.to_csv(file_name)
        
    def plot(self, chemical = None):
        ''' Plot the concentration '''
        if not self.flag_run:
            raise RuntimeError('You have to run before plot')
        n_chem = self.data.shape[1]
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        if chemical == None:
            ind = range(n_chem) # Plot all the concentration
        else:
            ind = [self.chemicals.index(x) for x in chemical] # Plot the concentration given by chemical argument
        color_map = cm.rainbow(np.linspace(0,1,len(ind))) # Continuous rainbow color
        for i,j in enumerate(ind):
            temp, = ax.plot(self.times, self.data[:,j], c=color_map[i],label = self.chemicals[j])
        ax.legend()
        fig.show()
        self.fig += [fig] # Save the figure to delete
        
    
    
