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
including 4 basic steps:
+ add chemicals
+ add reactions
+ init
+ run'''
    def __init__(self,delta = None, check = False, criterion = None):
        '''Variable declaration'''
        self.delta = delta # the time step, must be small enough
        self.check = check # Whether check the criterion or not
        self.criterion = criterion # the criterion of the ratio between the change and current concentration of any chemicals
        self.status = 0 # 1 for add chemical, 2 for add reaction, 3 for initialization, 4 for already run
        # Chemical properties
        self.chemicals = [] # Contain the name of reaction
        self.concentrations = [] # Contain the initial concentration
        self.stables = [] # True if the concentration does not change during the simulation
        # Reaction properties
        self.reaction_inputs = [] # List of reactants
        self.reaction_outputs = [] # List of products
        self.reaction_constants = [] # List of rate constants
        # Other properties
        self.fig = [] # List of figure object, for clean-up only
        
        
    def add_chemical(self,name,concentration,stable = False):
        '''Add chemical nomenclature to the model'''
        if self.status == 0:
            self.status = 1 # Update the status
        if self.status == 1:
            if name in self.chemicals:
                raise ValueError('The chemical name is already in the Kinetic object')
            self.chemicals += [name]
            self.concentrations += [concentration]
            self.stables += [stable]
        else: # Error if already add a reaction
            raise RuntimeError('Cannot add chemical after adding reaction')
    
    
    def add_reaction(self, inp, outp, constant):
        '''Add reactions to model, including reactants, products and rate constant'''
        if self.status == 1:
            self.status = 2
        if self.status == 2:
            # Ensure that all the names are valid
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
            raise RuntimeError('Cannot add reaction before adding chemical') # Check for right sequence
        else:
            raise RuntimeError('Cannot add reaction after running initialization step')
    
    def init(self, delta = None, check = None, criterion = None):
        '''Initialize the variable necessary in the simulation'''
        if self.status == 2:
            if delta != None:
                self.delta = delta
            if criterion != None:
                self.criterion = criterion
            if check != None:
                self.check = check
            if self.check == False:
                if self.delta == None:
                    raise ValueError('Missing step size')
            elif self.delta == None or self.criterion == None:
                raise ValueError('Missing step size or criterion') # Check for missing required parameters
            self.status = 3
            n_chem = len(self.chemicals)
            n_react = len(self.reaction_constants)
            self.k = np.array(self.reaction_constants) * self.delta # The constant with the delta multiplied in advance
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
            self.data = np.array(self.concentrations).reshape((1,n_chem)) # Matrix with each row is run, each column is a chemical
            self.stables = np.array(self.stables) # Convert to numpy object
        else:
            if self.status < 2: # Check the right sequence
                raise RuntimeError('Cannot initialize before adding reaction')
            else:
                raise RuntimeError('Already initialized before')
    
    def run(self,times):
        if self.status == 3: # If run for the first time
            n_react, n_chem = self.inp_marker.shape
            self.status = 4
            self.data = np.zeros((times+1, n_chem)) # Create empty data
            self.data[0,:] = np.array(self.concentrations)
            for i in range(times):
                rate = ( self.data[i,:] ** self.inp_marker ).prod(axis=1) * self.k # Rate formula in numpy format
                                                                                   # Using numpy broadcasting rules
                change_outp = (self.outp_marker.T * rate).sum(axis = 1)
                change_inp = (self.inp_marker.T * rate).sum(axis = 1)
                change = np.where(self.stables, 0.0, change_outp - change_inp) # Do not change if the concentration is stable
                # Check if the criterion is violated
                # Only check for consumed amount, not forming amount
                if self.check and np.where(self.data[i,:] != 0, np.abs(change_inp/self.data[i,:]) , 0.0).max() > self.criterion:
                    raise ValueError('The change in one step is too large, decrease the step size')
                self.data[i+1,:] = self.data[i,:] + change
        elif self.status == 4: # If already run
            n_react, n_chem = self.inp_marker.shape
            data = np.zeros((times+1,n_chem)) # Create different data
            data[0,:] = self.data[-1,:] # Copy the last simulation concentration
            for i in range(times):
                # The same as the case of status=3
                rate = ( data[i,:] ** self.inp_marker ).prod(axis=1) * self.k
                change_outp = (self.outp_marker.T * rate).sum(axis = 1)
                change_inp = (self.inp_marker.T * rate).sum(axis = 1)
                change = np.where(self.stables, 0.0, change_outp - change_inp)
                if self.check and np.where(data[i,:] !=0 , np.abs(change_inp/data[i,:]) , 0.0).max() > self.criterion:
                    raise ValueError('The change in one step is too large, decrease the step size')
                data[i+1,:] = data[i,:] + change
            self.data = np.concatenate((self.data,data[1:,:]),axis = 0) # Append the data to self.data
        if self.status < 3:
            raise RuntimeError('Cannot run before initialization')
    
    def reset(self):
        '''Free the data and figure, keep only the last concentration'''
        self.data = self.data[-1,:]
        self.status = 3
        for fig in self.fig:
            plt.close(fig)
            
    def save(self,file_name):
        '''Save the data to a csv file'''
        data = pd.DataFrame(self.data,columns = self.chemicals)
        data.to_csv(file_name)
        
    def plot(self, chemical = None):
        ''' Plot the concentration '''
        if self.status < 4:
            raise RuntimeError('Cannot plot before run')
        n_run, n_chem = self.data.shape
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        if chemical == None:
            ind = range(n_chem) # Plot all the concentration
        else:
            ind = [self.chemicals.index(x) for x in chemical] # Plot the concentration given by chemical argument
        color_map = cm.rainbow(np.linspace(0,1,len(ind))) # Continuous rainbow color
        for i,j in enumerate(ind):
            temp, = ax.plot(self.delta*np.arange(n_run), self.data[:,j], c=color_map[i],label = self.chemicals[j])
        ax.legend()
        fig.show()
        self.fig += [fig] # Save the figure to delete
        
    
    
