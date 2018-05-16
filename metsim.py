import pandas as pd

verbose = True

class metabolite_pool(object):
    
    def __init__(self, metabolite_name, number_of_carbons, pool_size):
        self.metabolite_name = metabolite_name
        self.number_of_carbons = number_of_carbons
        self.pool_size = pool_size

        self.pool = None
        self.tmp = None

    def initialize_pool(self):
        global verbose
        if verbose:
            print('running %s.initialize_pool' % self.metabolite_name)
        #Create dataframe according to number of carbons and pool size with unlabeled molecules
        columns = {}
        for i in range(self.number_of_carbons): # create dictionary from C1 to CX (number of carbons with unlabeled carbons (0)
            columns['C' + str(i+1)] = self.pool_size*[0]
        self.pool = pd.DataFrame(columns)
        if verbose:
            print(self.pool)
            print('')

    def calculate_enrichment(self, target):
        global verbose
        if verbose:
            print('running %s.calculate_enrichment' % self.metabolite_name)
        # Total excess
        number_13C = self.pool.values.sum() # Sum all values
        #print('Number of 13C-atoms: ' + str(number_13C))
        number_totalC = self.pool.shape[0]*self.pool.shape[1]
        #print('Total number of carbons: ' + str(number_totalC))
        self.total_excess = number_13C/number_totalC

        # Calculate isotopologue distribution
        # Create dataframe where every row only contains the number of 13C
        ID = pd.DataFrame({'Number_of_13C': self.pool.sum(axis=1)})
        # M+1
        metabolites_M1 = ID[ID.Number_of_13C == 1] # create new dataframe only according to condition
        count_M1 = metabolites_M1.shape[0] # count rows in dataframe
        M1 = count_M1/ID.shape[0]
        # M+2
        metabolites_M2 = ID[ID.Number_of_13C == 2]
        count_M2 = metabolites_M2.shape[0]
        M2 = count_M2/ID.shape[0]
        # M+3
        metabolites_M3 = ID[ID.Number_of_13C == 3]
        count_M3 = metabolites_M3.shape[0]
        M3 = count_M3/ID.shape[0]
        # M+4
        metabolites_M4 = ID[ID.Number_of_13C == 4]
        count_M4 = metabolites_M4.shape[0]
        M4 = count_M4/ID.shape[0]
        # M+5
        metabolites_M5 = ID[ID.Number_of_13C == 5]
        count_M5 = metabolites_M5.shape[0]
        M5 = count_M5/ID.shape[0]
        # M+6
        metabolites_M6 = ID[ID.Number_of_13C == 6]
        count_M6 = metabolites_M6.shape[0]
        M6 = count_M6/ID.shape[0]

        if target == 'terminal':
            print('Total Excess: ' + str(self.total_excess))
            print('M+1: ' + str(M1))
            print('M+2: ' + str(M2))
            print('M+3: ' + str(M3))
            print('M+4: ' + str(M4))
            print('M+5: ' + str(M5))
            print('M+6: ' + str(M6))
            print('')

        if target == 'csv':
            #TODO: export to csv
            pass

    def to_tmp(self, number_of_molecules):
        global verbose
        if verbose:
            print('running %s.to_tmp' % self.metabolite_name)
        self.tmp = self.pool.sample(n = number_of_molecules)
        self.pool = self.pool.loc[~self.pool.index.isin(self.tmp.index)]
        if verbose:
            print(self.tmp)
            print('')
    
    def from_tmp(self, number_of_molecules, source):
        global verbose
        if verbose:
            print('running %s.from_tmp' % self.metabolite_name)
        tmp = source.sample(n = number_of_molecules)
        
        #test
        source.tmp = source.loc[~source.index.isin(tmp.index)] # Consume molecules from source pool

        self.pool = pd.concat([self.pool, tmp])
        # Sanity check
        if self.pool.shape[0] != self.pool_size:
            print('ERROR in function "%s.from_tmp": Pool size changed' % self.metabolite_name)
            quit()
        self.pool = self.pool.reset_index(drop=True)
        if verbose:
            print(self.pool)
            print('')

    def introduce_molecule(self, number_of_molecules, molecule):
        global verbose
        if verbose:
            print('running %s.introduce_molecule' % self.metabolite_name)
        if type(molecule) is not dict:
            print('ERROR: variable "molecule" in function "introduce_molecule" is not a dictionary.')
            quit()
        new_molecules = pd.DataFrame(number_of_molecules*[molecule])
        self.pool = pd.concat([self.pool, new_molecules])
        self.pool = self.pool.reset_index(drop=True)
        # Sanity check
        if self.pool.shape[0] != self.pool_size:
            print('ERROR in function "introduce_molecule": Pool size of %s changed' % self.metabolite_name)
            quit()
        if verbose:
            print(self.pool)
        

    def mirror_symmetry(self):
        #TODO: generalize mirror symmetry
        pass
        
# Metabolite conversions

def oxaloacetate_to_citrate(number_of_molecules):
    if verbose:
        print('running oxaloacetate_to_citrate')
    tmp1 = oxaloacetate.tmp.sample(n = number_of_molecules)
    tmp1 = tmp1.rename(columns={'C1': 'C6', 'C2': 'C3', 'C3': 'C4', 'C4': 'C5'}) # Rename carbons according to new position in citrate, side chain is assigned C6
    tmp1 = tmp1.reset_index(drop = True)

    tmp2 = pyruvate.tmp.sample(n = number_of_molecules)
    # TODO: check which carbon of pyruvate is removed
    tmp2 = tmp2.drop(['C1'], axis = 1)
    tmp2 = tmp2.rename(columns = {'C2': 'C1', 'C3': 'C2'})
    tmp2 = tmp2.reset_index(drop = True)

    oxaloacetate_to_citrate.tmp = tmp1.join(tmp2)

    if verbose:
        print(oxaloacetate_to_citrate.tmp)
        print('')

def citrate_to_glutamate(number_of_molecules):
    if verbose:
        print('running citrate_to_glutamate')
    tmp = citrate.pool.sample(n = number_of_molecules)
    citrate_to_glutamate.tmp = tmp.drop(['C6'], axis = 1)
    if verbose:
        print(citrate_to_glutamate.tmp)
        print('')

def glutamate_to_succinate(number_of_molecules):
    if verbose:
        print('running glutamate_to_succinate')
    tmp = glutamate.pool.sample(n = number_of_molecules)
    glutamate_to_succinate.tmp = tmp.drop(['C5'], axis = 1)
    if verbose:
        print(glutamate_to_succinate.tmp)
        print('')
    



# Testing

pool1 = metabolite_pool('pool1', 4, 10)
pool1.initialize_pool()

pool2 = metabolite_pool('pool2', 4, 10)
pool2.initialize_pool()

pool1.to_tmp(3)
pool2.to_tmp(6)

pool1.from_tmp(3, pool2.tmp)
print(pool2.tmp)


# TCA cycle

# Initialte pools
#pyruvate = metabolite_pool('pyruvate', 3, 10) # will be replenished infinitely with labeled carbons
#pyruvate.initialize_pool()
#
#citrate = metabolite_pool('citrate', 6, 10)
#citrate.initialize_pool()
#
#glutamate = metabolite_pool('glutamate', 5, 10)
#glutamate.initialize_pool()
#
#aspartate = metabolite_pool('aspartate', 4, 10)
#aspartate.initialize_pool()


# move molecules to tmp for every metabolite pool
#TODO: rerwite from_tmp and metabolite specific functions in a way that they 'consume' molecules from tmp, i.e. they delete the ones they used
# Implementation: see to_tmp, make class.tmp variables writable from outside?
