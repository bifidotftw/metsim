import pandas as pd

verbose = True

class metabolite_pool(object):
    


    def __init__(self, metabolite_name, number_of_carbons, pool_size):
        self.metabolite_name = metabolite_name
        self.number_of_carbons = number_of_carbons
        self.pool_size = pool_size

        # Initialize empty variables for later use
        self.pool = None
        self.tmp = None

        # Initialize dataframes with first value(s) of 0 for exporting to csv
        self.collection_total_excess = pd.DataFrame({'total_excess': [0]})
        # Creates dataframe according to number_of_carbons of the metabolite
        first_row = {}
        for i in range(self.number_of_carbons):
            i = i+1
            first_row['M+' + str(i)] = 0
            self.collection_isotop_distr = pd.DataFrame.from_dict([first_row])



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
    def to_tmp(self, number_of_molecules):
        global verbose
        if verbose:
            print('running %s.to_tmp' % self.metabolite_name)

        self.tmp = self.pool.sample(n = number_of_molecules)
        self.pool = self.pool.loc[~self.pool.index.isin(self.tmp.index)] # rows which are not in self.tmp

        if verbose:
            print(self.tmp)
            print('')


    
    def from_tmp(self, number_of_molecules, source):
        global verbose
        if verbose:
            print('running %s.from_tmp' % self.metabolite_name)

        tmp = source.sample(n = number_of_molecules)
        # TODO: Debug
        # source == pool2.tmp (specified in line 95)
        source = source.loc[~source.index.isin(tmp.index)] # Consume molecules from source pool
        #pool2.tmp = source.loc[~source.index.isin(tmp.index)] # Produces desired output!

        self.pool = pd.concat([self.pool, tmp])

        # Sanity check
        #if self.pool.shape[0] != self.pool_size:
        #    print('ERROR in function "%s.from_tmp": Pool size changed' % self.metabolite_name)
        #    quit()

        self.pool = self.pool.reset_index(drop=True)

        if verbose:
            print(self.pool)
            print('')


# Sequence

pool1 = metabolite_pool('pool1', 4, 10)
pool1.initialize_pool()

pool2 = metabolite_pool('pool2', 4, 10)
pool2.initialize_pool()

pool1.to_tmp(3)
pool2.to_tmp(6)

pool1.from_tmp(3, pool2.tmp) # see line 64
print(pool2.tmp)
# print: pool2.tmp
#   C1  C2  C3  C4
#7   0   0   0   0
#8   0   0   0   0
#3   0   0   0   0
#5   0   0   0   0
#6   0   0   0   0
#4   0   0   0   0
#(Indices are random)
#
# Desired output:
#   C1  C2  C3  C4
#7   0   0   0   0
#8   0   0   0   0
#3   0   0   0   0
#(Indices are random)
