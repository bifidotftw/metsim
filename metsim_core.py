import pandas as pd
import os.path
import csv
import matplotlib.pyplot as plt



def verbosity(flag):
    global verbose
    verbose = flag



class metabolite_pool(object):
    


    def __init__(self, metabolite_name, number_of_carbons, pool_size):
        self.metabolite_name = metabolite_name
        self.number_of_carbons = number_of_carbons
        self.pool_size = pool_size

        # Initialize empty variables for later use
        self.pool = None
        self.tmp = None
        self.collection_total_excess = pd.DataFrame()
        self.collection_isotop_distr = pd.DataFrame()



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



    def calculate_enrichment(self):
        global verbose
        if verbose:
            print('running %s.calculate_enrichment' % self.metabolite_name)

        # Total excess
        number_13C = self.pool.values.sum()
        number_totalC = self.pool.shape[0]*self.pool.shape[1]
        self.total_excess = number_13C/number_totalC

        # Calculate isotopologue distribution
        
        # Sums up all rows individually, sum represents number of 13C in the molecule
        self.row_sum = pd.DataFrame({'Number_of_13C': self.pool.sum(axis = 1)})

        # Calculates for every M+# the enrichment and stores it in a dict
        self.isotop_distr = {}
        for i in range(self.number_of_carbons):
            i = i+1
            tmp = self.row_sum[self.row_sum.Number_of_13C == i] # Sum of row == i -> M+i
            count = tmp.shape[0] # count rows where sum of rows == i -> number of molecules with M+i
            enrichment = count/self.row_sum.shape[0] # relative to the total amount of rows
            self.isotop_distr['M+' + str(i)] = enrichment

        if verbose:
            print('Total excess: ' + str(self.total_excess))
            print('Isotop. Distr.: ' + str(self.isotop_distr))
            print('')

        # Collect dicts in dataframes
        self.collection_total_excess = self.collection_total_excess.append({'total_excess': self.total_excess}, ignore_index = True)
        self.collection_isotop_distr = self.collection_isotop_distr.append(self.isotop_distr, ignore_index = True)



    def export_csv(self):
        if verbose:
            print('running %s export_csv' % self.metabolite_name)

        directory = os.path.join(os.getcwd(), 'csv')
        if not os.path.exists(directory):
            os.makedirs(directory)

        self.collection_total_excess.to_csv(os.path.join(directory, str(self.metabolite_name) + '_total_excess.csv'))
        self.collection_isotop_distr.to_csv(os.path.join(directory, str(self.metabolite_name) + '_isotop_distr.csv'))

        # Calculate relative isotopologue distribution
        self.collection_relative_isotop_distr = pd.DataFrame()
        for i in range(len(self.collection_isotop_distr.columns)):
            self.collection_relative_isotop_distr['M+' + str(i+1)] = self.collection_isotop_distr.loc[:, 'M+' + str(i+1)]/self.collection_isotop_distr.sum(axis = 1)
        
        # Export relative isotopologue distribution
        self.collection_relative_isotop_distr.to_csv(os.path.join(directory, str(self.metabolite_name) + '_relative_isotop_distr.csv'))


        if verbose:
            #print(self.collection_total_excess)
            #print(self.collection_isotop_distr)
            print('exported to folder "csv"')
            print('')
            print('drawing total excess and relative isotopologue distribution')
            print('')


        # Draw figures
        self.collection_total_excess.plot()
        #self.collection_isotop_distr.plot()

        # relative isotop distr
        # TODO: Debug
        fig, ax = plt.subplots()
        x = range(self.collection_relative_isotop_distr.shape[0])
        ax.stackplot(x, self.collection_relative_isotop_distr.T)
        ax.margins(0, 0) # Set margins to avoid "whitespace"

        plt.show()



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

        tmp = source.tmp.sample(n = number_of_molecules)
        source.tmp = source.tmp.loc[~source.tmp.index.isin(tmp.index)] # Consume molecules from source pool
        self.pool = pd.concat([self.pool, tmp])

        # Sanity check
        #if self.pool.shape[0] != self.pool_size:
        #    print('ERROR in function "%s.from_tmp": Pool size changed' % self.metabolite_name)
        #    quit()

        self.pool = self.pool.reset_index(drop=True)

        if verbose:
            print(self.pool)
            print('')



    def introduce_molecules(self, number_of_molecules, molecule):
        global verbose
        if verbose:
            print('running %s.introduce_molecules' % self.metabolite_name)

        # Sanity check
        if type(molecule) is not str:
            print('ERROR: variable "molecule" in function "introduce_molecules" is not a string.')
            quit()

        # Convert input string to dict
        molecule_tuple = tuple(molecule)

        # Sanity check
        if len(molecule_tuple) != self.number_of_carbons:
            print('ERROR in function introduce_molecules: Molecule specified has the wrong number of carbons')
            quit()

        molecule_dict = {}
        for i in range(len(molecule)):
            molecule_dict['C' + str(i+1)] =  int(molecule_tuple[i])


        new_molecules = pd.DataFrame(number_of_molecules*[molecule_dict])
        self.pool = pd.concat([self.pool, new_molecules])
        self.pool = self.pool.reset_index(drop = True)

        # Sanity check
        #if self.pool.shape[0] != self.pool_size:
        #    print('ERROR in function "introduce_molecules": Pool size of %s changed' % self.metabolite_name)
        #    quit()

        if verbose:
            print(self.pool)
            print('')
        


    def mirror_symmetry(self):
        if verbose:
            print('running mirror_symmetry')
            print('')

        # Randomly choose 50% of the rows
        to_be_rotated = self.pool.sample(frac=0.5)
        # Selects all rows which are NOT in rotate and updates pool
        self.pool = self.pool.loc[~self.pool.index.isin(to_be_rotated.index)]

        have_been_rotated = to_be_rotated.rename(columns={a: b for a, b in zip(to_be_rotated.columns, reversed(to_be_rotated.columns))})
        print(have_been_rotated)

        # Concatenate dataframes
        self.pool = pd.concat([have_been_rotated, self.pool])

        # Sanity check
        #if succinate.pool.shape[0] != succinate.pool_size:
        #    print('ERROR in function "mirror_symmetry_succinate": Pool size changed')
        #    quit()
        



# Metabolite conversions

def oxaloacetate_to_citrate(number_of_molecules):
    if verbose:
        print('running oxaloacetate_to_citrate')

    tmp1 = oxaloacetate.tmp.sample(n = number_of_molecules)
    tmp1 = tmp1.rename(columns={'C1': 'C6', 'C2': 'C3', 'C3': 'C4', 'C4': 'C5'}) # Rename carbons according to new position in citrate, side chain is assigned C6
    tmp1 = tmp1.reset_index(drop = True)

    tmp2 = pyruvate.tmp.sample(n = number_of_molecules)
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
#pool1 = metabolite_pool('pool1', 4, 10)
#pool1.initialize_pool()
#
#pool2 = metabolite_pool('pool2', 4, 10)
#pool2.initialize_pool()
#
#pool1.to_tmp(3)
#pool2.to_tmp(6)
#
#pool1.from_tmp(3, pool2)
#print(pool2.tmp)
#print('')



# Initialize pools

verbose = True

pyruvate = metabolite_pool('pyruvate', 3, 5) # will be replenished infinitely with labeled carbons
pyruvate.initialize_pool()

citrate = metabolite_pool('citrate', 6, 100)
citrate.initialize_pool()

glutamate = metabolite_pool('glutamate', 5, 100)
glutamate.initialize_pool()

succinate = metabolite_pool('succinate', 4, 100)
succinate.initialize_pool()

oxaloacetate = metabolite_pool('oxaloacetate', 4, 100)
oxaloacetate.initialize_pool()


# Sequence

# Start with labeled pyruvate pool
pyruvate.to_tmp(5)
pyruvate.introduce_molecules(5, '111')

for i in range(1000):

    succinate.mirror_symmetry()

    # Move to tmp
    pyruvate.to_tmp(5)
    citrate.to_tmp(5)
    glutamate.to_tmp(10)
    succinate.to_tmp(5)
    oxaloacetate.to_tmp(5)

    # Move from tmp
    pyruvate.introduce_molecules(5, '111')

    oxaloacetate_to_citrate(5)
    citrate.from_tmp(5, oxaloacetate_to_citrate)

    citrate_to_glutamate(5)
    glutamate.introduce_molecules(5, '00000')
    glutamate.from_tmp(5, citrate_to_glutamate)

    glutamate_to_succinate(5)
    succinate.from_tmp(5, glutamate_to_succinate)

    oxaloacetate.from_tmp(5, succinate)

    # Calculate enrichment
    citrate.calculate_enrichment()
    glutamate.calculate_enrichment()
    oxaloacetate.calculate_enrichment()

citrate.export_csv()
