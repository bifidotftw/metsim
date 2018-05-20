Scope
=====

Simulate distribution of 13C for a given metabolic network at dynamic equilibrium (= constant metabolite pool sizes).


Programming Design
==================

MetSim is employing object oriented programming where each metabolite pool is an object.
Metabolite pools are stored in pandas dataframes where each column represents a carbon and each row a molecule.
Unlabeled (12C) carbons are represented by a 0, while labeled (13C) carbons are represented by a 1.
Consequently a molecule is a sequence of 0s and 1s, for example {'C1': 0, 'C2': 0, 'C3': 1} represents a three-carbon metabolite, in this case lactate, which is labeled once.

Calculation is done step-wise and every step consists of two sub-steps.
First, from every metabolite pool random molecules are moved to an intermediate pool called tmp and subsequently distributed to the specified neighboring metabolite pools.

All necessary sub-steps are defined as class methods or functions and ultimately strung together to form a sequence which makes up a step which is iterated to simulate the metabolism.
Enrichment is calculated, stored in a dataframe and can afterwards be saved to csv-files and plotted.
