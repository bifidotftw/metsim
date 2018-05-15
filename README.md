Scope
=====

Simulate enrichment kinetic for a given metabolic network at steady state (= constant metabolite pool sizes).


Functionality
=============

In constant (nutrient) conditions influx and outflux of every metabolite pool and consequently metabolite pool sizes are constant.
These values are used to simulate the flux in the metabolic network.

Like real metabolic pathways, flux is controlled by a few rate-limiting, and often irreversible, reactions while most other reactions are close to equilibrium and therefore exhibit fast back- and forth reactions, possibly even behaving like a single pool.
Check $diagram for a detailed behavioral diagram.

Programming Design
==================

MetSim is employing object oriented programming where each metabolite pool is an object.
Metabolite pools are stored in pandas dataframes where each column represents a carbon and each row a molecule.
Unlabeled (12C) carbons are represented by a 0, while labeled (13C) carbons are represented by a 1.
Consequently a molecule is a sequence of 0s and 1s, for example {'C1': 0, 'C2': 0, 'C3': 1} represents a three-carbon metabolite, in this case lactate, which is labeled once.

Calculation is done step-wise and every step consists of two sub-steps.
First, from every metabolite pool random molecules are moved to a intermediate pool called temp and subsequently distributed to the specified neighboring metabolite pools.

All necessary steps are defined as class methods or functions and ultimately strung together to form a sequence which makes up the step which is iterated to simulate the metabolism.
Enrichment is calculated, stored in a dataframe and can afterwards be saved to csv-files and plotted.
Check $diagram for a detailed class diagram.


Planned Features
================

- GUI
- Use machine learning to generate desired enrichment patterns by changing parameters


ToDo
====

- Write function which calculates enrichment after enrichment has been equilibrated between two pools
- Calculate enrichment for every carbon
