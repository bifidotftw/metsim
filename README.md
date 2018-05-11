Scope
=====

Simulate enrichment for a given metabolic network with constant pool sizes.
Considered variables:
- Pool sizes
- Flux
Omitted variables:
- Enzyme kinetics
- Changing pool sizes
- All processes are instantaneous, no delays e.g. for transport across membranes


Functionality
=============

Calculation is done step-wise.
Flux is driven by the "suction" created by removing (i.e. consuming) metabolites in contrast to pushing metabolites through the network.
Consequently, each step starts with removing endproducts.
Subsequently, pools consisting of less molecules than their constant concentration are filled up from upstream metabolite pools.
Although this behavior is counterintuitive when reading an writing code it is deemed a more accurate representation of reality in conditions where the conditions (nutrient concentrations) are approximately constant.
A step is finished once all pools are replenished.
Check $diagram for a detailed behavioral diagram.


Programming Design
==================

MetSim is employing object oriented programming where each metabolite pool is an object.
Metabolite pools are stored in pandas dataframes where each column represents a carbon and each row a molecule.
Unlabeled (12C) carbons are represented by a 0, while labeled (13C) carbons are represented by a 1.
Consequently a molecule is a sequence of 0s and 1s, for example {'C1': 0, 'C2': 0, 'C3': 1} represents a three-carbons metabolite, e.g. lactate, which is labeled once.
All necessary functions are defined as methods and ultimately strung together to form a sequence which makes up a step which is iterated to simulate the metabolism.
Enrichment is calculated, stored in a dataframe and can afterwards be saved to csv-files and plotted.
Check $diagram for a detailed class diagram.


Planned Features
================




ToDo
====


