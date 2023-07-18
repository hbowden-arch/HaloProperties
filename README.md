# HaloProperties
Code for predicting halo properties from observable environment information.

## Overview

The goal of this project is to use observable environmental information (kNN distances or counts in cylinders) to predict the halo properties of a given galaxy.
As of 6/15/23 this repository contains work on halo masses, including code for training neural networks, an array of pre-trained models, and code for analyzing model performances. Future work is planned for exploring halo properties beyond mass.

For more information, see https://arxiv.org/abs/2307.07549 or contact me at hbowden@arizona.edu


## Requirements

Tensorflow version >= 2.0



## Simulation data

Bolshoi-Planck/SMDPL + UniverseMachine

http://hipacc.ucsc.edu/Bolshoi/index.html

https://bitbucket.org/pbehroozi/universemachine/src/main/
 

## Nearest Neighbors

"find_neighbors.c" uses existing UniverseMachine code to:

1. Search for an object's nearest neighbors (up to 50) within the simulation box with periodic boundary conditions
2. Calculate the number of neighbors within cylindrical bins around the target object

## Predicting Halo Masses

See HaloMassNetwork.ipynb

Models made with tensorflow. Several pre-trained models are available in the "Models" directory.

Sample statistics for our training data (SMDPL) are shown in "SampleStats.ipynb"
