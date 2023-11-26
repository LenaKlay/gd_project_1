# Understanding the temporal and spatial spread of gene drive alleles through modeling

## Introduction

Artificial gene drive is a genetic engineering technology that could be used for the control of natural populations. Gene drive alleles bias their own transmission and can therefore spread in a population within a relatively small number of generations, even if they are deleterious. Understanding the potential outcomes of this technology, including the modification and/or the eradication of a natural population, is essential before real-world applications are considered. Here is the code we used to simulate the spatial and temporal spread of gene drive alleles in the population.

## Authors

I (Léna Kläy) wrote the code but this project was carried out in collaboration with Florence Débarre, Vincent Calvez and Léo Girardin.

## Contents

This repository contains the pdf of my first article, my PhD final presentation and two main folders:   

`Determinist` to simulate the deterministic models (PDE, reaction-diffusion), 

`Stochastic` to simulate the stochastic models (birth–death process with migration).  

They always follow the same organisation inside:

1) `Functions` contains the code to run the simulations (.py), as well as a README.rmd file for the explanations,

2) `Outputs` stores the results of the simulations,

3) `Illustrations` contains the important figures, sometimes improved with Inkscape.

4) `Migale` contains the code to run the heaviest simulations on the cluster Migale (INRAE, doi: 10.15454/1.5572390655343293E12) as well as some outputs of previous simulations.

(only in the `Determinist` folder:)

5) `Mathematica` contains the mathematica files used in preliminary mathematical analyses.

6) `Poster` contains some posters I presented during my PhD.

