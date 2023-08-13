Readme file for the paper

Bagged Pretested Portfolio Selection

Kazak and Pohlmeier, 2022

https://doi.org/10.1080/07350015.2022.2110880

The following repository contains MATLAB codes

-replication.m is the main file, which runs using the data in 100.mat to replicate the main simulaiton results of the paper

- f_in.m, f_out.m and f_bag.m are the functions which implement the investment strategy based on the in-sample CE pretesting, out-of-sample and bagged pretest estimator respectively

- deltalw_ce.m - funciton which computes standard errors for CE difference based on Delta Method from Ledoit and Wolf (2004)

-covCor.m is the funciton from https://www.econ.uzh.ch/en/people/faculty/wolf/publications.html#Programming_Code from the Ledoit Wolf (2004) (c) Olivier Ledoit and Michael Wolf
