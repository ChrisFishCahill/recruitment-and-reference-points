# Course material for the Quantitative Fisheries Center's recruitment and reference points workshop

This repository contains all scripts and lecture materials for Cahill's portion of the course on reference points and feedback control held in **Duluth, Minnesota**, from **23–25 June 2025**. 

# Fmsy folder

This folder contains scripts for equilibrium Fmsy calculations based on Botsford incidence functions and simulation based Fmsy calculations where simulations account for variable recruitment.

# Depensation folder

There are examples in this folder of two different ways to specify depensatory stock-recruitment dynamics, including an attempt to fit those models to the Escanaba Lake Walleye dataset using \`optim()\`

# HCR folder

Harvest control rule folder. This folder contains an omniscient manager script that was constructed using RTMB (`om.R`) and contains four additional approximation in policy space scripts. Those additional scripts estimate linear or rectlinear feedback policies according to output vs. input (i.e., direct exploitation management or U-control) control for fisheries management systems. These breakpoint models are approximated using softplus and softmin approximations as they are non-differentiable at the breakpoints.

# EIV folder

We did not explicitly cover the code in this folder during class, but it simulates the effect of errors-in-variables (i.e., error in the measure of S) on the estimation of stock-recruitment parameters.

# Extended Ricker folder

We also did not cover this code in class, but did discuss these topics. This folder just contains another way to model environmental effects when using the Ricker stock-recruitment function.
