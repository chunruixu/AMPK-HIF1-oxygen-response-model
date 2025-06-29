main.m:  run main.m for the simulation results of temporal plots for normoxia and different conditions (including hypoxia, hyperoxia, AMPK activating, and SIRT1 activating).

getParam.m and getMutantParam.m: load parameters for normal or specific situations (like hypoxia, AMPK activator)

Sim.m: get simulation results; input simulation time, parameters and initial values of variables; out time course simulated results

ODE_ROSAMPK.m: ordinary differential equations

PlotResult.m, PlotResult2.m, etc: plot temporal dynamics

PlotDrug.m: plot relationships between molecules when varying activity of AMPK or SIRT1

findSquares.m: scale the experimental data to the simulated results



%%%%%%%%%%%%%%%%Parameter Estimation%%%%%%%%%%%%%%%%%
RunParam.m: a main file for parameterization
getCost.m: compute cost
ParameterEstimationGA.m: use the genetic algorithm to get a new set of parameters