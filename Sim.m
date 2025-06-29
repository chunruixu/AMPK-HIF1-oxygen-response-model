
function [tout, yout] = Sim(TimeOfRun,param, y0) %input simulation time, parameters and y0

tspan = [0 TimeOfRun];
tstart = tspan(1);
tfinal = tspan(2);

[tout, yout] = ode15s(@(t,y) ODE_ROSAMPK(t,y,param), [tstart tfinal], y0);
