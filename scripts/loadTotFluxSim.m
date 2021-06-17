function [fluxsim,tvals] = loadTotFluxSim(filename)
% load flux from simulation results
% that output *total* flux only

data = dlmread(filename);

% is this flux out of permeable nodes
isperm = data(1,4);
% number of fields
nfield = data(1,2);

    %% read in flux from permeable or fixed nodes
    step = 1+nfield;
    tvals = data(1:step:end,1);
    fluxsim = data(2:step:end,2);    

end