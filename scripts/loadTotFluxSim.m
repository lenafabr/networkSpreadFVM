function [fluxsim,tvals,cumflux,tavg] = loadTotFluxSim(filename)
% load flux from simulation results
% that output *total* flux only
% also returns cumulative flux and average times corresponding to
% cumulative flux evaluation

data = dlmread(filename);

% is this flux out of permeable nodes
isperm = data(1,4);
% number of fields
nfield = data(1,2);

%% read in flux from permeable or fixed nodes
step = 1+nfield;
tvals = data(1:step:end,1);
fluxsim = data(2:step:end,2);

tvals = [0; tvals];
fluxsim = [fluxsim(1); fluxsim];
vals = sum(fluxsim,2);
integ = (vals(2:end)+vals(1:end-1))/2;
dt = diff(tvals);
cumflux = cumsum(integ.*dt);
%tavg = (tvals(2:end)+tvals(1:end-1))/2;
tavg = tvals(2:end);

end