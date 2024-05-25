addpath('../../networktools')

%% load network data
NT = NetworkObj('../testing/linear2.net',struct('dim',2));
NT.interpolateEdgePaths(2);
NT.setCumEdgeLen();

%% load snapshot data
resdir = '../testing/'; 
runname = 'testspreadLinSin';
filename = [resdir runname '.mesh.txt'];
MSH = MeshObj(filename);

snapfile = [resdir runname '.snap.txt'];
[field,snaptimes,vels] = loadSnapshotFVM(snapfile);
nsnap = length(snaptimes);

%% show the snapshot
showtime = 2; % time at which to show in sec
 [~,showsnap] = min(abs(snaptimes - showtime));
 Rabs = 0.25;
 tubeR = MSH.rad*1;
 clf
for sc =100%:nsnap
    [sc nsnap]
    ind = find(MSH.resvind==0);    
    %intpos = interpolateMeshPos(MSH,NT);    
    %scatter(intpos(ind,1),intpos(ind,2),15,field(ind,1,sc)/0.007854,'filled')   
    %scatter(MSH.pos(ind,1),MSH.pos(ind,2),15,field(ind,1,sc)/0.007854,'filled')   
       
    [X,Y,patchind] = meshPatches(MSH,NT,tubeR,0);
    C = field(patchind,1,sc); % colors
    %C = field(patchind,1,sc)/(pi*tubeR^2); % colors
    patch(X,Y,C);
    shading flat    
   
    colormap copper
        
    caxis([0,0.05])

    %title(sprintf('Snap %d time %f', sc, snaptimes(sc)))
    % title('WT region 1')    
   
    %set(gca,'Visible','off')    
    colorbar; 
    plot_cleanup(gca,'FontSize',14,'pcolor',true)
    
    axis equal       
    drawnow
    
end
hold off

%% plot (3D) concentration profiles
% sort mesh cells by position
[xplot,sortind] = sort(MSH.pos(:,1));
meshvol = MSH.len.*MSH.rad'.^2*pi;

fieldplot = squeeze(field(sortind,1,2:2:50));
cmap = jet(size(fieldplot,2));

% 3D concentrations, do not need to be scaled by area to be flat at
% steady-state
plot(xplot,fieldplot,'--')
%plot(xplot,fieldplot./MSH.rad(sortind)'.^2/pi,'.-')
colororder(cmap)

xlabel('position')
ylabel('3D conc, mM')
%% load concs and check total concentration over time
snapfile = [resdir 'testreleaseLinNoBufEq.snap.txt'];
[field,snaptimes,vels] = loadSnapshotFVM(snapfile);
nsnap = length(snaptimes)

Kd = 0.001583D0; % dissociation constant
S = 0.021342D0; % 

% plot avg conc over time
scl = 1/(pi*0.05^2);
Ufield = squeeze(field(1:end,1,:));
Sfield = squeeze(field(1:end,2,:));
Cfield = Ufield.*(1 + Sfield./(Ufield+Kd));
avgU = sum(Ufield.*MSH.len(1:end),1)/sum(MSH.len);

% get total calcium
% total binding sites and total calcium
avgC = sum(Cfield.*MSH.len(1:end),1)/sum(MSH.len);
avgS = sum(Sfield.*MSH.len)/sum(MSH.len);

a = 0.05;
Vtot = sum(MSH.len(1:end-1));
Vg = MSH.len(end);
Atot = sum(MSH.len(1:end-1));

ind = find(MSH.nodeind>0 | MSH.edgeind(:,1)'>0);
Anuc = 0;
Atot = sum(MSH.len(ind))+Anuc;
Vnuc = 0;
Vtot = sum(MSH.len(ind))+Vnuc;

plot(snaptimes, avgC*scl,'.-')
xlabel('time (s)')
ylabel('ER total calcium')
legend off

% -------------
%% Test release from a permeable node

%% load snapshot data
resdir = '../testing/';
runname = 'testreleaseLinNoBuf'
filename = [resdir runname '.mesh.txt'];
MSH = MeshObj(filename);

snapfile = [resdir runname '.snap.txt'];
[field,snaptimes,vels] = loadSnapshotFVM(snapfile);
nsnap = length(snaptimes);

%% plot concentrations over time
[xplot,sortind] = sort(MSH.pos(:,1));
fieldplot = squeeze(field(sortind,1,2:2:50));
cmap = jet(size(fieldplot,2));

plot(xplot,fieldplot*scl,'.-')
colororder(cmap)
xlabel('position (um)')
ylabel('free calcium (mM)')
%% load and plot flux over time
[flux,tvals,cumflux,tavg] = loadTotFluxSim([resdir runname '.out']);    

% convert from mM to #/um
mMum = 1e-3*6e23*1e-3*1e-12;
plot(tvals,flux*mMum,'.-')
xlabel('time (sec)')
ylabel('flux (ions per sec)')

%%
plot(tavg,cumflux*mMum,'.-')
xlabel('time (sec)')
ylabel('cumulative released (ions)')

%% explicitly 3D concentrations
runname = 'testreleaseLin3DnoBufEq'
%% load and plot flux over time

[flux,tvals,cumflux,tavg] = loadTotFluxSim([resdir runname '.out']);    

plot(tvals,flux*mMum,'.-')
xlabel('time (sec)')
ylabel('flux (ions per sec)')

%%
plot(tavg,cumflux,'.-')
xlabel('time (sec)')
ylabel('cumulative released (ions)')
