addpath('../../networktools')

%% load network data
NT = NetworkObj('../testing/linear2.net',struct('dim',2));
NT.interpolateEdgePaths(2);
NT.setCumEdgeLen();

%% load snapshot data
resdir = '../testing/'; 
runname = 'testfptCone';
filename = [resdir runname '.mesh.txt'];
MSH = MeshObj(filename);

snapfile = [resdir runname '.snap.txt'];
[field,snaptimes,vels] = loadSnapshotFVM(snapfile);
nsnap = length(snaptimes);
% total starting mass
totstart = sum(field(:,1,1).*MSH.len)
%% show the snapshot
showtime = 2; % time at which to show in sec
 [~,showsnap] = min(abs(snaptimes - showtime));
 Rabs = 0.25;
 tubeR = 0.05;
 clf
for sc =4%:nsnap
    [sc nsnap]
    ind = find(MSH.resvind==0);    
    %intpos = interpolateMeshPos(MSH,NT);    
    %scatter(intpos(ind,1),intpos(ind,2),15,field(ind,1,sc)/0.007854,'filled')   
    %scatter(MSH.pos(ind,1),MSH.pos(ind,2),15,field(ind,1,sc)/0.007854,'filled')   
       
    [X,Y,patchind] = meshPatches(MSH,NT,tubeR,tubeR);
    C = field(patchind,1,sc); % colors
    patch(X,Y,C);
    shading flat    
   
    colormap copper
        
    caxis([0,0.3])

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
meshvol = MSH.len

fieldplot = squeeze(field(sortind,1,1));
cmap = jet(size(fieldplot,2));

% 3D concentrations, do not need to be scaled by area to be flat at
% steady-state
plot(xplot,fieldplot,'.-')
%plot(xplot,fieldplot./MSH.rad(sortind)'.^2/pi,'.-')
colororder(cmap)

xlabel('position')
ylabel('3D conc, mM')

%%
runname = 'testfptCone'
%% load and plot flux over time

[flux,tvals,cumflux,tavg] = loadTotFluxSim([resdir runname '.out']);    

%plot(tvals,flux*mMum,'.-')
xlabel('time (sec)')
ylabel('flux (ions per sec)')

%%
% going from narrow to wid
[flux,tvals,cumflux,tavg] = loadTotFluxSim([resdir 'testfptCone.out']);    
% going from wide to narrow
[flux2,tvals2,cumflux2,tavg2] = loadTotFluxSim([resdir 'testfptConeW2N.out']);    

plot(tavg,cumflux,'.-', tavg2, cumflux2, '.-')
xlabel('time (sec)')
ylabel('cumulative released (ions)')


% estimate MFPT
tot = sum(MSH.len.*field(:,1,1))
dt = tavg(2)-tavg(1);
avghit = sum((2-cumflux(2:end)-cumflux(1:end-1))/2)*dt
avghitW2N = sum((2-cumflux2(2:end)-cumflux2(1:end-1))/2)*dt

%% compare to berezhkovsii solution
L = 2-0.025; 
lam = 0.5;
D = 1/sqrt(1+lam^2);
taunw = L^2/6/D*(3+lam*L)/(1+lam*L)
tauwn = L^2/6/D*(3+2*lam*L)