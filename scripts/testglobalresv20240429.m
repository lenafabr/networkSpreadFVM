% load in and visualize snapshots

NT = NetworkObj('../test/circlenuchexresv.net',struct('dim',2));
NT.setCumEdgeLen(1:NT.nedge);
%% load snapshot data

resdir = '../test/';
filename = [resdir 'test.mesh.txt'];
MSH = MeshObj(filename);

snapfile = [resdir 'test.snap.txt'];
[field,snaptimes,vels] = loadSnapshotFVM(snapfile);
nsnap = length(snaptimes)

%[~,centcell] = min(field(:,1,2))
% TODO:
% FIGURE out why different mesh cells are getting different concentrations
% with global recovery.
% figure out why total ER concentration is not asymptoting correctly

%% show the snapshot
showtime = 2; % time at which to show in sec
 [~,showsnap] = min(abs(snaptimes - showtime));
 Rabs = 0.25;
 tubeR = 0.05;
 clf
for sc = nsnap
    [sc nsnap]
    ind = find(MSH.resvind==0);
    ind(ind==MSH.globalresvind) = [];
    %intpos = interpolateMeshPos(MSH,NT);    
    scatter(MSH.pos(ind,1),MSH.pos(ind,2),15,field(ind,1,sc)/0.007854,'filled')   
    %[X,Y,patchind] = meshPatches(MSH,NT,tubeR,tubeR);
    %C = field(patchind,1,sc)/(pi*tubeR^2); % colors
    %patch(X,Y,C);
    %shading flat    
    hold all    
%     cent =MSH.pos(centcell,:);
%     th=linspace(0,2*pi,100)';
%     circ = cent+Rabs*[cos(th) sin(th)];
%     plot(circ(:,1),circ(:,2),'g-','LineWidth',2)
    colormap copper
    
    % show nucleus
    nucind = find(MSH.resvind>0); 
    nuccent = MSH.pos(nucind,:);
    th = linspace(0,2*pi)';
    nucpts = nuccent + 7.5*[cos(th) sin(th)];           
    cmap = colormap(gca);
    nuccol = interp1(linspace(0,1,size(cmap,1)),cmap,field(nucind,1,sc)/(pi*0.05^2),'linear','extrap');
    fill(nucpts(:,1),nucpts(:,2),nuccol,'LineStyle','none','FaceAlpha',1)

    %title(sprintf('Snap %d time %f', sc, snaptimes(sc)))
    % title('WT region 1')    
   
    set(gca,'Visible','off')
    
    caxis([0,0.5])
    colorbar 
    plot_cleanup(gca,'FontSize',14,'pcolor',true)
    
    axis equal       
    
end
hold off

%% plot nucleus conc over time
scl = 1/(pi*0.05^2);
plot(snaptimes, squeeze(field(nucind,1,:))*scl,'.-')