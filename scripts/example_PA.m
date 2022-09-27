% example code to load and visualize results
% of bolus spreading from FVM simulations

% load network object for visualizing (optional)
%NT = NetworkObj('../examples/WT_COS7_KDEL_Cell3_um.net');

%% Load in mesh data
dirname = '../examples/';
name = 'example_PA';
filename = sprintf([dirname name '.mesh']);
MSH = MeshObj(filename);


%% load in snapshot data
[field,snaptimes] = loadSnapshotFVM([dirname name '.snap.txt']);
nsnap = length(snaptimes);

%% Show network snapshots over time
cmap = copper(200);

clear M

ct = 0;
figure
colormap(cmap)    
for sc = 1%:nsnap
    ct = ct+1;
    sc
       
    scatter(MSH.pos(:,1),MSH.pos(:,2),20,field(:,1,sc),'filled')    
    
    caxis([0,0.1])   
        
    set(gca,'Visible','off','Position',[0 0 1 1],'YDir','reverse')
    axis equal
        
    cb = colorbar;                    
    drawnow
    
%    M(ct) = getframe(gcf);

end