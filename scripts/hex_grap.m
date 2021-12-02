dir_name ='~/ER_Ca_rot/networkSpreadFVM-main/examples/density_op';
dirname ='~/ER_Ca_rot/networkSpreadFVM-main/examples/';
fil_n=1:6;
fil_arr = dir_name+"/hex_"+fil_n+"/example1_hex_"+fil_n+"_rd1.out";
cumflux = cell(6,1);
tavg_arr = cell(6,1);
for i = 1:6
    
    [fluxsim,tvals] = loadTotFluxSim(fil_arr(i));
    if (tvals(1)~=0)
        tvals = [0; tvals];
        fluxsim = [fluxsim(1); fluxsim];
    end
    vals = sum(fluxsim,2);
    integ = (vals(2:end)+vals(1:end-1))/2;
    dt = diff(tvals);
    cumflux{i,1} = size(cumsum(integ.*dt))
    tavg_arr{i,1} = size((tvals(2:end)+tvals(1:end-1))/2)
    fluxsim0 = fluxsim; tvals0 = tvals;
end

scl = 1e-3*6e23/(1e3*1e12)/1e7;
for i =1:6
    cumflux{i,1} =cumflux{i,1}*scl;
end

%plot(tavg_arr{1,1},cumflux{1,1},tavg_arr{2,1},cumflux{2,1},tavg_arr{3,1},cumflux{3,1},tavg_arr{4,1},cumflux{4,1},tavg_arr{5,1},cumflux{5,1},tavg_arr{6,1},cumflux{6,1});

%plot(tavg_arr,cumflux);

%% look at network snapshots over time
% pull https://github.com/lenafabr/networktools to get the tools for
% working with networks
% ----------
%% load network
NT = NetworkObj([dirname 'hex_4_rd1.net'])
NT.edgelens = NT.edgevals;
NT.dim = 2;
NT.nodepos = NT.nodepos(:,1:2);
%NT.plotNetwork()

% find permeable region center
permcentnode = find(cellfun(@(x) strcmp(x,'P1'),NT.nodelabels));

%% load network mesh
filename = 'example1_hex_4_rd1.mesh.txt';
filename = [dirname filename];
MSH = MeshObj(filename);

% get total length of mesh cells within the permeable region
diffs = MSH.pos - NT.nodepos(permcentnode,:);
dists = sqrt(sum(diffs.^2,2));
permind = find(dists<3);
permlen = sum(MSH.len(permind))

%% load snapshot data
[field,snaptimes,vels] = loadSnapshotFVM([dirname 'example1_hex_4_rd1.snap.txt']);
nsnap = length(snaptimes);

%% show snapshot evolution
clear M
for sc = 1:nsnap
    [sc nsnap]
    scatter(MSH.pos(:,1),MSH.pos(:,2),30,field(:,1,sc)/0.007854,'filled')   
    caxis([0,1])
    hold all
    colormap copper
    set(gcf,'Color','w')
    %colorbar
    hold off
    title(sprintf('Time %f',snaptimes(sc)))
    set(gca,'FontSize',14,'Visible','off','Position',[0.03,0.02,0.92,0.9])
    
    % draw release region
    pt = NT.nodepos(permcentnode,:);
    hold all
    %plot(pt(1),pt(2),'r*')
    th = linspace(0,2*pi,50);
    circ = pt + 3*[sin(th') cos(th')];
    plot(circ(:,1),circ(:,2),'g','LineWidth',3)
    hold off
    
    text(16,19,sprintf('Time $%0.1f$ sec',snaptimes(sc)),'FontSize',16,'Interpreter','latex')   
    cb=colorbar;
    cb.TickLabelInterpreter='latex';
    cb.FontSize = 14;    
    axis equal
   
    %text(-22,19,sprintf('asterB37dpt4Cpt6.1, Ltot %g, Lperm %g',sum(NT.edgelens),permlen), 'FontSize',12,'Interpreter','latex')
    drawnow       
    
   M(sc) = getframe(gcf);
end