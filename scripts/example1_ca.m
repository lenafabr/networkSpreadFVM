dirname = '~/proj/networkFVM/examples/';

%% total flux out of holes
filename = [dirname 'example1_ca.out'];

[fluxsim,tvals] = loadTotFluxSim(filename);
if (tvals(1)~=0)
    tvals = [0; tvals];
    fluxsim = [fluxsim(1); fluxsim];
end
vals = sum(fluxsim,2);
integ = (vals(2:end)+vals(1:end-1))/2;
dt = diff(tvals);
cumflux0 = cumsum(integ.*dt);
tavg0 = (tvals(2:end)+tvals(1:end-1))/2;
fluxsim0 = fluxsim; tvals0 = tvals;

%% plot cumulative flux
scl = 1e-3*6e23/(1e3*1e12)/1e7;

h0 = plot(tavg0,cumflux0*scl,'r-','LineWidth',3)
xlabel('time (s)')
ylabel('cumulative Ca ions released')
xlim([0,20])


% ---------
%% look at network snapshots over time
% pull https://github.com/lenafabr/networktools to get the tools for
% working with networks
% ----------
%% load network
NT = NetworkObj([dirname 'circlenuchexresv.net'])
NT.edgelens = NT.edgevals;
NT.dim = 2;
NT.nodepos = NT.nodepos(:,1:2);
NT.plotNetwork()

% find permeable region center
permcentnode = find(cellfun(@(x) strcmp(x,'P1'),NT.nodelabels));

%% load network mesh
filename = 'example1_ca.mesh.txt';
filename = [dirname filename];
MSH = MeshObj(filename);

% get total length of mesh cells within the permeable region
diffs = MSH.pos - NT.nodepos(permcentnode,:);
dists = sqrt(sum(diffs.^2,2));
permind = find(dists<3);
permlen = sum(MSH.len(permind))

%% load snapshot data
[field,snaptimes,vels] = loadSnapshotFVM([dirname 'example1_ca.snap.txt']);
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
