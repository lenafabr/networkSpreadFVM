%% load in mesh output by fortran file

meshfile = '../../networkFVMsims/test/circpolyresv.mesh.txt'
MSH = MeshObj(meshfile);


%% load snapshots
snapfile = '../../networkFVMsims/test/circpolyresv.snap.txt';    
[fields,snaptimes,vels] = loadSnapshotFVM(snapfile);

%% plot field on triangular mesh
% for this to work, need to have matlab mesh object
% matched up to the one used by fortran
nm = size(mesh.Elements,2);

X = zeros(3,nm); Y = zeros(3,nm);
for cc = 1:nm
    X(:,cc) = mesh.Nodes(1,mesh.Elements(:,cc));
    Y(:,cc) = mesh.Nodes(2,mesh.Elements(:,cc));
end

delete(findobj(gca,'Type','patch'))
clear h
for sc = length(snaptimes)
    C = fields(1:nm,1,sc)/(pi*0.05^2);    
    % clear prior plots
    if (exist('h','var'))
        if (isvalid(h)); delete(h); end
    end
    
    % make new plots
    h = patch(X,Y,C)
    hold all
    % also draw scatterdots for network mesh elements (should turn these to
    % rectangular patches really
    ind = find(MSH.edgeind(:,1)>0);
    scatter(MSH.pos(ind,1),MSH.pos(ind,2),MSH.len(ind)*500,fields(ind,1,sc)','filled')
    hold off
    caxis([0,1])
    title(sprintf('Time %0.4f', snaptimes(sc)))
    drawnow
    
end
%trisurf(T,x,y,z)

%% look at calcium released over time
% original sims
dirname = '../results/netleak/struct2D/';
filename = [dirname 'sheetnuchex_Rpt29_mpt1_lowca_P20.out'];
[flux0,tvals0,cumflux0,tavg0] = loadTotFluxSim(filename); 
%% new sims with meshed reservoir
dirname = '../results/netleak/sheets/';
filename = [dirname 'circlenuchexmesh_lowca_P20.out'];
[flux,tvals,cumflux,tavg] = loadTotFluxSim(filename);   
%%
sclmM = 1e-3*6e23/1000*1e-12;% scale to convert mM to per um^3

plot(tavg0,cumflux0*sclmM,tavg,cumflux*sclmM)
xlabel('time')
ylabel('release')

%% read in a snapshot and get total amount of calcium in system
[field0,snaptimes0] = loadSnapshotFVM([dirname '../struct2D/sheetnuchex_Rpt29_mpt1_lowca_P20.snap.txt']);
MSH0 = MeshObj([dirname '../struct2D/sheetnuchex_Rpt29_mpt1_lowca_P20.mesh.txt']);
%%
tot0 = MSH0.len'*field0(:,1,1)

%%
[field,snaptimes] = loadSnapshotFVM([dirname '../sheets/circlenuchexmesh_lowca_P20.snap.txt']);
MSH = MeshObj([dirname '../sheets/circlenuchexmesh_lowca_P20.mesh.txt']);

%% compare volumes (excluding reservoir)
a=0.05;
ind0 = find(MSH0.resvind==0);
V0 = sum(MSH0.len(ind0))*pi*a^2

ind = find(MSH.resvind==0);
V = sum(MSH.len(ind))

%% compare total concentration (excluding reservoir)
tot0 = MSH0.len(ind0)'*field0(ind0,1,end)
tot = MSH.len(ind)'*field(ind,1,1)


%% plot field concentrations over time
% for this to work, need to have matlab mesh object
% matched up to the one used by fortran
nm = size(mesh.Elements,2);

X = zeros(3,nm); Y = zeros(3,nm);
for cc = 1:nm
    X(:,cc) = mesh.Nodes(1,mesh.Elements(:,cc));
    Y(:,cc) = mesh.Nodes(2,mesh.Elements(:,cc));
end

delete(findobj(gca,'Type','patch'))
clear h
for sc = length(snaptimes)
    C = field(1:nm,1,sc);
    % clear prior plots
    if (exist('h','var'))
        if (isvalid(h)); delete(h); end
    end
    
    % make new plots
    h = patch(X,Y,C); shading flat
    hold all
    % also draw scatterdots for network mesh elements (should turn these to
    % rectangular patches really
    ind = find(MSH.edgeind(:,1)>0);
    scatter(MSH.pos(ind,1),MSH.pos(ind,2),MSH.len(ind)*5000,field(ind,1,sc)','filled')
    hold off
    caxis([0,1])
    title(sprintf('Time %0.4f', snaptimes(sc)))
    drawnow
    
end
%trisurf(T,x,y,z)
