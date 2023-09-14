% load original network

NT = NetworkObj('../networks/circlenuchexresv_N20.net',struct('dim',2))

%%
NT.plotNetwork()
hold all
%% set up and mesh single flat perinuclear sheet
% pick out reservoir nodes
resvnodes = [];
for nc = 1:NT.nnode
    if (~isempty(NT.nodelabels{nc}))
        if (NT.nodelabels{nc}(1)=='R')
            resvnodes(end+1) = nc;
        end
    end
end

resvnodesnuc = resvnodes;
%%
resvnodepos = NT.nodepos(resvnodes,:)
cent = mean(resvnodepos);

rnuc = 7; % desired radius

rvec = resvnodepos - cent;
% order nodes by angle
nodeang = atan2(rvec(:,2),rvec(:,1));
[~,ind] = sort(nodeang,'descend');
resvnodes = resvnodes(ind);
resvnodepos = resvnodepos(ind,:);
rvec = rvec(ind,:);


% adjust each connecting node to the right radius
for nc = 1:length(resvnodes)
     nr = norm(rvec(nc,:));
     rr = rvec(nc,:)/nr;
     dr = rnuc-nr;
    newpos(nc,:) = NT.nodepos(resvnodes(nc),:) + rr*dr;
end
%newpos = NT.nodepos(resvnodes,:);

NT.nodepos(resvnodes,:) = newpos;
resvnodepos = newpos;
NT.interpolateEdgePaths(2)

hold all
plot(newpos(:,1),newpos(:,2),'go-')
%%
rtube =0.05;

% set up vertices, including short edges at node connections
dv = resvnodepos(1,:) - resvnodepos(end,:);
dvn = dv/norm(dv);
verts = resvnodepos(1,:) - dvn*rtube;
for ic = 2:size(resvnodepos,1)
    dv = resvnodepos(ic,:) - resvnodepos(ic-1,:);
    dvn = dv/norm(dv);
    verts(end+1,:) = resvnodepos(ic-1,:) + dvn*rtube;
    verts(end+1,:) = resvnodepos(ic,:) - dvn*rtube;
end
verts(end+1,:) = resvnodepos(end,:) + dvn*rtube;

hold off
NT.plotNetwork()
hold all

plot(verts(:,1),verts(:,2),'.-')

vertsnuc = verts;
%% make the polygon in pde modeler
pdepoly(verts(:,1),verts(:,2))

%% place randomly positioned peripheral sheets
Rcell = 15;
avgRsheet = 2; % average sheet radius
sigRsheet = avgRsheet*0.1; % standard dev
nsheet = 10;

allsheetcent = [];
for sc = 1:nsheet
    while 1
        r = sqrt(rand())*Rcell;
        th = rand()*2*pi;

        Rsheet = randn()*sigRsheet + avgRsheet;
        % try again if too small
        if (Rsheet<0.2); continue; end

        sheetcent = [r*cos(th),r*sin(th)];

        % check for overlap with nucleus
        d2 = sum((verts - sheetcent).^2,2);
        nucoverlap = min(d2)<Rsheet^2;
        inp = inpolygon(sheetcent(1), sheetcent(2), verts(:,1),verts(:,2));

        % check for overlap with other sheets
        if (isempty(allsheetcent))
            overlap = false;
        else
            d2 = sum((sheetcent - allsheetcent(:,1:2)).^2,2);      
            overlap = any(d2<(Rsheet + allsheetcent(:,3)).^2);
        end

        if (~nucoverlap & ~inp & ~overlap); break; end
    end

    allsheetcent(end+1,:) = [sheetcent Rsheet];
    pdecirc(sheetcent(1), sheetcent(2), Rsheet);
end

%%
% hit export from graphical pdemodeler
% plot to check
[dl,bt] = decsg(gd);


model = createpde('thermal','transient');
%g = decsg([gd(:,1:2) gsheet],sf,ns)
geometryFromEdges(model,dl)
pdegplot(dl,'FaceLabels','on','EdgeLabels','on')

%%
save('../networks/circlenuchexmesh_sheets10.mat')


%% plot network and mesh
NT.plotNetwork()
hold all
%pdeplot(mesh)
pdegplot(dl,'FaceLabels','on')
hold off

%%
NT0 = copy(NT);
%% Restructure network around the sheets

ct = 0;
newnodes = []; newcon = {};
dokeepnode = true(1,NT.nnode);

for sc = 1:size(allsheetcent,1)
    %%
    sheet2node{sc} = [];
    
    sheetcent = allsheetcent(sc,1:2);
    Rsheet = allsheetcent(sc,3);

    % find nodes within this sheet
    dd = sqrt(sum((NT.nodepos - sheetcent).^2,2));
    nodeisin = (dd<Rsheet);

    nodelist = find(nodeisin);
    for ic = 1:length(nodelist)
        nc = nodelist(ic);
        
        dokeepnode(nc) = false; % will drop nodes inside the sheet
        
        p1 = NT.nodepos(nc,:);
        neighb = NT.nodenodes(nc,1:NT.degrees(nc));

        for ic2 = 1:length(neighb)
            n2 = neighb(ic2);
            % only look at neighbors outside the circle
            if (nodeisin(n2)); continue; end
            p2 = NT.nodepos(n2,:);

            % edge going to this neighbor
            ec = NT.nodeedges(nc,ic2);
            
            % get intersection point
            [z,pint] = segcircintersect(p1,p2,sheetcent,Rsheet);
            
            % do not make connections that would lead to very short edges
            if (z > 0.9); continue; end
            
            hold all
            plot(pint(1),pint(2),'co')
            
            % create a new node and new connection to it
            ct = ct+1;
            newnodes(ct,:) = pint;           
            newcon{ct} = [n2];       
            node2sheet(ct+NT.nnode) = sc;
            sheet2node{sc} = [sheet2node{sc} ct+NT.nnode];
        end
    end
end
%% add new nodes to the network

nnprev = NT.nnode;
NT.addNodes(newnodes,newcon);
NT.interpolateEdgePaths(2);
nn = NT.nnode;

allkeepnode = true(1,NT.nnode);
allkeepnode(find(~dokeepnode)) = false;
%% remove internal nodes
mapold2new = NT.keepNodes(find(allkeepnode));

resvnodesnuc = mapold2new(resvnodesnuc);
%%
for sc = 1:length(sheet2node)
    sheet2node{sc} = mapold2new(sheet2node{sc});
end

%% plot network and sheet circles
NT.plotNetwork()
hold all
pdegplot(dl,'FaceLabels','on')
hold off

%% redefine sheet polygons around network nodes
%%
rtube =0.05;

for sc = 1:size(allsheetcent,1)
    resvnodes = sheet2node{sc};
    resvnodepos = NT.nodepos(resvnodes,:);
    
    
    rvec = resvnodepos - allsheetcent(sc,1:2);
    % order nodes by angle
    nodeang = atan2(rvec(:,2),rvec(:,1));
    [nodeangsort,ind] = sort(nodeang,'descend');
    resvnodes = resvnodes(ind);
    resvnodepos = resvnodepos(ind,:);
    rvec = rvec(ind,:);
    
    
    % set up vertices, including short edges at node connections
    dv = resvnodepos(1,:) - resvnodepos(end,:);
    dvn = dv/norm(dv);
    verts = resvnodepos(1,:) - dvn*rtube;
    for ic = 2:size(resvnodepos,1)
        dv = resvnodepos(ic,:) - resvnodepos(ic-1,:);
        dvn = dv/norm(dv);
        
        newpt = resvnodepos(ic-1,:) + dvn*rtube;
        newrvec = newpt- allsheetcent(sc,1:2);
        newang = atan2(newrvec(2),newrvec(1));
        if (newang<nodeangsort(ic))
            verts(end+1,:) = resvnodepos(ic-1,:)*2/3+resvnodepos(ic,:)*1/3;
            verts(end+1,:) = resvnodepos(ic-1,:)*1/3+resvnodepos(ic,:)*2/3;
        else        
            verts(end+1,:) = resvnodepos(ic-1,:) + dvn*rtube;
            verts(end+1,:) = resvnodepos(ic,:) - dvn*rtube;
        end
    end
    verts(end+1,:) = resvnodepos(end,:) + dvn*rtube;

    sheetverts{sc} = verts;
end

%%
pdepoly(vertsnuc(:,1),vertsnuc(:,2))

%
for sc = 1:length(sheetverts)
    pdepoly(sheetverts{sc}(:,1),sheetverts{sc}(:,2))
end

%%
% hit export from graphical pdemodeler
% plot to check
[dl,bt] = decsg(gd);


model = createpde('thermal','transient');
%g = decsg([gd(:,1:2) gsheet],sf,ns)
geometryFromEdges(model,dl)

NT.plotNetwork()
hold all
pdegplot(dl,'FaceLabels','on','EdgeLabels','on')
hold off


%% make the mesh
mesh = generateMesh(model,'Hmax',0.3,'Hmin',0.2,'GeometricOrder','linear');
NT.plotNetwork()
hold all
pdeplot(mesh)
hold off

%% get list of all edges in the triangular mesh object
nvert = size(mesh.Nodes,2); % number of vertices in the mesh
% connectivity matrix for vertices
% stores 0 if two vertices are not connected
% and the index of the edge if they are
conmat = zeros(nvert,nvert); 

% thickness of the reservoir, for defining boundary areas
resvheight = 0.06;

% info for mesh edges:
% for each edge, this stores the 2 vertices of the edge, and the 2 elements
% on either side of it

% areas of mesh elements:
[A,AE] =area(mesh);

clear edgelist
edgect = 0;
for mc = 1:size(mesh.Elements,2)
    ind = mesh.Elements(:,mc);
    ind = [ind; ind(1)];

    for vc = 1:3
        i1 = ind(vc); i2 = ind(vc+1);
        if conmat(i1,i2) == 0
            % create a new edge
            edgect = edgect+1;
            edge = struct();
            % vertices for this edge
            edge.vert = [i1,i2];

            % elements surrounding this edge
            edge.elem = [mc,0];

            % area of the edge
            edge.area = resvheight*norm(mesh.Nodes(:,i1)-mesh.Nodes(:,i2));

            conmat(i1,i2) = edgect;
            conmat(i2,i1) = edgect;

            edgelist(edgect) = edge;

            % keep track of edges around each reservoir
            resvedges(mc,vc) = edgect;
        else
            % this edge already exists, save the other element around it
            ec = conmat(i1,i2);
            edgelist(ec).elem(2) = mc;
            resvedges(mc,vc) = ec;
        end    
    end
    
    % compute volume of the reservoir
    resvvol(mc) = resvheight*AE(mc);

    % surface area (both faces of triangle)
    resvSA(mc) = 2*AE(mc);

    % compute effective length used to calculate escape into tubules
    % for sheets: L = r^2/h*ln(R/r) (h=sheet thickness, pi*R^2 = A = sheet area)
    reff = sqrt(AE(mc)/pi);
    rtube = 0.05; % tubule radius
    resvleneff(mc) = rtube^2/resvheight*log(reff/rtube);
end

%% make list of reservoir indices to be connected to tubules

% remove old node labels
for nc = 1:NT.nnode
    NT.nodelabels{nc} = '';
end

% tubecon resv lists the reservoir and mesh edge index for each connection
tubeconresv = []; nodeconresv = []; nodeconedge = [];
ct=0;
for sc= 0:size(allsheetcent,1)
    if (sc==0)
        resvnodes = resvnodesnuc;
    else
        resvnodes = sheet2node{sc};
    end
    resvnodepos = NT.nodepos(resvnodes,:);
    
    for ic = 1:size(resvnodepos,1)
        % mesh nodes associated with this connection
        nid = findNodes(mesh,'radius',resvnodepos(ic,:),rtube*1.001);
        
        % mesh edge connected to these nodes
        eid = conmat(nid(1),nid(2));
        rid = edgelist(eid).elem;
        % reservoir id adjacent to this edge
        rid = max(rid);
        
        if (~isempty(tubeconresv))
            if (any(tubeconresv(2,:)==eid))
                % avoid duplicate connections through a single edge
                continue
            end
        end
        
        ct=ct+1;
        tubeconresv(1,ct) = rid;
        tubeconresv(2,ct) = eid;
        
        
        
        nc = resvnodes(ic); % network node index
        
        % for this node, list connecting reservoir and edge
        nodeconresv(nc) = rid;
        nodeconedge(nc) = eid;
        
        % update reservoir label on network node
        NT.nodelabels{nc} = sprintf('R%d',rid);
    end
end

pdemesh(model)
%NT.plotNetwork()
hold all
pdemesh(mesh.Nodes,mesh.Elements(:,tubeconresv(1,:)),"EdgeColor","green")
hold off


%% Output mesh information to a text file
% each line lists:
% type of object, info
% type V = vertex, info = coordinates x,y,z
% type E = edge, info = vertex i1, i2, mesh elements 1,2, area
% type R = mesh element, info = vertex i1,i2,i3, edges e1,e2,e3, volume,SA

fname = '../networks/circlenuchexmesh_sheets10_resv.txt';
OF = fopen(fname,'w');

fprintf(OF,'%s\n%s\n\n','# meshed pentagon reservoirs','# made with polygon_reservoir_mesh')

% output vertices
for vc = 1:nvert
    fprintf(OF,'%s %d %15.10f %15.10f\n','VERT', vc, mesh.Nodes(:,vc)');
end

% output edges (faces if we were in 3D)
for ec = 1:length(edgelist)
    edge = edgelist(ec);
    fprintf(OF,'%s %d %d %d %d %d %15.10f\n','EDGE', ec, edge.vert,edge.elem,edge.area);
end

% output reservoirs
% index, vertices, edges, volume, SA, effective len (not used)
for rc = 1:size(mesh.Elements,2)
    fprintf(OF,'%s %d %d %d %d %d %d %d %15.10f %15.10f %15.10f\n','RESV', ...
        rc, mesh.Elements(:,rc)',resvedges(rc,:),resvvol(mc),resvSA(mc),resvleneff(mc));
end


% output connections to network nodes
ind = find(nodeconresv>0); % nodes that have a reservoir connected
for ic = 1:length(ind)
    nc = ind(ic); % which node
    rc = nodeconresv(nc); % which reservoir
    ec = nodeconedge(nc); % which edge

    fprintf(OF,'%s %d %d %d\n','NODECON', ...
        nc, rc, ec)
end
%
fclose(OF)


%% output network
NT.outputNetwork('../networks/circlenuchexmesh_sheets10.net')

