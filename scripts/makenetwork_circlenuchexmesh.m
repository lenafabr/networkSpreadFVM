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

%% mesh the perinuclear region
gd = zeros(7,height(verts));
gd(1,:) = 2; % all elements are line segments
% starting x, y coordinates
gd([2,4],:) = verts(:,1:2)';
% ending x, y coordinates
gd([3,5],:) = [verts(2:end,1:2); verts(1,1:2)]';
% region to left of segment (outside, traverse polygon clockwise)
gd(6,:) = 0; % external
% region to right of segment
gd(7,:) = 1; % internal

model = createpde('thermal','transient');
g = geometryFromEdges(model,gd)
%pdegplot(g,'FaceLabels','on','EdgeLabels','on')

mesh = generateMesh(model,'Hmax',0.3,'Hmin',0.2,'GeometricOrder','linear');
pdeplot(mesh)
hold all
plot(verts(:,1),verts(:,2),'k.-','MarkerSize',20)
hold off

%% get list of all edges in the triangular mesh object
nvert = width(mesh.Nodes); % number of vertices in the mesh
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
for mc = 1:width(mesh.Elements)
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
% tubecon resv lists the reservoir and mesh edge index for each connection
tubeconresv = []; nodeconresv = []; nodeconedge = [];
ct=0;
for ic = 1:size(resvnodepos,1)
    % mesh nodes associated with this connection
    nid = findNodes(mesh,'radius',resvnodepos(ic,:),rtube*1.001);

    % mesh edge connected to these nodes
    eid = conmat(nid(1),nid(2));
    rid = edgelist(eid).elem;
    % reservoir id adjacent to this edge
    rid = max(rid);

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

pdemesh(model)
hold all
pdemesh(mesh.Nodes,mesh.Elements(:,tubeconresv(1,:)),"EdgeColor","green")
hold off

%% Output mesh information to a text file
% each line lists:
% type of object, info
% type V = vertex, info = coordinates x,y,z
% type E = edge, info = vertex i1, i2, mesh elements 1,2, area
% type R = mesh element, info = vertex i1,i2,i3, edges e1,e2,e3, volume,SA

fname = '../networks/circlenuchexmesh_resv.txt';
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
for rc = 1:width(mesh.Elements)
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
NT.outputNetwork('../networks/circlenuchexmesh.net')

