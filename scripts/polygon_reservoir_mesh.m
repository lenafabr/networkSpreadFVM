% this script builds a bunch of connected reservoirs
% representing a meshing of a polygonal region


% ----------------
%% Create a decomposed geometry matrix, with line segments giving boundary of polygon
% documented here: 
% https://www.mathworks.com/help/pde/ug/create-geometry-at-the-command-line.html

% set up the polygon
nside=100;
Rpoly = 3; % polygon radius
pgon = nsidedpoly(nside);
pgon.Vertices = pgon.Vertices*Rpoly;

% break up some edges of the polygon
rad = 0.01; 
verts = pgon.Vertices(1,:);
for sc = 1:nside
    if (sc==nside) 
        vertend = pgon.Vertices(1,:);
    else
        vertend = pgon.Vertices(sc+1,:);
    end
    cent = (pgon.Vertices(sc,:)+vertend)/2
    evec = vertend-pgon.Vertices(sc,:); evec = evec/norm(evec);
    coords = [cent-evec*rad; cent+evec*rad];
    if (mod(sc,20)==0)
        verts = [verts; coords; vertend];
    else
        verts = [verts; vertend];
    end
end
verts = verts(1:end-1,:);
plot(verts(:,1),verts(:,2),'.-')
%%

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
pdegplot(g,'FaceLabels','on','EdgeLabels','on')

%% create a triangular mesh object
mesh = generateMesh(model,'Hmax',0.1,'Hmin',0.1,'GeometricOrder','linear');
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

%% Output mesh information to a text file
% each line lists:
% type of object, info
% type V = vertex, info = coordinates x,y,z
% type E = edge, info = vertex i1, i2, mesh elements 1,2, area
% type R = mesh element, info = vertex i1,i2,i3, edges e1,e2,e3, volume,SA

fname = '../networks/circpolyresv.txt';
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

%% make list of reservoir indices to be connected to tubules
% tubecon resv lists the reservoir and mesh edge index for each connection
tubeconresv = [];
ct=0;

% get the shortest edges
d = diff(g.Vertices);
nd = sqrt(sum(d.^2,2));
minlen =  min(nd);
ind = find(nd <= minlen*1.0000001);
%
for sc = ind'%2:3:g.NumEdges
    % get mesh nodes associated with the short boundary edges
    nid = findNodes(mesh,'region',Edge=sc);
    % mesh edge for the short segment
    eid = conmat(nid(1),nid(2));
    rid = edgelist(eid).elem;
    % reservoir id adjacent to this edge
    rid = max(rid);

    ct=ct+1;
    tubeconresv(1,ct) = rid;
    tubeconresv(2,ct) = eid;
end

pdemesh(model)
hold all
pdemesh(mesh.Nodes,mesh.Elements(:,tubeconresv(1,:)),"EdgeColor","green")
hold off

%% make a minimal 'network' of independent short edges
% linked up to the mesh reservoirs
edgelen = 0.2;
nc=1; ec=1;
nodepos = []; edgenodes = []; noceconresv = []; nodeconedge = [];
for rc = 1:length(tubeconresv)
    eid = tubeconresv(2,rc);
    i1 = edgelist(eid).vert(1); i2 = edgelist(eid).vert(2);
    p1 = (mesh.Nodes(:,i1)+mesh.Nodes(:,i2))/2;

    % get perpendicular vector
    v = mesh.Nodes(:,i2)-mesh.Nodes(:,i1);
    v = v/norm(v); v = [v(2); -v(1)];
    p2 = p1+v*edgelen;
    
    nodepos(nc,:) = [p1'];     
    nodepos(nc+1,:) = [p2'];
    edgenodes(ec,:) = [nc nc+1];

    % label which reservoir the nodes belong to
    nodelabels{nc} = sprintf('R%d',tubeconresv(1,rc));
    nodelabels{nc+1} = '';

    % which reservoir, through which edge is this node connected
    nodeconresv(nc) = tubeconresv(1,rc);
    nodeconedge(nc) = eid;

    nc = nc+2;
    ec=ec+1;
end

NT = NetworkObj();
NT.nodepos = nodepos;
NT.edgenodes = edgenodes;
NT.setupNetwork();
%%
pdemesh(model)
hold all
pdemesh(mesh.Nodes,mesh.Elements(:,tubeconresv(1,:)),"EdgeColor","green")
NT.plotNetwork()
hold all

%% output network
options = struct();
options.nodelabels = nodelabels;
NT.outputNetwork('../networks/circlepolyresv_tubes.net',options);