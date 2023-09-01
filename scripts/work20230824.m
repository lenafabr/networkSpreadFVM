%% draw a regular polygon to approximate a circle:
nside=5;
pgon = nsidedpoly(nside);

%%
model = createpde("thermal","transient")
%%
model = createpde;
g = geometryFromEdges(model,@circleg);
pdegplot(model)

%%
mesh = generateMesh(model)

% ----------------
%% Create a decomposed geometry matrix, with line segments giving boundary of polygon
% documented here: 
% https://www.mathworks.com/help/pde/ug/create-geometry-at-the-command-line.html

% set up the polygon
nside=5;
pgon = nsidedpoly(nside);

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
    verts = [verts; coords; vertend];
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
pdegplot(g)
%% create a triangular mesh object
mesh = generateMesh(model,'Hmax',1,'Hmin',0.05,'GeometricOrder','linear');
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

            % area of the element
            edge.area = resvheight*norm(mesh.Nodes(:,i1)-mesh.Nodes(:,i2));

            conmat(i1,i2) = edgect;
            conmat(i2,i1) = edgect;

            edgelist(edgect) = edge;
        else
            % this edge already exists, save the other element around it
            ec = conmat(i1,i2);
            edgelist(ec).elem(2) = mc;
        end    
    end
end

%% Output mesh information to a text file
% each line lists:
% type of object, info
% type 0 = vertex, info = coordinates x,y,z
% type 1 = edge, info = vertex i1, i2, mesh elements 1,2, area
% type 2 = mesh element, info = vertex i1,i2,i3, edges e1,e2,e3, volume,SA

fname = ''
OF = open('')
for vc = 1:nvert

end