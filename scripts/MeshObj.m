classdef MeshObj < handle
    % object defining a FVM mesh
properties                               
    ncell
    dim
    bounds    
    len
    pos
    deg
    maxdeg
    resvind % reservoir index
    nodeind % which node index
    edgeind % which edge index
end

methods
    function MSH = MeshObj(fname)
        % create a network object
        % optionally, load data from file
        
        MSH.ncell = 0; % number of nodes        
        MSH.dim = 2; % spatial dimension        
        MSH.maxdeg = 10; % maximum allowed degree              
        
        if (nargin>0)
            MSH.loadMesh(fname);
        end
    end        
    
    function loadMesh(MSH, filename)
        % load in mesh object from file        
        data = dlmread(filename);
        ncell = data(1,1);
        dim = data(1,2);
        
        MSH.ncell = ncell;
        MSH.dim = dim;
        MSH.pos = data(2:1+ncell,2:dim+1);
        MSH.len = data(2:1+ncell,dim+2);
        MSH.deg = data(2:ncell+1,dim+3);
        MSH.maxdeg = max(MSH.deg);
        
        % Load in connectivity info and lengths
        MSH.bounds = zeros(ncell,MSH.maxdeg);   
        MSH.edgeind = zeros(ncell,2);
        for cc = 1:MSH.ncell
            deg = MSH.deg(cc);
            MSH.bounds(cc,1:deg) = data(cc+1,dim+4:dim+3+deg);  
            if (size(data,2)>dim+deg+3)
                if (size(data,2)>dim+deg+4) % nodeind and edgeind included
                    MSH.nodeind(cc) = data(cc+1,dim+deg+4);
                    MSH.edgeind(cc,:) = data(cc+1,dim+deg+(5:6));
                    MSH.resvind(cc) = data(cc+1,dim+deg+7);
                else
                    % no node/edge info
                    MSH.nodeind(cc) = 0; MSH.edgeind(cc,:) = 0;
                    MSH.resvind(cc) = data(cc+1,dim+deg+4);
                end
            else
                MSH.nodeind(cc) = 0; MSH.edgeind(cc,:) = 0;
                MSH.resvind(cc) = 0; % no reservoir info
            end
        end
    end
    
    function plotField(MSH,field,climits)
        % plot field value on the mesh as colored lines
        
        colormap(copper)
        caxis(climits)
        cmap = copper;
        
        allX1 = []; allX2 = []; allY1 = []; allY2 = []; allcvals=[];
        for bct = 1:MSH.maxdeg
            
            % direction to the bct-th boundary for each node
            blist = MSH.bounds(:,bct);
            ind = find(blist>0);
            bdir = MSH.pos(blist(ind),:)-MSH.pos(ind,:);
            bdir = bdir./sqrt(sum(bdir.^2,2));
            %
            X1 = MSH.pos(ind,1);
            X2 = MSH.pos(ind,1)+MSH.len(ind)/2.*bdir(:,1);
            Y1 = MSH.pos(ind,2);
            Y2 = MSH.pos(ind,2)+MSH.len(ind)/2.*bdir(:,2);
            
            allX1 = [allX1; X1];
            allX2 = [allX2; X2];
            allY1 = [allY1; Y1];
            allY2 = [allY2; Y2];
            
            
            clamped = field(ind);
            clamped(clamped<climits(1)) = climits(1);
            clamped(clamped>climits(2)) = climits(2);
            cvals = interp1(linspace(climits(1),climits(2),size(cmap,1)), cmap, clamped);
            allcvals = [allcvals; cvals];
        end
           
        set(gca,'ColorOrder',allcvals,'NextPlot', 'replacechildren');
           
        plot([allX1 allX2]', [allY1 allY2]','LineWidth',3)
            
        axis equal
    end    
end
end