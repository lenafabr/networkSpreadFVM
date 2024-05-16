% create rectangular patches for displaying a meshed network

function [X,Y,patchind] = meshPatches(MSH,NT,radii,extend)
% convert a mesh of a network into coordinates for plotting
% patch rectangles along the edges
% patchind gives the order of mesh indices for the patches
% to get the right field values for the color

% tubeR = is the half-width of the patch objects (tube radius)
% extend = extend all edges outward by this length

if (~exist('extend','var'))
    extend = 0;
end

%%
ct = 0; X = []; Y = []; patchind = []; 
for ec = 1:NT.nedge
    
    n1 = NT.edgenodes(ec,1); n2 = NT.edgenodes(ec,2);    
    edgepath = NT.edgepath{ec};
    
    % get ordered mesh cells for this edge
    mind = find(MSH.edgeind(:,1)==ec);
    [a,b] = sort(MSH.edgeind(mind,2));
    mind = mind(b);
    % linear displacement along linearized edge
    dv = NT.nodepos(n2,:) - NT.nodepos(n1,:);
    dn = dv/norm(dv)*NT.edgelens(ec);
    %MSH.pos(mind(end),:)-NT.nodepos(n1,:);
    %if (isempty(mind)); continue; end
    
    % get fraction along edge for each mesh cell    
    if (abs(dn(1))<eps)
        frac = (MSH.pos(mind,2)-NT.nodepos(n1,2))./dn(2);
    else
        frac = (MSH.pos(mind,1)-NT.nodepos(n1,1))./dn(1);
    end
    
    startfrac = [0; (frac(2:end)+frac(1:end-1))/2];
    endfrac = [(frac(2:end)+frac(1:end-1))/2; 1];
    
    [param,lens] = arclenparam(edgepath');
    
    startpos = interp1(NT.cumedgelen{ec},NT.edgepath{ec},startfrac*NT.cumedgelen{ec}(end));
    endpos = interp1(NT.cumedgelen{ec},NT.edgepath{ec},endfrac*NT.cumedgelen{ec}(end));
    
    if (length(radii)>1)
        tubeR = radii(ec);  
    else
        tubeR = radii;
    end

    % plot the edge mesh segments    
    for ic = 1:size(startpos,1)
        % make rectangle
        dv = endpos(ic,:)-startpos(ic,:); dv = dv/norm(dv);
        pv = [-dv(2) dv(1)]; pv = pv/norm(pv);
        
        vert = [startpos(ic,:)+pv*tubeR-dv*extend; startpos(ic,:)-pv*tubeR-dv*extend; ...
            endpos(ic,:)-pv*tubeR+dv*extend; endpos(ic,:)+pv*tubeR+dv*extend];
        
        X(:,end+1) = vert(:,1);
        Y(:,end+1) = vert(:,2);    
        if (any(isnan(X(:))))
            error('NaN in mesh patches')
        end
        %plot([startpos(ic,1) endpos(ic,1)],[startpos(ic,2) endpos(ic,2)],'go-')                
    end
    if (size(startpos,1)~=length(mind))
        error('something wrong')
    end
    patchind(end+1:end+length(mind)) = mind;    
end
