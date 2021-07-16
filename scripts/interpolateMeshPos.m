function intpos = interpolateMeshPos(MSH,NT)
% get mesh cell positions interpolating over network edge paths

intpos = zeros(MSH.ncell,NT.dim);
for mc = 1:MSH.ncell
    %%
    ec = MSH.edgeind(mc,1);
    nc = MSH.nodeind(mc);
    if (nc>0) % nodal mesh cell
        intpos(mc,:) = NT.nodepos(MSH.nodeind(mc),:);    
    elseif (ec>0) % edge mesh cell
        pos = MSH.pos(mc,:);
        n1 = NT.edgenodes(ec,1); n2 = NT.edgenodes(ec,2);
        dn = NT.nodepos(n2,:)-NT.nodepos(n1,:);
        
        % what fraction along the edge is this mesh cell?
        frac = (pos(1)-NT.nodepos(n1,1))./dn(1);
        
        intpos(mc,:) = interp1(NT.cumedgelen{ec},NT.edgepath{ec},frac*NT.cumedgelen{ec}(end));
    else % reservoir cell
        intpos(mc,:) = MSH.pos(mc,:);
    end
end