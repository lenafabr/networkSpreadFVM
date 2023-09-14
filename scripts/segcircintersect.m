function [z,pint] = segcircintersect(p1,p2,cent,R)
% intersection of line segment between p0 and p1, with circle 
% with center cent and radius R
% p0 and p1 should be row vectors for the segment endpoint
% returns z = fraction along segment where intersection occurs
% and pint = intersection coords

%%

dx = p2-p1;
p0 = p1 - cent;

% find point where it crossed the boundary
A = sum(dx.^2,2);
B = 2*sum(p0.*dx,2);
C = sum(p0.^2,2)-R^2;

disc = sqrt(B.^2 - 4*A.*C);
z1 = (-B - disc)./(2*A);
z2 = (-B+disc)./(2*A);
z= z1;
for cc = 1:size(p1,1)
    if (z1(cc)> 0 & z1(cc)<1)
        z(cc) = z1(cc);
    elseif (z2(cc)>0 & z2(cc)<1)
        z(cc) = z2(cc);
    else
        error('no crossing found!')
    end
end

pint = p1+dx.*z; % intersection point

end