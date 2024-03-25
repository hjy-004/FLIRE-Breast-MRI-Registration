function [gradvols1,gradvols2,gradvols3] = volGradient_mc(vols,indmat1,indmat2,indmat3)

dx = 0.1; % Is this the optimal value -- seems arbitrary

if ~iscell(vols), vols = {vols}; end

nvols = length(vols);
gradvols1 = cell(1,nvols); gradvols2 = cell(1,nvols); gradvols3 = cell(1,nvols);

for vi = 1:nvols
  vol = vols{vi};
  gradvol1 = (volsamp_trilin(vol,indmat1+dx,indmat2,indmat3)-volsamp_trilin(vol,indmat1-dx,indmat2,indmat3))/(2*dx);
  gradvol2 = (volsamp_trilin(vol,indmat1,indmat2+dx,indmat3)-volsamp_trilin(vol,indmat1,indmat2-dx,indmat3))/(2*dx);
  gradvol3 = (volsamp_trilin(vol,indmat1,indmat2,indmat3+dx)-volsamp_trilin(vol,indmat1,indmat2,indmat3-dx))/(2*dx);
  gradvols1{vi} = gradvol1;
  gradvols2{vi} = gradvol2;
  gradvols3{vi} = gradvol3;
end

% ToDo
%   Use image gradient calculations from other code (e.g., PROMO), instead of volsamp_trilin?
