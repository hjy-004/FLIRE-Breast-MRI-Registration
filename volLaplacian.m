function lapvol = volLaplacian(vol)

dims = size(vol);
lapvol = 0;
D = sparse(diff(eye(dims(1))));
DtD = D'*D;
lapvol = lapvol + reshape(DtD*reshape(double(vol),[dims(1) dims(2)*dims(3)]),dims);
D = sparse(diff(eye(dims(2))));
DtD = D'*D;
lapvol = lapvol + permute(reshape(DtD*reshape(permute(double(vol),[2 1 3]),[dims(2) dims(1)*dims(3)]),[dims(2) dims(1) dims(3)]),[2 1 3]);
D = sparse(diff(eye(dims(3))));
DtD = D'*D;
lapvol = lapvol + permute(reshape(DtD*reshape(permute(double(vol),[3 2 1]),[dims(3) dims(2)*dims(1)]),[dims(3) dims(2) dims(1)]),[3 2 1]);
% Shoudl have opposite sign??
