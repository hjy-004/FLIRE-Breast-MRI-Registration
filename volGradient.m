function [gradvol1,gradvol2,gradvol3] = volGradient(vol,indmat1,indmat2,indmat3)

dx = 0.1;
gradvol1 = (volsamp_trilin(vol,indmat1+dx,indmat2,indmat3)-volsamp_trilin(vol,indmat1-dx,indmat2,indmat3))/(2*dx);
gradvol2 = (volsamp_trilin(vol,indmat1,indmat2+dx,indmat3)-volsamp_trilin(vol,indmat1,indmat2-dx,indmat3))/(2*dx);
gradvol3 = (volsamp_trilin(vol,indmat1,indmat2,indmat3+dx)-volsamp_trilin(vol,indmat1,indmat2,indmat3-dx))/(2*dx);
return
dims = size(vol);
tmpvol = diff(vol,1,1);
gradvol1 = (cat(1,tmpvol,zeros(1,dims(2),dims(3))) + cat(1,zeros(1,dims(2),dims(3)),tmpvol))/2;
tmpvol = diff(vol,1,2);
gradvol2 = (cat(2,tmpvol,zeros(dims(1),1,dims(3))) + cat(2,zeros(dims(1),1,dims(3)),tmpvol))/2;
tmpvol = diff(vol,1,3);
gradvol3 = (cat(3,tmpvol,zeros(dims(1),dims(2),1)) + cat(3,zeros(dims(1),dims(2),1),tmpvol))/2;

