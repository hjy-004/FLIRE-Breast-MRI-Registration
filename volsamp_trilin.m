function  volr = volsamp_trilin(vol,indmat1,indmat2,indmat3)

volr = interp3(vol,indmat2,indmat1,indmat3,'linear',0);
%volr = ba_interp3(vol,indmat2,indmat1,indmat3,'linear',0); % What the hell -- this is actually like interp2, not interp3!!

return


% Old version

dims = size(vol);
volr = zeros(size(indmat1));;
indmat = double([indmat1(:) indmat2(:) indmat3(:) ones(length(indmat1(:)),1)]);
volr(:) = single(volgetvxlsvalMEX(indmat, dims, double(vol), eye(4),  1, 1));

