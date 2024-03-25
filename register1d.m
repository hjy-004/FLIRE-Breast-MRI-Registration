function [tx v1r costvec dxvec] = register1d(v1,v2,wvec)

if ~exist('wvec','var'),  wvec = ones(size(v2)); end

xvec = v1; xvec(1:end) = 1:length(v1);
dxvec = [-length(v1)/3:length(v1)/3];
costvec = NaN(size(dxvec)); ivec = find(wvec);
for dxi = 1:length(dxvec)
  v1_res = interp1(v1,xvec+dxvec(dxi));
  costvec(dxi) = 1-nancorr(v1_res(ivec),v2(ivec));;
end
[mv mi] = min(costvec);
tx = dxvec(mi);
v1r = interp1(v1,xvec+tx);

% ToDo
%   Use local quadratic fit to achieve super-resolution

