function y = Hfun_mc(x,lambda0,lambda1,lambda2,gs1,gs2,gs3,ws)

nchan = length(gs1);
dims = size(gs1{1});
x = reshape(x,[dims 3]);
x1 = squeeze(x(:,:,:,1));
x2 = squeeze(x(:,:,:,2));
x3 = squeeze(x(:,:,:,3));
y1 = 0; y2 = 0; y3 = 0;
for i = 1:nchan
  y1 = y1 + lambda0*ws{i}.*((gs1{i}.*gs1{i}).*x1 + (gs1{i}.*gs2{i}).*x2 + (gs1{i}.*gs3{i}).*x3); % Not clear if ws / wvol should be included here, since already in gs  -- should check with non-binary wvol
  y2 = y2 + lambda0*ws{i}.*((gs2{i}.*gs1{i}).*x1 + (gs2{i}.*gs2{i}).*x2 + (gs2{i}.*gs3{i}).*x3);
  y3 = y3 + lambda0*ws{i}.*((gs3{i}.*gs1{i}).*x1 + (gs3{i}.*gs2{i}).*x2 + (gs3{i}.*gs3{i}).*x3);
end
y1 = y1 + lambda1*x1;
y2 = y2 + lambda1*x2;
y3 = y3 + lambda1*x3;
y1 = y1 + lambda2*volLaplacian(x1);
y2 = y2 + lambda2*volLaplacian(x2);
y3 = y3 + lambda2*volLaplacian(x3);
y = double([colvec(y1); colvec(y2); colvec(y3)]);

% Should check why multi-channel reg is not consistent with single-channel, even if channels identical
%   Compare gs* and H*x across channels (looks like H scales with nchan^2, while g scales with nchan)

% Should check if gradient & Hessian are correctly computed
%   Plot cost as function of x = x0 + s*dx 

