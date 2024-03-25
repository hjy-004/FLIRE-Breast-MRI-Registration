function [di1vol,di2vol,di3vol,vols1r] = NonlinRegisterVolumes_mc_MT(vols1_orig0,vols2_orig0,lambdas0,lambdas1,lambdas2,smfacts1,smfacts2,spacings,wvols_orig0,di1vol,di2vol,di3vol,plotflag,niter)
%% Part 1: Initalize variables --------------------------------------------
if ~exist('plotflag','var') | isempty(plotflag)
  plotflag = false;
end

if ~exist('niter','var') | isempty(niter)
  niter = 1;
end

vols1_orig={}; vols1_orig{1}=vols1_orig0;   % (HY - 2024-3-19)
vols2_orig={}; vols2_orig{1}=vols2_orig0;   % (HY - 2024-3-19)
wvols_orig={}; wvols_orig{1}=wvols_orig0;   % (HY - 2024-3-19)

dims_orig = size(vols1_orig{1});
% dims_orig = size(vols1_orig);    % HY (2024-3-17)
nvols = length(vols1_orig);
vols1r = cell(1,nvols); vols1 = cell(1,nvols); vols2 = cell(1,nvols); wvols = cell(1,nvols);

% remove nonfinite numbers from input volumes
for vi = 1:nvols
  vols1_orig{vi}(find(~isfinite(vols1_orig{vi}))) = 0;
  vols2_orig{vi}(find(~isfinite(vols2_orig{vi}))) = 0;
end

if ~exist('wvols_orig','var') | isempty(wvols_orig)
  for vi = 1:nvols
    wvols_orig{vi} = ones(dims_orig);
  end
end

if ~exist('di3vol','var') | isempty(di3vol) % Initialize displacements to zero, if none specified
  di1vol = zeros(dims_orig);
  di2vol = zeros(dims_orig);
  di3vol = zeros(dims_orig);
end
%% Part 2: Register volumes -----------------------------------------------
nsteps = max([length(smfacts1),length(smfacts2),length(spacings),length(lambdas0),length(lambdas1),length(lambdas2),length(spacings)]);

for stepi = 1:nsteps
  % Intialize params from the input. Use next val in vector or last val
  lambda0 = lambdas0(min(length(lambdas0),stepi));
  lambda1 = lambdas1(min(length(lambdas1),stepi));
  lambda2 = lambdas2(min(length(lambdas2),stepi));
  smfact1 = smfacts1(min(length(smfacts1),stepi));
  smfact2 = smfacts2(min(length(smfacts2),stepi));
  spacing = spacings(min(length(spacings),stepi));
  clear vol vol2
  % guassian smoothing
  for vi = 1:nvols
    vols1{vi} = real(smooth3d(vols1_orig{vi},max(smfact1,2*spacing),max(smfact1,2*spacing),max(smfact1,2*spacing)));
    vols2{vi} = real(smooth3d(vols2_orig{vi},max(smfact2,2*spacing),max(smfact2,2*spacing),max(smfact2,2*spacing)));
  end
  % specify voxels used in reg
  i1vec = 1:spacing:dims_orig(1);
  i2vec = 1:spacing:dims_orig(2);
  i3vec = 1:spacing:dims_orig(3);
  for vi = 1:nvols
    vols1{vi} = vols1{vi}(i1vec,i2vec,i3vec);
    vols2{vi} = vols2{vi}(i1vec,i2vec,i3vec);
    wvols{vi} = wvols_orig{vi}(i1vec,i2vec,i3vec);
  end
  % 
  dims = size(vols1{1}); 
  [indvol1,indvol2,indvol3] = ndgrid(1:dims(1),1:dims(2),1:dims(3));
  indvol1_orig = 1 + spacing*(indvol1-1); indvol2_orig = 1 + spacing*(indvol2-1); indvol3_orig = 1 + spacing*(indvol3-1);
  x1 = volsamp_trilin(di1vol/spacing,indvol1_orig,indvol2_orig,indvol3_orig); x1(find(~isfinite(x1))) = 0;
  x2 = volsamp_trilin(di2vol/spacing,indvol1_orig,indvol2_orig,indvol3_orig); x2(find(~isfinite(x2))) = 0;
  x3 = volsamp_trilin(di3vol/spacing,indvol1_orig,indvol2_orig,indvol3_orig); x3(find(~isfinite(x3))) = 0;
  xmat = cat(4,x1,x2,x3);
  for vi = 1:nvols
    vol1r = volsamp_trilin(vols1{vi},indvol1+x1,indvol2+x2,indvol3+x3); vol1r(find(~isfinite(vol1r))) = 0;
    vols1r{vi} = vol1r;
  end
  for iter = 0:niter % Not clear we need support for multiple iterations
    tic
    niters = NaN;
    [dvols1dx1,dvols1dx2,dvols1dx3] = volGradient_mc(vols1r,indvol1,indvol2,indvol3);
    if iter > 0
      g1 = lambda1*x1; g2 = lambda1*x2; g3 = lambda1*x3;
      g1 = g1 + lambda2*volLaplacian(x1);
      g2 = g2 + lambda2*volLaplacian(x2);
      g3 = g3 + lambda2*volLaplacian(x3);
      for vi = 1:nvols
        vol1r = vols1r{vi};
        vol2 = vols2{vi};
        wvol = wvols{vi};
        dvol1dx1 = dvols1dx1{vi};
        dvol1dx2 = dvols1dx2{vi};
        dvol1dx3 = dvols1dx3{vi};
        g1 = g1 + lambda0*wvol.*(vol1r-vol2).*dvol1dx1;
        g2 = g2 + lambda0*wvol.*(vol1r-vol2).*dvol1dx2;
        g3 = g3 + lambda0*wvol.*(vol1r-vol2).*dvol1dx3;
      end
      g = double([colvec(g1); colvec(g2); colvec(g3)]);
      H = @ (x) Hfun_mc(x,lambda0,lambda1,lambda2,dvols1dx1,dvols1dx2,dvols1dx3,wvols); % Should check if gradient and Hessian are correct (test with tiny perturbations)
      [dx,flag,relres,niters,resvec] = bicgstab(H,-g,1e-2,1000);
      dx = reshape(dx,[dims 3]);
      xmat = xmat+dx;
    end
    x1 = xmat(:,:,:,1);
    x2 = xmat(:,:,:,2);
    x3 = xmat(:,:,:,3);
    for vi = 1:nvols
      vol1r = volsamp_trilin(vols1{vi},indvol1+x1,indvol2+x2,indvol3+x3); 
      vol1r(find(~isfinite(vol1r))) = 0;
      vols1r{vi} = vol1r;
    end
    cost0 = 0; sum0 = 0;
    for vi = 1:nvols
      cost0 = cost0 + sum(1/2*lambda0*colvec(wvols{vi}.*(vols1r{vi}-vols2{vi}).^2));
      sum0 = sum0 + sum(1/2*lambda0*colvec(wvols{vi}.*(vols2{vi}).^2)); 
    end 
    mserr = cost0/sum0;  
    costvol1 = x1.^2+x2.^2+x3.^2;
    costvol2 = 0;
    costvol2 = costvol2 + (cat(1,diff(x1,1,1),zeros(1,dims(2),dims(3))).^2) + (cat(1,diff(x2,1,1),zeros(1,dims(2),dims(3))).^2) + (cat(1,diff(x3,1,1),zeros(1,dims(2),dims(3))).^2);
    costvol2 = costvol2 + (cat(2,diff(x1,1,2),zeros(dims(1),1,dims(3))).^2) + (cat(2,diff(x2,1,2),zeros(dims(1),1,dims(3))).^2) + (cat(2,diff(x3,1,2),zeros(dims(1),1,dims(3))).^2);
    costvol2 = costvol2 + (cat(3,diff(x1,1,3),zeros(dims(1),dims(2),1)).^2) + (cat(3,diff(x2,1,3),zeros(dims(1),dims(2),1)).^2) + (cat(3,diff(x3,1,3),zeros(dims(1),dims(2),1)).^2);
    cost1 = 1/2*lambda1*sum(costvol1(:)); 
    cost2 = 1/2*lambda2*sum(costvol2(:));   
    cost = cost0 + cost1 + cost2;
    fprintf(1,'iter=%d: cost=%f (%f %f %f), mserr=%f, niters=%.1f\n',iter,cost,cost0,cost1,cost2,mserr,niters)
    toc
    if plotflag
      for vi = 1:nvols
        showVol(ctx_mgh2ctx(vols2{1}),ctx_mgh2ctx(vols1r{1}),ctx_mgh2ctx(vols1r{1}-vols2{1}),ctx_mgh2ctx(dvols1dx1{1}),ctx_mgh2ctx(dvols1dx2{1}),ctx_mgh2ctx(dvols1dx3{1})); drawnow;
      end
    end
  end
  % Upsample displacement fields
  [indvol1,indvol2,indvol3] = ndgrid(1:dims_orig(1),1:dims_orig(2),1:dims_orig(3));
  indvol1 = (1-1/spacing)/2 + 1/spacing*indvol1;
  indvol2 = (1-1/spacing)/2 + 1/spacing*indvol2;
  indvol3 = (1-1/spacing)/2 + 1/spacing*indvol3;
  di1vol = volsamp_trilin(x1*spacing,indvol1,indvol2,indvol3); di1mat(find(~isfinite(di1vol))) = 0;
  di2vol = volsamp_trilin(x2*spacing,indvol1,indvol2,indvol3); di2mat(find(~isfinite(di2vol))) = 0;
  di3vol = volsamp_trilin(x3*spacing,indvol1,indvol2,indvol3); di3mat(find(~isfinite(di3vol))) = 0;
end

clear xmat costvol0 costvol1 costvol2 diffvol g1 g2 g3 dx vol1 vol2 wvol vol1r vols1r vols2 wvols
[indvol1,indvol2,indvol3] = ndgrid(1:dims_orig(1),1:dims_orig(2),1:dims_orig(3));
for vi = 1:nvols
  vol1r = volsamp_trilin(vols1_orig{vi},indvol1+di1vol,indvol2+di2vol,indvol3+di3vol); 
  vol1r(find(~isfinite(vol1r))) = 0;
  vols1r{vi} = vol1r;
end

% ToDo
%   Test w. non-binary wvol
%   Test w. different wvols per channel