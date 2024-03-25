function [di1vol,di2vol,di3vol,vol1r] = InitialRegisterVolumes_new(vol1,vol2)

% Should perform 2d affine reg per axial slice per breast, or 3d affine reg per breast

rvec1 = colvec(mean(mean(vol1,3),2)); rvec2 = colvec(mean(mean(vol2,3),2));
cvec1 = colvec(mean(mean(vol1,3),1)); cvec2 = colvec(mean(mean(vol2,3),1));
svec1 = colvec(mean(mean(vol1,2),1)); svec2 = colvec(mean(mean(vol2,2),1));

[tr rvec1r costvecr dxvecr] = register1d(rvec1,rvec2);
[tc cvec1r costvecc dxvecc] = register1d(cvec1,cvec2);
[ts svec1r costvecs dxvecs] = register1d(svec1,svec2);

if 0
  figure(666); 
  subplot(2,2,1); plot([rvec1 rvec2 rvec1r],'LineWidth',2);
  subplot(2,2,2); plot([cvec1 cvec2 cvec1r],'LineWidth',2);
  subplot(2,2,3); plot([svec1 svec2 svec1r],'LineWidth',2);
  figure(667); 
  subplot(2,2,1); plot(dxvecr,costvecr,'LineWidth',2);
  subplot(2,2,2); plot(dxvecc,costvecc,'LineWidth',2);
  subplot(2,2,3); plot(dxvecs,costvecs,'LineWidth',2);
  drawnow;
end

%tc = 0;

di1vol = tr*ones(size(vol1),'single');
di2vol = tc*ones(size(vol1),'single');
di3vol = ts*ones(size(vol1),'single');

vol1r = movepixels(vol1,di1vol,di2vol,di3vol,1);

%keyboard

return

% ToDos
%   Modify register1d to take weighting vector
%   Find midpoint for cvec2,  box for each breast

% Try piecewise linear: edges and middles of breasts, chest wall and "nipples", top and bottom of breasts
smf = 10;
tmp1 = real(smooth1d_amd(cvec1,smf)); tmp2 = real(smooth1d_amd(cvec2,smf));
xvecc = interp1(cumsum(tmp2)/sum(tmp2),[1:length(tmp2)],cumsum(tmp1)/sum(tmp1),'linear','extrap');
dxvecc = xvecc-[1:length(xvecc)]';
figure; plot([tmp2 interp1(tmp2,xvecc) tmp1],'LineWidth',2) 

di2vol = repmat(-dxvecc,[size(vol1,1) 1 size(vol1,3)]);

showVol(vol1,vol1r,vol2)

rvec1r = colvec(mean(mean(vol1r,3),2));
cvec1r = colvec(mean(mean(vol1r,3),1)); 
svec1r = colvec(mean(mean(vol1r,2),1));

if 1
  figure(668); 
  subplot(2,2,1); plot([rvec1 rvec2 rvec1r],'LineWidth',2);
  subplot(2,2,2); plot([cvec1 cvec2 cvec1r],'LineWidth',2);
  subplot(2,2,3); plot([svec1 svec2 svec1r],'LineWidth',2);
  drawnow;
end

% Old crap, using Matlab built-in imregtform

[optimizer, metric] = imregconfig('monomodal');
%transformType = 'affine';
transformType = 'rigid';

fixed = squeeze(max(vol2,[],3)); moving = squeeze(max(vol1,[],3));
tform = imregtform(moving,fixed,transformType,optimizer,metric);
movingRegistered = imwarp(moving,tform,'OutputView',imref2d(size(fixed)));
%figure; imshowpair(fixed, moving,'Scaling','joint')
%figure; imshowpair(fixed, movingRegistered,'Scaling','joint')
tr1 = tform.T(3,2); tc = tform.T(3,1);

fixed = squeeze(max(vol2,[],2)); moving = squeeze(max(vol1,[],2));
tform = imregtform(moving,fixed,transformType,optimizer,metric);
movingRegistered = imwarp(moving,tform,'OutputView',imref2d(size(fixed)));
%figure; imshowpair(fixed, moving,'Scaling','joint')
%figure; imshowpair(fixed, movingRegistered,'Scaling','joint')
tr2 = tform.T(3,2); ts = tform.T(3,1);

tr = (tr1+tr2)/2;

