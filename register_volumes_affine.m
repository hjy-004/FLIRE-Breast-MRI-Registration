function [M_af regStruct vol_af di1vol_af di2vol_af di3vol_af] = register_volumes_affine(vol, volm, volstd)

if ~exist('volstd','var')
   volstd = [];
end

volm_sub = subsample_volume(subsample_volume(volm));
volstd_sub = subsample_volume(subsample_volume(volstd));
vol_sub = subsample_volume(subsample_volume(vol));

if ~isempty(volstd)
  wvol = max(1,volstd_sub.imgs).^-2;
else
  wvol = 1;
end

costfun = @(x)nanmean(colvec((getfield(vol_resample_amd(vol_sub,volm_sub,pars2M_af(x),1),'imgs')-volm_sub.imgs).^2.*wvol)); % Should perhaps use different cost function for binary volumes?

initp = [0 0 0 1 1 1 0 0 0 0 0 0];

vol_sub_res0 = vol_resample_amd(vol_sub,volm_sub,pars2M_af(initp),1);

parvec_opt = initp;

options.randscale = [3 3 3 0.0 0.0 0.0 5*pi/180 5*pi/180 5*pi/180 0 0 0]; % Start with coarse rb-reg
[parvec_opt cost_opt] = fminsearch_stochastic(costfun,parvec_opt,options);

if 0
  vol_sub_res1 = vol_resample_amd(vol_sub,volm_sub,pars2M_af(parvec_opt),1);
  showVol(volm_sub.imgs,vol_sub_res0.imgs,vol_sub_res1.imgs,log10(max(eps,(vol_sub_res0.imgs-volm_sub.imgs).^2.*wvol)),log10(max(eps,(vol_sub_res1.imgs-volm_sub.imgs).^2.*wvol)))
end

[parvec_opt cost_opt] = fminsearch(costfun,parvec_opt,statset('Display','final','MaxIter',2000,'MaxFunEvals',2000));

Mreg = eye(4);
[M_af Mss Ms Mrl Mrp Mrh Mt] = pars2M_af(parvec_opt);
M_rb =  Mrl*Mrp*Mrh*Mt*Mreg;
vol_af = vol_resample_amd(vol,volm,M_af,1);

regStruct.min_cost_af = cost_opt;
regStruct.parvec_af = parvec_opt;
regStruct.M_atl_to_vol_rb = M_rb;
regStruct.M_atl_to_vol_af = M_af;

%showVol(vol_af,volm)

if nargout>=6 % Note that this only works if vol and volm have same dimensions
%  ctr = (size(volm.imgs)+1)/2; [I2 I1 I3] = meshgrid([1:size(volm.imgs,2)]-ctr(2),[1:size(volm.imgs,1)]-ctr(1),[1:size(volm.imgs,3)]-ctr(3)); % To work with the Matlab affine reg, which is around center of volume
  M_reg = M_af;
  [I2 I1 I3] = meshgrid([1:size(volm.imgs,2)],[1:size(volm.imgs,1)],[1:size(volm.imgs,3)]);
  coords = cat(2,I1(:),I2(:),I3(:),ones(size(I1(:))))';
  Mvxl2vxl = inv(volm.Mvxl2lph)*M_reg*vol.Mvxl2lph; % Assume M_reg is in LPH coordinates 
  coords_trans = Mvxl2vxl*coords;
  di1vol_af = reshape((coords_trans(1,:)-coords(1,:))',size(volm.imgs));
  di2vol_af = reshape((coords_trans(2,:)-coords(2,:))',size(volm.imgs));
  di3vol_af = reshape((coords_trans(3,:)-coords(3,:))',size(volm.imgs));
  if 0
    vol_reg = vol_resample_amd(vol,volm,M_reg,1); % Same as vol_af
    vol_reg2 = movepixels(vol.imgs,di1vol_af,di2vol_af,di3vol_af,0);
    showVol(volm.imgs,vol_reg.imgs,vol_reg2)
  end
end

