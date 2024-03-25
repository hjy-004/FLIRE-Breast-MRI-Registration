function vol_out = conform_vol_ctx(vol_in,vol_ref)

% Force same dimensions and slice direction as vol_ref

vol_out = vol_ref;
if dot(vol_in.Mvxl2lph(:,3),vol_ref.Mvxl2lph(:,3))<0
  vol_out.imgs = vol_conform(flipdim(vol_in.imgs,3),vol_ref.imgs);
else
  vol_out.imgs = vol_conform(vol_in.imgs,vol_ref.imgs);
end

