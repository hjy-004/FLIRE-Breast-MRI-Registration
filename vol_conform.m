function vol_out = vol_conform(vol_in,vol_target)

dims_target = size(vol_target);
dims_in = size(vol_in);
if isequal(dims_in,dims_target)
  vol_out = vol_in; % Nothing to do
else
  vol_out = zeros(dims_target,class(vol_in));
  if dims_target(3)<dims_in(3)
    svec_in = round([1:dims_target(3)]-(1+dims_target(3))/2+(1+dims_in(3))/2); ivec = find(svec_in>=1&svec_in<=dims_in(3));
    vol_out(:,:,ivec) = vol_in(:,:,svec_in(ivec));
  else
    svec_target = round([1:dims_in(3)]+(1+dims_in(3))/2-(1+dims_target(3))/2); ivec = find(svec_target>=1&svec_target<=dims_target(3));
    vol_out(:,:,svec_target(ivec)) = vol_in(:,:,ivec);
  end
end

