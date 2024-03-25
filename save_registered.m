
function save_registered(fname, vol_reg, version)
    
    % save as .mat
    save(fname,'vol_reg',version);
    
    % save as .mgz
    %vol_reg.imgs = squeeze(permute(vol_reg.imgs,[2 1 3 4]));
    %QD_ctx_save_mgh(vol_reg,strrep(fname,'.mat','.mgz'));
end
