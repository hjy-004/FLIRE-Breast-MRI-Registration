function FLIRE_for_T1_Breast_MRI(dirlist,dirname_proc,outdir,method,veclist,inputType,preregType,forceflag)
% formerly RSI_breast_UCSD_long_new2.m

%% Part 1: Initalize Data ---------------------------------------------------
dirname_annot = dirname_proc;
dirname_annot2 = outdir;
version = '-v7.3';

%% Part 2: Analysis -----------------------------------------
% Read in data and save to cache dir -- should change to keep geometry info for each volume
%for each patient
for veci = 1:length(veclist)
    diri = veclist(veci);
    dirname = sprintf('%s/%s',dirname_proc,dirlist{diri});
    subdirlist = dir(sprintf('%s/20*',dirname)); subdirlist = {subdirlist.name};
    fprintf(1,'diri=%d/%d %s subdirs=%d (now=%s)\n',diri,length(dirlist),dirlist{diri},length(subdirlist),datestr(now));
    
    %check if longitudinal study
    if length(subdirlist)<2, continue; end
    
    %initialize variables
    datastructs = cell(1,length(subdirlist)); 
    failedflag = false;
    fname_out = sprintf('%s/%s/%s/datastructs_%s.mat',outdir,dirlist{diri},method,method);
    
    if ~(exist(fullfile(outdir,dirlist{diri}))); mkdir(fullfile(outdir,dirlist{diri})); end
    if ~(exist(fullfile(outdir,dirlist{diri},method))); mkdir(fullfile(outdir,dirlist{diri},method)); end
    
    if ~exist(fname_out,'file') || forceflag
        try
        try
            %% Part 2.1 Load Data
            t0_total_vec = [];
            %for each date
            for sdiri = 1:length(subdirlist)
                subdirname = sprintf('%s/%s',dirname,subdirlist{sdiri});
                fprintf(1,'Reading sdiri=%d %s\n',sdiri,subdirname);
                datastructs{sdiri} = struct();
                try
                    %load as imgs
                    t0_start=clock;
                    tmp = load(sprintf('%s/%s/%s/T1_noFS_PRE.mat',dirname_annot,dirlist{diri},subdirlist{sdiri}));
                    vol_T1_noFS_PRE = getfielddefault(tmp,'vol');   % (HY - 2024-3-17)
%                     vol_T1_noFS_PRE = getfield_regexp(tmp,'^vol.*');
                    vol_T1_noFS_PRE.imgs = single(vol_T1_noFS_PRE.imgs);
                    t0_end=clock;
                    t0_total=etime(t0_end,t0_start);
                    t0_total=t0_total./60;
                    t0_total_vec(sdiri) = t0_total;
                    voxsz = sqrt(sum(tmp.vol.Mvxl2lph(1:3,1:3).^2,1));
                    
                    if sdiri > 1
                        if ~isequal(size(vol_T1_noFS_PRE),size(datastructs{1}.vol_T1_noFS_PRE))
                            vol_T1_noFS_PRE = vol_conform(vol_T1_noFS_PRE,datastructs{1}.vol_T1_noFS_PRE);
                        end
                    end
                    datastructs{sdiri}.vol_T1_noFS_PRE = vol_T1_noFS_PRE;
                    datastructs{sdiri}.voxsz = voxsz;
                    datastructs{sdiri}.mins_load = t0_total;
                    
                catch ME
                    fprintf(1,'  **Failed*** ID:%s Message:%s\n',ME.identifier,ME.message);
                    disp(ME.stack)
                end
            end
            
            %% Part 2.2 Register data
            
            prctileval = 98;
            if 1 % Register whatever data exists %  ~failedflag % If no failure occurred
                
                % load volume and normalize by percentile
                if strcmp(inputType,'T1_noFS_PRE')
                    vol1_ctx = datastructs{1}.vol_T1_noFS_PRE;
                    [vol1 M_vol1] = ctx_ctx2mgh(vol1_ctx); 
                    vol1 = vol1/prctile(vol1(:),prctileval); 
                    voxsz1 = sqrt(sum(M_vol1(1:3,1:3).^2,1)); 
                    %exclude vol1, get matrix info
                end
                
                %initialize variables
                win = [5 5 2]; vol_mask = zeros(size(vol1));
                vol_mask(win(1):(end-win(1)+1),win(2):(end-win(2)+1),win(3):(end-win(3)+1)) = 1;
                wvol = max(0.0,vol_mask); 
                
                ctr = round((size(vol1)+1)/2);
                wvol1 = zeros(size(vol1));
                wvol1(round(0.3*size(vol1,1)):end,1:(ctr(2)-1),:) = 1;
                wvol2 = zeros(size(vol1)); wvol2(round(0.3*size(vol1,1)):end,(ctr(2)+1):end,:) = 1;
                
                for sdiri = 2:length(datastructs)
                    fprintf(1,'\nRegistering diri=%d/%d sdiri=%d/%d (now=%s)\n',diri,length(dirlist),sdiri,length(datastructs),datestr(now));
                    
                    % load volume and normalize by percentile
                    if strcmp(inputType,'T1_noFS_PRE')
                        t1_start=clock;
                        vol2_ctx = datastructs{sdiri}.vol_T1_noFS_PRE;
                        [vol2 M_vol2] = ctx_ctx2mgh(vol_conform_ctx(vol2_ctx,vol1_ctx)); 
                        vol2 = vol2/prctile(vol2(:),prctileval); 
                        voxsz2 = sqrt(sum(M_vol2(1:3,1:3).^2,1));
                        
                    end
                    
                    %Pre-registration step
                    [di1vol0_fwd di2vol0_fwd di3vol0_fwd vol2_reg0] = InitialRegisterVolumes(vol2,vol1); % showVol(vol2,vol2_reg0,vol1)
                    vol1_ctx = ctx_mgh2ctx(vol1,M_vol1); 
                    vol2_ctx = ctx_mgh2ctx(vol2,M_vol2); 
                    vol_std1 = ctx_mgh2ctx(1./max(1e-1,wvol1),M_vol1); 
                    vol_std2 = ctx_mgh2ctx(1./max(1e-1,wvol2),M_vol1);
                    curr_dir = pwd

                    [M_af1 regStruct1 vol2_af1_ctx di1vol_af1 di2vol_af1 di3vol_af1] = register_volumes_affine(vol2_ctx,vol1_ctx,vol_std1); 
                    vol2_af1 = vol2_af1_ctx.imgs;
                    [M_af2 regStruct2 vol2_af2_ctx di1vol_af2 di2vol_af2 di3vol_af2] = register_volumes_affine(vol2_ctx,vol1_ctx,vol_std2); 
                    vol2_af2 = vol2_af2_ctx.imgs;
                    cd(curr_dir)
                    if 0
                        coordmat_tmp = RSI_find_centers(wvol1); showVol(vol2,vol2_af1,vol1,wvol1,struct('RCS',round(coordmat_tmp(1,:))))
                        coordmat_tmp = RSI_find_centers(wvol2); showVol(vol2,vol2_af2,vol1,wvol2,struct('RCS',round(coordmat_tmp(1,:))))
                    end
                    wvol0 = 1-(wvol1+wvol2);
                    
                    di1vol1_fwd = (wvol0.*di1vol0_fwd + wvol1.*di1vol_af1 + wvol2.*di1vol_af2)./(wvol0 + wvol1 + wvol2); % Should perhaps smooth wvol1 amd wvol2?
                    di2vol1_fwd = (wvol0.*di2vol0_fwd + wvol1.*di2vol_af1 + wvol2.*di2vol_af2)./(wvol0 + wvol1 + wvol2); % Should perhaps smooth wvol1 amd wvol2?
                    di3vol1_fwd = (wvol0.*di3vol0_fwd + wvol1.*di3vol_af1 + wvol2.*di3vol_af2)./(wvol0 + wvol1 + wvol2); % Should perhaps smooth wvol1 amd wvol2?
                    di1vol1_pre = di1vol1_fwd;
                    di2vol1_pre = di2vol1_fwd;
                    di3vol1_pre = di3vol1_fwd;
                    vol2_reg1 = (wvol0.*vol2_reg0 + wvol1.*vol2_af1 + wvol2.*vol2_af2)./(wvol0 + wvol1 + wvol2);
                    datastructs{sdiri}.vol_T1_noFS_PRE_reg0 = movepixels(vol2,di1vol0_fwd,di2vol0_fwd,di3vol0_fwd,0);
                    datastructs{sdiri}.vol_T1_noFS_PRE_reg1 = movepixels(vol2,di1vol1_fwd,di2vol1_fwd,di3vol1_fwd,0);
                    % showVol(datastructs{sdiri}.vol_T1_noFS_PRE,datastructs{sdiri}.vol_T1_noFS_PRE_reg0,datastructs{sdiri}.vol_T1_noFS_PRE_reg1,datastructs{1}.vol_T1_noFS_PRE,struct('RCS',round(coordmat(1,:))))
                    
                    %Nonlinear reg
                    %initialize variables
                    lambda0 = 1;
                    lambda1 = 0;
                    lambda2 = [10 1 0.1 0.03 0.01 0.003 0.001];
                    smfact1 = [32 32 32 32 16 8 4];
                    smfact2 = smfact1;
                    spacing = 4*ones(size(smfact1)); spacing(end) = 2; % Should run with much greater smoothing / spacing

                    [di1vol1_fwd,di2vol1_fwd,di3vol1_fwd,vol2_reg] = NonlinRegisterVolumes_mc(vol2,vol1,lambda0,lambda1,lambda2,[smfact1],[smfact2],[spacing],wvol,di1vol1_fwd,di2vol1_fwd,di3vol1_fwd,false);
                    
                    datastructs{sdiri}.vol_T1_noFS_PRE_reg = movepixels(getfield(vol_conform_ctx(datastructs{sdiri}.vol_T1_noFS_PRE,vol1_ctx),'imgs'),di1vol1_fwd,di2vol1_fwd,di3vol1_fwd,0);
                    
                    t1_end=clock;
                    t1_total=etime(t1_end,t1_start);
                    t1_total=t1_total./60;
                    fprintf('--- %.2f mins\n',t1_total)
                    
                    datastructs{sdiri}.di1vol1_fwd = single(di1vol1_fwd);
                    datastructs{sdiri}.di2vol1_fwd = single(di2vol1_fwd);
                    datastructs{sdiri}.di3vol1_fwd = single(di3vol1_fwd);
%                     if ~strcmp(preregType,'none')
%                         datastructs{sdiri}.di1vol1_pre = single(di1vol1_pre);
%                         datastructs{sdiri}.di2vol1_pre = single(di2vol1_pre);
%                         datastructs{sdiri}.di3vol1_pre = single(di3vol1_pre);
%                     end
                    datastructs{sdiri}.mins_load = t0_total_vec(1)+t0_total_vec(sdiri);
                    datastructs{sdiri}.mins_reg = t1_total;
                    datastructs{sdiri}.mins = t0_total_vec(1)+t0_total_vec(sdiri)+t1_total;
                    
                end     % end registering volumes to reference
            end
        catch ME
            fprintf(1,'  **Failed*** ID:%s Message:%s\n',ME.identifier,ME.message);
            disp(ME.stack)
        end
        %save(fname_out,'datastructs',version);
        %fprintf(1,'--- file %s written\n',fname_out);
        
        % save displacement_field because it is a smaller file for faster reading of imgs
        fname_out2 = sprintf('%s/%s/%s/displacement_field.mat',outdir,dirlist{diri},method);
        if ~exist(fname_out2) || forceflag
            displacement_field = {};
            for datei = 1:length(datastructs)
                displacement_field{datei}.voxsz = datastructs{datei}.voxsz;
                try
                    displacement_field{datei}.di1vol1_pre = datastructs{datei}.di1vol1_pre;
                    displacement_field{datei}.di2vol1_pre = datastructs{datei}.di2vol1_pre;
                    displacement_field{datei}.di3vol1_pre = datastructs{datei}.di3vol1_pre;
                end
                try
                    displacement_field{datei}.di1vol1_fwd = datastructs{datei}.di1vol1_fwd;
                    displacement_field{datei}.di2vol1_fwd = datastructs{datei}.di2vol1_fwd;
                    displacement_field{datei}.di3vol1_fwd = datastructs{datei}.di3vol1_fwd;
                    displacement_field{datei}.mins = datastructs{datei}.mins;
                end
                try
                    displacement_field{datei}.di1vol1_rev = datastructs{datei}.di1vol1_rev;
                    displacement_field{datei}.di2vol1_rev = datastructs{datei}.di2vol1_rev;
                    displacement_field{datei}.di3vol1_rev = datastructs{datei}.di3vol1_rev;
                end
            end
            save(fname_out2,'displacement_field','-v7.3');
            fprintf('diri=%d method=%s %s saved\n',diri,method,dirlist{diri})
        end
        
        %% Save images
        %% Part 3.2 Add registered volumes to the datastructs
        % initialize registered volumes
        vols_T1_noFS_PRE_reg = []; 

        % for each date
        for sdiri = 1:length(datastructs)
            try
                subdirname = sprintf('%s/%s',dirname,subdirlist{sdiri});
                fprintf(1,'Part 3: sdiri=%d Concatenating %s\n',sdiri,subdirname);
                % first volume is already registered (?)
                if sdiri == 1
                    vol_T1_noFS_PRE_reg = datastructs{sdiri}.vol_T1_noFS_PRE.imgs;
                else
                    vol_T1_noFS_PRE_reg = datastructs{sdiri}.vol_T1_noFS_PRE_reg;

                end
                
                % concatenate vols into a 4D matrix with the 4th dimension is time
                vols_T1_noFS_PRE_reg = cat(4,vols_T1_noFS_PRE_reg,vol_T1_noFS_PRE_reg);
            catch ME
                fprintf(1,'**Failed*** ID:%s Message:%s\n',ME.identifier,ME.message);
                disp(ME.stack)
            end
        end % end date
        
        %% Part 3.3: Save registered volumes
        if ~(exist(fullfile(outdir,dirlist{diri}))); mkdir(fullfile(outdir,dirlist{diri})); end
        if ~(exist(fullfile(outdir,dirlist{diri},method,'analysis'))); mkdir(fullfile(outdir,dirlist{diri},method,'analysis')); end
        
        % files name
        fname_T1_noFS_PRE_reg = sprintf('%s/%s/%s/T1_noFS_PRE_reg.mat',outdir,dirlist{diri},method);
        
        % add vol to ctx struct and save
        struct_T1_noFS_PRE_reg = vol1_ctx; struct_T1_noFS_PRE_reg.imgs = vols_T1_noFS_PRE_reg; struct_T1_noFS_PRE_reg.visits = subdirlist;
       
        % save registered volumes as .mat and .mgz
        fprintf(1,'Part 3: diri=%d Saving registered volumes %s\n',diri,subdirname);
        t_start = clock;
        save_registered(fname_T1_noFS_PRE_reg,struct_T1_noFS_PRE_reg,version);
        t_end=clock;
        t_total=etime(t_end,t_start);
        t_total=t_total./60;
        for sdiri = 2:length(datastructs)
            t_min_save = t_total/(length(datastructs)-1);
            datastructs{sdiri}.mins_save = t_min_save;
            datastructs{sdiri}.mins = datastructs{sdiri}.mins_load + datastructs{sdiri}.mins_reg + datastructs{sdiri}.mins_save;
            fprintf(1,'--- %d mins\n',datastructs{sdiri}.mins);
        end
        save(fname_out,'datastructs',version);
        fprintf(1,'--- file %s written\n',fname_out);
        
        catch ME
            fprintf(1,'**Failed*** ID:%s Message:%s\n',ME.identifier,ME.message);
            disp(ME.stack)
        end
    end
end
    