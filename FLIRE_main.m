%% Part 0: Set Up Matlab ---------------------------------------------------
% Set the path to where all the Matlab dependencies are located
addpath(genpath('/directory_to_where_the_Matlab_scripts_are_located'));

%% Part 1: Prepare Data ---------------------------------------------------
% Part 1.1: Initialize Data
% set paths to where the source images and the output are located
dirname_proc = '/directory_to_where_the_source_images_are_located';
outdir = '/directory_to_where_the_output_should_be_directed';

patientList = {...
%     'Cti00229D2Abc5802Ba4C6',...
%     'Cti056E752A6A1356C2Db5', ... 
    'Cti24A313973399Bcce8A9', ...
%     'Cti32Ae96288Efe442E2Df',...
    };

dirlist = cat(2,patientList);

%TO EDIT
veclist=[1:length(dirlist)]; %idx of patients in dirlist for analysis
dirlist(veclist)
%FLIRE1508 for figures and analysis, FLIRE1509 for runtime
methodlist ={'FLIRE1508'}; 
inputType = 'T1_noFS_PRE'; %Options: 'T1_noFS_PRE','T1_FS_PRE','T1_noFS_PRE,T1_FS_PRE'
preregType = 'affine';     %Options: 'none','translation','rigid','affine'
forceflag = true;          %Options: true, false, 0, 1

%% Part 2: Registration
method = 'FLIRE_for_github';

if contains(method,'FLIRE')
    FLIRE_for_T1_Breast_MRI(dirlist,dirname_proc,outdir,method,veclist,inputType,preregType,forceflag)
end