function Conn_ROI2ROI_Rest

% the cdirectory with this script becomes the current directory
WD = pwd;

degreeOfSmoothing = 6;
ScriptsDir = '/Users/mohamed/Documents/GitHub/CatRestScripts/CatRest_scripts';

% batch filename to save the conn mat file
BATCH.filename = '~/Data/CatRest_BIDS/CONN_output/Conn_ROI2ROI_Rest.mat';

% we add all the subfunctions that are in the sub directories
addpath(genpath(ScriptsDir))
opt = getOption();

[group, opt, BIDS] = getData(opt);

% Get TR from metadata
TR = opt.metadata.RepetitionTime;  %TR = 2.25  ;

% Get the number of sessions from the first subject (assuming all subjects have the same sessions)
subNumber = group(1).subNumber{1} ;   % Get the Subject ID
[~, numSessions] = getSessions(BIDS, subNumber, opt);  %numSessions = 1;

% creates prefix to look for
[prefix, ~] = getPrefix('FFX', opt, degreeOfSmoothing);

%% ROIs : separate Nifti files
ROI_folder = fullfile(pwd,'ROI') ;
ROI_names = {'V1_diff_L',...
    'V1_diff_R',...
    'V5_conj_L',...
    'V5_conj_R',...
    'LGN_L',...
    'LGN_R'
    } ;
%     'Pulvinar_L',...
%     'Pulvinar_R' } ;

ROI_fileName = {
    'V1_diff_L_ROI_x=-14_y=-98_z=2_515voxels_Sphere10',...
    'V1_diff_R_ROI_x=16_y=-92_z=2_515voxels_Sphere10',...
    'V5_conj_L_ROI_x=-38_y=-76_z=4_515voxels_Sphere10',...
    'V5_conj_R_ROI_x=44_y=-66_z=2_515voxels_Sphere10',...
    'L_LGN_mask_MNI_thr25',...
    'R_LGN_mask_MNI_thr25'};
%     'rL_pulvinar',...
%     'rR_pulvinar'};

% For testing - test 1 ROI
% ROI_fileName = ROI_fileName(1);
% ROI_names=ROI_names(1);

% Check what to run
doSetup = 1;
doDenoising = 1;
doFirstLevel = 1;

% Get Data
% [data, greyMatter, whiteMatter, csf, structural] = getAnatomicalStruct ;

% data = data(1:4);
% structural = structural(1:4);
% greyMatter = greyMatter(1:4);
% whiteMatter = whiteMatter(1:4);
% csf = csf(1:4);
numSubj = 0;
for i=1:length(group) 
    numSubj = numSubj + group(i).numSub;
end

SubCounter = 0;
for iGroup= 1:length(group)              % For each group
    groupName = group(iGroup).name ;     % Get the Group name
    
    for iSub = 1:group(iGroup).numSub    % For each Subject in the group
        
        SubCounter = SubCounter + 1;
        % clear previous matlabbatch
        matlabbatch = [];
        
        subNumber = group(iGroup).subNumber{iSub} ;   % Get the Subject ID
        
        % identify sessions for this subject
        [sessions, numSessions] = getSessions(BIDS, subNumber, opt);
        
        for iSes = 1:numSessions % For each session
            
            % get all runs for that subject across all sessions
            [runs, numRuns] = getRuns(BIDS, subNumber, sessions{iSes}, opt);
            
            for iRun = 1:numRuns % For each run
                
                fprintf(1,'PROCESSING GROUP: %s SUBJECT No.: %i SUBJECT ID : %s SESSION: %i RUN:  %i \n',...
                    groupName,iSub,subNumber,iSes,iRun)
                

                [fileName, subFuncDataDir]= getBoldFilename(...
                    BIDS, ...
                    subNumber, sessions{iSes}, runs{iRun}, opt);
                
                files = inputFileValidation(subFuncDataDir, prefix, fileName);
                assert(length(files)==1)
                [~,filename,extention]=fileparts(files{1});
                
                functionalDir{SubCounter} = subFuncDataDir;
                functionalFilename{SubCounter} = filename;
                structSession = 1; % assuming the strucural is obtained in session 1
                struct = spm_BIDS(BIDS, 'data', ...
                    'sub', subNumber, ...
                    'ses', sessions{structSession}, ...
                    'type', 'T1w');
                % we assume that the first T1w is the correct one (could be an
                % issue for data set with more than one
                struct = struct{1};
                
                [filepath,~,~] = fileparts(struct);
                
                anat = dir(fullfile(filepath,'wm*'));
                maskGrey=dir(fullfile(filepath,'wc1*'));
                maskWhite=dir(fullfile(filepath,'wc2*'));
                maskCSF=dir(fullfile(filepath,'wc3*'));
                
                assert(length(anat)==1);
                assert(length(maskGrey)==1);
                assert(length(maskWhite)==1);
                assert(length(maskCSF)==1);
                
                BATCH.Setup.masks.Grey{SubCounter} = fullfile(maskGrey.folder,maskGrey.name);
                BATCH.Setup.masks.White{SubCounter} =  fullfile(maskWhite.folder,maskWhite.name) ;
                BATCH.Setup.masks.CSF{SubCounter}   = fullfile(maskCSF.folder,maskCSF.name);
                BATCH.Setup.structurals{SubCounter} = fullfile(anat.folder,anat.name);
                
                ffxDir = getFFXdir(subNumber, degreeOfSmoothing, opt);
                BATCH.Setup.spmfiles{SubCounter}{1} = fullfile(ffxDir ,'SPM.mat') ;
                
                
            end
            
        end
    end
end

                        
%BATCH.Setup.Preprocessing.steps = {'art_thresholds'};

                        %%%%%%%%%%%%%%%%%%%%%%%

% parallel processing
BATCH.parallel.N = 0;


BATCH.Setup.isnew = 0;
BATCH.Setup.done = doSetup;
BATCH.Setup.overwrite = 0;

BATCH.Setup.nsubjects = numSubj ;
BATCH.Setup.RT = TR ;
BATCH.Setup.acquisitiontype = 1; % 1-Continous  2-Sparse
BATCH.Setup.voxelmask = 1 ;
BATCH.Setup.voxelmaskfile = fullfile(pwd,'ROI','wholeBrainMask.nii') ;
BATCH.Setup.voxelresolution = 1 ;
BATCH.Setup.analysisunits = 1;
BATCH.Setup.analyses = 1; %[1 2] ; % 1-ROI to ROI      2-Seed-to-voxel
                               % 3-Voxel-to-voxel  4-Dynamic FC

                               
% Optional output files (outputfiles(1): 1/0 confound beta-maps; outputfiles(2): 1/0 confound-corrected timeseries
% Outputfiles(3): 1/0 seed-to-voxel r-maps)  Outputfiles(4): 1/0 seed-to-voxel p-maps) 
% Outputfiles(5): 1/0 seed-to-voxel FDR-p-maps); outputfiles(6): 1/0 ROI-extraction REX files; [0,0,0,0,0,0] 
%BATCH.Setup.outputfiles = [1 1 1 1 1 0] ;
BATCH.Setup.outputfiles = [0 0 0 0 0 0] ;

% Source of functional data for ROI timeseries extraction (typically unsmoothed BOLD signal volumes); 1: same as 'functionals' field; 
% 2: same as 'functionals' field after removing leading 's' from filename; 
% 3: other (define rule programmatically; see help conn_rulebasedfilename);
% 4: other (different set of functional volume files) [2] 
BATCH.Setup.roiextract = 2 ;  

% for iSub = 1: BATCH.Setup.nsubjects
%     
%     BATCH.Setup.spmfiles{iSub}{1} = fullfile(ffxDir{iSub} ,'SPM.mat') ;
%     
%     BATCH.Setup.masks.Grey{iSub}  = greyMatter{iSub} ;
%     BATCH.Setup.masks.White{iSub} =  whiteMatter{iSub} ;
%     BATCH.Setup.masks.CSF{iSub}   = csf{iSub} ;
%     BATCH.Setup.structurals{iSub} = structural{iSub} ;
% 
% end

BATCH.Setup.add = 0 ;

%% ROIs
% load ROIs
BATCH.Setup.rois.names = ROI_names ;
% 
for SubCounter = 1:BATCH.Setup.nsubjects
    for iSess = 1:numSessions
        for iROI = 1:length(ROI_names)
            BATCH.Setup.rois.files{iROI}{SubCounter}{iSess} = fullfile(ROI_folder,[ROI_fileName{iROI},'.nii']);            
            BATCH.Setup.rois.dimensions{iROI} = 1 ;
        end
    end
end

dimensions = ones(1,length(BATCH.Setup.rois.names));

regresscovariates = zeros(1,length(BATCH.Setup.rois.names)) ;
regresscovariates (2:3) = 1 ;  % only regress covariates of the white matter and csf

BATCH.Setup.rois.mask = zeros(1,length(BATCH.Setup.rois.names)) ;
%BATCH.Setup.rois.dimensions = {dimensions} ;
BATCH.Setup.rois.regresscovariates = regresscovariates ;
BATCH.Setup.rois.roiextract = 1 ;


%% Covariates  (Add the ART Outliers)
BATCH.Setup.covariates.names{1}= 'ART_Outliers';  
for iSub = 1:BATCH.Setup.nsubjects
    for iSess = 1:numSessions
        artFile = dir(fullfile(functionalDir{SubCounter},['art_regression_outliers_w*']));
        assert(length(artFile)==1);
        artFile = artFile.name;
        %BATCH.Setup.covariates.files{1}{iSub}{iSess} = fullfile(functionalDir{SubCounter},['art_regression_outliers_',functionalFilename{SubCounter},'.mat']) ;
        BATCH.Setup.covariates.files{1}{iSub}{iSess} = fullfile(functionalDir{SubCounter},artFile) ;
    end
end
BATCH.Setup.covariates.add = 1 ;   % Add this covariate to any previous covariates without deleting the old ones.


%% CONDITIONS
% conditions
%        names         : conditions.names{ncondition} char array of condition name
%        onsets        : conditions.onsets{ncondition}{nsub}{nses} vector of condition onsets (in seconds)
%        durations     : conditions.durations{ncondition}{nsub}{nses} vector of condition durations (in seconds)

BATCH.Setup.conditions.names= [];
BATCH.Setup.conditions.names{1}='Rest';

for SubCounter = 1:BATCH.Setup.nsubjects
        % Onsets
        BATCH.Setup.conditions.onsets{1}{SubCounter}{1}= [0];
        % durations
        BATCH.Setup.conditions.durations{1}{SubCounter}{1}= [Inf];

end

%% 2nd level groups 
groupID=[];
for iGroup=1:length(group)
    BATCH.Setup.subjects.group_names{iGroup} = group(iGroup).name ; % subjects.group_names{ngroup} char array of second-level group name
    groupID = [groupID, ones(1,group(iGroup).numSub)*iGroup];
end

% 
% numberControls = 15 ; 
% numberPatients = 11 ;
% numberControls = (ones(1,numberControls)*1);
% numberPatients = (ones(1,numberPatients)*2);


BATCH.Setup.subjects.groups = [groupID] ;


% Age as 2nd level covariate
BATCH.Setup.subjects.effect_names{1}= 'Age';
BATCH.Setup.subjects.effects{1}= [ 
    31    25    28    36    24    23    24    26    17    26    23    29    24    24    22 ... % CONTROLS
    22    24    23    20    34    29    23    35    34    19    32  ];                         % PATIENTS


%% DENOISING

BATCH.Denoising.done = doDenoising ;
BATCH.Denoising.overwrite = 0;
BATCH.Denoising.filter = [0.008 0.09] ;
BATCH.Denoising.detrending = 1 ;  %0/1/2/3: BOLD times-series polynomial detrending order (0: no detrending; 1: linear detrending; ... 3: cubic detrending) 
BATCH.Denoising.despiking = 0 ;
BATCH.Denoising.regbp = 1 ;

BATCH.Denoising.confounds.names = {'White Matter','CSF','SPM covariates','ART_Outliers','Effect of rest'} ;
BATCH.Denoising.confounds.dimensions = {5,5,Inf,Inf,Inf} ;
BATCH.Denoising.confounds.deriv = {0,0,1,0,1} ;


%% FIRST LEVEL  ROI-to-ROI & ROI-to-Voxel

 % BATCH.Analysis PERFORMS FIRST-LEVEL ANALYSES (ROI-to-ROI and seed-to-voxel) %!
 % Analysis            
%  
    BATCH.Analysis.done = doFirstLevel ;
    BATCH.Analysis.overwrite = 0 ;
    BATCH.Analysis.analysis_number = 1 ;
    BATCH.Analysis.measure = 1 ;        % Connectivity measure used, 1 = 'correlation (bivariate)', 2 = 'correlation (semipartial)', 3 = 'regression (bivariate)', 4 = 'regression (multivariate)'; [1] 
    BATCH.Analysis.weight  = 2 ;       % Within-condition weight, 1 = 'none', 2 = 'hrf', 3 = 'hanning'; [2] 
    BATCH.Analysis.modulation = 0 ;    % temporal modulation, 0 = standard weighted GLM analyses; 1 = gPPI analyses of condition-specific temporal modulation factor, or a string for PPI analyses of other temporal modulation factor (same for all conditions; valid strings are ROI names and 1st-level covariate names)'; [0] 
    %conditions      : (for modulation==1 only) list of condition names to be simultaneously entered in gPPI model (leave empty for default 'all existing conditions') [] 
    BATCH.Analysis.type = 1 ;          % analysis type, 1 = 'ROI-to-ROI', 2 = 'Seed-to-Voxel', 3 = 'all'; [3] 
    

%% BATCH.Results PERFORMS SECOND-LEVEL ANALYSES (ROI-to-ROI and Seed-to-Voxel analyses) %!
  %Results             
 
%     BATCH.Results.done            = 1 ;
%     BATCH.Results.overwrite       = 1 ;
%     BATCH.Results.analysis_number = 1 ;
%     BATCH.Results.foldername      = '' ;
%  
% %    BATCH.Results.between_subjects
%        BATCH.Results.between_subjects.effect_names  = {'AllSubjs'} ;
%        BATCH.Results.between_subjects.contrast  = [1] ;   % contrast vector (same size as effect_names)
% %  
%     between_conditions [defaults to multiple analyses, one per condition]
%       effect_names  : cell array of condition names (as in Setup.conditions.names)
%       contrast      : contrast vector (same size as effect_names)
%  
%     between_sources    [defaults to multiple analyses, one per source]
%       effect_names  : cell array of source names (as in Analysis.regressors, typically appended with _1_1; generally they are appended with _N_M -where N is an index ranging from 1 to 1+derivative order, and M is an index ranging from 1 to the number of dimensions specified for each ROI; for example ROINAME_2_3 corresponds to the first derivative of the third PCA component extracted from the roi ROINAME) 
%       contrast      : contrast vector (same size as effect_names)

 
%%

save('CONN_BATCH.mat')
conn_batch(BATCH)