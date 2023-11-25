#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
%% This script is for automated downloading of data from the '100 Synaesthetic Brains' HCP database
% Racey, C., Kampoureli, C., Bowen-Hill, O., Simpson, I., Rae, C., del Rio,
% m., Simner, J., & Ward, J. An Open Science MRI Database of over 100
% Synaesthetic Brains and Accompanying Deep Phenotypic Information. Sci
% Data 10, 766 (2023). https://doi.org/10.1038/s41597-023-02664-4
%
% https://www.nature.com/articles/s41597-023-02664-4 
%
% NITRC data repository page: https://syn_hcp.projects.nitrc.org
%
% This data release is in HCP format. The dataset can be downloaded in full
% or subsets of participants and datafiles can be selected. The current example
% configuration downloads just the first 3 subjects, and selects only
% packages of motion and tnsr figures. These files are small for easy test
% downloads. 
%
% If downloaded in full, each subject's data is ~37GB. There are 127
% subjects, the full size of the dataset is ~4.7TB. For most practical
% analysis purposes, only a subset of the files are neededfor each subject. 
% For example `/MNINonLinear/Results/rfMRI_REST.tar` (12GB) contains the
% concattenated resting state data, cleaned, in MSMAll standard space.
% `/unprocessed.tar` (2.3GB) contains the raw data.
"""
import os
import sys  
import requests

# specify output folder
outputFolder = os.getcwd() + '/wcs_download/'
baseURL = 'https://syn_hcp.projects.nitrc.org/data/'
#full list of subject id's (n=127)
subjectList = ['644', '13877', '13930', '19478', '26062', '27599', '27600', '27616', '27618', '27619', '27681', '27834', '28696', '28807', '28808', '28809', '28865', '28866', '29067', '29105', '29106', '29107', '29108', '29187', '29257', '29258', '29267', '29378', '29397', '29398', '29407', '29797', '29866', '29867', '29906', '29908', '29909', '29910', '29911', '30017', '30054', '30130', '30273', '30274', '30275', '30283', '30344', '30353', '30459', '30500', '30501', '30502', '30503', '30571', '30572', '30614', '30615', '30616', '30625', '30701', '30737', '30795', '30796', '30797', '30843', '30864', '30873', '30903', '30913', '30914', '30922', '30942', '30946', '30958', '30966', '30983', '31013', '31039', '31040', '31042', '31043', '31063', '31215', '31217', '31287', '31331', '31336', '31374', '31382', '31422', '31466', '31467', '31469', '31565', '31567', '31568', '31569', '31570', '31637', '31692', '31693', '31694', '31734', '31741', '31742', '31743', '31744', '31745', '31746', '31817', '31818', '31829', '31838', '31840', '31871', '31873', '31874', '31875', '31877', '31965', '31967', '32049', '32050', '32051', '32052', '32115', '32116']
#list of all files available for download for each subject (21 entries)
fileList = ['/MNINonLinear/fsaverage_LR32k.tar', '/MNINonLinear/fsaverage.tar', '/MNINonLinear/Native.tar', '/MNINonLinear/other.tar', '/MNINonLinear/Results/motion_figures.tar', '/MNINonLinear/Results/rfMRI_REST1_AP.tar', '/MNINonLinear/Results/rfMRI_REST1_PA.tar', '/MNINonLinear/Results/rfMRI_REST2_AP.tar', '/MNINonLinear/Results/rfMRI_REST2_PA.tar', '/MNINonLinear/Results/rfMRI_REST.tar', '/MNINonLinear/Results/tsnr_figures.tar', '/MNINonLinear/ROIs.tar', '/MNINonLinear/StructuralQC.tar', '/MNINonLinear/xfms.tar', '/rfMRI_REST1_AP.tar', '/rfMRI_REST1_PA.tar', '/rfMRI_REST2_AP.tar', '/rfMRI_REST2_PA.tar', '/T1w.tar', '/T2w.tar', '/unprocessed.tar']
#list of approx file sizes for each available file (21 entries), useful for decidng what to download/where to save etc.
approxSize = {'220M', '17M', '857M', '423M', '90K', '2.3G', '2.3G', '2.3G', '2.3G', '12G', '1010K', '220K', '68M', '142M', '2.3G', '2.3G', '2.3G', '2.3G', '3.0G', '669M', '2.3G'}

# Select desired subjects and files
# specify subjects to download
desiredSubjects = list(range(0, 3))  # for this example case, just download the first 3 subjects
downloadSubjects = [subjectList[i] for i in desiredSubjects]

# specify files to download. see the table variable 'fileList' for all availble files and approximate sizes
desiredFiles = [4, 10]  # for this example case, just download the small files motion and tsnr figures. Change this to [9] to download the rfMRI_REST.tar for selected subjects
downloadFiles = [fileList[i] for i in desiredFiles]

# Loop over subject selection downloading files
for s in range(len(downloadSubjects)):
    subject = downloadSubjects[s]
    for f in range(len(downloadFiles)):
        os.makedirs(outputFolder + subject + os.path.dirname(downloadFiles[f]), exist_ok=True)
        url = baseURL + subject + downloadFiles[f]
        response = requests.get(url, stream=True)
        if response.status_code == 200:
            with open(outputFolder + subject + downloadFiles[f], 'wb') as file:
                for chunk in response.iter_content(chunk_size=8192):
                    file.write(chunk)
        else:
            print(f"Failed to download: {url}")