% 3PAD study MATLAB source adapted from the preprocessing workflow supplied by Marcus Grueschow.
% Repository copy with local inputs and tools supplied at run time. See README.md.
homePath = (getenv('THREEPAD_EST_ROOT'));
cd(homePath);
addpath(fileparts(mfilename('fullpath')));
allSubsEyeStruct = []; allSubsBevMat = [];
mySubDirs = dir(homePath);
currSubsIDs = [];
for k = 1:numel(mySubDirs)
    itemName = mySubDirs(k).name;
    if mySubDirs(k).isdir && ~ismember(itemName, {'.', '..'})
        numericPart = regexp(itemName, '\d+', 'match');
        if ~isempty(numericPart)
            subID = str2double(numericPart{1});
            if ~isnan(subID)
                currSubsIDs = [currSubsIDs; subID];
            end
        end
    end
end
disp('Participant IDs:');
disp(currSubsIDs);
subType = 'R1';
subs2Anal = currSubsIDs';
for i = 1:numel(subs2Anal)
    currentSubID = subs2Anal(i);
    disp(['Processing participant: ', num2str(currentSubID)]);
    [currSubEyeStruct, currSubBevMat] = AE_EST_RAW_R1(currentSubID, 1);
    if isempty(currSubEyeStruct) || isempty(currSubBevMat)
        warning(['Skipping participant ', num2str(currentSubID), ' due to missing R1 data.']);
        continue;
    end
    allSubsEyeStruct(currentSubID).sub = currSubEyeStruct;
    allSubsBevMat = cat(1, allSubsBevMat, currSubBevMat);
    close all;
    disp(['Finished processing participant: ', num2str(currentSubID)]);
end
homePath = (fullfile(getenv('THREEPAD_EST_ROOT'), 'Results_EST'));
cd(homePath);
if strcmp(subType , 'CONTROL')
    AE_EYE_STRUCT_EST_CONTROL = allSubsEyeStruct;
    save('AE_EYE_STRUCT_EST_CONTROL', 'AE_EYE_STRUCT_EST_CONTROL', '-v7.3')
    AE_BEV_MAT_EST_CONTROL = allSubsBevMat;
    save AE_BEV_MAT_EST_CONTROL AE_BEV_MAT_EST_CONTROL
elseif strcmp(subType , 'R1')
    AE_EYE_STRUCT_EST_R1 = allSubsEyeStruct;
    save('AE_EYE_STRUCT_EST_R1', 'AE_EYE_STRUCT_EST_R1', '-v7.3')
    AE_BEV_MAT_EST_R1 = allSubsBevMat;
    save AE_BEV_MAT_EST_R1 AE_BEV_MAT_EST_R1
elseif strcmp(subType , 'OFC')
    AE_EYE_STRUCT_EST_OFC = allSubsEyeStruct;
    save('AE_EYE_STRUCT_EST_OFC', 'AE_EYE_STRUCT_EST_OFC', '-v7.3')
    AE_BEV_MAT_EST_OFC = allSubsBevMat;
    save AE_BEV_MAT_EST_OFC AE_BEV_MAT_EST_OFC
elseif strcmp(subType , 'INS')
    AE_EYE_STRUCT_EST_INS = allSubsEyeStruct;
    save('AE_EYE_STRUCT_EST_INS', 'AE_EYE_STRUCT_EST_INS', '-v7.3')
    AE_BEV_MAT_EST_INS = allSubsBevMat;
    save AE_BEV_MAT_EST_INS AE_BEV_MAT_EST_INS
end
close all; clear all;
homePath = (getenv('THREEPAD_EST_ROOT'));
cd (homePath);
 AE_EST_EYE_ANAL_R1(1, 'R1');
