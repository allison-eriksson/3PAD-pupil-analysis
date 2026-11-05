% 3PAD study MATLAB source adapted from the preprocessing workflow supplied by Marcus Grueschow.
% Repository copy with local inputs and tools supplied at run time. See README.md.
function  [currSubEyeStruct, currSubBevMat] = AE_EST_RAW_R2(subs2Anal, task2Anal)
addpath(getenv('THREEPAD_MATLAB_TOOLS'));
addpath(fullfile(getenv('THREEPAD_MATLAB_TOOLS'), 'analysisScriptEyetrackingHelenaMarcB'));
addpath(genpath(fullfile(getenv('THREEPAD_MATLAB_TOOLS'), 'eeglab2024.0')));
addpath(genpath(fullfile(getenv('THREEPAD_MATLAB_TOOLS'), 'edf-converter-master')));
if nargin < 1
    error('Supply participant codes from the approved local study inputs.');
end
if nargin < 2
task2Anal = 1;
end
if task2Anal == 1
taskName = 'EST';
elseif task2Anal == 2
taskName = 'Brightness';
elseif task2Anal == 3
taskName = 'Oddball';
end
plotFigureSingleRunsBlinks = 1;
plotFigureSingleRunsXYPupVsConditions = 1;
plotFiguresHeatMap = 1;
plotGaussX = 1;
plotGaussY = 1;
saveRunFigs = 1;
saveGazeFigs = 1;
saveHeatFigs = 1;
motherpath = getenv('THREEPAD_EST_ROOT');
mkdir(fullfile(motherpath, 'Results_EST', 'Eye_Figures_Single_Subs2'));
myRunNames = struct();
runCount = 1;
for subIdx = 1:numel(subs2Anal)
    subID = sprintf('%03d', subs2Anal(subIdx));
    currSubFolderName = [filesep, subID];
    potentialRunNames = dir([motherpath, currSubFolderName, filesep, 'r*']);
    isR2 = strcmp({potentialRunNames.name}, 'r2');
    if any(isR2)
        currRunNames = potentialRunNames(isR2);
        disp(['Run 2 Directory found for Participant ID: ' subID]);
        disp(currRunNames);
        currFolderName = [motherpath, currSubFolderName, filesep, currRunNames(1).name, filesep];
        currFieldID = runCount;
        myRunNames(currFieldID).participantID = subID;
        myRunNames(currFieldID).runFolder = currRunNames(1).name;
        myRunNames(currFieldID).EST_EXP_PATH = currFolderName;
        currEDFFolderName = dir([currFolderName, 'EST*']);
        if ~isempty(currEDFFolderName)
            currEDFFolderName = currEDFFolderName(1).name;
            a = dir(fullfile(currFolderName, currEDFFolderName, '*.edf'));
            if ~isempty(a)
                myRunNames(currFieldID).EST_EXP_PATH_with_EDF = fullfile(currFolderName, currEDFFolderName, a(1).name);
            else
                myRunNames(currFieldID).EST_EXP_PATH_with_EDF = '';
            end
            a = dir(fullfile(currFolderName, currEDFFolderName, '*ResultsMat*.mat'));
            if ~isempty(a)
                myRunNames(currFieldID).EST_EXP_PATH_with_MAT = fullfile(currFolderName, currEDFFolderName, a(1).name);
            else
                myRunNames(currFieldID).EST_EXP_PATH_with_MAT = '';
            end
        else
            myRunNames(currFieldID).EST_EXP_PATH_with_EDF = '';
            myRunNames(currFieldID).EST_EXP_PATH_with_MAT = '';
        end
        runCount = runCount + 1;
    else
        warning(['Run 2 directory not found for participant ', subID]);
        currSubEyeStruct = [];
        currSubBevMat = [];
        return;
    end
    disp(['Finished processing Participant ', subID]);
    disp(myRunNames(runCount - 1));
end
disp('All runs data:');
disp(struct2table(myRunNames));
heatGaussSize  = 80;
heatGaussSigma = 12;
centerGaussSize = 0.1;
params.res = [1024 768];
params.sz = [35.5 26.8];
params.vdist =  [78];
params.objSize = [5.2 6.8];
[pixperdeg, degperpix, objectInVisAngle] = VisAngMarc(params);
centralCuePos = params.res/2;
coverInVisAng = objectInVisAngle;
coverInPixel = [pixperdeg(1)*coverInVisAng(1), pixperdeg(2)*coverInVisAng(2)];
coverDiagEuclDistInPixel = sqrt([coverInPixel(1)].^2 + [coverInPixel(2)].^2);
coverDiagEuclDistInVisAng = coverDiagEuclDistInPixel * mean(degperpix);
screenInVisAng = [26.0769   19.6862];
screenInPixel = [pixperdeg(1)*screenInVisAng(1), pixperdeg(2)*screenInVisAng(2)];
screenDiagEuclDistInPixel = sqrt([0 - params.res(1)].^2 + [0 - params.res(2)].^2);
screenDiagEuclDistInVisAng = screenDiagEuclDistInPixel * mean(degperpix);
runs2Anal = 1:size(myRunNames, 2);
edf = cell(numel(subs2Anal), numel(runs2Anal));
bev = cell(numel(subs2Anal), numel(runs2Anal));
heatMapCell = cell(numel(subs2Anal), numel(runs2Anal));
subjCount = 1;
task2Anal = 1;
taskName = 'EST';
for subjIdx = 1:numel(subs2Anal)
    currentSubj = subs2Anal(subjIdx);
    disp(['Working on subject: ', num2str(currentSubj)]);
    runCount = 1;
    for runIdx = runs2Anal
        participantID = myRunNames(runIdx).participantID;
        if str2double(participantID) == currentSubj
            if task2Anal == 1
                currentEdfFilePathAndName = myRunNames(runIdx).EST_EXP_PATH_with_EDF;
                currentMatFilePathAndName = myRunNames(runIdx).EST_EXP_PATH_with_MAT;
            elseif task2Anal == 2
                currentEdfFilePathAndName = myRunNames(runIdx).EST_EXP_Brightness_with_EDF;
                currentMatFilePathAndName = myRunNames(runIdx).EST_EXP_Brightness_with_MAT;
            elseif task2Anal == 3
                currentEdfFilePathAndName = myRunNames(runIdx).EST_EXP_OddBall_with_EDF;
                currentMatFilePathAndName = myRunNames(runIdx).EST_EXP_OddBall_with_MAT;
            end
            if exist(currentEdfFilePathAndName, 'file')
                [~, name, ~] = fileparts(currentEdfFilePathAndName);
                saveMatFilePath = fullfile(fileparts(currentEdfFilePathAndName), [name, '.mat']);
                disp(['Processing EDF file: ', currentEdfFilePathAndName]);
                try
                    edf{subjCount, runCount} = Edf2Mat(currentEdfFilePathAndName, false, true);
                    save(saveMatFilePath, 'edf', '-v7.3');
                    if exist(currentMatFilePathAndName, 'file')
                        a = load(currentMatFilePathAndName);
                        bev{subjCount, runCount} = a.myResultsMat;
                    else
                        disp(['Behavior file not found: ', currentMatFilePathAndName]);
                    end
                    disp(['Processed ', name, '.edf and saved as ', name, '.mat']);
                catch ME
                    disp(['Error with Edf2Mat on file: ', name, '.edf']);
                    disp(ME.message);
                end
            else
                disp(['EDF file not found: ', currentEdfFilePathAndName]);
            end
            runCount = runCount + 1;
        end
    end
    subjCount = subjCount + 1;
end
disp('All participants processed successfully!');
   overallEdfSavePath = fullfile(motherpath , filesep, 'Results_EST',  filesep, 'OverallEdfData.mat');
save(overallEdfSavePath, 'edf', '-v7.3');
disp('Processed successfully and EDF data saved!');
edfData  = reshape(edf.', [], 1);
runs     = repmat(1:numel(runs2Anal), 1, numel(subs2Anal));
subjects = repmat(1:numel(subs2Anal), numel(runs2Anal), 1);
fullGazeData = cell(numel(edfData), 1);
gazeDataChoiceOnlyInterp = cell(numel(edfData), 1);
gazeDataRatingOnlyInterp = cell(numel(edfData), 1);
gazeDataFixationOnlyInterp = cell(numel(edfData), 1);
screenWidth  = cell(numel(edfData), 1);
screenHeight = cell(numel(edfData), 1);
center       = cell(numel(edfData), 1);
centerC      = cell(numel(edfData), 1);
centerR      = cell(numel(edfData), 1);
centerF      = cell(numel(edfData), 1);
heatMapCell  = cell(numel(edfData), 1);
currRunStartTimeCell = cell(numel(edfData), 1);
myCurrRun_eyeT_markers = cell(numel(edfData), 1);
BevStruct = cell(numel(edfData), 1);
currSubPupilStimLocked_m1p9 = [];
currSubPupilRespLocked_m2p8 = [];
currSubBevMat = [];
currSubEyeStruct = [];
participantIndex = 1;
for currSamp = 1:numel(edfData)
    if ~isempty(edfData{currSamp})
        currSubj   = subs2Anal(subjects(currSamp));
        currRun    = runs2Anal(runs(currSamp));
        currData   = edfData{currSamp};
        currBevDataCellName = sprintf('bev{currSamp}');
        currBev.runInfo    = eval(currBevDataCellName);
        tmpRunBlinks2SPM = [];
        tmpRunMicroSacc2SPM = [];
        tmpRunRealSacc2SPM = [];
        try
            infoName = 'GAZE_COORDS';
            screenSize = cellfun(@(x) strncmp(x, infoName, numel(infoName)), currData.Events.Messages.info(:));
            screenSize = strsplit(currData.Events.Messages.info{screenSize}, " ");
            screenSize = str2double(screenSize(end-3:end));
            screenWidth{currSamp}  = screenSize(3) + 1 - screenSize(1);
            screenHeight{currSamp} = screenSize(4) + 1 - screenSize(2);
            currExpStart = 'EXPERIMENT START';
            timestampExpStart = cellfun(@(x) strncmp(x, currExpStart, numel(currExpStart)), currData.Events.Messages.info(:));
            timestampExpStart = currData.Events.Messages.time(timestampExpStart);
        catch
            disp(['CAUGHT: currSub:', num2str(currSubj), ' run: ', num2str(currSamp)]);
            screenWidth{currSamp}  = 1280;
            screenHeight{currSamp} = 1024;
            timestampExpStart = 0;
        end
        currWidth  = screenWidth{currSamp};
        currHeight = screenHeight{currSamp};
        if isempty(currData) || isempty(currSubj) || isempty(currRun)
            error('Data missing for sample index: %d', currSamp);
        end
        ExpStartTime = timestampExpStart(1);
        currDataTimeStepRes = currData.Samples.time(2) - currData.Samples.time(1);
        saccID2Use = currData.Events.Esacc.start > ExpStartTime;
        currSaccOns = currData.Events.Esacc.start(saccID2Use) - ExpStartTime;
        currSaccEnds = currData.Events.Esacc.end(saccID2Use) - ExpStartTime;
        currSaccDur  = currSaccEnds - currSaccOns;
        currSaccEuclDist = sqrt([currData.Events.Esacc.posX(saccID2Use) - currData.Events.Esacc.posXend(saccID2Use)].^2 + [currData.Events.Esacc.posY(saccID2Use) - currData.Events.Esacc.posYend(saccID2Use)].^2);
        currSaccEuclDistInVisAng = currSaccEuclDist * mean(degperpix);
        currSaccOnsInSec = currSaccOns/1000;
        currMicroSaccID = currSaccEuclDistInVisAng < 0.5;
        currMacroSaccID = currSaccEuclDistInVisAng > screenDiagEuclDistInVisAng;
        L =  0.5;
        U = screenDiagEuclDistInVisAng;
        data = currSaccEuclDistInVisAng;
        data(isnan(data)) = 0;
        currRealSaccID = not(abs(sign(sign(L - data) + sign(U - data))));
        currMicroSaccOnsInSec  = currSaccOns(currMicroSaccID)/1000;
        currMicroSaccDur  = currSaccDur(currMicroSaccID)/1000;
        currMicroSaccDist = currSaccEuclDistInVisAng(currMicroSaccID);
        tmpRunMicroSacc2SPM = [currMicroSaccOnsInSec', currMicroSaccDur', currMicroSaccDist'];
        currRealSaccOnsInSec  = currSaccOns(currRealSaccID)/1000;
        currRealSaccDur  = currSaccDur(currRealSaccID)/1000;
        currRealSaccDist = currSaccEuclDistInVisAng(currRealSaccID);
        tmpRunRealSacc2SPM = [currRealSaccOnsInSec', currRealSaccDur', currRealSaccDist'];
        blinkID2Use = currData.Events.Eblink.start > ExpStartTime;
        currBlinkOns = currData.Events.Eblink.start(blinkID2Use) - ExpStartTime;
        currBlinkEnds = currData.Events.Eblink.end(blinkID2Use) - ExpStartTime;
        currBlinkDur  = currBlinkEnds - currBlinkOns;
        blinkBelow80msID = currBlinkDur < 80;
        blinkAbove2000msID = currBlinkDur > 2000;
        currBlinkOns(logical(sum([blinkBelow80msID', blinkAbove2000msID'], 2))) = [];
        currBlinkEnds(logical(sum([blinkBelow80msID', blinkAbove2000msID'], 2))) = [];
        currBlinkDur(logical(sum([blinkBelow80msID', blinkAbove2000msID'], 2))) = [];
        currBlinkOnsInSec = currBlinkOns/1000;
        myBlinkRemovedLinIntPupilSize = currData.Samples.pupilSize(:,1);
        myBlinkRemovedLinIntXPosSize = currData.Samples.posX(:,1);
        myBlinkRemovedLinIntYPosSize = currData.Samples.posY(:,1);
        myBlinkOnsID = [];
        myBlinkEndsID = [];
        currPupilSize = currData.Samples.pupilSize(:,1);
        blinkIDVect = currPupilSize <= 100;
        blinkXAxis = 1:length(currPupilSize);
        vector = blinkIDVect;
        windowWidth = int16(1000);
        halfWidth = windowWidth / 2;
        gaussFilter = gausswin(200);
        gaussFilter = gaussFilter / sum(gaussFilter);
        smoothedVector = conv(double(vector), gaussFilter, 'same');
        smoothedIDVector = logical(ceil(smoothedVector));
        currPosX = currData.Samples.posX(:,1);
        currPosY = currData.Samples.posY(:,1);
        if plotFigureSingleRunsBlinks
            fh = figure(100+currSamp);
            subplot(3,1,1)
            plot(blinkXAxis/1000*currDataTimeStepRes,currPupilSize, 'linewidth', 1); hold on
            plot(blinkXAxis(smoothedIDVector)/1000*currDataTimeStepRes, currPupilSize(smoothedIDVector), 'r-' ); hold on
            title([taskName, ' PUPIL DIAMETER - IDENTIFIED BLINKS IN RED:   SUB: ', num2str(currSubj),'   RUN: ', num2str(currRun) ]);
            subplot(3,1,2)
            plot(blinkXAxis/1000*currDataTimeStepRes, currPosX, 'linewidth', 1); hold on
            plot(blinkXAxis(smoothedIDVector)/1000*currDataTimeStepRes, currPosX(smoothedIDVector), 'r-' ); hold on
            title([taskName, ' X-GAZE - IDENTIFIED BLINKS IN RED:   SUB: ', num2str(currSubj),'   RUN: ', num2str(currRun) ]);
            subplot(3,1,3)
            plot(blinkXAxis/1000*currDataTimeStepRes, currPosY, 'linewidth', 1); hold on
            plot(blinkXAxis(smoothedIDVector)/1000*currDataTimeStepRes, currPosY(smoothedIDVector), 'r-' ); hold on
            title([taskName, ' Y-GAZE - IDENTIFIED BLINKS IN RED:   SUB: ', num2str(currSubj),'   RUN: ', num2str(currRun) ]);
        end
        smoothedIDVectorDiff = diff(smoothedIDVector);
        myBlinkOnsIDVect = logical(cat(1, 0, smoothedIDVectorDiff == 1));
        myBlinkEndsIDVect = logical(cat(1, 0,  smoothedIDVectorDiff == -1));
        if plotFigureSingleRunsBlinks
            subplot(3,1,1)
            plot(blinkXAxis(myBlinkOnsIDVect)/1000*currDataTimeStepRes, currPupilSize(myBlinkOnsIDVect), 'go', 'linewidth', 3 ); hold on
            plot(blinkXAxis(myBlinkEndsIDVect)/1000*currDataTimeStepRes, currPupilSize(myBlinkEndsIDVect), 'ko', 'linewidth', 3); hold on
            subplot(3,1,2)
            plot(blinkXAxis(myBlinkOnsIDVect)/1000*currDataTimeStepRes, currPosX(myBlinkOnsIDVect), 'go', 'linewidth', 3 ); hold on
            plot(blinkXAxis(myBlinkEndsIDVect)/1000*currDataTimeStepRes, currPosX(myBlinkEndsIDVect), 'ko', 'linewidth', 3); hold on
            subplot(3,1,3)
            plot(blinkXAxis(myBlinkOnsIDVect)/1000*currDataTimeStepRes, currPosY(myBlinkOnsIDVect), 'go', 'linewidth', 3 ); hold on
            plot(blinkXAxis(myBlinkEndsIDVect)/1000*currDataTimeStepRes, currPosY(myBlinkEndsIDVect), 'ko', 'linewidth', 3); hold on
        end
        myBlinkOnsID = find(myBlinkOnsIDVect == 1);
        myBlinkEndsID = find(myBlinkEndsIDVect == 1);
        if ~isempty(myBlinkOnsID)
            blinkStartFlag = 0;
            blinkEndFlag = 0;
            if myBlinkOnsID(1) > myBlinkEndsID(1)
                myBlinkOnsID = cat(1, 1, myBlinkOnsID);
                blinkStartFlag = 1;
            end
            if myBlinkOnsID(end) > myBlinkEndsID(end)
                myBlinkEndsID = cat(1, myBlinkEndsID, length(currPupilSize));
                blinkEndFlag = 1;
            end
            for currBlink = 1 : length(myBlinkEndsID)
                myBlinkRemovedLinIntPupilSize(myBlinkOnsID(currBlink) : myBlinkEndsID(currBlink)) = linspace(currPupilSize(myBlinkOnsID(currBlink)), currPupilSize(myBlinkEndsID(currBlink)), (myBlinkEndsID(currBlink) - myBlinkOnsID(currBlink))+1 ) ;
                myBlinkRemovedLinIntXPosSize(myBlinkOnsID(currBlink) : myBlinkEndsID(currBlink)) = linspace(currPosX(myBlinkOnsID(currBlink)), currPosX(myBlinkEndsID(currBlink)), (myBlinkEndsID(currBlink) - myBlinkOnsID(currBlink))+1 ) ;
                myBlinkRemovedLinIntYPosSize(myBlinkOnsID(currBlink) : myBlinkEndsID(currBlink)) = linspace(currPosY(myBlinkOnsID(currBlink)), currPosY(myBlinkEndsID(currBlink)), (myBlinkEndsID(currBlink) - myBlinkOnsID(currBlink))+1 ) ;
                if  myBlinkOnsID(currBlink)/1000*currDataTimeStepRes
                    tmpRunBlinks2SPM = cat(1, tmpRunBlinks2SPM, [myBlinkOnsID(currBlink)/1000*currDataTimeStepRes (myBlinkEndsID(currBlink) - myBlinkOnsID(currBlink))/1000*currDataTimeStepRes]);
                end
            end
            if blinkStartFlag == 1
                myBlinkRemovedLinIntPupilSize(myBlinkOnsID(1) : myBlinkEndsID(1)) = linspace(currPupilSize(myBlinkEndsID(1)), currPupilSize(myBlinkEndsID(1)), (myBlinkEndsID(1) - myBlinkOnsID(1))+1 ) ;
                myBlinkRemovedLinIntXPosSize(myBlinkOnsID(1) : myBlinkEndsID(1)) = linspace(currPosX(myBlinkEndsID(1)), currPosX(myBlinkEndsID(1)), (myBlinkEndsID(1) - myBlinkOnsID(1))+1 ) ;
                myBlinkRemovedLinIntYPosSize(myBlinkOnsID(1) : myBlinkEndsID(1)) = linspace(currPosY(myBlinkEndsID(1)), currPosY(myBlinkEndsID(1)), (myBlinkEndsID(1) - myBlinkOnsID(1))+1 ) ;
                blinkStartFlag = 0;
            end
            if blinkEndFlag == 1
                myBlinkRemovedLinIntPupilSize(myBlinkOnsID(end) : myBlinkEndsID(end)) = linspace(currPupilSize(myBlinkOnsID(end)), currPupilSize(myBlinkOnsID(end)), (myBlinkEndsID(end) - myBlinkOnsID(end))+1 ) ;
                myBlinkRemovedLinIntXPosSize(myBlinkOnsID(end) : myBlinkEndsID(end)) = linspace(currPosX(myBlinkOnsID(end)), currPosX(myBlinkOnsID(end)), (myBlinkEndsID(end) - myBlinkOnsID(end))+1 ) ;
                myBlinkRemovedLinIntYPosSize(myBlinkOnsID(end) : myBlinkEndsID(end)) = linspace(currPosY(myBlinkOnsID(end)), currPosY(myBlinkOnsID(end)), (myBlinkEndsID(end) - myBlinkOnsID(end))+1 ) ;
                blinkEndFlag = 0;
            end
            if plotFigureSingleRunsBlinks
                subplot(3,1,1)
                plot([1:length(currData.Samples.pupilSize)]/1000*currDataTimeStepRes, myBlinkRemovedLinIntPupilSize, 'g-', 'linewidth', 2); hold on
                subplot(3,1,2)
                plot([1:length(currData.Samples.pupilSize)]/1000*currDataTimeStepRes, myBlinkRemovedLinIntXPosSize, 'g-', 'linewidth', 2); hold on
                subplot(3,1,3)
                plot([1:length(currData.Samples.pupilSize)]/1000*currDataTimeStepRes, myBlinkRemovedLinIntYPosSize, 'g-', 'linewidth', 2); hold on
            end
        end
        if saveRunFigs==1
            if plotFigureSingleRunsBlinks
                myHome = pwd;
                cd([motherpath, filesep, 'Results_EST', filesep, 'Eye_Figures_Single_Subs2']);
                set(fh, 'position', [100 100 3700 900])
                figure(100+currSamp); set(gcf, 'paperunits', 'centimeters', 'paperposition', [0 0 60 15])
                print([taskName, '_BlinkDetectRemove_Sub',num2str(currSubj),'_Runs', num2str(currRun) ,'.tif'], '-dtiff', '-r300')
                cd(myHome)
            end
        end
        gazeData = [currPosY, currPosX, currData.Samples.time - timestampExpStart(1), currPupilSize];
        gazeDataInterp = [myBlinkRemovedLinIntYPosSize, myBlinkRemovedLinIntXPosSize, currData.Samples.time - timestampExpStart(1), myBlinkRemovedLinIntPupilSize];
        dataselectorStart = gazeData(:,3) >= 0;
        gazeDataExpOnly = gazeData(dataselectorStart,:);
        gazeDataExpOnlyInterp = gazeDataInterp(dataselectorStart,:);
        expTimeBeam = gazeDataExpOnly(:,3);
        currRunMat = currBev.runInfo;
        currRunIRE = 0;
        for currT = 2 : size(currRunMat,1)
            if task2Anal == 1
                if currRunMat(currT-1, 3) == 1 && currRunMat(currT, 3) == 1
                    currRunIRE(currT) = 1;
                elseif currRunMat(currT-1, 3) == 2 && currRunMat(currT, 3) == 1
                    currRunIRE(currT) = 2;
                elseif currRunMat(currT-1, 3) == 1 && currRunMat(currT, 3) == 2
                    currRunIRE(currT) = 3;
                elseif currRunMat(currT-1, 3) == 2 && currRunMat(currT, 3) == 2
                    currRunIRE(currT) = 4;
                end
            elseif task2Anal == 2
                if currRunMat(currT, 4) == 1
                    currRunIRE(currT) = 1;
                elseif currRunMat(currT, 4) == 2
                    currRunIRE(currT) = 2;
                elseif currRunMat(currT, 4) == 3
                    currRunIRE(currT) = 3;
                elseif currRunMat(currT, 4) == 4
                    currRunIRE(currT) = 4;
                end
            else
                if currRunMat(currT, 4) == 0
                    currRunIRE(currT) = 3;
                elseif currRunMat(currT, 4) == 1
                    currRunIRE(currT) = 1;
                elseif currRunMat(currT, 4) == 2
                    currRunIRE(currT) = 2;
                end
            end
        end
        currRunMat = cat(2, currRunMat, currRunIRE');
        currBev.runInfo = currRunMat;
        if task2Anal == 1
            currBevE.runInfo = currBev.runInfo(currBev.runInfo(:,end) <= 2, :);
            currBevP.runInfo = currBev.runInfo(currBev.runInfo(:,end) >= 3, :);
        elseif task2Anal == 2
            currBevE.runInfo = currBev.runInfo(currBev.runInfo(:,end) <= 2, :);
            currBevP.runInfo = currBev.runInfo(currBev.runInfo(:,end) >= 3, :);
        else
            currBevE.runInfo = currBev.runInfo(currBev.runInfo(:,end) == 0, :);
            currBevP.runInfo = currBev.runInfo(currBev.runInfo(:,end) >= 1, :);
        end
        Ons_E = currBevE.runInfo(:,8);
        if task2Anal == 3
            RT_E = currBevE.runInfo(:,5)+100;
        else
            RT_E = currBevE.runInfo(:,5);
        end
        Ons_P = currBevP.runInfo(:,8);
         if task2Anal == 3
             RT_P = currBevP.runInfo(:,5)+100;
         else
             RT_P = currBevP.runInfo(:,5);
         end
        logicalVectChoiceE = false(size(expTimeBeam));
        for currTrial = 1 : size(currBevE.runInfo, 1)
            L =  Ons_E(currTrial);
            U = Ons_E(currTrial) + RT_E(currTrial);
            data = expTimeBeam;
            if ~isnan(U)
                inds = not(abs(sign(sign(L - data) + sign(U - data))));
                logicalVectChoiceE(inds) = 1;
            end
        end
        logicalVectChoiceP = false(size(expTimeBeam));
        for currTrial = 1 : size(currBevP.runInfo, 1)
            L =  Ons_P(currTrial);
            U = Ons_P(currTrial) + RT_P(currTrial);
            data = expTimeBeam;
            if ~isnan(U)
                inds = not(abs(sign(sign(L - data) + sign(U - data))));
                logicalVectChoiceP(inds) = 1;
            end
        end
        currExpRunPupil = gazeDataExpOnly(:,4);
        currExpRunX = gazeDataExpOnly(:,2);
        currExpRunY = gazeDataExpOnly(:,1);
        if plotFigureSingleRunsXYPupVsConditions
            fh = figure(200+currSamp);
            subplot(3,1,1)
            plot(expTimeBeam/1000, logicalVectChoiceE*max(currExpRunPupil), 'c'); hold on;
            plot(expTimeBeam/1000, logicalVectChoiceP*max(currExpRunPupil), 'm'); hold on;
            plot(expTimeBeam/1000, gazeDataExpOnlyInterp(:,4), 'g-', 'linewidth', 2); hold on
            title([taskName, ' Pupil dilation during experiment (cyan = congruent trials; magenta = incongruent trials) SUB: ', num2str(currSubj),'   RUN: ', num2str(currRun) ]);
            xlim([0 round(expTimeBeam(end)/1000+5)]);
            subplot(3,1,2)
            plot(expTimeBeam/1000, logicalVectChoiceE*max(gazeDataExpOnlyInterp(:,2)), 'c'); hold on;
            plot(expTimeBeam/1000, logicalVectChoiceP*max(gazeDataExpOnlyInterp(:,2)), 'm'); hold on;
            plot(expTimeBeam/1000, gazeDataExpOnlyInterp(:,2), 'g-', 'linewidth', 2); hold on
            title([taskName, ' X-position during experiment (cyan = congruent trials; magenta = incongruent trials) SUB: ', num2str(currSubj),'   RUN: ', num2str(currRun) ]);
            xlim([0 round(expTimeBeam(end)/1000+5)]);
            subplot(3,1,3)
            plot(expTimeBeam/1000, logicalVectChoiceE*max(gazeDataExpOnlyInterp(:,1)), 'c'); hold on;
            plot(expTimeBeam/1000, logicalVectChoiceP*max(gazeDataExpOnlyInterp(:,1)), 'm'); hold on;
            plot(expTimeBeam/1000, gazeDataExpOnlyInterp(:,1), 'g-', 'linewidth', 2); hold on
            title([taskName, ' Y-position during experiment (cyan = congruent trials; magenta = incongruent trials) SUB: ', num2str(currSubj),'   RUN: ', num2str(currRun) ]);
            xlim([0 round(expTimeBeam(end)/1000+5)]);
            set(fh, 'position', [100 100 3700 900])
            if saveRunFigs == 1
                myHome = pwd;
                cd([motherpath, filesep, 'Results_EST', filesep, 'Eye_Figures_Single_Subs2']);
                set(fh, 'position', [100 100 3700 900])
                figure(200+currSamp); set(gcf, 'paperunits', 'centimeters', 'paperposition', [0 0 60 15])
                print([taskName, ' PupilAndGazeVsChoiceAndRating_Sub',num2str(currSubj),'_Runs', num2str(currRun) ,'.tif'], '-dtiff', '-r300')
                cd(myHome)
            end
        end
        logicalVectFixation = ~logical(logicalVectChoiceE + logicalVectChoiceP);
        gazeDataChoiceOnlyInterp{currSamp} = gazeDataExpOnlyInterp(logicalVectChoiceE,:);
        gazeDataRatingOnlyInterp{currSamp} = gazeDataExpOnlyInterp(logicalVectChoiceP,:);
        gazeDataFixationOnlyInterp{currSamp} = gazeDataExpOnlyInterp(logicalVectFixation,:);
        fullGazeData{currSamp} = gazeData(any(~isnan(gazeData(:, 1:2)), 2) & dataselectorStart, :);
        histcounts2_fulldata = histcounts2(fullGazeData{currSamp}(:,1), fullGazeData{currSamp}(:,2), 1:currHeight+1, 1:currWidth+1);
        histcounts2_Choice = histcounts2(gazeDataChoiceOnlyInterp{currSamp}(:,1), gazeDataChoiceOnlyInterp{currSamp}(:,2), 1:currHeight+1, 1:currWidth+1);
        histcounts2_Rating = histcounts2(gazeDataRatingOnlyInterp{currSamp}(:,1), gazeDataRatingOnlyInterp{currSamp}(:,2), 1:currHeight+1, 1:currWidth+1);
        histcounts2_Fixation = histcounts2(gazeDataFixationOnlyInterp{currSamp}(:,1), gazeDataFixationOnlyInterp{currSamp}(:,2), 1:currHeight+1, 1:currWidth+1);
        kernel = createGauss(heatGaussSize, heatGaussSigma);
        heatMapData = conv2(histcounts2_fulldata, kernel, 'same');
        heatMapCell{currSamp} = heatMapData;
        heatMapDataChoice = conv2(histcounts2_Choice, kernel, 'same');
        heatMapCellC{currSamp} = heatMapDataChoice;
        heatMapDataRating = conv2(histcounts2_Rating, kernel, 'same');
        heatMapCellR{currSamp} = heatMapDataRating;
        heatMapDataFixation = conv2(histcounts2_Fixation, kernel, 'same');
        heatMapCellF{currSamp} = heatMapDataFixation;
        gaussSizeX = centerGaussSize .* currWidth;
        gaussSizeY = centerGaussSize .* currHeight;
        gaussFuncX = gaussmf(1:currWidth, [gaussSizeX, currWidth/2]).';
        gaussFuncY = gaussmf(1:currHeight, [gaussSizeY, (currHeight/2)]).';
        profileHorizontal = max(heatMapData, [], 1).';
        profileVertical   = max(heatMapData, [], 2);
        profileHorizontalChoice = max(heatMapDataChoice, [], 1).';
        profileVerticalChoice   = max(heatMapDataChoice, [], 2);
        profileHorizontalRating = max(heatMapDataRating, [], 1).';
        profileVerticalRating   = max(heatMapDataRating, [], 2);
        profileHorizontalFixation = max(heatMapDataFixation, [], 1).';
        profileVerticalFixation   = max(heatMapDataFixation, [], 2);
        [peakValueX, peakIdxX] = findpeaks(profileHorizontal);
        [peakValueY, peakIdxY] = findpeaks(profileVertical);
        [peakValueXC, peakIdxXC] = findpeaks(profileHorizontalChoice);
        [peakValueYC, peakIdxYC] = findpeaks(profileVerticalChoice);
        [peakValueXR, peakIdxXR] = findpeaks(profileHorizontalRating);
        [peakValueYR, peakIdxYR] = findpeaks(profileVerticalRating);
        [peakValueXF, peakIdxXF] = findpeaks(profileHorizontalFixation);
        [peakValueYF, peakIdxYF] = findpeaks(profileVerticalFixation);
        [gaussPeakValueX, gaussPeakIdxX] = findpeaks(profileHorizontal .* gaussFuncX);
        [~, gaussMaxPeakX]  = max(gaussPeakValueX);
        [~, centerPeakIdxX] = min(abs(peakIdxX - gaussPeakIdxX(gaussMaxPeakX)));
        [gaussPeakValueXC, gaussPeakIdxXC] = findpeaks(profileHorizontalChoice .* gaussFuncX);
        [~, gaussMaxPeakXC]  = max(gaussPeakValueXC);
        [~, centerPeakIdxXC] = min(abs(peakIdxXC - gaussPeakIdxXC(gaussMaxPeakXC)));
        [gaussPeakValueXR, gaussPeakIdxXR] = findpeaks(profileHorizontalRating .* gaussFuncX);
        [~, gaussMaxPeakXR]  = max(gaussPeakValueXR);
        [~, centerPeakIdxXR] = min(abs(peakIdxXR - gaussPeakIdxXR(gaussMaxPeakXR)));
        [gaussPeakValueXF, gaussPeakIdxXF] = findpeaks(profileHorizontalFixation .* gaussFuncX);
        [~, gaussMaxPeakXF]  = max(gaussPeakValueXF);
        [~, centerPeakIdxXF] = min(abs(peakIdxXF - gaussPeakIdxXF(gaussMaxPeakXF)));
        [gaussPeakValueY, gaussPeakIdxY] = findpeaks(profileVertical .* gaussFuncY);
        [~, gaussMaxPeakY]  = max(gaussPeakValueY);
        [~, centerPeakIdxY] = min(abs(peakIdxY - gaussPeakIdxY(gaussMaxPeakY)));
        [gaussPeakValueYC, gaussPeakIdxYC] = findpeaks(profileVerticalChoice .* gaussFuncY);
        [~, gaussMaxPeakYC]  = max(gaussPeakValueYC);
        [~, centerPeakIdxYC] = min(abs(peakIdxYC - gaussPeakIdxYC(gaussMaxPeakYC)));
        [gaussPeakValueYR, gaussPeakIdxYR] = findpeaks(profileVerticalRating .* gaussFuncY);
        [~, gaussMaxPeakYR]  = max(gaussPeakValueYR);
        [~, centerPeakIdxYR] = min(abs(peakIdxYR - gaussPeakIdxYR(gaussMaxPeakYR)));
        [gaussPeakValueYF, gaussPeakIdxYF] = findpeaks(profileVerticalFixation .* gaussFuncY);
        [~, gaussMaxPeakYF]  = max(gaussPeakValueYF);
        [~, centerPeakIdxYF] = min(abs(peakIdxYF - gaussPeakIdxYF(gaussMaxPeakYF)));
        center{currSamp} = [peakIdxX(centerPeakIdxX), peakIdxY(centerPeakIdxY)];
        centerC{currSamp} = [peakIdxXC(centerPeakIdxXC), peakIdxYC(centerPeakIdxYC)];
        centerR{currSamp} = [peakIdxXR(centerPeakIdxXR), peakIdxYR(centerPeakIdxYR)];
        centerF{currSamp} = [peakIdxXF(centerPeakIdxXF), peakIdxYF(centerPeakIdxYF)];
        minTot = min(min(heatMapData));
        maxTot = max(max(heatMapData));
        minTotC = min(min(heatMapDataChoice));
        maxTotC = max(max(heatMapDataChoice));
        minTotR = min(min(heatMapDataRating));
        maxTotR = max(max(heatMapDataRating));
        minTotF = min(min(heatMapDataFixation));
        maxTotF = max(max(heatMapDataFixation));
        if plotGaussX
            fh1 = figure(1);
            subplot(4, numel(runs2Anal), currSamp);
            hold on;
            plotGaussXFunction(currWidth, profileHorizontal, gaussFuncX, center, minTot, maxTot, [taskName, ' All Data'], currSamp, currSubj, currRun)
            figure(1); subplot(4, numel(runs2Anal), currSamp+numel(runs2Anal)); hold on;
            plotGaussXFunction(currWidth, profileHorizontalChoice, gaussFuncX, centerC, minTotC, maxTotC, 'CONGRUENT TRIALS', currSamp, currSubj, currRun)
            figure(1); subplot(4, numel(runs2Anal), currSamp+numel(runs2Anal)*2); hold on;
            plotGaussXFunction(currWidth, profileHorizontalRating, gaussFuncX, centerR, minTotR, maxTotR, 'INCONGRUENT TRIALS', currSamp, currSubj, currRun)
            figure(1); subplot(4, numel(runs2Anal), currSamp+numel(runs2Anal)*3); hold on;
            plotGaussXFunction(currWidth, profileHorizontalFixation, gaussFuncX, centerF, minTotF, maxTotF, 'FIXATION', currSamp, currSubj, currRun)
            set(fh1, 'position', [20 100 1000 900])
        end
        if plotGaussY
            fh2 = figure(2);
            subplot(4, numel(runs2Anal), currSamp);
            hold on;
            plotGaussYFunction(currHeight, profileVertical, gaussFuncY, center, minTot, maxTot, [taskName, ' All Data'], currSamp, currSubj, currRun)
            figure(2); subplot(4, numel(runs2Anal), currSamp+numel(runs2Anal)); hold on;
            plotGaussYFunction(currHeight, profileVerticalChoice, gaussFuncY, centerC, minTotC, maxTotC, 'CONGRUENT TRIALS', currSamp, currSubj, currRun)
            figure(2); subplot(4, numel(runs2Anal), currSamp+numel(runs2Anal)*2); hold on;
            plotGaussYFunction(currHeight, profileVerticalRating, gaussFuncY, centerR, minTotR, maxTotR, 'INCONGRUENT TRIALS', currSamp, currSubj, currRun)
            figure(2); subplot(4, numel(runs2Anal), currSamp+numel(runs2Anal)*3); hold on;
            plotGaussYFunction(currHeight, profileVerticalFixation, gaussFuncY, centerF, minTotF, maxTotF, 'FIXATION', currSamp, currSubj, currRun)
            set(fh2, 'position', [920 100 1000 900])
        end
        if currData.Samples.time(2) - currData.Samples.time(1) == 2
            tmp = [gazeDataExpOnlyInterp(:,2)'; gazeDataExpOnlyInterp(:,2)']; currX = tmp(:);
            tmp = [gazeDataExpOnlyInterp(:,1)'; gazeDataExpOnlyInterp(:,1)']; currY = tmp(:);
            tmp = [gazeDataExpOnlyInterp(:,4)'; gazeDataExpOnlyInterp(:,4)']; currP = tmp(:);
        elseif currData.Samples.time(2) - currData.Samples.time(1) == 4
            tmp = [gazeDataExpOnlyInterp(:,2)'; gazeDataExpOnlyInterp(:,2)'; gazeDataExpOnlyInterp(:,2)'; gazeDataExpOnlyInterp(:,2)']; currX = tmp(:);
            tmp = [gazeDataExpOnlyInterp(:,1)'; gazeDataExpOnlyInterp(:,1)'; gazeDataExpOnlyInterp(:,1)'; gazeDataExpOnlyInterp(:,1)']; currY = tmp(:);
            tmp = [gazeDataExpOnlyInterp(:,4)'; gazeDataExpOnlyInterp(:,4)'; gazeDataExpOnlyInterp(:,4)'; gazeDataExpOnlyInterp(:,4)']; currP = tmp(:);
        else
            tmp = [gazeDataExpOnlyInterp(:,2)']; currX = tmp(:);
            tmp = [gazeDataExpOnlyInterp(:,1)']; currY = tmp(:);
            tmp = [gazeDataExpOnlyInterp(:,4)']; currP = tmp(:);
        end
        currT = 1:length(currX);
        currGazeFixCenterCoords = centerF{currSamp};
        currSubBevMat = cat(1, currSubBevMat, [currBevE.runInfo,  ones(size(currBevE.runInfo,1),1)*currSubj, ones(size(currBevE.runInfo,1),1)*currRun]);
        currSubBevMat = cat(1, currSubBevMat, [currBevP.runInfo,  ones(size(currBevP.runInfo,1),1)*currSubj, ones(size(currBevP.runInfo,1),1)*currRun]);
        currSubEyeStruct(participantIndex).subjectID = currSubj;
        currSubEyeStruct(participantIndex).run(currRun).currX = currX;
        currSubEyeStruct(participantIndex).run(currRun).currY = currY;
        currSubEyeStruct(participantIndex).run(currRun).currP = currP;
        currSubEyeStruct(participantIndex).run(currRun).currT = currT';
        currSubEyeStruct(participantIndex).run(currRun).center = currGazeFixCenterCoords;
        participantIndex = participantIndex + 1;
    end
end
if plotFiguresHeatMap == 1
    for currSamp = 1:numel(edfData)
        fh3 = figure(3);
        subplot(4, numel(runs2Anal), currSamp);
        hold on
        if  ~isempty(edfData{currSamp})
            plotHeatMap(subs2Anal, runs2Anal, subjects, runs, edfData, screenWidth, screenHeight, currSamp, heatMapCell, center, [taskName, ' All Data'])
            figure(3); subplot(4, numel(runs2Anal), currSamp+numel(runs2Anal)); hold on;
            plotHeatMap(subs2Anal, runs2Anal, subjects, runs, edfData, screenWidth, screenHeight, currSamp, heatMapCellC, centerC, 'CONGRUENT TRIALS')
            figure(3); subplot(4, numel(runs2Anal), currSamp+numel(runs2Anal)*2); hold on;
            plotHeatMap(subs2Anal, runs2Anal, subjects, runs, edfData, screenWidth, screenHeight, currSamp, heatMapCellR, centerR, 'INCONGRUENT TRIALS')
            figure(3); subplot(4, numel(runs2Anal), currSamp++numel(runs2Anal)*3); hold on;
            plotHeatMap(subs2Anal, runs2Anal, subjects, runs, edfData, screenWidth, screenHeight, currSamp, heatMapCellF, centerF, 'FIXATION')
            set(fh3, 'position', [1920 100 1000 900])
        end
    end
end
if saveGazeFigs==1
    myHome = pwd;
    cd([motherpath, filesep, 'Results_EST', filesep, 'Eye_Figures_Single_Subs2']);
    figure(1); set(gcf, 'paperunits', 'centimeters', 'paperposition', [0 0 30 15])
    print([taskName, ' GazeX_Sub',num2str(currSubj),'_AllChoiceRuns.tif'], '-dtiff', '-r300')
    figure(2); set(gcf, 'paperunits', 'centimeters', 'paperposition', [0 0 30 15])
    print([taskName, ' GazeY_Sub',num2str(currSubj),'_AllChoiceRuns.tif'], '-dtiff', '-r300')
    cd(myHome)
end
if saveHeatFigs==1
    myHome = pwd;
    cd([motherpath, filesep, 'Results_EST', filesep, 'Eye_Figures_Single_Subs2']);
    figure(3); set(gcf, 'paperunits', 'centimeters', 'paperposition', [0 0 30 15])
    print([taskName, ' HeatMaps_Sub',num2str(currSubj),'_AllChoiceRuns.tif'], '-dtiff', '-r300')
    close all;
    cd(myHome)
end
    function plotGaussXFunction(currWidth, profileHorizontal, gaussFuncX, center, minTot, maxTot, titleTextString, currSamp, currSubj, currRun)
        plot(1:currWidth, profileHorizontal,'b');
        plot(1:currWidth, profileHorizontal .* gaussFuncX,'r--');
        plot(ones(1, 2) * center{currSamp}(1), [minTot, maxTot * 1.05], 'g-.');
        plot(1:currWidth, gaussFuncX .* maxTot, 'k:');
        plot(ones(1, 2) * currWidth/2, [minTot, maxTot], 'k:');
        hold off;
        axis tight
        title([ titleTextString, ': ', sprintf('Sub %d, Run %d', currSubj, currRun)])
        xlabel('X')
        ylabel('frequency')
    end
    function plotGaussYFunction(currHeight, profileVertical, gaussFuncY, center, minTot, maxTot, titleTextString, currSamp, currSubj, currRun)
        plot(1:currHeight, profileVertical,'b');
        plot(1:currHeight, profileVertical .* gaussFuncY,'r--');
        plot(ones(1, 2) * center{currSamp}(2), [minTot, maxTot * 1.05], 'g-.');
        plot(1:currHeight, gaussFuncY .* maxTot, 'k:');
        plot(ones(1, 2) * currHeight/2, [minTot, maxTot], 'k:');
        hold off;
        axis tight
        title([ titleTextString, ': ', sprintf('Sub %d, Run %d', currSubj, currRun)])
        xlabel('Y')
        ylabel('frequency')
    end
    function plotHeatMap(subs2Anal, runs2Anal, subjects, runs, edfData, screenWidth, screenHeight, currSamp, heatMapData, center, titleTextString)
        currSubj   = subs2Anal(subjects(currSamp));
        currRun    = runs2Anal(runs(currSamp));
        currData   = edfData{currSamp};
        currWidth  = screenWidth{currSamp};
        currHeight = screenHeight{currSamp};
        xRange = [-currWidth/2, currWidth/2];
        yRange = [-currHeight/2, currHeight/2];
        imgPtr = imagesc(heatMapData{runs(currSamp)});
        hold on;
        plot(xRange, yRange(1) + ones(1,2) * center{currSamp}(2), 'k--');
        plot(xRange(1) + ones(1,2) * center{currSamp}(1), yRange, 'k--');
        hold off;
        axis xy;
        axis equal;
        set(gca,'Ydir','reverse')
        title([ titleTextString, ': ', sprintf('Sub %d, Run %d', currSubj, currRun)])
        imgPtr.XData       = xRange;
        imgPtr.YData       = yRange;
        imgPtr.Parent.XLim = xRange;
        imgPtr.Parent.YLim = yRange;
        colormap hsv;
        drawnow;
    end
 end
