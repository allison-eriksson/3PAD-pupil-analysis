% 3PAD study MATLAB source adapted from the preprocessing workflow supplied by Marcus Grueschow.
% Repository copy with local inputs and tools supplied at run time. See README.md.
function AE_EST_EYE_ANAL_R1(task2Anal,subType )
task2Anal = 1;
if task2Anal == 1
taskName = 'EST';
elseif task2Anal == 2
taskName = 'Brightness';
elseif task2Anal == 3
taskName = 'Oddball';
end
subType = 'R1';
saveFigures = 1;
myMotherPath = getenv('THREEPAD_EST_ROOT');
addpath(getenv('THREEPAD_MATLAB_TOOLS'));
addpath(fullfile(getenv('THREEPAD_MATLAB_TOOLS'), 'analysisScriptEyetrackingHelenaMarcB'));
addpath(genpath(fullfile(getenv('THREEPAD_MATLAB_TOOLS'), 'eeglab2024.0')));
addpath(genpath(fullfile(getenv('THREEPAD_MATLAB_TOOLS'), 'edf-converter-master')));
tic
doExtraction = 1;
doLowPassFilteringOnly = 0;
doBandPassFiltering = 1;
saveCoditionResultsMats = 1;
saveMatForDDM = 0;
doMoveAvg = 0;
moveAvgKernel = 500;
doLTSO_BestAIC = 1;
doLTSO_PupDistOnly = 1;
doLTSO_FullMo_LC_Pup_PPI = 1;
doLTSO_FullMo_LC_Pup = 1;
doLTSO_FullMo_PupOnly = 1;
doLTSO_FullMo_LCOnly = 1;
doLTSO_FullMo_Only = 1;
doLTSO_FullMo_noBevCA = 1;
nI = 1000;
whichPup = 1;
whichROI = 0;
withEye = 0;
VOIradius = 5;
myPThresh = 0.05;
tic
zscorePupil = 1;
plotSingleSubsMeanPupil = 1;
plotGranMeanPupil = 1;
plotSingleTrialYGazeChoice = 0;
myLineWidth = 2;
singleHue5ShadesGrey =  [247,247,247; 204,204,204; 150,150,150; 99,99,99;    37,37,37];
singleHue5ShadesGreen = [237,248,233; 186,228,179; 116,196,118; 49,163,84;   0,109,44];
singleHue5ShadesBlue =  [239,243,255; 189,215,231; 107,174,214; 49,130,189;  8,81,156];
singleHue5ShadesRed =   [254,229,217; 252,174,145; 251,106,74;  222,45,38;   165,15,21];
singleHue5ShadesLilac = [242,240,247; 203,201,226; 158,154,200; 117,107,177; 84,39,143];
singleHue5ShadesBrown = [254,237,222; 253,190,133; 253,141,60;  230,85,13;   166,54,3];
singleHue5ShadesPurple =[254,235,226; 251,180,185; 247,104,161; 197,27,138;  122,1,119];
myEconColor = [singleHue5ShadesGreen(3,:)]/255;
myPercColor = [singleHue5ShadesPurple(3,:)]/255;
myEconColor2 = [singleHue5ShadesGreen(4,:)]/255;
myPercColor2 = [singleHue5ShadesPurple(4,:)]/255;
myRunColor = [0 0 0; 0.2 0.2 0.2; 0.4 0.4 0.4; 0.6 0.6 0.6];
myEconStimOnColor = [.2 .2 .2];
myEconRTColor = [.4 .4 .4];
myPercStimOnColor = [.2 .2 .2];
myPercRTColor = [.4 .4 .4];
myCongruentColor = [0 255 0]/255;
myIncongruentColor = [255 0 0]/255;
myICColor = [0;139;0]/255;
myCCColor = [0;205;0]/255;
myCIColor = [255;36;0]/255;
myIIColor = [255;127;0]/255;
myTtestStarColor = [0 127 255]/255;
myCorrelationStarColor = [108 166 205]/255;
myRed = [255, 0, 0]/255;
myOrange = [255, 165, 0]/255;
myBlue = [0, 165, 255]/255;
myDarkRed = [0.5430 0 0];
myIndianRed = [0.8,0.3,0.3];
myGrey = singleHue5ShadesGrey(3,:)/255;
myPupilPlotXlimVals = [6000 15000];
myPupilPlotYlimVals = [-3 2];
myAllTrialsColor = [0.7 0.7 0.7 ];
myBarWidth = 0.5;
myErrBarLineWidth = 1;
myXYFontSize = 18;
myFontSize = 20;
if doExtraction
    if task2Anal == 1
        if strcmp(subType , 'CONTROL')
            load([myMotherPath ,filesep, 'AE_EYE_STRUCT_EST_CONTROL.mat']);
            EYE_MAT_R1 = AE_EYE_STRUCT_EST_CONTROL;
            load([myMotherPath , filesep, 'AE_BEV_MAT_EST_CONTROL.mat']);
            BEV_MAT_RAW = AE_BEV_MAT_EST_CONTROL;
        elseif strcmp(subType , 'R1')
            load([myMotherPath , filesep, 'Results_EST', filesep, 'AE_EYE_STRUCT_EST_R1.mat']);
            EYE_MAT_RAW = AE_EYE_STRUCT_EST_R1;
            load([myMotherPath ,filesep,'Results_EST', filesep,'AE_BEV_MAT_EST_R1.mat']);
            BEV_MAT_RAW = AE_BEV_MAT_EST_R1;
        elseif strcmp(subType , 'OFC')
            load([myMotherPath ,'AE_EYE_STRUCT_EST_OFC.mat']);
            EYE_MAT_R1 = AE_EYE_STRUCT_EST_OFC;
            load([myMotherPath ,'AE_BEV_MAT_EST_OFC.mat']);
            BEV_MAT_R1 = AE_BEV_MAT_EST_OFC;
        elseif strcmp(subType , 'DLPFC')
            load([myMotherPath ,'AE_EYE_STRUCT_EST_DLPFC.mat']);
            EYE_MAT_R1 = AE_EYE_STRUCT_EST_DLPFC;
            load([myMotherPath ,'AE_BEV_MAT_EST_DLPFC.mat']);
            BEV_MAT_R1 = AE_BEV_MAT_EST_DLPFC;
        elseif strcmp(subType , 'INS')
            load([myMotherPath ,'AE_EYE_STRUCT_EST_INS.mat']);
            EYE_MAT_R1 = AE_EYE_STRUCT_EST_INS;
            load([myMotherPath ,'AE_BEV_MAT_EST_INS.mat']);
            BEV_MAT_R1 = AE_BEV_MAT_EST_INS;
        end
    elseif task2Anal == 2
        if strcmp(subType , 'CONTROL')
            load([myMotherPath ,'IREM_EYE_STRUCT_BRIGHTNESS_CONTROL.mat']);
            EYE_MAT_R1 = IREM_EYE_STRUCT_BRIGHTNESS_CONTROL;
            load([myMotherPath ,'IREM_BEV_MAT_BRIGHTNESS_CONTROL.mat']);
            BEV_MAT_R1 = IREM_BEV_MAT_BRIGHTNESS_CONTROL;
        elseif strcmp(subType , 'TREAT')
            load([myMotherPath ,'IREM_EYE_STRUCT_BRIGHTNESS_TREAT.mat']);
            EYE_MAT_R1 = IREM_EYE_STRUCT_BRIGHTNESS_TREAT;
            load([myMotherPath ,'IREM_BEV_MAT_BRIGHTNESS_TREAT.mat']);
            BEV_MAT_R1 = IREM_BEV_MAT_BRIGHTNESS_TREAT;
        elseif strcmp(subType , 'OFC')
            load([myMotherPath ,'IREM_EYE_STRUCT_BRIGHTNESS_OFC.mat']);
            EYE_MAT_R1 = IREM_EYE_STRUCT_BRIGHTNESS_OFC;
            load([myMotherPath ,'IREM_BEV_MAT_BRIGHTNESS_OFC.mat']);
            BEV_MAT_R1 = IREM_BEV_MAT_BRIGHTNESS_OFC;
        elseif strcmp(subType , 'DLPFC')
            load([myMotherPath ,'IREM_EYE_STRUCT_BRIGHTNESS_DLPFC.mat']);
            EYE_MAT_R1 = IREM_EYE_STRUCT_BRIGHTNESS_DLPFC;
            load([myMotherPath ,'IREM_BEV_MAT_BRIGHTNESS_DLPFC.mat']);
            BEV_MAT_R1 = IREM_BEV_MAT_BRIGHTNESS_DLPFC;
        elseif strcmp(subType , 'INS')
            load([myMotherPath ,'IREM_EYE_STRUCT_BRIGHTNESS_INS.mat']);
            EYE_MAT_R1 = IREM_EYE_STRUCT_BRIGHTNESS_INS;
            load([myMotherPath ,'IREM_BEV_MAT_BRIGHTNESS_INS.mat']);
            BEV_MAT_R1 = IREM_BEV_MAT_BRIGHTNESS_INS;
        end
    elseif task2Anal == 3
        if strcmp(subType , 'CONTROL')
            load([myMotherPath ,'IREM_EYE_STRUCT_ODDBALL_CONTROL.mat']);
            EYE_MAT_R1 = IREM_EYE_STRUCT_ODDBALL_CONTROL;
            load([myMotherPath ,'IREM_BEV_MAT_ODDBALL_CONTROL.mat']);
            BEV_MAT_R1 = IREM_BEV_MAT_ODDBALL_CONTROL;
        elseif strcmp(subType , 'TREAT')
            load([myMotherPath ,'IREM_EYE_STRUCT_ODDBALL_TREAT.mat']);
            EYE_MAT_R1 = IREM_EYE_STRUCT_ODDBALL_TREAT;
            load([myMotherPath ,'IREM_BEV_MAT_ODDBALL_TREAT.mat']);
            BEV_MAT_R1 = IREM_BEV_MAT_ODDBALL_TREAT;
        elseif strcmp(subType , 'OFC')
            load([myMotherPath ,'IREM_EYE_STRUCT_ODDBALL_OFC.mat']);
            EYE_MAT_R1 = IREM_EYE_STRUCT_ODDBALL_OFC;
            load([myMotherPath ,'IREM_BEV_MAT_ODDBALL_OFC.mat']);
            BEV_MAT_R1 = IREM_BEV_MAT_ODDBALL_OFC;
        elseif strcmp(subType , 'DLPFC')
            load([myMotherPath ,'IREM_EYE_STRUCT_ODDBALL_DLPFC.mat']);
            EYE_MAT_R1 = IREM_EYE_STRUCT_ODDBALL_DLPFC;
            load([myMotherPath ,'IREM_BEV_MAT_ODDBALL_DLPFC.mat']);
            BEV_MAT_R1 = IREM_BEV_MAT_ODDBALL_DLPFC;
        elseif strcmp(subType , 'INS')
            load([myMotherPath ,'IREM_EYE_STRUCT_ODDBALL_INS.mat']);
            EYE_MAT_R1 = IREM_EYE_STRUCT_ODDBALL_INS;
            load([myMotherPath ,'IREM_BEV_MAT_ODDBALL_INS.mat']);
            BEV_MAT_R1 = IREM_BEV_MAT_ODDBALL_INS;
        end
    else
        disp('Task not Specified. Please define task2Anal variable --> 1=EST; 2=Brightness; 3=Oddball;')
    end
 BEV_MAT_R1 = BEV_MAT_RAW;
 EYE_MAT_R1 = EYE_MAT_RAW;
motherpath = getenv('THREEPAD_EST_ROOT');
mkdir([motherpath, filesep, 'Results_EST', filesep, 'EYE_ANAL', filesep, 'RUN1']);
    if task2Anal == 1
        allSubRTs_R1 = [];
        allSubACCs_R1 = [];
        theSubs = unique(BEV_MAT_R1(:,17));
        for currSubID = 1 : length(theSubs)
            currSub = theSubs(currSubID);
            currSubBevMat_R1 = BEV_MAT_R1(BEV_MAT_R1(:,17) == currSub, :);
            currSubBevMat_R1(currSubBevMat_R1(:,2)==99,:) = [];
            currSubBevMatCAOnly_R1 = currSubBevMat_R1;
            currSubBevMatCAOnly_R1(currSubBevMatCAOnly_R1(:,2)==0,:) = [];
            C_ID = currSubBevMatCAOnly_R1(:,3) == 1;
            C_RTs = currSubBevMatCAOnly_R1(C_ID,5);
            C_ID = currSubBevMat_R1(:,3) == 1;
            C_CAs = currSubBevMat_R1(C_ID,2);
            I_ID = currSubBevMatCAOnly_R1(:,3) == 2;
            I_RTs = currSubBevMatCAOnly_R1(I_ID,5);
            I_ID = currSubBevMat_R1(:,3) == 2;
            I_CAs = currSubBevMat_R1(I_ID,2);
            CC_ID = currSubBevMatCAOnly_R1(:,16) == 1;
            CC_RTs = currSubBevMatCAOnly_R1(CC_ID,5);
            CC_ID = currSubBevMat_R1(:,16) == 1;
            CC_CAs = currSubBevMat_R1(CC_ID,2);
            IC_ID = currSubBevMatCAOnly_R1(:,16) == 2;
            IC_RTs = currSubBevMatCAOnly_R1(IC_ID,5);
            IC_ID = currSubBevMat_R1(:,16) == 2;
            IC_CAs = currSubBevMat_R1(IC_ID,2);
            CI_ID = currSubBevMatCAOnly_R1(:,16) == 3;
            CI_RTs = currSubBevMatCAOnly_R1(CI_ID,5);
            CI_ID = currSubBevMat_R1(:,16) == 3;
            CI_CAs = currSubBevMat_R1(CI_ID,2);
            II_ID = currSubBevMatCAOnly_R1(:,16) == 4;
            II_RTs = currSubBevMatCAOnly_R1(II_ID,5);
            II_ID = currSubBevMat_R1(:,16) == 4;
            II_CAs = currSubBevMat_R1(II_ID,2);
            wordH_ID = currSubBevMatCAOnly_R1(:,6) == 1;
            wordH_RTs = currSubBevMatCAOnly_R1(wordH_ID,5);
            wordH_ID = currSubBevMat_R1(:,6) == 1;
            wordH_CAs = currSubBevMat_R1(wordH_ID,2);
            wordF_ID = currSubBevMatCAOnly_R1(:,6) == 2;
            wordF_RTs = currSubBevMatCAOnly_R1(wordF_ID,5);
            wordF_ID = currSubBevMat_R1(:,6) == 2;
            wordF_CAs = currSubBevMat_R1(wordF_ID,2);
            faceH_ID = currSubBevMatCAOnly_R1(:,7) == 1;
            faceH_RTs = currSubBevMatCAOnly_R1(faceH_ID,5);
            faceH_ID = currSubBevMat_R1(:,7) == 1;
            faceH_CAs = currSubBevMat_R1(faceH_ID,2);
            faceF_ID = currSubBevMatCAOnly_R1(:,7) == 2;
            faceF_RTs = currSubBevMatCAOnly_R1(faceF_ID,5);
            faceF_ID = currSubBevMat_R1(:,7) == 2;
            faceF_CAs = currSubBevMat_R1(faceF_ID,2);
            allSubRTs_R1 = cat(1, allSubRTs_R1, [mean(C_RTs), mean(I_RTs), mean(CC_RTs), mean(IC_RTs), mean(CI_RTs), mean(II_RTs), mean(wordH_RTs), mean(wordF_RTs), mean(faceH_RTs), mean(faceF_RTs)]);
            allSubACCs_R1 = cat(1, allSubACCs_R1, [sum(C_CAs)/length(C_CAs), sum(I_CAs)/length(I_CAs), sum(CC_CAs)/length(CC_CAs), sum(IC_CAs)/length(IC_CAs), sum(CI_CAs)/length(CI_CAs), sum(II_CAs)/length(II_CAs), sum(wordH_CAs)/length(wordH_CAs), sum(wordF_CAs)/length(wordF_CAs), sum(faceH_CAs)/length(faceH_CAs), sum(faceF_CAs)/length(faceF_CAs)]);
        end
        allSubsRTs_MeanR1 = mean(allSubRTs_R1, 1);
        allSubsRTs_SEMR1 = std(allSubRTs_R1, 1)/sqrt(size(allSubRTs_R1, 1));
        allSubsACCs_MeanR1 = mean(allSubACCs_R1, 1);
        allSubsACCs_SEMR1 = std(allSubACCs_R1, 1)/sqrt(size(allSubACCs_R1, 1));
        fh = figure;
        bar(1, allSubsRTs_MeanR1(1),'FaceColor', myCongruentColor,'EdgeColor',[1 1 1],'LineWidth',1, 'BarWidth', 0.5); hold on;
        bar(2, allSubsRTs_MeanR1(2),'FaceColor',myIncongruentColor,'EdgeColor',[1 1 1],'LineWidth',1, 'BarWidth', 0.5); hold on;
        errorbar([1,2], allSubsRTs_MeanR1(1:2), allSubsRTs_SEMR1(1:2), 'linestyle', 'none', 'color', 'k'); hold on;
        set(gca, 'xtick', [1,2] ,'xticklabel', {'No-Conflict', 'Conflict'}, 'fontsize', 17)
        xlim([.5, 2.5]);
        ylim([500, 1000]);
        title( 'Reaction time');
        box off;
        axis square
        if saveFigures
            myHome = pwd;
            cd([myMotherPath ,filesep, fullfile('Results_EST', 'EYE_ANAL')]);
            set(gcf, 'paperunits', 'centimeters', 'paperposition', [0 0 8 8])
            print(['BEV_Results_RT_Task_',num2str(task2Anal),'_subType_', subType, '_IvsC' ,'.tif'], '-dtiff', '-r300')
            cd(myHome)
        end
        fh = figure;
        bar(1, allSubsRTs_MeanR1(3),'FaceColor', myCCColor,'EdgeColor',[1 1 1],'LineWidth',1, 'BarWidth', 0.5); hold on;
        bar(2, allSubsRTs_MeanR1(4),'FaceColor',myICColor,'EdgeColor',[1 1 1],'LineWidth',1, 'BarWidth', 0.5); hold on;
        bar(3, allSubsRTs_MeanR1(5),'FaceColor', myCIColor,'EdgeColor',[1 1 1],'LineWidth',1, 'BarWidth', 0.5); hold on;
        bar(4, allSubsRTs_MeanR1(6),'FaceColor',myIIColor,'EdgeColor',[1 1 1],'LineWidth',1, 'BarWidth', 0.5); hold on;
        errorbar([1:4], allSubsRTs_MeanR1(3:6), allSubsRTs_SEMR1(3:6), 'linestyle', 'none', 'color', 'k'); hold on;
        set(gca, 'xtick', [1:4] ,'xticklabel', {'CC', 'IC', 'CI', 'II'}, 'fontsize', 17)
        xlim([.5, 4.5]);
        ylim([500, 1000]);
        title( 'Reaction time');
        box off;
        axis square
        if saveFigures
            myHome = pwd;
            cd([myMotherPath ,filesep, fullfile('Results_EST', 'EYE_ANAL', 'RUN1')]);
            set(gcf, 'paperunits', 'centimeters', 'paperposition', [0 0 8 8])
            print(['BEV_Results_RT_Task_',num2str(task2Anal),'_subType_', subType, '_CIvsII' ,'.tif'], '-dtiff', '-r300')
            cd(myHome)
        end
        fh = figure;
        bar(1, allSubsACCs_MeanR1(1),'FaceColor', myCongruentColor,'EdgeColor',[1 1 1],'LineWidth',1, 'BarWidth', 0.5); hold on;
        bar(2, allSubsACCs_MeanR1(2),'FaceColor',myIncongruentColor,'EdgeColor',[1 1 1],'LineWidth',1, 'BarWidth', 0.5); hold on;
        errorbar([1,2], allSubsACCs_MeanR1(1:2), allSubsACCs_SEMR1(1:2), 'linestyle', 'none', 'color', 'k'); hold on;
        set(gca, 'xtick', [1,2] ,'xticklabel', {'No-Conflict', 'Conflict'}, 'fontsize', 17)
        xlim([.5, 2.5]);
        ylim([0.9, 1.1]);
        title( 'Accuracy');
        box off;
        axis square
        if saveFigures
            myHome = pwd;
            cd([myMotherPath ,filesep, fullfile('Results_EST', 'EYE_ANAL', 'RUN1')]);
            set(gcf, 'paperunits', 'centimeters', 'paperposition', [0 0 8 8])
            print(['BEV_Results_ACC_Task_',num2str(task2Anal),'_subType_', subType, '_IvsC' ,'.tif'], '-dtiff', '-r300')
            cd(myHome)
        end
        fh = figure ;
        bar(1, allSubsACCs_MeanR1(3),'FaceColor', myCCColor,'EdgeColor',[1 1 1],'LineWidth',1, 'BarWidth', 0.5); hold on;
        bar(2, allSubsACCs_MeanR1(4),'FaceColor',myICColor,'EdgeColor',[1 1 1],'LineWidth',1, 'BarWidth', 0.5); hold on;
        bar(3, allSubsACCs_MeanR1(5),'FaceColor', myCIColor,'EdgeColor',[1 1 1],'LineWidth',1, 'BarWidth', 0.5); hold on;
        bar(4, allSubsACCs_MeanR1(6),'FaceColor',myIIColor,'EdgeColor',[1 1 1],'LineWidth',1, 'BarWidth', 0.5); hold on;
        errorbar([1:4], allSubsACCs_MeanR1(3:6), allSubsACCs_SEMR1(3:6), 'linestyle', 'none', 'color', 'k'); hold on;
        set(gca, 'xtick', [1:4] ,'xticklabel', {'CC', 'IC', 'CI', 'II'}, 'fontsize', 17)
        xlim([.5, 4.5]);
        ylim([0.9, 1.1]);
        title( 'Accuracy');
        box off;
        axis square
        if saveFigures
            myHome = pwd;
            cd([myMotherPath , filesep, fullfile('Results_EST', 'EYE_ANAL', 'RUN1')]);
            set(gcf, 'paperunits', 'centimeters', 'paperposition', [0 0 8 8])
            print(['BEV_Results_ACC_Task_',num2str(task2Anal),'_subType_', subType, '_CIvsII' ,'.tif'], '-dtiff', '-r300')
            cd(myHome)
        end
    end
    msBelow = 10000;
    msAbove = 10000;
    theSubs = unique(BEV_MAT_R1(:,17));
    allSubsPupilStimLocked_AL = nan(length(theSubs), length(-msBelow:msAbove));
    allSubsPupilRespLocked_AL = nan(length(theSubs), length(-msBelow:msAbove));
    allSubsPupilStimLocked_CC = nan(length(theSubs), length(-msBelow:msAbove));
    allSubsPupilRespLocked_CC = nan(length(theSubs), length(-msBelow:msAbove));
    allSubsPupilStimLocked_IC = nan(length(theSubs), length(-msBelow:msAbove));
    allSubsPupilRespLocked_IC = nan(length(theSubs), length(-msBelow:msAbove));
    allSubsPupilStimLocked_CI = nan(length(theSubs), length(-msBelow:msAbove));
    allSubsPupilRespLocked_CI = nan(length(theSubs), length(-msBelow:msAbove));
    allSubsPupilStimLocked_II = nan(length(theSubs), length(-msBelow:msAbove));
    allSubsPupilRespLocked_II = nan(length(theSubs), length(-msBelow:msAbove));
    allSubsPupilStimLocked_Cong = nan(length(theSubs), length(-msBelow:msAbove));
    allSubsPupilRespLocked_Cong = nan(length(theSubs), length(-msBelow:msAbove));
    allSubsPupilStimLocked_Incong = nan(length(theSubs), length(-msBelow:msAbove));
    allSubsPupilRespLocked_Incong = nan(length(theSubs), length(-msBelow:msAbove));
    allSubsGazeBev = [];
    allSubsAllBev = [];
    allSubsAllPup = [];
    allSubUpDownSequenceMat = [];
    allSubUpDownDurMat      = [];
    allSubUpDownOnsetMat    = [];
    for currSubID = 1 : length(theSubs)
        currSub = theSubs(currSubID);
        currSubBevMat_R1 = BEV_MAT_R1(BEV_MAT_R1(:,17) == currSub, :);
        theRuns = unique(currSubBevMat_R1(:,18));
        currSubPupilStimLocked_AL = nan(length(theRuns), length(-msBelow:msAbove));
        currSubPupilRespLocked_AL = nan(length(theRuns), length(-msBelow:msAbove));
        currSubPupilStimLocked_CC = nan(length(theRuns), length(-msBelow:msAbove));
        currSubPupilRespLocked_CC = nan(length(theRuns), length(-msBelow:msAbove));
        currSubPupilStimLocked_IC = nan(length(theRuns), length(-msBelow:msAbove));
        currSubPupilRespLocked_IC = nan(length(theRuns), length(-msBelow:msAbove));
        currSubPupilStimLocked_CI = nan(length(theRuns), length(-msBelow:msAbove));
        currSubPupilRespLocked_CI = nan(length(theRuns), length(-msBelow:msAbove));
        currSubPupilStimLocked_II = nan(length(theRuns), length(-msBelow:msAbove));
        currSubPupilRespLocked_II = nan(length(theRuns), length(-msBelow:msAbove));
        currSubGazeBev = [];
        currSubAllBev = [];
        currSubAllPup = [];
        currSubUpDownSequenceMat= [];
        currSubUpDownDurMat     = [];
        currSubUpDownOnsetMat   = [];
        for currRunID = 1:length(theRuns)
            currRun = theRuns(currRunID);
            disp(['Working on SUB: ', num2str(currSub), ' Run: ', num2str(currRun)])
            currRunPupilSizeVect = EYE_MAT_R1(currSub).sub(currRun).run(currRun).currP;
            currRunPupilTimeVect = EYE_MAT_R1(currSub).sub(currRun).run(currRun).currT;
            currRunPupilSizeVectMean = mean(currRunPupilSizeVect);
            srate = 1/mean(diff(currRunPupilTimeVect/1000));
            if doLowPassFilteringOnly
                disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
                disp(['LOW-PASS-Filtering Sub: ', num2str(currSub), ' Run: ', num2str(currRun) ]);
                disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
                [currRunPupilSizeVectLowPassFiltered] = eegfilt(currRunPupilSizeVect', srate, 0, 4);
                currP = [nan(msBelow,1); currRunPupilSizeVectLowPassFiltered'; nan(msAbove,1)];
            elseif doBandPassFiltering
                disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
                disp(['BAND-PASS-Filtering Sub: ', num2str(currSub), ' Run: ', num2str(currRun) ]);
                disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
                [currRunPupilSizeVectHighPassFiltered] = eegfilt(currRunPupilSizeVect', srate, 0.05, 0);
                [currRunPupilSizeVectBandPassFiltered] = eegfilt(currRunPupilSizeVectHighPassFiltered, srate, 0, 4);
                currP = [nan(msBelow,1); currRunPupilSizeVectBandPassFiltered'; nan(msAbove,1)];
            else
                currP = [nan(msBelow,1); EYE_MAT_R1(currSub).sub(currRun).run(currRun).currP; nan(msAbove,1)];
            end
            if zscorePupil; currP = (currP - nanmean(currP))/nanstd(currP); end
            currPupilT = [nan(msBelow,1); EYE_MAT_R1(currSub).sub(currRun).run(currRun).currT; nan(msAbove,1)];
            currX = EYE_MAT_R1(currSub).sub(currRun).run(currRun).currX;
            currY = EYE_MAT_R1(currSub).sub(currRun).run(currRun).currY;
            currT = EYE_MAT_R1(currSub).sub(currRun).run(currRun).currT;
            currRunCenter = EYE_MAT_R1(currSub).sub(currRun).run(currRun).center;
            trialOnsetColumn = 8; RTColumn = 5;
            currBev = currSubBevMat_R1(currSubBevMat_R1(:,18) == currRun, :);
            currBev(:,trialOnsetColumn) = (round(currBev(:,trialOnsetColumn)/10))*10;
            currBev(:,RTColumn) = (round(currBev(:,RTColumn)/10))*10 ;
            currRunPupilStimLocked_AL = nan(size(currBev, 1),  length(-msBelow:msAbove));
            currRunPupilRespLocked_AL = nan(size(currBev, 1),  length(-msBelow:msAbove));
            currRunPupilStimLocked_CC = nan(size(currBev, 1),  length(-msBelow:msAbove));
            currRunPupilRespLocked_CC = nan(size(currBev, 1),  length(-msBelow:msAbove));
            currRunPupilStimLocked_IC = nan(size(currBev, 1),  length(-msBelow:msAbove));
            currRunPupilRespLocked_IC = nan(size(currBev, 1),  length(-msBelow:msAbove));
            currRunPupilStimLocked_CI = nan(size(currBev, 1),  length(-msBelow:msAbove));
            currRunPupilRespLocked_CI = nan(size(currBev, 1),  length(-msBelow:msAbove));
            currRunPupilStimLocked_II = nan(size(currBev, 1),  length(-msBelow:msAbove));
            currRunPupilRespLocked_II = nan(size(currBev, 1),  length(-msBelow:msAbove));
            currGazeBev = nan(size(currBev, 1), 8);
            currPupBev =  nan(size(currBev, 1), 3);
            allInfoCurrRun =  nan(size(currBev, 1), 52);
            currUpDownSequenceMat   = nan(size(currBev, 1), 9);
            currUpDownDurMat        = nan(size(currBev, 1), 9);
            currUpDownOnsetMat      = nan(size(currBev, 1), 9);
            for currTrial = 1 : size(currBev, 1)
                curOn = currBev(currTrial,trialOnsetColumn) - currBev(currTrial,15);
                curRT = currBev(currTrial,RTColumn);
                if currBev(currTrial,16) == 1
                    ccOn = currBev(currTrial,trialOnsetColumn) - currBev(currTrial,15);
                    ccRT = currBev(currTrial,RTColumn);
                elseif currBev(currTrial,16) == 2
                    icOn = currBev(currTrial,trialOnsetColumn) - currBev(currTrial,15);
                    icRT = currBev(currTrial,RTColumn);
                elseif currBev(currTrial,16) == 3
                    ciOn = currBev(currTrial,trialOnsetColumn) - currBev(currTrial,15);
                    ciRT = currBev(currTrial,RTColumn);
                elseif currBev(currTrial,16) == 4
                    iiOn = currBev(currTrial,trialOnsetColumn) - currBev(currTrial,15);
                    iiRT = currBev(currTrial,4);
                else
                    if currBev(currTrial,3) == 1
                        firstTrial_conOn = currBev(currTrial,trialOnsetColumn) - currBev(currTrial,15);
                        firstTrial_conRT = currBev(currTrial,RTColumn);
                    elseif currBev(currTrial,3) == 2
                        firstTrial_inconOn = currBev(currTrial,trialOnsetColumn) - currBev(currTrial,15);
                        firstTrial_inconRT = currBev(currTrial,RTColumn);
                    else
                        firstTrial_conOn = currBev(currTrial,trialOnsetColumn) - currBev(currTrial,15);
                        firstTrial_conRT = currBev(currTrial,RTColumn);
                    end
                end
                stimOnID = find(currPupilT == curOn);
                currRunPupilStimLocked_AL(currTrial,:) = currP(stimOnID-msBelow : stimOnID+msAbove)';
                if ~isnan(curRT)
                    curRTID = find(currPupilT == curOn + curRT);
                    currRunPupilRespLocked_AL(currTrial,:) = currP(curRTID-msBelow : curRTID+msAbove)';
                end
                if currBev(currTrial,16) == 1
                    stimOnID = find(currPupilT == ccOn);
                    currRunPupilStimLocked_CC(currTrial,:) = currP(stimOnID-msBelow : stimOnID+msAbove)';
                    PBL = mean(currP(stimOnID-500 : stimOnID));
                    currPupBev(currTrial,1) = PBL;
                    if ~isnan(ccRT)
                        ccRTID = find(currPupilT == ccOn + ccRT);
                        currRunPupilRespLocked_CC(currTrial,:) = currP(ccRTID-msBelow : ccRTID+msAbove)';
                        PAR = mean(currP(ccRTID-250 : ccRTID+250));
                        currPupBev(currTrial,2) = PAR;
                        ROD = sum(diff(currP(stimOnID : ccRTID)));
                        currPupBev(currTrial,3) = ROD;
                    end
                    PCU = mean(currP(stimOnID+500 : stimOnID+2000));
                    currPupBev(currTrial,4) = PCU;
                    PTR = mean(currP(stimOnID+500 : stimOnID+850));
                    currPupBev(currTrial,5) = PTR;
                    PAA = mean(currP(stimOnID+2128 : stimOnID+4138));
                    currPupBev(currTrial,6) = PAA;
                    stimOn = ccOn;
                    choiceRT = ccRT;
                elseif currBev(currTrial,16) == 2
                    stimOnID = find(currPupilT == icOn);
                    currRunPupilStimLocked_IC(currTrial,:) = currP(stimOnID-msBelow : stimOnID+msAbove)';
                    PBL = mean(currP(stimOnID-500 : stimOnID));
                    currPupBev(currTrial,1) = PBL;
                    if ~isnan(icRT)
                        icRTID = find(currPupilT == icOn + icRT);
                        currRunPupilRespLocked_IC(currTrial,:) = currP(icRTID-msBelow : icRTID+msAbove)';
                        PAR = mean(currP(icRTID-250 : icRTID+250));
                        currPupBev(currTrial,2) = PAR;
                        ROD = sum(diff(currP(stimOnID : icRTID)));
                        currPupBev(currTrial,3) = ROD;
                    end
                    PCU = mean(currP(stimOnID+500 : stimOnID+2000));
                    currPupBev(currTrial,4) = PCU;
                    PTR = mean(currP(stimOnID+500 : stimOnID+850));
                    currPupBev(currTrial,5) = PTR;
                    PAA = mean(currP(stimOnID+2128 : stimOnID+4138));
                    stimOn = icOn;
                    choiceRT = icRT;
                elseif currBev(currTrial,16) == 3
                    stimOnID = find(currPupilT == ciOn);
                    currRunPupilStimLocked_CI(currTrial,:) = currP(stimOnID-msBelow : stimOnID+msAbove)';
                    PBL = mean(currP(stimOnID-500 : stimOnID));
                    currPupBev(currTrial,1) = PBL;
                    if ~isnan(ciRT)
                        ciRTID = find(currPupilT == ciOn + ciRT);
                        currRunPupilRespLocked_CI(currTrial,:) = currP(ciRTID-msBelow : ciRTID+msAbove)';
                        PAR = mean(currP(ciRTID-250 : ciRTID+250));
                        currPupBev(currTrial,2) = PAR;
                        ROD = sum(diff(currP(stimOnID : ciRTID)));
                        currPupBev(currTrial,3) = ROD;
                    end
                    PCU = mean(currP(stimOnID+500 : stimOnID+2000));
                    currPupBev(currTrial,4) = PCU;
                    PTR = mean(currP(stimOnID+500 : stimOnID+850));
                    currPupBev(currTrial,5) = PTR;
                    PAA = mean(currP(stimOnID+2128 : stimOnID+4138));
                    currPupBev(currTrial,6) = PAA;
                    stimOn = ciOn;
                    choiceRT = ciRT;
                elseif currBev(currTrial,16) == 4
                    stimOnID = find(currPupilT == iiOn);
                    currRunPupilStimLocked_II(currTrial,:) = currP(stimOnID-msBelow : stimOnID+msAbove)';
                    PBL = mean(currP(stimOnID-500 : stimOnID));
                    currPupBev(currTrial,1) = PBL;
                    if ~isnan(iiRT)
                        iiRTID = find(currPupilT == iiOn + iiRT);
                        currRunPupilRespLocked_II(currTrial,:) = currP(iiRTID-msBelow : iiRTID+msAbove)';
                        PAR = mean(currP(iiRTID-250 : iiRTID+250));
                        currPupBev(currTrial,2) = PAR;
                        ROD = sum(diff(currP(stimOnID : iiRTID)));
                        currPupBev(currTrial,3) = ROD;
                    end
                    PCU = mean(currP(stimOnID+500 : stimOnID+2000));
                    currPupBev(currTrial,4) = PCU;
                    PTR = mean(currP(stimOnID+500 : stimOnID+850));
                    currPupBev(currTrial,5) = PTR;
                    PAA = mean(currP(stimOnID+2128 : stimOnID+4138));
                    currPupBev(currTrial,6) = PAA;
                    stimOn = iiOn;
                    choiceRT = iiRT;
                elseif currBev(currTrial,16) == 0
                    if currBev(currTrial,3) == 1
                        stimOnID = find(currPupilT == firstTrial_conOn);
                        currRunPupilStimLocked_FTCON(currTrial,:) = currP(stimOnID-msBelow : stimOnID+msAbove)';
                        PBL = mean(currP(stimOnID-500 : stimOnID));
                        currPupBev(currTrial,1) = PBL;
                        if ~isnan(firstTrial_conRT)
                            firstTrial_conRTID = find(currPupilT == firstTrial_conOn + firstTrial_conRT);
                            currRunPupilRespLocked_FTCON(currTrial,:) = currP(firstTrial_conRTID-msBelow : firstTrial_conRTID+msAbove)';
                            PAR = mean(currP(firstTrial_conRTID-250 : firstTrial_conRTID+250));
                            currPupBev(currTrial,2) = PAR;
                            ROD = sum(diff(currP(stimOnID : firstTrial_conRTID)));
                            currPupBev(currTrial,3) = ROD;
                        end
                        PCU = mean(currP(stimOnID+500 : stimOnID+2000));
                        currPupBev(currTrial,4) = PCU;
                        PTR = mean(currP(stimOnID+500 : stimOnID+850));
                        currPupBev(currTrial,5) = PTR;
                        PAA = mean(currP(stimOnID+2128 : stimOnID+4138));
                        currPupBev(currTrial,6) = PAA;
                        stimOn = firstTrial_conOn;
                        choiceRT = firstTrial_conRT;
                    elseif currBev(currTrial,3) == 2
                        stimOnID = find(currPupilT == firstTrial_inconOn);
                        currRunPupilStimLocked_FTCON(currTrial,:) = currP(stimOnID-msBelow : stimOnID+msAbove)';
                        PBL = mean(currP(stimOnID-500 : stimOnID));
                        currPupBev(currTrial,1) = PBL;
                        if ~isnan(firstTrial_inconRT)
                            firstTrial_inconRTID = find(currPupilT == firstTrial_inconOn + firstTrial_inconRT);
                            currRunPupilRespLocked_FTCON(currTrial,:) = currP(firstTrial_inconRTID-msBelow : firstTrial_inconRTID+msAbove)';
                            PAR = mean(currP(firstTrial_inconRTID-250 : firstTrial_inconRTID+250));
                            currPupBev(currTrial,2) = PAR;
                            ROD = sum(diff(currP(stimOnID : firstTrial_inconRTID)));
                            currPupBev(currTrial,3) = ROD;
                        end
                        PCU = mean(currP(stimOnID+500 : stimOnID+2000));
                        currPupBev(currTrial,4) = PCU;
                        PTR = mean(currP(stimOnID+500 : stimOnID+850));
                        currPupBev(currTrial,5) = PTR;
                        PAA = mean(currP(stimOnID+2128 : stimOnID+4138));
                        currPupBev(currTrial,6) = PAA;
                        stimOn = firstTrial_inconOn;
                        choiceRT = firstTrial_inconRT;
                    else
                        stimOnID = find(currPupilT == firstTrial_conOn);
                        currRunPupilStimLocked_FTCON(currTrial,:) = currP(stimOnID-msBelow : stimOnID+msAbove)';
                        PBL = mean(currP(stimOnID-500 : stimOnID));
                        currPupBev(currTrial,1) = PBL;
                        if ~isnan(firstTrial_conRT)
                            firstTrial_conRTID = find(currPupilT == firstTrial_conOn + firstTrial_conRT);
                            currRunPupilRespLocked_FTCON(currTrial,:) = currP(firstTrial_conRTID-msBelow : firstTrial_conRTID+msAbove)';
                            PAR = mean(currP(firstTrial_conRTID-250 : firstTrial_conRTID+250));
                            currPupBev(currTrial,2) = PAR;
                            ROD = sum(diff(currP(stimOnID : firstTrial_conRTID)));
                            currPupBev(currTrial,3) = ROD;
                        end
                        PCU = mean(currP(stimOnID : stimOnID+1500));
                        currPupBev(currTrial,4) = PCU;
                        PTR = mean(currP(stimOnID+500 : stimOnID+850));
                        currPupBev(currTrial,5) = PTR;
                        PAA = mean(currP(stimOnID+2128 : stimOnID+4138));
                        currPupBev(currTrial,6) = PAA;
                        stimOn = firstTrial_conOn;
                        choiceRT = firstTrial_conRT;
                    end
                end
                fixZoneHalfWidth = 60;
                if ~isnan(choiceRT)
                    data = currY(stimOn:stimOn+choiceRT) - currRunCenter(2);
                    aboveBoundID = data > fixZoneHalfWidth;
                    belowBoundID = data < -fixZoneHalfWidth;
                    if plotSingleTrialYGazeChoice
                        fh = figure;
                        plot(data); hold on;
                        line([0, length(stimOn:stimOn+choiceRT)], [fixZoneHalfWidth fixZoneHalfWidth], 'color', 'k', 'linestyle', ':')
                        line([0, length(stimOn:stimOn+choiceRT)], [0 0], 'color', 'k', 'linestyle', '-.')
                        line([0, length(stimOn:stimOn+choiceRT)], [-fixZoneHalfWidth -fixZoneHalfWidth], 'color', 'k', 'linestyle', ':')
                        line([100, 100], [min(data) max(data)], 'color', 'k', 'linestyle', ':');
                        line([200, 200], [min(data) max(data)], 'color', 'k', 'linestyle', ':');
                        line([300, 300], [min(data) max(data)], 'color', 'k', 'linestyle', ':');
                        title(['currSubID: ', num2str(currSub), '; currRun: ', num2str(currRun), '; currTrial: ', num2str(currTrial)]);
                        set(fh, 'position', [1950 500 1000 500], 'name', 'current Trial - Y-Fixation traces');
                    end
                    currStartsID = []; currEndsID = [];
                    diffAboveBoundID = diff(aboveBoundID);
                    currStartsID = find(diffAboveBoundID==1); currEndsID = find(diffAboveBoundID==-1);
                    if ~isempty(currStartsID) || ~isempty(currEndsID)
                        if isempty(currEndsID); currEndsID(1) = nan; end
                        if isreal(currEndsID(1)) && isempty(currStartsID); currStartsID(1) = 0; end
                        if currEndsID(1) < currStartsID(1)
                            SandE = nan(size(currStartsID,1)+1,2);
                            SandE(1:length(currStartsID)+1,1) = [nan; currStartsID]; SandE(1:length(currEndsID),2) = currEndsID;
                        else
                            SandE = nan(size(currStartsID,1),2);
                            SandE(1:length(currStartsID),1) = currStartsID; SandE(1:length(currEndsID),2) = currEndsID;
                        end
                        if isnan(SandE(1,1)); SandE(1,1) = 0; end
                        if isnan(SandE(end,end)); SandE(end,end) = choiceRT; end
                        properAboveFixTraceID = diff(SandE,1,2) > 50;
                        SandE_properAboveFixTrace = SandE(properAboveFixTraceID,:);
                    elseif isempty(currStartsID) && isempty(currEndsID)
                        SandE_properAboveFixTrace = [nan nan];
                    end
                    currStartsID = []; currEndsID = [];
                    diffBelowBoundID = diff(belowBoundID);
                    currStartsID = find(diffBelowBoundID==1); currEndsID = find(diffBelowBoundID==-1);
                    if ~isempty(currStartsID) || ~isempty(currEndsID)
                        if isempty(currEndsID); currEndsID(1) = nan; end
                        if isreal(currEndsID(1)) && isempty(currStartsID); currStartsID(1) = 0; end
                        if currEndsID(1) < currStartsID(1)
                            SandE = nan(size(currStartsID,1)+1,2);
                            SandE(1:length(currStartsID)+1,1) = [nan; currStartsID]; SandE(1:length(currEndsID),2) = currEndsID;
                        else
                            SandE = nan(size(currStartsID,1),2);
                            SandE(1:length(currStartsID),1) = currStartsID; SandE(1:length(currEndsID),2) = currEndsID;
                        end
                        if isnan(SandE(1,1)); SandE(1,1) = 0; end
                        if isnan(SandE(end,end)); SandE(end,end) = choiceRT; end
                        properBelowFixTraceID = diff(SandE,1,2) > 50;
                        SandE_properBelowFixTrace = SandE(properBelowFixTraceID,:);
                    elseif isempty(currStartsID) && isempty(currEndsID)
                        SandE_properBelowFixTrace = [nan nan];
                    end
                    if plotSingleTrialYGazeChoice
                        for currFix = 1 : size(SandE_properAboveFixTrace,1); line([SandE_properAboveFixTrace(currFix,1), SandE_properAboveFixTrace(currFix,2)], [fixZoneHalfWidth fixZoneHalfWidth], 'color', 'r', 'linewidth', 5); hold on; end
                        for currFix = 1 : size(SandE_properBelowFixTrace,1); line([SandE_properBelowFixTrace(currFix,1), SandE_properBelowFixTrace(currFix,2)], -[fixZoneHalfWidth fixZoneHalfWidth], 'color', 'g', 'linewidth', 5); hold on; end
                        fh = figure;
                        subplot(2,1,1); plot(diffAboveBoundID); hold on; title('Lower Item')
                        subplot(2,1,2); plot(diffBelowBoundID); hold on; title('Upper Item')
                        set(fh, 'position', [2950 500 800 500], 'name', 'current Trial - Y-Fixation start and end of lower and upper item fixation');
                    end
                    DTU = []; DTL = []; NS = []; FF = []; FFT = []; FL = []; FLT = [];
                    DTU = sum(SandE_properBelowFixTrace(:,2) - SandE_properBelowFixTrace(:,1));
                    DTL = sum(SandE_properAboveFixTrace(:,2) - SandE_properAboveFixTrace(:,1));
                    NS = size(SandE_properBelowFixTrace,1) + size(SandE_properAboveFixTrace,1);
                    currTrialFixSequenceData    = [];
                    upDownSequenceVect          = [];
                    upDownDurVect               = [];
                    upDownOnsetVect             = [];
                    if ~isempty(SandE_properAboveFixTrace) && ~isempty(SandE_properBelowFixTrace)
                        if SandE_properAboveFixTrace(1,1) < SandE_properBelowFixTrace(1,1)
                            FF = -1;
                            FFT = SandE_properAboveFixTrace(1, 2)-SandE_properAboveFixTrace(1, 1);
                        else
                            FF = 1;
                            FFT = SandE_properBelowFixTrace(1, 2)-SandE_properBelowFixTrace(1, 1);
                        end
                        currTrialFixSequenceData = sortrows([SandE_properBelowFixTrace(:,1), SandE_properBelowFixTrace(:,2), ones(size(SandE_properBelowFixTrace(:,1))); SandE_properAboveFixTrace(:,1), SandE_properAboveFixTrace(:,2), ones(size(SandE_properAboveFixTrace(:,1)))*-1], 1);
                        upDownSequenceVect = currTrialFixSequenceData(:,3)';
                        upDownDurVect = (currTrialFixSequenceData(:,2) - currTrialFixSequenceData(:,1))';
                        upDownOnsetVect = currTrialFixSequenceData(:,1)';
                    elseif isempty(SandE_properAboveFixTrace) && isempty(SandE_properBelowFixTrace)
                        FF = nan;
                        FFT = nan;
                        upDownSequenceVect = nan;
                        upDownDurVect = nan;
                        upDownOnsetVect = nan;
                    elseif isempty(SandE_properAboveFixTrace)
                        FF = 1;
                        FFT = SandE_properBelowFixTrace(1, 2)-SandE_properBelowFixTrace(1, 1);
                        currTrialFixSequenceData = sortrows([SandE_properBelowFixTrace(:,1), SandE_properBelowFixTrace(:,2), ones(size(SandE_properBelowFixTrace(:,1)))], 1);
                        upDownSequenceVect = currTrialFixSequenceData(:,3)';
                        upDownDurVect = (currTrialFixSequenceData(:,2) - currTrialFixSequenceData(:,1))';
                        upDownOnsetVect = currTrialFixSequenceData(:,1)';
                    elseif isempty(SandE_properBelowFixTrace)
                        FF = -1;
                        FFT = SandE_properAboveFixTrace(1, 2)-SandE_properAboveFixTrace(1, 1);
                        currTrialFixSequenceData = sortrows([SandE_properAboveFixTrace(:,1), SandE_properAboveFixTrace(:,2), ones(size(SandE_properAboveFixTrace(:,1)))], 1);
                        upDownSequenceVect = currTrialFixSequenceData(:,3)';
                        upDownDurVect = (currTrialFixSequenceData(:,2) - currTrialFixSequenceData(:,1))';
                        upDownOnsetVect = currTrialFixSequenceData(:,1)';
                    end
                    if ~isempty(SandE_properAboveFixTrace) && ~isempty(SandE_properBelowFixTrace)
                        if SandE_properAboveFixTrace(end,end) > SandE_properBelowFixTrace(end,end)
                            FL = -1;
                            FLT = SandE_properAboveFixTrace(end, 2)-SandE_properAboveFixTrace(end, 1);
                        else
                            FL = 1;
                            FLT = SandE_properBelowFixTrace(end, 2)-SandE_properBelowFixTrace(end, 1);
                        end
                    elseif isempty(SandE_properAboveFixTrace) && isempty(SandE_properBelowFixTrace)
                        FL = nan;
                        FLT = nan;
                    elseif isempty(SandE_properAboveFixTrace)
                        FL = 1;
                        FLT = SandE_properBelowFixTrace(end, 2)-SandE_properBelowFixTrace(end, 1);
                    elseif isempty(SandE_properBelowFixTrace)
                        FL = -1;
                        FLT = SandE_properAboveFixTrace(end, 2)-SandE_properAboveFixTrace(end, 1);
                    end
                    currTrialGazeResults = [DTU, DTL, NS, FF, FFT, FL, FLT, currSub];
                else
                    currTrialGazeResults= [nan, nan, nan, nan, nan, nan, nan, nan];
                    upDownSequenceVect  = nan;
                    upDownDurVect       = nan;
                    upDownOnsetVect     = nan;
                end
                currGazeBev(currTrial, :) = currTrialGazeResults;
                if ~isnan(upDownSequenceVect)
                    currUpDownSequenceMat(currTrial, 1:size(upDownSequenceVect,2))  = upDownSequenceVect;
                    currUpDownDurMat(currTrial, 1:size(upDownSequenceVect,2))       = upDownDurVect;
                    currUpDownOnsetMat(currTrial, 1:size(upDownSequenceVect,2))     = upDownOnsetVect;
                end
                allInfoCurrRun(currTrial, 1:32) = [currBev(currTrial,:), currPupBev(currTrial,:), currTrialGazeResults];
            end
            currSubPupilStimLocked_AL(currRun,:)  = nanmean(currRunPupilStimLocked_AL);
            currSubPupilStimLocked_CC(currRun,:)  = nanmean(currRunPupilStimLocked_CC);
            currSubPupilStimLocked_IC(currRun,:)  = nanmean(currRunPupilStimLocked_IC);
            currSubPupilStimLocked_CI(currRun,:)  = nanmean(currRunPupilStimLocked_CI);
            currSubPupilStimLocked_II(currRun,:)  = nanmean(currRunPupilStimLocked_II);
            currSubPupilRespLocked_AL(currRun,:)  = nanmean(currRunPupilRespLocked_AL);
            currSubPupilRespLocked_CC(currRun,:)  = nanmean(currRunPupilRespLocked_CC);
            currSubPupilRespLocked_IC(currRun,:)  = nanmean(currRunPupilRespLocked_IC);
            currSubPupilRespLocked_CI(currRun,:)  = nanmean(currRunPupilRespLocked_CI);
            currSubPupilRespLocked_II(currRun,:)  = nanmean(currRunPupilRespLocked_II);
            currSubGazeBev = cat(1, currSubGazeBev, currGazeBev);
            currSubAllBev = cat(1, currSubAllBev, allInfoCurrRun);
            currSubAllPup = cat(1, currSubAllPup, currRunPupilStimLocked_AL);
            try
                currSubUpDownSequenceMat = cat(1, currSubUpDownSequenceMat, currUpDownSequenceMat);
                currSubUpDownDurMat      = cat(1, currSubUpDownDurMat, currUpDownDurMat);
                currSubUpDownOnsetMat    = cat(1, currSubUpDownOnsetMat, currUpDownOnsetMat);
            catch
            end
        end
        if plotSingleSubsMeanPupil
            fh1 = figure(1); subplot(10,10,currSubID);
            data = currSubPupilStimLocked_CC;
            data(sum(isnan(data),2) >= length(data), :) = [];
            if size(data,1) > 1
                shadedErrorBar(1:length(nanmean(data)), nanmean(data), nansem(data), {'-', 'color', myPercColor, 'lineWidth', myLineWidth}, 1); hold on;
            else
                plot(1:length(data), data, '-', 'color', myPercColor, 'lineWidth', myLineWidth); hold on;
            end
            line([round(length(data)/2), round(length(data)/2)], [min(nanmean(data))-(max(nanmean(data))-min(nanmean(data)))*0.3, max(nanmean(data))+(max(nanmean(data))-min(nanmean(data)))*0.3], 'color', myPercStimOnColor )
            title(['currSubID: ', num2str(currSub)]);
            xlim(myPupilPlotXlimVals); ylim(myPupilPlotYlimVals); set(gca, 'XTickLabel', []);
            set(fh1, 'position', [50 100 1600 900], 'name', 'currSubPupilStimLocked_CC');
            fh1 = figure(1); subplot(10,10,currSubID);
            data = currSubPupilStimLocked_IC;
            data(sum(isnan(data),2) >= length(data), :) = [];
            if size(data,1) > 1
                shadedErrorBar(1:length(nanmean(data)), nanmean(data), nansem(data), {'-', 'color', myEconColor, 'lineWidth', myLineWidth}, 1); hold on;
            else
                plot(1:length(data), data, '-', 'color', myEconColor, 'lineWidth', myLineWidth); hold on;
            end
            line([round(length(data)/2), round(length(data)/2)], [min(nanmean(data))-(max(nanmean(data))-min(nanmean(data)))*0.3, max(nanmean(data))+(max(nanmean(data))-min(nanmean(data)))*0.3], 'color', myEconStimOnColor )
            title(['currSubID: ', num2str(currSub)]);
            set(fh1, 'position', [50 100 1600 900], 'name', 'currSubPupilStimLocked_ CC vs IC');
            fh3 = figure(3); subplot(10,10,currSubID);
            data = currSubPupilStimLocked_CI;
            data(sum(isnan(data),2) >= length(data), :) = [];
            if size(data,1) > 1
                shadedErrorBar(1:length(nanmean(data)), nanmean(data), nansem(data), {'-', 'color', myPercColor, 'lineWidth', myLineWidth}, 1); hold on;
            else
                plot(1:length(data), data, '-', 'color', myPercColor, 'lineWidth', myLineWidth); hold on;
            end
            line([round(length(data)/2), round(length(data)/2)], [min(nanmean(data))-(max(nanmean(data))-min(nanmean(data)))*0.3, max(nanmean(data))+(max(nanmean(data))-min(nanmean(data)))*0.3], 'color', myPercStimOnColor )
            title(['currSubID: ', num2str(currSub)]);
            xlim(myPupilPlotXlimVals); ylim(myPupilPlotYlimVals); set(gca, 'XTickLabel', []);
            set(fh3, 'position', [1950 100 1600 900], 'name', 'currSubPupilStimLocked_CI');
            fh3 = figure(3); subplot(10,10,currSubID);
            data = currSubPupilStimLocked_II;
            data(sum(isnan(data),2) >= length(data), :) = [];
            if size(data,1) > 1
                shadedErrorBar(1:length(nanmean(data)), nanmean(data), nansem(data), {'-', 'color', myEconColor, 'lineWidth', myLineWidth}, 1); hold on;
            else
                plot(1:length(data), data, '-', 'color', myEconColor, 'lineWidth', myLineWidth); hold on;
            end
            line([round(length(data)/2), round(length(data)/2)], [min(nanmean(data))-(max(nanmean(data))-min(nanmean(data)))*0.3, max(nanmean(data))+(max(nanmean(data))-min(nanmean(data)))*0.3], 'color', myEconStimOnColor )
            title(['currSubID: ', num2str(currSub)]);
            xlim(myPupilPlotXlimVals); ylim(myPupilPlotYlimVals); set(gca, 'XTickLabel', []);
            set(fh3, 'position', [1950 100 1600 900], 'name', 'currSubPupilStimLocked_ CI vs II');
            fh5 = figure(5); subplot(10,10,currSubID);
            data = [currSubPupilStimLocked_CC;currSubPupilStimLocked_IC];
            data(sum(isnan(data),2) >= length(data), :) = [];
            if size(data,1) > 1
                shadedErrorBar(1:length(nanmean(data)), nanmean(data), nansem(data), {'-', 'color', myPercColor, 'lineWidth', myLineWidth}, 1); hold on;
            else
                plot(1:length(data), data, '-', 'color', myPercColor, 'lineWidth', myLineWidth); hold on;
            end
            line([round(length(data)/2), round(length(data)/2)], [min(nanmean(data))-(max(nanmean(data))-min(nanmean(data)))*0.3, max(nanmean(data))+(max(nanmean(data))-min(nanmean(data)))*0.3], 'color', myPercStimOnColor )
            title(['currSubID: ', num2str(currSub)]);
            xlim(myPupilPlotXlimVals); ylim(myPupilPlotYlimVals); set(gca, 'XTickLabel', []);
            fh5 = figure(5); subplot(10,10,currSubID);
            if task2Anal == 3
                data = [currSubPupilStimLocked_CI];
            else
                data = [currSubPupilStimLocked_CI;currSubPupilStimLocked_II];
            end
            data(sum(isnan(data),2) >= length(data), :) = [];
            if size(data,1) > 1
                shadedErrorBar(1:length(nanmean(data)), nanmean(data), nansem(data), {'-', 'color', myEconColor, 'lineWidth', myLineWidth}, 1); hold on;
            else
                plot(1:length(data), data, '-', 'color', myEconColor, 'lineWidth', myLineWidth); hold on;
            end
            line([round(length(data)/2), round(length(data)/2)], [min(nanmean(data))-(max(nanmean(data))-min(nanmean(data)))*0.3, max(nanmean(data))+(max(nanmean(data))-min(nanmean(data)))*0.3], 'color', myEconStimOnColor )
            title(['currSubID: ', num2str(currSub)]);
            xlim(myPupilPlotXlimVals); ylim(myPupilPlotYlimVals); set(gca, 'XTickLabel', []);
            set(fh5, 'position', [1950 100 1600 900], 'name', 'currSubPupilStimLocked_ Congruent vs Incongruent');
        end
        if size(currSubPupilStimLocked_CC, 1) > 1
            allSubsPupilStimLocked_AL(currSubID, :) = nanmean(currSubPupilStimLocked_AL);
            allSubsPupilRespLocked_AL(currSubID, :) = nanmean(currSubPupilRespLocked_AL);
            allSubsPupilStimLocked_CC(currSubID, :) = nanmean(currSubPupilStimLocked_CC);
            allSubsPupilRespLocked_CC(currSubID, :) = nanmean(currSubPupilRespLocked_CC);
            allSubsPupilStimLocked_IC(currSubID, :) = nanmean(currSubPupilStimLocked_IC);
            allSubsPupilRespLocked_IC(currSubID, :) = nanmean(currSubPupilRespLocked_IC);
            allSubsPupilStimLocked_CI(currSubID, :) = nanmean(currSubPupilStimLocked_CI);
            allSubsPupilRespLocked_CI(currSubID, :) = nanmean(currSubPupilRespLocked_CI);
            allSubsPupilStimLocked_II(currSubID, :) = nanmean(currSubPupilStimLocked_II);
            allSubsPupilRespLocked_II(currSubID, :) = nanmean(currSubPupilRespLocked_II);
            allSubsPupilStimLocked_Cong(currSubID, :) = nanmean([currSubPupilStimLocked_CC; currSubPupilStimLocked_IC]);
            allSubsPupilRespLocked_Cong(currSubID, :) = nanmean([currSubPupilRespLocked_CC; currSubPupilRespLocked_IC]);
            allSubsPupilStimLocked_Incong(currSubID, :) = nanmean([currSubPupilStimLocked_CI; currSubPupilStimLocked_II]);
            allSubsPupilRespLocked_Incong(currSubID, :) = nanmean([currSubPupilRespLocked_CI; currSubPupilRespLocked_II]);
        elseif size(currSubPupilStimLocked_CC, 1) == 1
            allSubsPupilStimLocked_AL(currSubID, :) = (currSubPupilStimLocked_AL);
            allSubsPupilRespLocked_AL(currSubID, :) = (currSubPupilRespLocked_AL);
            allSubsPupilStimLocked_CC(currSubID, :) = (currSubPupilStimLocked_CC);
            allSubsPupilRespLocked_CC(currSubID, :) = (currSubPupilRespLocked_CC);
            allSubsPupilStimLocked_IC(currSubID, :) = (currSubPupilStimLocked_IC);
            allSubsPupilRespLocked_IC(currSubID, :) = (currSubPupilRespLocked_IC);
            allSubsPupilStimLocked_CI(currSubID, :) = (currSubPupilStimLocked_CI);
            allSubsPupilRespLocked_CI(currSubID, :) = (currSubPupilRespLocked_CI);
            allSubsPupilStimLocked_II(currSubID, :) = (currSubPupilStimLocked_II);
            allSubsPupilRespLocked_II(currSubID, :) = (currSubPupilRespLocked_II);
            allSubsPupilStimLocked_Cong(currSubID, :) = nanmean([currSubPupilStimLocked_CC; currSubPupilStimLocked_IC]);
            allSubsPupilRespLocked_Cong(currSubID, :) = nanmean([currSubPupilRespLocked_CC; currSubPupilRespLocked_IC]);
            allSubsPupilStimLocked_Incong(currSubID, :) = (currSubPupilStimLocked_CI);
            allSubsPupilRespLocked_Incong(currSubID, :) = (currSubPupilRespLocked_CI);
        else
        end
        allSubsGazeBev = cat(1, allSubsGazeBev, currSubGazeBev);
        allSubsAllBev = cat(1, allSubsAllBev, currSubAllBev);
        allSubsAllPup = cat(1, allSubsAllPup, currSubAllPup);
        allSubUpDownSequenceMat = cat(1, allSubUpDownSequenceMat, currSubUpDownSequenceMat);
        allSubUpDownDurMat      = cat(1, allSubUpDownDurMat, currSubUpDownDurMat);
        allSubUpDownOnsetMat    = cat(1, allSubUpDownOnsetMat, currSubUpDownOnsetMat);
    end
    if saveCoditionResultsMats
        myHome = pwd;
        if task2Anal == 1
            if doLowPassFilteringOnly
                cd([myMotherPath , filesep, fullfile('Results_EST', 'EYE_ANAL', 'RUN1')]);
                eval(['save allSubsAllPup_LowPass_', subType,' allSubsAllPup']);
                eval(['save allSubsAllBev_LowPass_', subType,' allSubsAllBev']);
                eval(['save allSubsGazeBev_LowPass_', subType,' allSubsGazeBev']);
                eval(['save allSubsPupilStimLocked_AL_LowPass_', subType,' allSubsPupilStimLocked_AL']);
                eval(['save allSubsPupilRespLocked_AL_LowPass_', subType,' allSubsPupilRespLocked_AL']);
                eval(['save allSubsPupilStimLocked_CC_LowPass_', subType,' allSubsPupilStimLocked_CC']);
                eval(['save allSubsPupilRespLocked_CC_LowPass_', subType,' allSubsPupilRespLocked_CC']);
                eval(['save allSubsPupilStimLocked_IC_LowPass_', subType,' allSubsPupilStimLocked_IC']);
                eval(['save allSubsPupilRespLocked_IC_LowPass_', subType,' allSubsPupilRespLocked_IC']);
                eval(['save allSubsPupilStimLocked_CI_LowPass_', subType,' allSubsPupilStimLocked_CI']);
                eval(['save allSubsPupilRespLocked_CI_LowPass_', subType,' allSubsPupilRespLocked_CI']);
                eval(['save allSubsPupilStimLocked_II_LowPass_', subType,' allSubsPupilStimLocked_II']);
                eval(['save allSubsPupilRespLocked_II_LowPass_', subType,' allSubsPupilRespLocked_II']);
                eval(['save allSubUpDownSequenceMat_LowPass_', subType,' allSubUpDownSequenceMat']);
                eval(['save allSubUpDownDurMat_LowPass_', subType,' allSubUpDownDurMat']);
                eval(['save allSubUpDownOnsetMat_LowPass_', subType,' allSubUpDownOnsetMat']);
                cd(myHome);
            elseif doBandPassFiltering
                cd([myMotherPath , filesep, fullfile('Results_EST', 'EYE_ANAL', 'RUN1')]);
                eval(['save allSubsAllPup_BandPass_', subType,' allSubsAllPup']);
                eval(['save allSubsAllBev_BandPass_', subType,' allSubsAllBev']);
                eval(['save allSubsGazeBev_BandPass_', subType,' allSubsGazeBev']);
                eval(['save allSubsPupilStimLocked_AL_BandPass_', subType,' allSubsPupilStimLocked_AL']);
                eval(['save allSubsPupilRespLocked_AL_BandPass_', subType,' allSubsPupilRespLocked_AL']);
                eval(['save allSubsPupilStimLocked_CC_BandPass_', subType,' allSubsPupilStimLocked_CC']);
                eval(['save allSubsPupilRespLocked_CC_BandPass_', subType,' allSubsPupilRespLocked_CC']);
                eval(['save allSubsPupilStimLocked_IC_BandPass_', subType,' allSubsPupilStimLocked_IC']);
                eval(['save allSubsPupilRespLocked_IC_BandPass_', subType,' allSubsPupilRespLocked_IC']);
                eval(['save allSubsPupilStimLocked_CI_BandPass_', subType,' allSubsPupilStimLocked_CI']);
                eval(['save allSubsPupilRespLocked_CI_BandPass_', subType,' allSubsPupilRespLocked_CI']);
                eval(['save allSubsPupilStimLocked_II_BandPass_', subType,' allSubsPupilStimLocked_II']);
                eval(['save allSubsPupilRespLocked_II_BandPass_', subType,' allSubsPupilRespLocked_II']);
                eval(['save allSubUpDownSequenceMat_BandPass_', subType,' allSubUpDownSequenceMat']);
                eval(['save allSubUpDownDurMat_BandPass_', subType,' allSubUpDownDurMat']);
                eval(['save allSubUpDownOnsetMat_BandPass_', subType,' allSubUpDownOnsetMat']);
                cd(myHome);
            else
                cd([myMotherPath , filesep, fullfile('Results_EST', 'EYE_ANAL')]);
                eval(['save allSubsAllPup_', subType,' allSubsAllPup']);
                eval(['save allSubsAllBev_', subType,' allSubsAllBev']);
                eval(['save allSubsGazeBev_', subType,' allSubsGazeBev']);
                eval(['save allSubsPupilStimLocked_AL_', subType,' allSubsPupilStimLocked_AL']);
                eval(['save allSubsPupilRespLocked_AL_', subType,' allSubsPupilRespLocked_AL']);
                eval(['save allSubsPupilStimLocked_CC_', subType,' allSubsPupilStimLocked_CC']);
                eval(['save allSubsPupilRespLocked_CC_', subType,' allSubsPupilRespLocked_CC']);
                eval(['save allSubsPupilStimLocked_IC_', subType,' allSubsPupilStimLocked_IC']);
                eval(['save allSubsPupilRespLocked_IC_', subType,' allSubsPupilRespLocked_IC']);
                eval(['save allSubsPupilStimLocked_CI_', subType,' allSubsPupilStimLocked_CI']);
                eval(['save allSubsPupilRespLocked_CI_', subType,' allSubsPupilRespLocked_CI']);
                eval(['save allSubsPupilStimLocked_II_', subType,' allSubsPupilStimLocked_II']);
                eval(['save allSubsPupilRespLocked_II_', subType,' allSubsPupilRespLocked_II']);
                eval(['save allSubUpDownSequenceMat_', subType,' allSubUpDownSequenceMat']);
                eval(['save allSubUpDownDurMat_', subType,' allSubUpDownDurMat']);
                eval(['save allSubUpDownOnsetMat_', subType,' allSubUpDownOnsetMat']);
                cd(myHome);
            end
        elseif task2Anal == 2
            if doLowPassFilteringOnly
                cd([myMotherPath , filesep, fullfile('Results_EST', 'EYE_ANAL', 'RUN1')]);
                eval(['save allSubsAllPup_LowPass_BRIGHTNESS_', subType,' allSubsAllPup']);
                eval(['save allSubsAllBev_LowPass_BRIGHTNESS_', subType,' allSubsAllBev']);
                eval(['save allSubsGazeBev_LowPass_BRIGHTNESS_', subType,' allSubsGazeBev']);
                eval(['save allSubsPupilStimLocked_AL_LowPass_BRIGHTNESS_', subType,' allSubsPupilStimLocked_AL']);
                eval(['save allSubsPupilRespLocked_AL_LowPass_BRIGHTNESS_', subType,' allSubsPupilRespLocked_AL']);
                eval(['save allSubsPupilStimLocked_CC_LowPass_BRIGHTNESS_', subType,' allSubsPupilStimLocked_CC']);
                eval(['save allSubsPupilRespLocked_CC_LowPass_BRIGHTNESS_', subType,' allSubsPupilRespLocked_CC']);
                eval(['save allSubsPupilStimLocked_IC_LowPass_BRIGHTNESS_', subType,' allSubsPupilStimLocked_IC']);
                eval(['save allSubsPupilRespLocked_IC_LowPass_BRIGHTNESS_', subType,' allSubsPupilRespLocked_IC']);
                eval(['save allSubsPupilStimLocked_CI_LowPass_BRIGHTNESS_', subType,' allSubsPupilStimLocked_CI']);
                eval(['save allSubsPupilRespLocked_CI_LowPass_BRIGHTNESS_', subType,' allSubsPupilRespLocked_CI']);
                eval(['save allSubsPupilStimLocked_II_LowPass_BRIGHTNESS_', subType,' allSubsPupilStimLocked_II']);
                eval(['save allSubsPupilRespLocked_II_LowPass_BRIGHTNESS_', subType,' allSubsPupilRespLocked_II']);
                eval(['save allSubUpDownSequenceMat_LowPass_BRIGHTNESS_', subType,' allSubUpDownSequenceMat']);
                eval(['save allSubUpDownDurMat_LowPass_BRIGHTNESS_', subType,' allSubUpDownDurMat']);
                eval(['save allSubUpDownOnsetMat_LowPass_BRIGHTNESS_', subType,' allSubUpDownOnsetMat']);
                cd(myHome);
            elseif doBandPassFiltering
                cd([myMotherPath , filesep, fullfile('Results_EST', 'EYE_ANAL', 'RUN1')]);
                eval(['save allSubsAllPup_BandPass_BRIGHTNESS_', subType,' allSubsAllPup']);
                eval(['save allSubsAllBev_BandPass_BRIGHTNESS_', subType,' allSubsAllBev']);
                eval(['save allSubsGazeBev_BandPass_BRIGHTNESS_', subType,' allSubsGazeBev']);
                eval(['save allSubsPupilStimLocked_AL_BandPass_BRIGHTNESS_', subType,' allSubsPupilStimLocked_AL']);
                eval(['save allSubsPupilRespLocked_AL_BandPass_BRIGHTNESS_', subType,' allSubsPupilRespLocked_AL']);
                eval(['save allSubsPupilStimLocked_CC_BandPass_BRIGHTNESS_', subType,' allSubsPupilStimLocked_CC']);
                eval(['save allSubsPupilRespLocked_CC_BandPass_BRIGHTNESS_', subType,' allSubsPupilRespLocked_CC']);
                eval(['save allSubsPupilStimLocked_IC_BandPass_BRIGHTNESS_', subType,' allSubsPupilStimLocked_IC']);
                eval(['save allSubsPupilRespLocked_IC_BandPass_BRIGHTNESS_', subType,' allSubsPupilRespLocked_IC']);
                eval(['save allSubsPupilStimLocked_CI_BandPass_BRIGHTNESS_', subType,' allSubsPupilStimLocked_CI']);
                eval(['save allSubsPupilRespLocked_CI_BandPass_BRIGHTNESS_', subType,' allSubsPupilRespLocked_CI']);
                eval(['save allSubsPupilStimLocked_II_BandPass_BRIGHTNESS_', subType,' allSubsPupilStimLocked_II']);
                eval(['save allSubsPupilRespLocked_II_BandPass_BRIGHTNESS_', subType,' allSubsPupilRespLocked_II']);
                eval(['save allSubUpDownSequenceMat_BandPass_BRIGHTNESS_', subType,' allSubUpDownSequenceMat']);
                eval(['save allSubUpDownDurMat_BandPass_BRIGHTNESS_', subType,' allSubUpDownDurMat']);
                eval(['save allSubUpDownOnsetMat_BandPass_BRIGHTNESS_', subType,' allSubUpDownOnsetMat']);
                cd(myHome);
            else
                cd([myMotherPath , filesep, fullfile('Results_EST', 'EYE_ANAL', 'RUN1')]);
                eval(['save allSubsAllPup_BRIGHTNESS_', subType,' allSubsAllPup']);
                eval(['save allSubsAllBev_BRIGHTNESS_', subType,' allSubsAllBev']);
                eval(['save allSubsGazeBev_BRIGHTNESS_', subType,' allSubsGazeBev']);
                eval(['save allSubsPupilStimLocked_AL_BRIGHTNESS_', subType,' allSubsPupilStimLocked_AL']);
                eval(['save allSubsPupilRespLocked_AL_BRIGHTNESS_', subType,' allSubsPupilRespLocked_AL']);
                eval(['save allSubsPupilStimLocked_CC_BRIGHTNESS_', subType,' allSubsPupilStimLocked_CC']);
                eval(['save allSubsPupilRespLocked_CC_BRIGHTNESS_', subType,' allSubsPupilRespLocked_CC']);
                eval(['save allSubsPupilStimLocked_IC_BRIGHTNESS_', subType,' allSubsPupilStimLocked_IC']);
                eval(['save allSubsPupilRespLocked_IC_BRIGHTNESS_', subType,' allSubsPupilRespLocked_IC']);
                eval(['save allSubsPupilStimLocked_CI_BRIGHTNESS_', subType,' allSubsPupilStimLocked_CI']);
                eval(['save allSubsPupilRespLocked_CI_BRIGHTNESS_', subType,' allSubsPupilRespLocked_CI']);
                eval(['save allSubsPupilStimLocked_II_BRIGHTNESS_', subType,' allSubsPupilStimLocked_II']);
                eval(['save allSubsPupilRespLocked_II_BRIGHTNESS_', subType,' allSubsPupilRespLocked_II']);
                eval(['save allSubUpDownSequenceMat_BRIGHTNESS_', subType,' allSubUpDownSequenceMat']);
                eval(['save allSubUpDownDurMat_BRIGHTNESS_', subType,' allSubUpDownDurMat']);
                eval(['save allSubUpDownOnsetMat_BRIGHTNESS_', subType,' allSubUpDownOnsetMat']);
                cd(myHome);
            end
        elseif task2Anal == 3
            if doLowPassFilteringOnly
                cd([myMotherPath , filesep, fullfile('Results_EST', 'EYE_ANAL', 'RUN1')]);
                eval(['save allSubsAllPup_LowPass_ODDBALL_', subType,' allSubsAllPup']);
                eval(['save allSubsAllBev_LowPass_ODDBALL_', subType,' allSubsAllBev']);
                eval(['save allSubsGazeBev_LowPass_ODDBALL_', subType,' allSubsGazeBev']);
                eval(['save allSubsPupilStimLocked_AL_LowPass_ODDBALL_', subType,' allSubsPupilStimLocked_AL']);
                eval(['save allSubsPupilRespLocked_AL_LowPass_ODDBALL_', subType,' allSubsPupilRespLocked_AL']);
                eval(['save allSubsPupilStimLocked_CC_LowPass_ODDBALL_', subType,' allSubsPupilStimLocked_CC']);
                eval(['save allSubsPupilRespLocked_CC_LowPass_ODDBALL_', subType,' allSubsPupilRespLocked_CC']);
                eval(['save allSubsPupilStimLocked_IC_LowPass_ODDBALL_', subType,' allSubsPupilStimLocked_IC']);
                eval(['save allSubsPupilRespLocked_IC_LowPass_ODDBALL_', subType,' allSubsPupilRespLocked_IC']);
                eval(['save allSubsPupilStimLocked_CI_LowPass_ODDBALL_', subType,' allSubsPupilStimLocked_CI']);
                eval(['save allSubsPupilRespLocked_CI_LowPass_ODDBALL_', subType,' allSubsPupilRespLocked_CI']);
                eval(['save allSubsPupilStimLocked_II_LowPass_ODDBALL_', subType,' allSubsPupilStimLocked_II']);
                eval(['save allSubsPupilRespLocked_II_LowPass_ODDBALL_', subType,' allSubsPupilRespLocked_II']);
                eval(['save allSubUpDownSequenceMat_LowPass_ODDBALL_', subType,' allSubUpDownSequenceMat']);
                eval(['save allSubUpDownDurMat_LowPass_ODDBALL_', subType,' allSubUpDownDurMat']);
                eval(['save allSubUpDownOnsetMat_LowPass_ODDBALL_', subType,' allSubUpDownOnsetMat']);
                cd(myHome);
            elseif doBandPassFiltering
                cd([myMotherPath , filesep, fullfile('Results_EST', 'EYE_ANAL', 'RUN1')]);
                eval(['save allSubsAllBev_BandPass_ODDBALL_', subType,' allSubsAllBev']);
                eval(['save allSubsGazeBev_BandPass_ODDBALL_', subType,' allSubsGazeBev']);
                eval(['save allSubsPupilStimLocked_AL_BandPass_ODDBALL_', subType,' allSubsPupilStimLocked_AL']);
                eval(['save allSubsPupilRespLocked_AL_BandPass_ODDBALL_', subType,' allSubsPupilRespLocked_AL']);
                eval(['save allSubsPupilStimLocked_CC_BandPass_ODDBALL_', subType,' allSubsPupilStimLocked_CC']);
                eval(['save allSubsPupilRespLocked_CC_BandPass_ODDBALL_', subType,' allSubsPupilRespLocked_CC']);
                eval(['save allSubsPupilStimLocked_IC_BandPass_ODDBALL_', subType,' allSubsPupilStimLocked_IC']);
                eval(['save allSubsPupilRespLocked_IC_BandPass_ODDBALL_', subType,' allSubsPupilRespLocked_IC']);
                eval(['save allSubsPupilStimLocked_CI_BandPass_ODDBALL_', subType,' allSubsPupilStimLocked_CI']);
                eval(['save allSubsPupilRespLocked_CI_BandPass_ODDBALL_', subType,' allSubsPupilRespLocked_CI']);
                eval(['save allSubsPupilStimLocked_II_BandPass_ODDBALL_', subType,' allSubsPupilStimLocked_II']);
                eval(['save allSubsPupilRespLocked_II_BandPass_ODDBALL_', subType,' allSubsPupilRespLocked_II']);
                eval(['save allSubUpDownSequenceMat_BandPass_ODDBALL_', subType,' allSubUpDownSequenceMat']);
                eval(['save allSubUpDownDurMat_BandPass_ODDBALL_', subType,' allSubUpDownDurMat']);
                eval(['save allSubUpDownOnsetMat_BandPass_ODDBALL_', subType,' allSubUpDownOnsetMat']);
                cd(myHome);
            else
                cd([myMotherPath , filesep, fullfile('Results_EST', 'EYE_ANAL', 'RUN1')]);
                eval(['save allSubsAllPup_ODDBALL_', subType,' allSubsAllPup']);
                eval(['save allSubsAllBev_ODDBALL_', subType,' allSubsAllBev']);
                eval(['save allSubsGazeBev_ODDBALL_', subType,' allSubsGazeBev']);
                eval(['save allSubsPupilStimLocked_AL_ODDBALL_', subType,' allSubsPupilStimLocked_AL']);
                eval(['save allSubsPupilRespLocked_AL_ODDBALL_', subType,' allSubsPupilRespLocked_AL']);
                eval(['save allSubsPupilStimLocked_CC_ODDBALL_', subType,' allSubsPupilStimLocked_CC']);
                eval(['save allSubsPupilRespLocked_CC_ODDBALL_', subType,' allSubsPupilRespLocked_CC']);
                eval(['save allSubsPupilStimLocked_IC_ODDBALL_', subType,' allSubsPupilStimLocked_IC']);
                eval(['save allSubsPupilRespLocked_IC_ODDBALL_', subType,' allSubsPupilRespLocked_IC']);
                eval(['save allSubsPupilStimLocked_CI_ODDBALL_', subType,' allSubsPupilStimLocked_CI']);
                eval(['save allSubsPupilRespLocked_CI_ODDBALL_', subType,' allSubsPupilRespLocked_CI']);
                eval(['save allSubsPupilStimLocked_II_ODDBALL_', subType,' allSubsPupilStimLocked_II']);
                eval(['save allSubsPupilRespLocked_II_ODDBALL_', subType,' allSubsPupilRespLocked_II']);
                eval(['save allSubUpDownSequenceMat_ODDBALL_', subType,' allSubUpDownSequenceMat']);
                eval(['save allSubUpDownDurMat_ODDBALL_', subType,' allSubUpDownDurMat']);
                eval(['save allSubUpDownOnsetMat_ODDBALL_', subType,' allSubUpDownOnsetMat']);
                cd(myHome);
            end
        end
    end
end
myHome = pwd;
task2Anal = 1;
if task2Anal == 1
    if doLowPassFilteringOnly
        cd([myMotherPath , filesep, fullfile('Results_EST', 'EYE_ANAL', 'RUN1')]);
        eval(['load allSubsAllPup_LowPass_', subType]);
        eval(['load allSubsAllBev_LowPass_', subType]);
        eval(['load allSubsGazeBev_LowPass_', subType]);
        eval(['load allSubsPupilStimLocked_AL_LowPass_', subType]);
        eval(['load allSubsPupilRespLocked_AL_LowPass_', subType]);
        eval(['load allSubsPupilStimLocked_CC_LowPass_', subType]);
        eval(['load allSubsPupilRespLocked_CC_LowPass_', subType]);
        eval(['load allSubsPupilStimLocked_IC_LowPass_', subType]);
        eval(['load allSubsPupilRespLocked_IC_LowPass_', subType]);
        eval(['load allSubsPupilStimLocked_CI_LowPass_', subType]);
        eval(['load allSubsPupilRespLocked_CI_LowPass_', subType]);
        eval(['load allSubsPupilStimLocked_II_LowPass_', subType]);
        eval(['load allSubsPupilRespLocked_II_LowPass_', subType]);
        eval(['load allSubUpDownSequenceMat_LowPass_', subType]);
        eval(['load allSubUpDownDurMat_LowPass_', subType]);
        eval(['load allSubUpDownOnsetMat_LowPass_', subType]);
        cd(myHome);
    elseif doBandPassFiltering
        cd([myMotherPath , filesep, fullfile('Results_EST', 'EYE_ANAL', 'RUN1')]);
        eval(['load allSubsAllPup_BandPass_', subType]);
        eval(['load allSubsAllBev_BandPass_', subType]);
        eval(['load allSubsGazeBev_BandPass_', subType]);
        eval(['load allSubsPupilStimLocked_AL_BandPass_', subType]);
        eval(['load allSubsPupilRespLocked_AL_BandPass_', subType]);
        eval(['load allSubsPupilStimLocked_CC_BandPass_', subType]);
        eval(['load allSubsPupilRespLocked_CC_BandPass_', subType]);
        eval(['load allSubsPupilStimLocked_IC_BandPass_', subType]);
        eval(['load allSubsPupilRespLocked_IC_BandPass_', subType]);
        eval(['load allSubsPupilStimLocked_CI_BandPass_', subType]);
        eval(['load allSubsPupilRespLocked_CI_BandPass_', subType]);
        eval(['load allSubsPupilStimLocked_II_BandPass_', subType]);
        eval(['load allSubsPupilRespLocked_II_BandPass_', subType]);
        eval(['load allSubUpDownSequenceMat_BandPass_', subType]);
        eval(['load allSubUpDownDurMat_BandPass_', subType]);
        eval(['load allSubUpDownOnsetMat_BandPass_', subType]);
        cd(myHome);
    else
        cd([myMotherPath , filesep, fullfile('Results_EST', 'EYE_ANAL', 'RUN1')]);
        eval(['load allSubsAllPup_', subType]);
        eval(['load allSubsAllBev_', subType]);
        eval(['load allSubsGazeBev_', subType]);
        eval(['load allSubsPupilStimLocked_AL_', subType]);
        eval(['load allSubsPupilRespLocked_AL_', subType]);
        eval(['load allSubsPupilStimLocked_CC_', subType]);
        eval(['load allSubsPupilRespLocked_CC_', subType]);
        eval(['load allSubsPupilStimLocked_IC_', subType]);
        eval(['load allSubsPupilRespLocked_IC_', subType]);
        eval(['load allSubsPupilStimLocked_CI_', subType]);
        eval(['load allSubsPupilRespLocked_CI_', subType]);
        eval(['load allSubsPupilStimLocked_II_', subType]);
        eval(['load allSubsPupilRespLocked_II_', subType]);
        eval(['load allSubUpDownSequenceMat_', subType]);
        eval(['load allSubUpDownDurMat_', subType]);
        eval(['load allSubUpDownOnsetMat_', subType]);
        cd(myHome);
    end
elseif task2Anal == 2
    if doLowPassFilteringOnly
        cd([myMotherPath , filesep, fullfile('Results_EST', 'EYE_ANAL', 'RUN1')]);
        eval(['load allSubsAllPup_LowPass_BRIGHTNESS_', subType]);
        eval(['load allSubsAllBev_LowPass_BRIGHTNESS_', subType]);
        eval(['load allSubsGazeBev_LowPass_BRIGHTNESS_', subType]);
        eval(['load allSubsPupilStimLocked_AL_LowPass_BRIGHTNESS_', subType]);
        eval(['load allSubsPupilRespLocked_AL_LowPass_BRIGHTNESS_', subType]);
        eval(['load allSubsPupilStimLocked_CC_LowPass_BRIGHTNESS_', subType]);
        eval(['load allSubsPupilRespLocked_CC_LowPass_BRIGHTNESS_', subType]);
        eval(['load allSubsPupilStimLocked_IC_LowPass_BRIGHTNESS_', subType]);
        eval(['load allSubsPupilRespLocked_IC_LowPass_BRIGHTNESS_', subType]);
        eval(['load allSubsPupilStimLocked_CI_LowPass_BRIGHTNESS_', subType]);
        eval(['load allSubsPupilRespLocked_CI_LowPass_BRIGHTNESS_', subType]);
        eval(['load allSubsPupilStimLocked_II_LowPass_BRIGHTNESS_', subType]);
        eval(['load allSubsPupilRespLocked_II_LowPass_BRIGHTNESS_', subType]);
        eval(['load allSubUpDownSequenceMat_LowPass_BRIGHTNESS_', subType]);
        eval(['load allSubUpDownDurMat_LowPass_BRIGHTNESS_', subType]);
        eval(['load allSubUpDownOnsetMat_LowPass_BRIGHTNESS_', subType]);
        cd(myHome);
    elseif doBandPassFiltering
        cd([myMotherPath , filesep, fullfile('Results_EST', 'EYE_ANAL', 'RUN1')]);
        eval(['load allSubsAllPup_BandPass_BRIGHTNESS_', subType]);
        eval(['load allSubsAllBev_BandPass_BRIGHTNESS_', subType]);
        eval(['load allSubsGazeBev_BandPass_BRIGHTNESS_', subType]);
        eval(['load allSubsPupilStimLocked_AL_BandPass_BRIGHTNESS_', subType]);
        eval(['load allSubsPupilRespLocked_AL_BandPass_BRIGHTNESS_', subType]);
        eval(['load allSubsPupilStimLocked_CC_BandPass_BRIGHTNESS_', subType]);
        eval(['load allSubsPupilRespLocked_CC_BandPass_BRIGHTNESS_', subType]);
        eval(['load allSubsPupilStimLocked_IC_BandPass_BRIGHTNESS_', subType]);
        eval(['load allSubsPupilRespLocked_IC_BandPass_BRIGHTNESS_', subType]);
        eval(['load allSubsPupilStimLocked_CI_BandPass_BRIGHTNESS_', subType]);
        eval(['load allSubsPupilRespLocked_CI_BandPass_BRIGHTNESS_', subType]);
        eval(['load allSubsPupilStimLocked_II_BandPass_BRIGHTNESS_', subType]);
        eval(['load allSubsPupilRespLocked_II_BandPass_BRIGHTNESS_', subType]);
        eval(['load allSubUpDownSequenceMat_BandPass_BRIGHTNESS_', subType]);
        eval(['load allSubUpDownDurMat_BandPass_BRIGHTNESS_', subType]);
        eval(['load allSubUpDownOnsetMat_BandPass_BRIGHTNESS_', subType]);
        cd(myHome);
    else
        cd([myMotherPath , filesep, fullfile('Results_EST', 'EYE_ANAL', 'RUN1')]);
        eval(['load allSubsAllPup_BRIGHTNESS_', subType]);
        eval(['load allSubsAllBev_BRIGHTNESS_', subType]);
        eval(['load allSubsGazeBev_BRIGHTNESS_', subType]);
        eval(['load allSubsPupilStimLocked_AL_BRIGHTNESS_', subType]);
        eval(['load allSubsPupilRespLocked_AL_BRIGHTNESS_', subType]);
        eval(['load allSubsPupilStimLocked_CC_BRIGHTNESS_', subType]);
        eval(['load allSubsPupilRespLocked_CC_BRIGHTNESS_', subType]);
        eval(['load allSubsPupilStimLocked_IC_BRIGHTNESS_', subType]);
        eval(['load allSubsPupilRespLocked_IC_BRIGHTNESS_', subType]);
        eval(['load allSubsPupilStimLocked_CI_BRIGHTNESS_', subType]);
        eval(['load allSubsPupilRespLocked_CI_BRIGHTNESS_', subType]);
        eval(['load allSubsPupilStimLocked_II_BRIGHTNESS_', subType]);
        eval(['load allSubsPupilRespLocked_II_BRIGHTNESS_', subType]);
        eval(['load allSubUpDownSequenceMat_BRIGHTNESS_', subType]);
        eval(['load allSubUpDownDurMat_BRIGHTNESS_', subType]);
        eval(['load allSubUpDownOnsetMat_BRIGHTNESS_', subType]);
        cd(myHome);
    end
elseif task2Anal == 3
    if doLowPassFilteringOnly
        cd([myMotherPath , filesep, fullfile('Results_EST', 'EYE_ANAL', 'RUN1')]);
        eval(['load allSubsAllPup_LowPass_ODDBALL_', subType]);
        eval(['load allSubsAllBev_LowPass_ODDBALL_', subType]);
        eval(['load allSubsGazeBev_LowPass_ODDBALL_', subType]);
        eval(['load allSubsPupilStimLocked_AL_LowPass_ODDBALL_', subType]);
        eval(['load allSubsPupilRespLocked_AL_LowPass_ODDBALL_', subType]);
        eval(['load allSubsPupilStimLocked_CC_LowPass_ODDBALL_', subType]);
        eval(['load allSubsPupilRespLocked_CC_LowPass_ODDBALL_', subType]);
        eval(['load allSubsPupilStimLocked_IC_LowPass_ODDBALL_', subType]);
        eval(['load allSubsPupilRespLocked_IC_LowPass_ODDBALL_', subType]);
        eval(['load allSubsPupilStimLocked_CI_LowPass_ODDBALL_', subType]);
        eval(['load allSubsPupilRespLocked_CI_LowPass_ODDBALL_', subType]);
        eval(['load allSubsPupilStimLocked_II_LowPass_ODDBALL_', subType]);
        eval(['load allSubsPupilRespLocked_II_LowPass_ODDBALL_', subType]);
        eval(['load allSubUpDownSequenceMat_LowPass_ODDBALL_', subType]);
        eval(['load allSubUpDownDurMat_LowPass_ODDBALL_', subType]);
        eval(['load allSubUpDownOnsetMat_LowPass_ODDBALL_', subType]);
        cd(myHome);
    elseif doBandPassFiltering
        cd([myMotherPath , filesep, fullfile('Results_EST', 'EYE_ANAL', 'RUN1')]);
        eval(['load allSubsAllPup_BandPass_ODDBALL_', subType]);
        eval(['load allSubsAllBev_BandPass_ODDBALL_', subType]);
        eval(['load allSubsGazeBev_BandPass_ODDBALL_', subType]);
        eval(['load allSubsPupilStimLocked_AL_BandPass_ODDBALL_', subType]);
        eval(['load allSubsPupilRespLocked_AL_BandPass_ODDBALL_', subType]);
        eval(['load allSubsPupilStimLocked_CC_BandPass_ODDBALL_', subType]);
        eval(['load allSubsPupilRespLocked_CC_BandPass_ODDBALL_', subType]);
        eval(['load allSubsPupilStimLocked_IC_BandPass_ODDBALL_', subType]);
        eval(['load allSubsPupilRespLocked_IC_BandPass_ODDBALL_', subType]);
        eval(['load allSubsPupilStimLocked_CI_BandPass_ODDBALL_', subType]);
        eval(['load allSubsPupilRespLocked_CI_BandPass_ODDBALL_', subType]);
        eval(['load allSubsPupilStimLocked_II_BandPass_ODDBALL_', subType]);
        eval(['load allSubsPupilRespLocked_II_BandPass_ODDBALL_', subType]);
        eval(['load allSubUpDownSequenceMat_BandPass_ODDBALL_', subType]);
        eval(['load allSubUpDownDurMat_BandPass_ODDBALL_', subType]);
        eval(['load allSubUpDownOnsetMat_BandPass_ODDBALL_', subType]);
        cd(myHome);
    else
        cd([myMotherPath , filesep, fullfile('Results_EST', 'EYE_ANAL', 'RUN1')]);
        eval(['load allSubsAllPup_ODDBALL_', subType]);
        eval(['load allSubsAllBev_ODDBALL_', subType]);
        eval(['load allSubsGazeBev_ODDBALL_', subType]);
        eval(['load allSubsPupilStimLocked_AL_ODDBALL_', subType]);
        eval(['load allSubsPupilRespLocked_AL_ODDBALL_', subType]);
        eval(['load allSubsPupilStimLocked_CC_ODDBALL_', subType]);
        eval(['load allSubsPupilRespLocked_CC_ODDBALL_', subType]);
        eval(['load allSubsPupilStimLocked_IC_ODDBALL_', subType]);
        eval(['load allSubsPupilRespLocked_IC_ODDBALL_', subType]);
        eval(['load allSubsPupilStimLocked_CI_ODDBALL_', subType]);
        eval(['load allSubsPupilRespLocked_CI_ODDBALL_', subType]);
        eval(['load allSubsPupilStimLocked_II_ODDBALL_', subType]);
        eval(['load allSubsPupilRespLocked_II_ODDBALL_', subType]);
        eval(['load allSubUpDownSequenceMat_ODDBALL_', subType]);
        eval(['load allSubUpDownDurMat_ODDBALL_', subType]);
        eval(['load allSubUpDownOnsetMat_ODDBALL_', subType]);
        cd(myHome);
    end
end
if doMoveAvg
    allSubsPupilStimLocked_CC = movmean(allSubsPupilStimLocked_CC, moveAvgKernel, 2);
    allSubsPupilStimLocked_IC = movmean(allSubsPupilStimLocked_IC, moveAvgKernel, 2);
    allSubsPupilStimLocked_II = movmean(allSubsPupilStimLocked_II, moveAvgKernel, 2);
    allSubsPupilStimLocked_CI = movmean(allSubsPupilStimLocked_CI, moveAvgKernel, 2);
end
keyboard
if plotGranMeanPupil
    if task2Anal == 1
        fh = figure; set(gcf, 'color', [1 1 1])
        allSubsPupilStimLocked_CONG = (allSubsPupilStimLocked_CC + allSubsPupilStimLocked_IC)./2;
        data = allSubsPupilStimLocked_CONG;
        shadedErrorBar(1:length(nanmean(data)), nanmean(data), nansem(data), {'-', 'color', myCongruentColor, 'lineWidth', myLineWidth}, 1); hold on;
        title('INCONGRUENT (Red) CONGRUENT (Green)'); hold on;
        xlim(myPupilPlotXlimVals); ylim([-1.3 .6]);
        set(gca, 'xtick', [min(myPupilPlotXlimVals):1000:max(myPupilPlotXlimVals)], 'xticklabel', [(min(myPupilPlotXlimVals)- 10000)/1000     : (max(myPupilPlotXlimVals)-10000)/1000]);
        allSubsPupilStimLocked_INCONG = (allSubsPupilStimLocked_CI + allSubsPupilStimLocked_II)./2;
        data = allSubsPupilStimLocked_INCONG;
        shadedErrorBar(1:length(nanmean(data)), nanmean(data), nansem(data), {'-', 'color', myIncongruentColor, 'lineWidth', myLineWidth}, 1); hold on;
        grid on
        xlim(myPupilPlotXlimVals); ylim([-1.3 .6]);
        set(gca, 'xtick', [min(myPupilPlotXlimVals):1000:max(myPupilPlotXlimVals)], 'xticklabel', [(min(myPupilPlotXlimVals)- 10000)/1000     : (max(myPupilPlotXlimVals)-10000)/1000]);
        ylabel('pupil dilation (z-scored)', 'fontsize', 14)
        axis square
        set(gca, 'fontsize', 17);
        ylim([-.8 .5])
        if saveFigures
            myHome = pwd;
            cd([myMotherPath , filesep, fullfile('Results_EST', 'EYE_ANAL', 'RUN1')]);
            set(gcf, 'paperunits', 'centimeters', 'paperposition', [0 0 15 15])
            if doLowPassFilteringOnly
                print(['PupilAvgLP_Task_',num2str(task2Anal),'_subType_', subType ,'.tif'], '-dtiff', '-r300')
            elseif doBandPassFiltering
                print(['PupilAvgBP_Task_',num2str(task2Anal),'_subType_', subType ,'.tif'], '-dtiff', '-r300')
            else
                print(['PupilAvg_Task_',num2str(task2Anal),'_subType_', subType ,'.tif'], '-dtiff', '-r300')
            end
            cd(myHome)
        end
        figure; set(gcf, 'color', [1 1 1])
        data = allSubsPupilStimLocked_INCONG - allSubsPupilStimLocked_CONG;
        [H_T, P_T, CI_T, STATS_T] = ttest(data); H_T(1) = 0;
        diff_H_T = diff(H_T); [startID] = find(diff_H_T == 1); [endID] = find(diff_H_T == -1);
maxLength = max(length(startID), length(endID));
if length(startID) < maxLength
    startID(end+1:maxLength) = NaN;
elseif length(endID) < maxLength
    endID(end+1:maxLength) = NaN;
end
clustersStartEndIDs = [startID, endID];
disp(clustersStartEndIDs);
        clustersStartEndIDs = [startID', endID'];
        if ~isempty(clustersStartEndIDs)
            foundClustersID = startID > 6000 & startID < 10000;
            foundClustersID = startID > 10000 & startID < 15000;
            foundClustersStartEndIDs = clustersStartEndIDs(foundClustersID,:);
            clusterSizes = clustersStartEndIDs(foundClustersID,2) - clustersStartEndIDs(foundClustersID,1);
            clusterIntervalToUse = foundClustersStartEndIDs(clusterSizes == max(clusterSizes), : );
            clusterIntervalToUsePreStim = clusterIntervalToUse;
            myLightGrey = [.7 .7 .7];
            myQs = [-.25 .251]; myXs = clusterIntervalToUsePreStim;
            try
                patch([myXs(1), myXs(1), myXs(2), myXs(2), myXs(1)], [myQs(1), myQs(2), myQs(2), myQs(1), myQs(1)], myLightGrey, 'FaceAlpha',.5, 'linestyle', 'none');
                Pupil_IncongHigherCong_Vect = mean(data(:,clusterIntervalToUse(1):clusterIntervalToUse(2)),2);
            catch
            end
                hold on
        end
        shadedErrorBar(1:length(nanmean(data)), nanmean(data), nansem(data), {'-', 'color', myIncongruentColor, 'lineWidth', myLineWidth}, 1); hold on;
        title('INCONGRUENT > CONGRUENT '); hold on;
        xlim([6000 15000]);
        ylim([-.15 .21]);
        set(gca, 'xtick', [min(myPupilPlotXlimVals):1000:max(myPupilPlotXlimVals)], 'xticklabel', [(min(myPupilPlotXlimVals)- 10000)/1000     : (max(myPupilPlotXlimVals)-10000)/1000]);
        ylabel('pupil dilation (z-scored)', 'fontsize', 14)
        box on; axis square; set(gca, 'fontsize', 17);
        if saveFigures
            myHome = pwd;
            cd([myMotherPath , filesep, fullfile('Results_EST', 'EYE_ANAL', 'RUN1')]);
            set(gcf, 'paperunits', 'centimeters', 'paperposition', [0 0 15 15])
            if doLowPassFilteringOnly
                print(['PupilAvgLP_Task_',num2str(task2Anal),'_subType_', subType, '_Diff_IvsC', '.tif'], '-dtiff', '-r300')
            elseif doBandPassFiltering
                print(['PupilAvgBP_Task_',num2str(task2Anal),'_subType_', subType, '_Diff_IvsC', '.tif'], '-dtiff', '-r300')
            else
                print(['PupilAvg_Task_',num2str(task2Anal),'_subType_', subType, '_Diff_IvsC', '.tif'], '-dtiff', '-r300')
            end
            cd(myHome)
        end
        figure; set(gcf, 'color', [1 1 1])
        data = allSubsPupilStimLocked_CC;
        shadedErrorBar(1:length(nanmean(data)), nanmean(data), nansem(data), {'-', 'color', myCCColor, 'lineWidth', myLineWidth}, 1); hold on;
        title('IC (dark)  CC (light)'); hold on;
        xlim(myPupilPlotXlimVals); ylim([-1.3 .6]);
        ylabel('pupil dilation (z-scored)', 'fontsize', 14);
        set(gca, 'xtick', [min(myPupilPlotXlimVals):1000:max(myPupilPlotXlimVals)], 'xticklabel', [(min(myPupilPlotXlimVals)- 10000)/1000     : (max(myPupilPlotXlimVals)-10000)/1000]);
        data = allSubsPupilStimLocked_IC;
        shadedErrorBar(1:length(nanmean(data)), nanmean(data), nansem(data), {'-', 'color', myICColor, 'lineWidth', myLineWidth}, 1); hold on;
        xlim(myPupilPlotXlimVals);
        ylim([-.6 .5]);
        set(gca, 'xtick', [min(myPupilPlotXlimVals):1000:max(myPupilPlotXlimVals)], 'xticklabel', [(min(myPupilPlotXlimVals)- 10000)/1000     : (max(myPupilPlotXlimVals)-10000)/1000]);
        axis square
        set(gca, 'fontsize', 17);
        if saveFigures
            myHome = pwd;
            cd([myMotherPath , filesep, fullfile('Results_EST', 'EYE_ANAL', 'RUN1')]);
            set(gcf, 'paperunits', 'centimeters', 'paperposition', [0 0 15 15])
            if doLowPassFilteringOnly
                print(['PupilAvgLP_Task_',num2str(task2Anal),'_subType_', subType, '_Avg_IC_CC', '.tif'], '-dtiff', '-r300')
            elseif doBandPassFiltering
                print(['PupilAvgBP_Task_',num2str(task2Anal),'_subType_', subType, '_Avg_IC_CC', '.tif'], '-dtiff', '-r300')
            else
                print(['PupilAvg_Task_',num2str(task2Anal),'_subType_', subType, '_Avg_IC_CC', '.tif'], '-dtiff', '-r300')
            end
            cd(myHome)
        end
        figure; set(gcf, 'color', [1 1 1])
        data = allSubsPupilStimLocked_IC - allSubsPupilStimLocked_CC;
        dataICMinusCC = data;
        shadedErrorBar(1:length(nanmean(data)), nanmean(data), nansem(data), {'-', 'color', myICColor, 'lineWidth', myLineWidth}, 1); hold on;
         title('IC - CC'); hold on;
        xlim(myPupilPlotXlimVals);
        ylim([-.15 .25]);
        grid on
        ylabel('pupil dilation (z-scored)', 'fontsize', 14);
        set(gca, 'xtick', [min(myPupilPlotXlimVals):1000:max(myPupilPlotXlimVals)], 'xticklabel', [(min(myPupilPlotXlimVals)- 10000)/1000     : (max(myPupilPlotXlimVals)-10000)/1000]);
        box on
        axis square
        set(gca, 'fontsize', 17);
        if saveFigures
            myHome = pwd;
            cd([myMotherPath , filesep, fullfile('Results_EST', 'EYE_ANAL', 'RUN1')]);
            set(gcf, 'paperunits', 'centimeters', 'paperposition', [0 0 15 15])
            if doLowPassFilteringOnly
                print(['PupilAvgLP_Task_',num2str(task2Anal),'_subType_', subType, '_Diff_IC_CC', '.tif'], '-dtiff', '-r300')
            elseif doBandPassFiltering
                print(['PupilAvgBP_Task_',num2str(task2Anal),'_subType_', subType, '_Diff_IC_CC', '.tif'], '-dtiff', '-r300')
            else
                print(['PupilAvg_Task_',num2str(task2Anal),'_subType_', subType, '_Diff_IC_CC', '.tif'], '-dtiff', '-r300')
            end
        cd(myHome)
        end
        figure; set(gcf, 'color', [1 1 1])
        data = allSubsPupilStimLocked_II;
        shadedErrorBar(1:length(nanmean(data)), nanmean(data), nansem(data), {'-', 'color', myIIColor, 'lineWidth', myLineWidth}, 1); hold on;
        title('CI (dark) II (light)'); hold on;
        xlim(myPupilPlotXlimVals);
        set(gca, 'xtick', [min(myPupilPlotXlimVals):1000:max(myPupilPlotXlimVals)], 'xticklabel', [(min(myPupilPlotXlimVals)- 10000)/1000     : (max(myPupilPlotXlimVals)-10000)/1000]);
        data = allSubsPupilStimLocked_CI;
        shadedErrorBar(1:length(nanmean(data)), nanmean(data), nansem(data), {'-', 'color', myCIColor, 'lineWidth', myLineWidth}, 1); hold on;
        xlim(myPupilPlotXlimVals);
        ylim([-.8 .5]);
        ylabel('pupil dilation (z-scored)', 'fontsize', 14);
        set(gca, 'xtick', [min(myPupilPlotXlimVals):1000:max(myPupilPlotXlimVals)], 'xticklabel', [(min(myPupilPlotXlimVals)- 10000)/1000     : (max(myPupilPlotXlimVals)-10000)/1000]);
        axis square
        set(gca, 'fontsize', 17);
        if saveFigures
            myHome = pwd;
            cd([myMotherPath , filesep, fullfile('Results_EST', 'EYE_ANAL', 'RUN1')]);
            set(gcf, 'paperunits', 'centimeters', 'paperposition', [0 0 15 15])
            if doLowPassFilteringOnly
                print(['PupilAvgLP_Task_',num2str(task2Anal),'_subType_', subType, '_Avg_CI_II', '.tif'], '-dtiff', '-r300')
            elseif doBandPassFiltering
                print(['PupilAvgBP_Task_',num2str(task2Anal),'_subType_', subType, '_Avg_CI_II', '.tif'], '-dtiff', '-r300')
            else
                print(['PupilAvg_Task_',num2str(task2Anal),'_subType_', subType, '_Avg_CI_II', '.tif'], '-dtiff', '-r300')
            end
            cd(myHome)
        end
        figure; set(gcf, 'color', [1 1 1])
        data = allSubsPupilStimLocked_CI - allSubsPupilStimLocked_II;
        dataCIMinusII = data;
 [H_T, P_T, CI_T, STATS_T] = ttest(data);
diff_H_T = diff(H_T);
startID = find(diff_H_T == 1);
endID = find(diff_H_T == -1);
startID = startID(:);
endID = endID(:);
if H_T(1) == 1
    if isempty(startID) || (startID(1) > 1)
        startID = [1; startID];
    end
end
if H_T(end) == 1
    if isempty(endID) || (endID(end) < length(H_T))
        endID = [endID; length(H_T)];
    end
end
if length(endID) ~= length(startID)
    warning('Unequal number of start and end IDs. Verifying indices...');
    if length(endID) > length(startID)
        if endID(1) < startID(1)
            startID = [1; startID];
        end
    elseif length(startID) > length(endID)
        if startID(end) > endID(end)
            endID = [endID; length(H_T)];
        end
    end
end
disp('Start IDs:');
disp(startID');
disp('End IDs:');
disp(endID');
        try
            clustersStartEndIDs = [startID, endID];
            foundClustersID_Pre = startID > 6000 & startID < 10000;
            foundClustersID_Post = startID > 10000 & startID < 15000;
            foundClustersStartEndIDs_Pre = clustersStartEndIDs(foundClustersID_Pre,:);
            foundClustersStartEndIDs_Post = clustersStartEndIDs(foundClustersID_Post,:);
            clusterSizesPre = clustersStartEndIDs(foundClustersID_Pre,2) - clustersStartEndIDs(foundClustersID_Pre,1);
            clusterSizesPost = clustersStartEndIDs(foundClustersID_Post,2) - clustersStartEndIDs(foundClustersID_Post,1);
            clusterIntervalToUsePreStim = foundClustersStartEndIDs_Pre(clusterSizesPre == max(clusterSizesPre), : );
            clusterIntervalToUsePostStim = foundClustersStartEndIDs_Post(clusterSizesPost == max(clusterSizesPost), : );
            clusterIntervalToUsePreStim = [6581  8698];
            clusterIntervalToUsePostStim = [10222 11741];
            myLightGrey = [.7 .7 .7];
            myQs = [-.15 .21]; myXs = clusterIntervalToUsePreStim;
            patch([myXs(1), myXs(1), myXs(2), myXs(2), myXs(1)], [myQs(1), myQs(2), myQs(2), myQs(1), myQs(1)], myLightGrey, 'FaceAlpha',.5, 'linestyle', 'none');
            myQs = [-.15 .21]; myXs = clusterIntervalToUsePostStim;
            patch([myXs(1), myXs(1), myXs(2), myXs(2), myXs(1)], [myQs(1), myQs(2), myQs(2), myQs(1), myQs(1)], myLightGrey, 'FaceAlpha',.5, 'linestyle', 'none');
            disp('Plotting patch...');
        catch
        end
         hold on
        shadedErrorBar(1:length(nanmean(data)), nanmean(data), nansem(data), {'-', 'color', myCIColor, 'lineWidth', myLineWidth}, 1); hold on;
       title('CI - II'); hold on;
        xlim(myPupilPlotXlimVals);
        ylabel('pupil dilation (z-scored)', 'fontsize', 14);
        ylim([-.15 .21]);
         grid on
        set(gca, 'xtick', [min(myPupilPlotXlimVals):1000:max(myPupilPlotXlimVals)], 'xticklabel', [(min(myPupilPlotXlimVals)- 10000)/1000     : (max(myPupilPlotXlimVals)-10000)/1000]);
        box on; axis square
        set(gca, 'fontsize', 17);
        if saveFigures
            myHome = pwd;
            cd([myMotherPath , filesep, fullfile('Results_EST', 'EYE_ANAL', 'RUN1')]);
            set(gcf, 'paperunits', 'centimeters', 'paperposition', [0 0 15 15])
            if doLowPassFilteringOnly
                print(['PupilAvgLP_Task_',num2str(task2Anal),'_subType_', subType, '_Diff_CI_II', '.tif'], '-dtiff', '-r300')
            elseif doBandPassFiltering
                print(['PupilAvgBP_Task_',num2str(task2Anal),'_subType_', subType, '_Diff_CI_II', '.tif'], '-dtiff', '-r300')
            else
                print(['PupilAvg_Task_',num2str(task2Anal),'_subType_', subType, '_Diff_CI_II', '.tif'], '-dtiff', '-r300')
            end
            cd(myHome)
        end
         figure; subplot(2,1,1); plot(H_T); ylim([-.2 1.2]); xlim([0 length(H_T)]); subplot(2,1,2); plot(diff(H_T)); ylim([-1.2 1.2]); xlim([0 length(H_T)]);
        diff_H_T = diff(H_T); [startID] = find(diff_H_T == 1); [endID] = find(diff_H_T == -1);
        try
            clustersStartEndIDs = [startID', endID'];
            foundClustersID = startID > 6000 & startID < 10000;
            foundClustersStartEndIDs = clustersStartEndIDs(foundClustersID,:);
            clusterSizes = clustersStartEndIDs(foundClustersID,2) - clustersStartEndIDs(foundClustersID,1);
            clusterIntervalToUse = foundClustersStartEndIDs(clusterSizes == max(clusterSizes), : );
            clusterIntervalToUsePreStim = clusterIntervalToUse;
            PupilCARVect = mean(data(:,clusterIntervalToUse(1):clusterIntervalToUse(2)),2);
            disp(['Correlation beween pre-trial pupil CI>II (sign. cluster) vs. LC CI>II responsivity: p = ', num2str(P(1,2)), ' R = ', num2str(R(1,2))])
            foundClustersID = startID > 10000 & startID < 15000;
            foundClustersStartEndIDs = clustersStartEndIDs(foundClustersID,:);
            clusterSizes = clustersStartEndIDs(foundClustersID,2) - clustersStartEndIDs(foundClustersID,1);
            clusterIntervalToUse = foundClustersStartEndIDs(clusterSizes == max(clusterSizes), : );
            clusterIntervalToUsePostStim = clusterIntervalToUse;
            PupilCARVectPostStim = mean(data(:,clusterIntervalToUse(1):clusterIntervalToUse(2)),2);
            [R,P] = corrcoef(mean(data(:,clusterIntervalToUse(1):clusterIntervalToUse(2)),2), CA);
        catch
        end
 clusterIntervalToUsePreStim = [6581 8698 ];
clusterIntervalToUsePostStim = [10222 11741];
pupilPRE_CI = mean(allSubsPupilStimLocked_CI(:,clusterIntervalToUsePreStim(1):clusterIntervalToUsePreStim(2)),2);
pupilPRE_II = mean(allSubsPupilStimLocked_II(:,clusterIntervalToUsePreStim(1):clusterIntervalToUsePreStim(2)),2);
pupilPRE_IC = mean(allSubsPupilStimLocked_IC(:,clusterIntervalToUsePreStim(1):clusterIntervalToUsePreStim(2)),2);
pupilPRE_CC = mean(allSubsPupilStimLocked_CC(:,clusterIntervalToUsePreStim(1):clusterIntervalToUsePreStim(2)),2);
pupilPOST_CI = mean(allSubsPupilStimLocked_CI(:,clusterIntervalToUsePostStim(1):clusterIntervalToUsePostStim(2)),2);
pupilPOST_II = mean(allSubsPupilStimLocked_II(:,clusterIntervalToUsePostStim(1):clusterIntervalToUsePostStim(2)),2);
pupilPOST_IC = mean(allSubsPupilStimLocked_IC(:,clusterIntervalToUsePostStim(1):clusterIntervalToUsePostStim(2)),2);
pupilPOST_CC = mean(allSubsPupilStimLocked_CC(:,clusterIntervalToUsePostStim(1):clusterIntervalToUsePostStim(2)),2);
pupilPOST_CONG = mean(allSubsPupilStimLocked_CONG(:,clusterIntervalToUsePostStim(1):clusterIntervalToUsePostStim(2)),2);
pupilPOST_INCONG = mean(allSubsPupilStimLocked_INCONG(:,clusterIntervalToUsePostStim(1):clusterIntervalToUsePostStim(2)),2);
theSubs = unique(allSubsAllBev(:,17));
results = [];
LC_act1 = (pupilPOST_CI-pupilPOST_II) - (pupilPRE_CI-pupilPRE_II);
results = table(theSubs, LC_act1, 'VariableNames', {'Subj_num', 'LC_act1'});
 LC_results_path = [myMotherPath, filesep, 'Results_EST', filesep,'EYE_ANAL', filesep, 'RUN1', filesep,  'LC_results_R1_Bandpass.xlsx'];
writetable(results, LC_results_path);
keyboard
    elseif task2Anal == 2
        fh = figure; set(gcf, 'color', [1 1 1])
        allSubsPupilStimLocked_CONG = (allSubsPupilStimLocked_CC + allSubsPupilStimLocked_IC)./2;
        data = allSubsPupilStimLocked_CONG;
        shadedErrorBar(1:length(nanmean(data)), nanmean(data), nansem(data), {'-', 'color', myCongruentColor, 'lineWidth', myLineWidth}, 1); hold on;
       title('INCONGRUENT (Red) CONGRUENT (Green)'); hold on;
        set(gca, 'xtick', [min(myPupilPlotXlimVals):1000:max(myPupilPlotXlimVals)], 'xticklabel', [(min(myPupilPlotXlimVals)- 10000)/1000     : (max(myPupilPlotXlimVals)-10000)/1000]);
        allSubsPupilStimLocked_INCONG = (allSubsPupilStimLocked_CI + allSubsPupilStimLocked_II)./2;
        data = allSubsPupilStimLocked_INCONG;
        shadedErrorBar(1:length(nanmean(data)), nanmean(data), nansem(data), {'-', 'color', myIncongruentColor, 'lineWidth', myLineWidth}, 1); hold on;
        xlim(myPupilPlotXlimVals);
        ylim([-2.1 2.0]);
        set(gca, 'xtick', [min(myPupilPlotXlimVals):1000:max(myPupilPlotXlimVals)], 'xticklabel', [(min(myPupilPlotXlimVals)- 10000)/1000     : (max(myPupilPlotXlimVals)-10000)/1000]);
        axis square
        set(gca, 'fontsize', 17);
        if saveFigures
            myHome = pwd;
            cd([myMotherPath , filesep, fullfile('Results_EST', 'EYE_ANAL', 'RUN1')]);
            set(gcf, 'paperunits', 'centimeters', 'paperposition', [0 0 15 15])
            if doLowPassFilteringOnly
                print(['PupilAvgLP_Task_',num2str(task2Anal),'_subType_', subType, '_Avg_D_B', '.tif'], '-dtiff', '-r300')
            elseif doBandPassFiltering
                print(['PupilAvgBP_Task_',num2str(task2Anal),'_subType_', subType, '_Avg_D_B', '.tif'], '-dtiff', '-r300')
            else
                print(['PupilAvg_Task_',num2str(task2Anal),'_subType_', subType, '_Avg_D_B', '.tif'], '-dtiff', '-r300')
            end
            cd(myHome)
        end
        figure; set(gcf, 'color', [1 1 1])
        data = allSubsPupilStimLocked_CC;
        shadedErrorBar(1:length(nanmean(data)), nanmean(data), nansem(data), {'-', 'color', myCCColor, 'lineWidth', myLineWidth}, 1); hold on;
        set(gca, 'xtick', [min(myPupilPlotXlimVals):1000:max(myPupilPlotXlimVals)], 'xticklabel', [(min(myPupilPlotXlimVals)- 10000)/1000     : (max(myPupilPlotXlimVals)-10000)/1000]);
        data = allSubsPupilStimLocked_IC;
        shadedErrorBar(1:length(nanmean(data)), nanmean(data), nansem(data), {'-', 'color', myICColor, 'lineWidth', myLineWidth}, 1); hold on;
        xlim(myPupilPlotXlimVals);
        ylim([-1.5 1.5]);
        set(gca, 'xtick', [min(myPupilPlotXlimVals):1000:max(myPupilPlotXlimVals)], 'xticklabel', [(min(myPupilPlotXlimVals)- 10000)/1000     : (max(myPupilPlotXlimVals)-10000)/1000]);
        axis square
        set(gca, 'fontsize', 17);
        data = allSubsPupilStimLocked_II;
        shadedErrorBar(1:length(nanmean(data)), nanmean(data), nansem(data), {'-', 'color', myIIColor, 'lineWidth', myLineWidth}, 1); hold on;
        set(gca, 'xtick', [min(myPupilPlotXlimVals):1000:max(myPupilPlotXlimVals)], 'xticklabel', [(min(myPupilPlotXlimVals)- 10000)/1000     : (max(myPupilPlotXlimVals)-10000)/1000]);
        data = allSubsPupilStimLocked_CI;
        shadedErrorBar(1:length(nanmean(data)), nanmean(data), nansem(data), {'-', 'color', myCIColor, 'lineWidth', myLineWidth}, 1); hold on;
        xlim(myPupilPlotXlimVals);
        ylim([-2.1 2.0]);
        set(gca, 'xtick', [min(myPupilPlotXlimVals):1000:max(myPupilPlotXlimVals)], 'xticklabel', [(min(myPupilPlotXlimVals)- 10000)/1000     : (max(myPupilPlotXlimVals)-10000)/1000]);
        axis square
        set(gca, 'fontsize', 17);
        if saveFigures
            myHome = pwd;
            cd([myMotherPath , filesep, fullfile('Results_EST', 'EYE_ANAL', 'RUN1')]);
            set(gcf, 'paperunits', 'centimeters', 'paperposition', [0 0 15 15])
            if doLowPassFilteringOnly
                print(['PupilAvgLP_Task_',num2str(task2Anal),'_subType_', subType, '_Avg_D2_B2', '.tif'], '-dtiff', '-r300')
            elseif doBandPassFiltering
                print(['PupilAvgBP_Task_',num2str(task2Anal),'_subType_', subType, '_Avg_D2_B2', '.tif'], '-dtiff', '-r300')
            else
                print(['PupilAvg_Task_',num2str(task2Anal),'_subType_', subType, '_Avg_D2_B2', '.tif'], '-dtiff', '-r300')
            end
            cd(myHome)
        end
    elseif task2Anal == 3
        fh = figure; set(gcf, 'color', [1 1 1])
        data = allSubsPupilStimLocked_CC;
        shadedErrorBar(1:length(nanmean(data)), nanmean(data), nansem(data), {'-', 'color', myCCColor, 'lineWidth', myLineWidth}, 1); hold on;
        xlim(myPupilPlotXlimVals); ylim([-0.4 1]);
        set(gca, 'xtick', [min(myPupilPlotXlimVals):1000:max(myPupilPlotXlimVals)], 'xticklabel', [(min(myPupilPlotXlimVals)- 10000)/1000     : (max(myPupilPlotXlimVals)-10000)/1000]);
        data = allSubsPupilStimLocked_IC;
        shadedErrorBar(1:length(nanmean(data)), nanmean(data), nansem(data), {'-', 'color', myICColor, 'lineWidth', myLineWidth}, 1); hold on;
        xlim(myPupilPlotXlimVals); ylim([-0.4 1]);
        set(gca, 'xtick', [min(myPupilPlotXlimVals):1000:max(myPupilPlotXlimVals)], 'xticklabel', [(min(myPupilPlotXlimVals)- 10000)/1000     : (max(myPupilPlotXlimVals)-10000)/1000]);
        data = allSubsPupilStimLocked_CI;
        shadedErrorBar(1:length(nanmean(data)), nanmean(data), nansem(data), {'-', 'color', myCIColor, 'lineWidth', myLineWidth}, 1); hold on;
        grid on;
        xlim(myPupilPlotXlimVals); ylim([-.4 .9]);
        set(gca, 'xtick', [min(myPupilPlotXlimVals):1000:max(myPupilPlotXlimVals)], 'xticklabel', [(min(myPupilPlotXlimVals)- 10000)/1000     : (max(myPupilPlotXlimVals)-10000)/1000]);
        axis square
        set(gca, 'fontsize', 17);
        if saveFigures
            myHome = pwd;
            cd([myMotherPath , filesep, fullfile('Results_EST', 'EYE_ANAL', 'RUN1')]);
            set(gcf, 'paperunits', 'centimeters', 'paperposition', [0 0 15 15])
            if doLowPassFilteringOnly
                print(['PupilAvgLP_Task_',num2str(task2Anal),'_subType_', subType ,'.tif'], '-dtiff', '-r300')
            elseif doBandPassFiltering
                print(['PupilAvgBP_Task_',num2str(task2Anal),'_subType_', subType ,'.tif'], '-dtiff', '-r300')
            else
                print(['PupilAvg_Task_',num2str(task2Anal),'_subType_', subType ,'.tif'], '-dtiff', '-r300')
            end
            cd(myHome)
        end
    end
end
if task2Anal == 1
    allUniqueSubs = unique(allSubsAllBev(:,17));
    nrOfAllCurrSubs = length(unique(allSubsAllBev(:,17)));
    ICDiffsSubs = [];
    TWC_RandP = [];
    for currSub = 1 : nrOfAllCurrSubs
        currSubID = allUniqueSubs(currSub);
        currSubMat = allSubsAllBev(allSubsAllBev(:,17)==currSubID, :);
        currSubMat(currSubMat(:,2) == 99, : ) = [];
        currSubMat_I = currSubMat(currSubMat(:,3)==2,:);
        currSubMat_C = currSubMat(currSubMat(:,3)==1,:);
        RTDiff = nanmean(currSubMat_I(:,5)) - nanmean(currSubMat_C(:,5));
        PARDiff = nanmean(currSubMat_I(:,20)) - nanmean(currSubMat_C(:,20));
        RodDiff = nanmean(currSubMat_I(:,21)) - nanmean(currSubMat_C(:,21));
        AccDiff = (nansum(currSubMat_I(:,2)))/length(currSubMat_I(:,2)) - (nansum(currSubMat_C(:,2)))/length(currSubMat_C(:,2));
        PCUDiff = nanmean(currSubMat_I(:,22)) - nanmean(currSubMat_C(:,22));
        PTRDiff = nanmean(currSubMat_I(:,23)) - nanmean(currSubMat_C(:,23));
        PAADiff = nanmean(currSubMat_I(:,24)) - nanmean(currSubMat_C(:,24));
        RTMean = nanmean(currSubMat(:,5));
        PARMean = nanmean(currSubMat(:,20));
        RodMean =  nanmean(currSubMat(:,21));
        AccMean =  nansum(currSubMat(:,2))/length(currSubMat(:,2));
        PCUMean =  nanmean(currSubMat(:,22));
        PTRMean =  nanmean(currSubMat(:,23));
        PAAMean =  nanmean(currSubMat(:,24));
cd([myMotherPath, filesep, 'Results_EST', filesep, 'EYE_ANAL', filesep, 'RUN1']);
ICDiffsSubs = cat(1, ICDiffsSubs, [RTDiff, PARDiff, RodDiff, AccDiff, RTMean, PARMean, RodMean, AccMean, PCUDiff, PCUMean, currSubID, 1, PTRDiff, PTRMean, PAADiff, PAAMean]);
if strcmp(subType, 'CONTROL')
    ICDiffsSubs = cat(1, ICDiffsSubs, [RTDiff, PARDiff, RodDiff, AccDiff, RTMean, PARMean, RodMean, AccMean, PCUDiff, PCUMean, currSubID, 1, PTRDiff, PTRMean, PAADiff, PAAMean]);
elseif strcmp(subType, 'TREAT')
    ICDiffsSubs = cat(1, ICDiffsSubs, [RTDiff, PARDiff, RodDiff, AccDiff, RTMean, PARMean, RodMean, AccMean, PCUDiff, PCUMean, currSubID, 2, PTRDiff, PTRMean, PAADiff, PAAMean]);
elseif strcmp(subType, 'OFC')
    ICDiffsSubs = cat(1, ICDiffsSubs, [RTDiff, PARDiff, RodDiff, AccDiff, RTMean, PARMean, RodMean, AccMean, PCUDiff, PCUMean, currSubID, 3, PTRDiff, PTRMean, PAADiff, PAAMean]);
else
    disp('PLEASE NOTE: some other than the main 3 groups has been analyzed');
end
        cd([myMotherPath , filesep, 'Results_EST', filesep, 'EYE_ANAL', filesep, 'RUN1']);
        ICDiffsSubs = cat(1, ICDiffsSubs,[RTDiff, PARDiff, RodDiff, AccDiff, RTMean, PARMean, RodMean, AccMean, PCUDiff, PCUMean, currSubID, 1, PTRDiff, PTRMean, PAADiff, PAAMean]);
        if strcmp(subType, 'CONTROL')
            ICDiffsSubs = cat(1, ICDiffsSubs,[RTDiff, PARDiff, RodDiff, AccDiff, RTMean, PARMean, RodMean, AccMean, PCUDiff, PCUMean, currSubID, 1, PTRDiff, PTRMean, PAADiff, PAAMean]);
            eval(['save bevPupMat_', subType,' ICDiffsSubs']);
        elseif strcmp(subType, 'TREAT')
            ICDiffsSubs = cat(1, ICDiffsSubs,[RTDiff, PARDiff, RodDiff, AccDiff, RTMean, PARMean, RodMean, AccMean, PCUDiff, PCUMean, currSubID, 2, PTRDiff, PTRMean, PAADiff, PAAMean]);
            eval(['save bevPupMat_', subType,' ICDiffsSubs']);
        elseif strcmp(subType, 'OFC')
            ICDiffsSubs = cat(1, ICDiffsSubs,[RTDiff, PARDiff, RodDiff, AccDiff, RTMean, PARMean, RodMean, AccMean, PCUDiff, PCUMean, currSubID, 3, PTRDiff, PTRMean, PAADiff, PAAMean]);
            eval(['save bevPupMat_', subType,' ICDiffsSubs']);
        else
            disp('PLEASE NOTE: some other than the main 3 groups has been analyzed')
        end
        [R,P,RLO,RUP] = corrcoef(currSubMat(:,5), currSubMat(:,20)); TWC_RT_PAR_AllTrials_R = R(1,2); TWC_RT_PAR_AllTrials_P = P(1,2);
        [R,P,RLO,RUP] = corrcoef(currSubMat(:,5), currSubMat(:,21)); TWC_RT_ROD_AllTrials_R = R(1,2); TWC_RT_ROD_AllTrials_P = P(1,2);
        [R,P,RLO,RUP] = corrcoef(currSubMat(:,5), currSubMat(:,22)); TWC_RT_PCU_AllTrials_R = R(1,2); TWC_RT_PCU_AllTrials_P = P(1,2);
        [R,P,RLO,RUP] = corrcoef(currSubMat_I(:,5), currSubMat_I(:,20)); TWC_RT_PAR_ITrials_R = R(1,2); TWC_RT_PAR_ITrials_P = P(1,2);
        [R,P,RLO,RUP] = corrcoef(currSubMat_I(:,5), currSubMat_I(:,21)); TWC_RT_ROD_ITrials_R = R(1,2); TWC_RT_ROD_ITrials_P = P(1,2);
        [R,P,RLO,RUP] = corrcoef(currSubMat_I(:,5), currSubMat_I(:,22)); TWC_RT_PCU_ITrials_R = R(1,2); TWC_RT_PCU_ITrials_P = P(1,2);
        [R,P,RLO,RUP] = corrcoef(currSubMat_C(:,5), currSubMat_C(:,20)); TWC_RT_PAR_CTrials_R = R(1,2); TWC_RT_PAR_CTrials_P = P(1,2);
        [R,P,RLO,RUP] = corrcoef(currSubMat_C(:,5), currSubMat_C(:,21)); TWC_RT_ROD_CTrials_R = R(1,2); TWC_RT_ROD_CTrials_P = P(1,2);
        [R,P,RLO,RUP] = corrcoef(currSubMat_C(:,5), currSubMat_C(:,22)); TWC_RT_PCU_CTrials_R = R(1,2); TWC_RT_PCU_CTrials_P = P(1,2);
        TWC_RandP = cat(1, TWC_RandP,[TWC_RT_PAR_AllTrials_R, TWC_RT_PAR_AllTrials_P, TWC_RT_ROD_AllTrials_R, TWC_RT_ROD_AllTrials_P, TWC_RT_PCU_AllTrials_R, TWC_RT_PCU_AllTrials_P, ...
            TWC_RT_PAR_ITrials_R, TWC_RT_PAR_ITrials_P, TWC_RT_ROD_ITrials_R, TWC_RT_ROD_ITrials_P, TWC_RT_PCU_ITrials_R, TWC_RT_PCU_ITrials_P, ...
            TWC_RT_PAR_CTrials_R, TWC_RT_PAR_CTrials_P, TWC_RT_ROD_CTrials_R, TWC_RT_ROD_CTrials_P, TWC_RT_PCU_CTrials_R, TWC_RT_PCU_CTrials_P]);
    end
   cd([myMotherPath , filesep, 'Results_EST', filesep, 'EYE_ANAL', filesep, 'RUN1']);
     ICDiffsSubs = cat(1, ICDiffsSubs,[RTDiff, PARDiff, RodDiff, AccDiff, RTMean, PARMean, RodMean, AccMean, PCUDiff, PCUMean, currSubID, 1, PTRDiff, PTRMean, PAADiff, PAAMean]);
     save('ICDiffsSubs.mat',  'ICDiffsSubs' );
columnHeaders = {'RTDiff', 'PARDiff', 'RodDiff', 'AccDiff', 'RTMean', 'PARMean', 'RodMean', 'AccMean', ...
                 'PCUDiff', 'PCUMean', 'currSubID', 'Constant', 'PTRDiff', 'PTRMean', 'PAADiff', 'PAAMean'};
dataWithHeaders = [columnHeaders; num2cell(ICDiffsSubs)];
 cd([myMotherPath , filesep, fullfile('Results_EST', 'EYE_ANAL', 'RUN1')]);
  eval(['load allSubsAllPup_LowPass_', subType]);
[uniqueSubs, ~, idxUnique] = unique(ICDiffsSubs(:, 11), 'stable');
uniqueICDiffsSubs = zeros(length(uniqueSubs), size(ICDiffsSubs, 2));
for i = 1:length(uniqueSubs)
    firstOccurrenceIndex = find(idxUnique == i, 1);
    uniqueICDiffsSubs(i, :) = ICDiffsSubs(firstOccurrenceIndex, :);
end
columnHeaders = {'RTDiff', 'PARDiff', 'RodDiff', 'AccDiff', 'RTMean', 'PARMean', 'RodMean', 'AccMean', ...
                 'PCUDiff', 'PCUMean', 'currSubID', 'Constant', 'PTRDiff', 'PTRMean', 'PAADiff', 'PAAMean'};
dataWithHeaders = [columnHeaders; num2cell(uniqueICDiffsSubs)];
dataWithHeaders = dataWithHeaders(~all(cellfun(@isempty, dataWithHeaders), 2), :);
eval(['save bevPupMat_', subType, ' dataWithHeaders']);
dataFileName = ['bevPupMat_', subType, '.mat'];
loadedData = load(dataFileName);
disp(loadedData);
if isfield(loadedData, 'dataWithHeaders')
    dataToExport = loadedData.dataWithHeaders;
else
    error('dataWithHeaders not found in the .mat file.');
end
ICDiffRun1 = [myMotherPath, filesep, 'Results_EST', filesep,'EYE_ANAL', filesep, 'RUN1', filesep,  'ICDiffRun1.xlsx'];
writecell(dataToExport, ICDiffRun1);
disp(['Data saved to Excel file: ', ICDiffRun1]);
    fh77 = figure(77);
    [R,P,RLO,RUP] = corrcoef(ICDiffsSubs(:,1), ICDiffsSubs(:,2));
    disp(['Correlation beween RT and PAR (IC-Diff): p = ', num2str(P(1,2)), ' R = ', num2str(R(1,2))])
    subplot(5,3,1); scatter(ICDiffsSubs(:,1), ICDiffsSubs(:,2), 'markerfacecolor', 'k'); xlabel('RT-Diff (I>C)'); ylabel('PAR-Diff (I>C)');
    h1 = lsline; h1.Color = 'k'; h1.LineWidth = 2;
    title(['p=', num2str(round(P(1,2),3)), ' R=', num2str(round(R(1,2),3))]); axis square; ylim([-0.6, 0.6]);
    fh7 = figure(7);
    scatter(ICDiffsSubs(:,1), ICDiffsSubs(:,2), 'markerfacecolor', 'k');
    xlabel('RT-Diff (I>C)'); ylabel('PAR-Diff (I>C)');
    h1 = lsline; h1.Color = 'k'; h1.LineWidth = 2;
    title(['p=', num2str(round(P(1,2),3)), ' R=', num2str(round(R(1,2),3))]);
    [b,stats] = robustfit(ICDiffsSubs(:,1), ICDiffsSubs(:,2));
    title(['p=', num2str(round(stats.p(2),3)), '(robust)']);
    axis square; ylim([-0.6, 0.6]); xlim([0, 250]);
    set(gca, 'fontsize', 17);
    set(fh7, 'name', [' PAR-Diff (I>C): ', subType ])
    set(fh7, 'position', [0 350 400 400])
    cd([myMotherPath , filesep, fullfile('Results_EST', 'EYE_ANAL', 'RUN1')]);
    set(gcf, 'paperunits', 'centimeters', 'paperposition', [0 0 15 15])
    print(['PAR_Diff_I_h_C_',num2str(task2Anal),'_subType_', subType ,'.tif'], '-dtiff', '-r300')
    figure(77)
    [R,P,RLO,RUP] = corrcoef(ICDiffsSubs(:,1), ICDiffsSubs(:,3));
    disp(['Correlation beween RT and ROD (IC-Diff): p = ', num2str(P(1,2)), ' R = ', num2str(R(1,2))])
    subplot(5,3,2); scatter(ICDiffsSubs(:,1), ICDiffsSubs(:,3), 'markerfacecolor', 'k'); xlabel('RT-Diff (I>C)'); ylabel('ROD-Diff (I>C)');
    h1 = lsline; h1.Color = 'k'; h1.LineWidth = 2;
    title(['p=', num2str(round(P(1,2),3)), ' R=', num2str(round(R(1,2),3))]); axis square; ylim([-0.6, 0.6]);
    [R,P,RLO,RUP] = corrcoef(ICDiffsSubs(:,1), ICDiffsSubs(:,9));
    disp(['Correlation beween RT and PCU (IC-Diff): p = ', num2str(P(1,2)), ' R = ', num2str(R(1,2))])
    subplot(5,3,3); scatter(ICDiffsSubs(:,1), ICDiffsSubs(:,9), 'markerfacecolor', 'k'); xlabel('RT-Diff (I>C)'); ylabel('Pupil Curve-Diff (I>C)');
    h1 = lsline; h1.Color = 'k'; h1.LineWidth = 2;
    title(['p=', num2str(round(P(1,2),3)), ' R=', num2str(round(R(1,2),3))]); axis square; ylim([-0.6, 0.6]);
    fh8 = figure(8);
    scatter(ICDiffsSubs(:,1), ICDiffsSubs(:,9), 'markerfacecolor', 'k');
     xlabel('RT-Diff (I>C)'); ylabel('PCU-Diff (I>C)');
    h1 = lsline; h1.Color = 'k'; h1.LineWidth = 2;
    title(['p=', num2str(round(P(1,2),3)), ' R=', num2str(round(R(1,2),3))]);
    [b,stats] = robustfit(ICDiffsSubs(:,1), ICDiffsSubs(:,9));
    title(['p=', num2str(round(stats.p(2),3)), '(robust)']);
    axis square; ylim([-0.6, 0.6]); xlim([0, 250]);
    set(gca, 'fontsize', 17);
    set(fh8, 'name', [' PCU-Diff (I>C): ', subType ])
    set(fh8, 'position', [400 350 400 400])
    cd([myMotherPath , filesep, fullfile('Results_EST', 'EYE_ANAL', 'RUN1')]);
    set(gcf, 'paperunits', 'centimeters', 'paperposition', [0 0 15 15])
    print(['PCU_Diff_I_h_C_',num2str(task2Anal),'_subType_', subType ,'.tif'], '-dtiff', '-r300')
    fh9 = figure(9);
    scatter(ICDiffsSubs(:,1), ICDiffsSubs(:,13), 'markerfacecolor', 'k');
     xlabel('RT-Diff (I>C)'); ylabel('PAR-Diff (I>C)');
    h1 = lsline; h1.Color = 'k'; h1.LineWidth = 2;
    title(['p=', num2str(round(P(1,2),3)), ' R=', num2str(round(R(1,2),3))]);
    [b,stats] = robustfit(ICDiffsSubs(:,1), ICDiffsSubs(:,13));
    title(['p=', num2str(round(stats.p(2),3)), '(robust)']);
    axis square; ylim([-0.6, 0.6]); xlim([0, 250]);
    set(gca, 'fontsize', 17);
    set(fh9, 'name', [' PTR-Diff (I>C): ', subType ])
    set(fh9, 'position', [800 350 400 400])
    cd([myMotherPath , filesep, fullfile('Results_EST', 'EYE_ANAL', 'RUN1')]);
    set(gcf, 'paperunits', 'centimeters', 'paperposition', [0 0 15 15])
    print(['PTR_Diff_I_h_C_',num2str(task2Anal),'_subType_', subType ,'.tif'], '-dtiff', '-r300')
    fh10 = figure(10);
    scatter(ICDiffsSubs(:,1), ICDiffsSubs(:,15), 'markerfacecolor', 'k');
     xlabel('RT-Diff (I>C)'); ylabel('PAA-Diff (I>C)');
    h1 = lsline; h1.Color = 'k'; h1.LineWidth = 2;
    title(['p=', num2str(round(P(1,2),3)), ' R=', num2str(round(R(1,2),3))]);
    [b,stats] = robustfit(ICDiffsSubs(:,1), ICDiffsSubs(:,15));
    title(['p=', num2str(round(stats.p(2),3)), '(robust)']);
    axis square; ylim([-0.6, 0.6]); xlim([0, 250]);
    set(gca, 'fontsize', 17);
    set(fh10, 'name', [' PAA-Diff (I>C): ', subType ])
    set(fh10, 'position', [1200 350 400 400])
    cd([myMotherPath , filesep, fullfile('Results_EST', 'EYE_ANAL', 'RUN1')]);
    set(gcf, 'paperunits', 'centimeters', 'paperposition', [0 0 15 15])
    print(['PAA_Diff_I_h_C_',num2str(task2Anal),'_subType_', subType ,'.tif'], '-dtiff', '-r300')
    figure(77)
    [R,P,RLO,RUP] = corrcoef(ICDiffsSubs(:,4), ICDiffsSubs(:,2));
    disp(['Correlation beween Accuracy and PAR (IC-Diff): p = ', num2str(P(1,2)), ' R = ', num2str(R(1,2))])
    subplot(5,3,4); scatter(ICDiffsSubs(:,4), ICDiffsSubs(:,2), 'markerfacecolor', 'k'); xlabel('Accuracy-Diff (I>C)'); ylabel('PAR-Diff (I>C)');
    h1 = lsline; h1.Color = 'k'; h1.LineWidth = 2;
    title(['p=', num2str(round(P(1,2),3)), ' R=', num2str(round(R(1,2),3))]); axis square; ylim([-0.6, 0.6]);
    [R,P,RLO,RUP] = corrcoef(ICDiffsSubs(:,4), ICDiffsSubs(:,3));
    disp(['Correlation beween Accuracy and ROD (IC-Diff): p = ', num2str(P(1,2)), ' R = ', num2str(R(1,2))])
    subplot(5,3,5); scatter(ICDiffsSubs(:,4), ICDiffsSubs(:,3), 'markerfacecolor', 'k'); xlabel('Accuracy-Diff (I>C)'); ylabel('ROD-Diff (I>C)');
    h1 = lsline; h1.Color = 'k'; h1.LineWidth = 2;
    title(['p=', num2str(round(P(1,2),3)), ' R=', num2str(round(R(1,2),3))]); axis square; ylim([-0.6, 0.6]);
    [R,P,RLO,RUP] = corrcoef(ICDiffsSubs(:,4), ICDiffsSubs(:,9));
    disp(['Correlation beween Accuracy and PCU (IC-Diff): p = ', num2str(P(1,2)), ' R = ', num2str(R(1,2))])
    subplot(5,3,6); scatter(ICDiffsSubs(:,4), ICDiffsSubs(:,9), 'markerfacecolor', 'k'); xlabel('Accuracy-Diff (I>C)'); ylabel('PCU-Diff (I>C)');
    h1 = lsline; h1.Color = 'k'; h1.LineWidth = 2;
    title(['p=', num2str(round(P(1,2),3)), ' R=', num2str(round(R(1,2),3))]); axis square; ylim([-0.6, 0.6]);
    [R,P,RLO,RUP] = corrcoef(ICDiffsSubs(:,5), ICDiffsSubs(:,6));
    disp(['Correlation beween RT and PAR (Means): p = ', num2str(P(1,2)), ' R = ', num2str(R(1,2))])
    subplot(5,3,7);
    bar(1, mean(ICDiffsSubs(:,6))); hold on; errorbar(1, mean(ICDiffsSubs(:,6)), std(ICDiffsSubs(:,6))/sqrt(length(ICDiffsSubs(:,6))), 'color', 'k'); hold on
    plot(ones(size(ICDiffsSubs(:,6)))+.1, ICDiffsSubs(:,6), 'ko', 'markerfacecolor', 'k'); axis square; box off
    title('PAR (MEAN-all trials)'); box off; axis square; ylim([-0.1, 0.1]);
    [R,P,RLO,RUP] = corrcoef(ICDiffsSubs(:,5), ICDiffsSubs(:,7));
    disp(['Correlation beween RT and ROD (Means): p = ', num2str(P(1,2)), ' R = ', num2str(R(1,2))])
    subplot(5,3,8);
    bar(1, mean(ICDiffsSubs(:,7))); hold on; errorbar(1, mean(ICDiffsSubs(:,7)), std(ICDiffsSubs(:,7))/sqrt(length(ICDiffsSubs(:,7))), 'color', 'k'); hold on
     plot(ones(size(ICDiffsSubs(:,7)))+.1, ICDiffsSubs(:,7), 'ko', 'markerfacecolor', 'k'); axis square; box off
    title('ROD (MEAN-all trials)'); box off; axis square; ylim([0, 0.4]);
    [R,P,RLO,RUP] = corrcoef(ICDiffsSubs(:,5), ICDiffsSubs(:,10));
    disp(['Correlation beween RT and PCU (Means): p = ', num2str(P(1,2)), ' R = ', num2str(R(1,2))])
    subplot(5,3,9);
    bar(1, mean(ICDiffsSubs(:,10))); hold on; errorbar(1, mean(ICDiffsSubs(:,10)), std(ICDiffsSubs(:,10))/sqrt(length(ICDiffsSubs(:,10))), 'color', 'k'); hold on
    plot(ones(size(ICDiffsSubs(:,10)))+.1, ICDiffsSubs(:,10), 'ko', 'markerfacecolor', 'k'); axis square; box off
    title('Pupil Curve (MEAN-all trials)'); box off; axis square; ylim([-0.1, 0.1]);
    [R,P,RLO,RUP] = corrcoef(ICDiffsSubs(:,8), ICDiffsSubs(:,6));
    disp(['Correlation beween Accuracy and PAR (Means): p = ', num2str(P(1,2)), ' R = ', num2str(R(1,2))])
    [R,P,RLO,RUP] = corrcoef(ICDiffsSubs(:,8), ICDiffsSubs(:,7));
    disp(['Correlation beween Accuracy and ROD (Means): p = ', num2str(P(1,2)), ' R = ', num2str(R(1,2))])
    [R,P,RLO,RUP] = corrcoef(ICDiffsSubs(:,8), ICDiffsSubs(:,10));
    disp(['Correlation beween Accuracy and PCU (Means): p = ', num2str(P(1,2)), ' R = ', num2str(R(1,2))])
    TWC_RandP_Mean = mean(TWC_RandP);
    TWC_RandP_SEM = std(TWC_RandP,1)./sqrt(size(TWC_RandP,1));
    subplot(5,3,10);
    bar([1,2,3], TWC_RandP_Mean(1:6:end)); hold on;
    errorbar([1,2,3], TWC_RandP_Mean(1:6:end), TWC_RandP_SEM(1:6:end), 'color', 'k', 'linestyle', 'none'); hold on
    ylabel('CorrCoef(R)'); xlabel('Trialtype'); xticklabels({'A','I','C'})
    title('Trialwise RT vs. PAR'); box off; axis square; ylim([0, 0.4])
    subplot(5,3,11);
    bar([1,2,3], TWC_RandP_Mean(3:6:end)); hold on;
    errorbar([1,2,3], TWC_RandP_Mean(3:6:end), TWC_RandP_SEM(3:6:end), 'color', 'k', 'linestyle', 'none'); hold on
    ylabel('CorrCoef(R)'); xlabel('Trialtype'); xticklabels({'A','I','C'})
    title('Trialwise RT vs. ROD'); box off; axis square; ylim([0, 0.4])
    subplot(5,3,12);
    bar([1,2,3], TWC_RandP_Mean(5:6:end)); hold on;
    errorbar([1,2,3], TWC_RandP_Mean(5:6:end), TWC_RandP_SEM(5:6:end), 'color', 'k', 'linestyle', 'none'); hold on
    ylabel('CorrCoef(R)'); xlabel('Trialtype'); xticklabels({'A','I','C'})
    title('Trialwise RT vs. PCU'); box off; axis square; ylim([0, 0.4])
    subplot(5,3,13);
    bar([1,2,3], TWC_RandP_Mean(2:6:end)); hold on;
    errorbar([1,2,3], TWC_RandP_Mean(2:6:end), TWC_RandP_SEM(2:6:end), 'color', 'k', 'linestyle', 'none'); hold on
    ylabel('P-value'); xlabel('Trialtype'); xticklabels({'A','I','C'});refline(0, 0.05);
    title('Trialwise RT vs. PAR'); box off; axis square; ylim([0, 0.4])
    subplot(5,3,14);
    bar([1,2,3], TWC_RandP_Mean(4:6:end)); hold on;
    errorbar([1,2,3], TWC_RandP_Mean(4:6:end), TWC_RandP_SEM(4:6:end), 'color', 'k', 'linestyle', 'none'); hold on
    ylabel('P-value'); xlabel('Trialtype'); xticklabels({'A','I','C'}); refline(0, 0.05);
    title('Trialwise RT vs. ROD'); box off; axis square; ylim([0, 0.4])
    subplot(5,3,15);
    bar([1,2,3], TWC_RandP_Mean(6:6:end)); hold on;
    errorbar([1,2,3], TWC_RandP_Mean(6:6:end), TWC_RandP_SEM(6:6:end), 'color', 'k', 'linestyle', 'none'); hold on
    ylabel('P-value'); xlabel('Trialtype'); xticklabels({'A','I','C'}); refline(0, 0.05);
    title('Trialwise RT vs. PCU'); box off; axis square; ylim([0, 0.4])
    set(fh77, 'name', [' Particpant-group: ', subType ])
    set(fh77, 'position', [900 50 600 680])
end
if task2Anal == 1
    fh = figure; set(gcf, 'color', [1 1 1])
    data = allSubsPupilStimLocked_INCONG - allSubsPupilStimLocked_CONG;
    line([round(length(data)/2), round(length(data)/2)], [min(nanmean(data))-(max(nanmean(data))-min(nanmean(data)))*0.3, max(nanmean(data))+(max(nanmean(data))-min(nanmean(data)))*0.3], 'color', myEconStimOnColor ); hold on;
    shadedErrorBar(1:length(nanmean(data)), nanmean(data), nansem(data), {'-', 'color', myCIColor, 'lineWidth', myLineWidth}, 1); hold on;
    title(['Conflict > No-Conflict (', subType ')']); hold on;
    xlim(myPupilPlotXlimVals); ylim([-.2 .3]); grid on
    set(gca, 'xtick', [min(myPupilPlotXlimVals):1000:max(myPupilPlotXlimVals)], 'xticklabel', [(min(myPupilPlotXlimVals)- 10000)/1000     : (max(myPupilPlotXlimVals)-10000)/1000]);
    box on; axis square
    [H_T, P_T, CI_T, STATS_T] = ttest(data);
    try plot(find(P_T < myPThresh), 0.01, 'o', 'color', myTtestStarColor, 'markerfacecolor', myTtestStarColor); catch; end
    set(fh, 'name', [' Particpant-group: ', subType ])
    set(fh, 'position', [0 350 400 400])
    cd([myMotherPath , filesep, fullfile('Results_EST', 'EYE_ANAL'), filesep, 'RUN1']);
    set(gcf, 'paperunits', 'centimeters', 'paperposition', [0 0 15 15])
    print(['PupilDiff_IhC_Task_',num2str(task2Anal),'_subType_', subType ,'.tif'], '-dtiff', '-r300')
    fh = figure; set(gcf, 'color', [1 1 1])
    data = allSubsPupilStimLocked_IC - allSubsPupilStimLocked_CC;
    line([round(length(data)/2), round(length(data)/2)], [min(nanmean(data))-(max(nanmean(data))-min(nanmean(data)))*0.3, max(nanmean(data))+(max(nanmean(data))-min(nanmean(data)))*0.3], 'color', myEconStimOnColor ); hold on;
    shadedErrorBar(1:length(nanmean(data)), nanmean(data), nansem(data), {'-', 'color', myCCColor, 'lineWidth', myLineWidth}, 1); hold on;
    title(['IC > CC (', subType ')']); hold on;
    xlim(myPupilPlotXlimVals); ylim([-.2 .3]); grid on
    set(gca, 'xtick', [min(myPupilPlotXlimVals):1000:max(myPupilPlotXlimVals)], 'xticklabel', [(min(myPupilPlotXlimVals)- 10000)/1000     : (max(myPupilPlotXlimVals)-10000)/1000]);
    box on; axis square
    [H_T, P_T, CI_T, STATS_T] = ttest(data);
    try plot(find(P_T < myPThresh), 0.01, 'o', 'color', myTtestStarColor, 'markerfacecolor', myTtestStarColor); catch; end
    set(fh, 'name', [' Particpant-group: ', subType ])
    set(fh, 'position', [360 350 400 400])
    cd([myMotherPath , filesep, fullfile('Results_EST', 'EYE_ANAL', 'RUN1')]);
    set(gcf, 'paperunits', 'centimeters', 'paperposition', [0 0 15 15])
    print(['PupilDiff_IChCC_Task_',num2str(task2Anal),'_subType_', subType ,'.tif'], '-dtiff', '-r300')
    fh = figure; set(gcf, 'color', [1 1 1])
    data = allSubsPupilStimLocked_CI - allSubsPupilStimLocked_II;
    line([round(length(data)/2), round(length(data)/2)], [min(nanmean(data))-(max(nanmean(data))-min(nanmean(data)))*0.3, max(nanmean(data))+(max(nanmean(data))-min(nanmean(data)))*0.3], 'color', myEconStimOnColor ); hold on;
    shadedErrorBar(1:length(nanmean(data)), nanmean(data), nansem(data), {'-', 'color', myIIColor, 'lineWidth', myLineWidth}, 1); hold on;
    title(['Pupil Diff CI > II (', subType ')']); hold on;
    xlim(myPupilPlotXlimVals); ylim([-.2 .3]); grid on
    set(gca, 'xtick', [min(myPupilPlotXlimVals):1000:max(myPupilPlotXlimVals)], 'xticklabel', [(min(myPupilPlotXlimVals)- 10000)/1000     : (max(myPupilPlotXlimVals)-10000)/1000]);
    box on; axis square
    [H_T, P_T, CI_T, STATS_T] = ttest(data);
    try plot(find(P_T < myPThresh), 0.01, 'o', 'color', myTtestStarColor, 'markerfacecolor', myTtestStarColor); catch; end
    set(fh, 'name', [' Particpant-group: ', subType ])
    set(fh, 'position', [730 350 400 400])
    cd([myMotherPath , filesep, fullfile('Results_EST', 'EYE_ANAL', 'RUN1')]);
    set(gcf, 'paperunits', 'centimeters', 'paperposition', [0 0 15 15])
    print(['PupilDiff_CIhII_Task_',num2str(task2Anal),'_subType_', subType ,'.tif'], '-dtiff', '-r300')
    fh = figure; set(gcf, 'color', [1 1 1])
    data = (allSubsPupilStimLocked_CI - allSubsPupilStimLocked_CC) - (allSubsPupilStimLocked_II - allSubsPupilStimLocked_IC);
    line([round(length(data)/2), round(length(data)/2)], [min(nanmean(data))-(max(nanmean(data))-min(nanmean(data)))*0.3, max(nanmean(data))+(max(nanmean(data))-min(nanmean(data)))*0.3], 'color', myEconStimOnColor ); hold on;
    shadedErrorBar(1:length(nanmean(data)), nanmean(data), nansem(data), {'-', 'color', myIIColor, 'lineWidth', myLineWidth}, 1); hold on;
    title(['CI-CC > II-IC (', subType ')']); hold on;
    xlim(myPupilPlotXlimVals); ylim([-.2 .3]); grid on
    set(gca, 'xtick', [min(myPupilPlotXlimVals):1000:max(myPupilPlotXlimVals)], 'xticklabel', [(min(myPupilPlotXlimVals)- 10000)/1000     : (max(myPupilPlotXlimVals)-10000)/1000]);
    box on; axis square
    [H_T, P_T, CI_T, STATS_T] = ttest(data);
    try plot(find(P_T < myPThresh), 0.01, 'o', 'color', myTtestStarColor, 'markerfacecolor', myTtestStarColor); catch; end
    set(fh, 'name', [' Particpant-group: ', subType ])
    set(fh, 'position', [1110 350 400 400])
    cd([myMotherPath , filesep, fullfile('Results_EST', 'EYE_ANAL', 'RUN1')]);
    set(gcf, 'paperunits', 'centimeters', 'paperposition', [0 0 15 15])
    print(['PupilDiff_CIhCChIIhIC_Task_',num2str(task2Anal),'_subType_', subType ,'.tif'], '-dtiff', '-r300')
    cd(myMotherPath)
elseif task2Anal == 2
    fh = figure; set(gcf, 'color', [1 1 1])
    data = allSubsPupilStimLocked_CC - allSubsPupilStimLocked_IC;
    shadedErrorBar(1:length(nanmean(data)), nanmean(data), nansem(data), {'-', 'color', 'g', 'lineWidth', myLineWidth}, 1); hold on;
    title(['DD > D (', subType ')']); hold on;
    xlim(myPupilPlotXlimVals); ylim([-.4 .8]); grid on
    set(gca, 'xtick', [min(myPupilPlotXlimVals):1000:max(myPupilPlotXlimVals)], 'xticklabel', [(min(myPupilPlotXlimVals)- 10000)/1000     : (max(myPupilPlotXlimVals)-10000)/1000]);
    box on; axis square
    [H_T, P_T, CI_T, STATS_T] = ttest(data);
    try plot(find(P_T < myPThresh), 0.01, 'o', 'color', myTtestStarColor, 'markerfacecolor', myTtestStarColor); catch; end
    set(fh, 'name', [' Particpant-group: ', subType ])
    set(fh, 'position', [0 350 400 400])
    cd([myMotherPath , filesep, fullfile('Results_EST', 'EYE_ANAL', 'RUN1')]);
    set(gcf, 'paperunits', 'centimeters', 'paperposition', [0 0 15 15])
    print(['PupilDiff_DDhD_Task_',num2str(task2Anal),'_subType_', subType ,'.tif'], '-dtiff', '-r300')
    fh = figure; set(gcf, 'color', [1 1 1])
    data = allSubsPupilStimLocked_CI - allSubsPupilStimLocked_II;
    line([round(length(data)/2), round(length(data)/2)], [min(nanmean(data))-(max(nanmean(data))-min(nanmean(data)))*0.3, max(nanmean(data))+(max(nanmean(data))-min(nanmean(data)))*0.3], 'color', myEconStimOnColor ); hold on;
    shadedErrorBar(1:length(nanmean(data)), nanmean(data), nansem(data), {'-', 'color', myCIColor, 'lineWidth', myLineWidth}, 1); hold on;
    title(['B > BB (', subType ')']); hold on;
    xlim(myPupilPlotXlimVals); ylim([-.4 .8]); grid on
    set(gca, 'xtick', [min(myPupilPlotXlimVals):1000:max(myPupilPlotXlimVals)], 'xticklabel', [(min(myPupilPlotXlimVals)- 10000)/1000     : (max(myPupilPlotXlimVals)-10000)/1000]);
    box on; axis square
    [H_T, P_T, CI_T, STATS_T] = ttest(data);
    try plot(find(P_T < myPThresh), 0.01, 'o', 'color', myTtestStarColor, 'markerfacecolor', myTtestStarColor); catch; end
    set(fh, 'name', [' Particpant-group: ', subType ])
    set(fh, 'position', [360 350 400 400])
    cd([myMotherPath , filesep, fullfile('Results_EST', 'EYE_ANAL', 'RUN1')]);
    set(gcf, 'paperunits', 'centimeters', 'paperposition', [0 0 15 15])
    print(['PupilDiff_BhBB_',num2str(task2Anal),'_subType_', subType ,'.tif'], '-dtiff', '-r300')
    fh = figure; set(gcf, 'color', [1 1 1])
    data = allSubsPupilStimLocked_CC - allSubsPupilStimLocked_II;
    line([round(length(data)/2), round(length(data)/2)], [min(nanmean(data))-(max(nanmean(data))-min(nanmean(data)))*0.3, max(nanmean(data))+(max(nanmean(data))-min(nanmean(data)))*0.3], 'color', myEconStimOnColor ); hold on;
    shadedErrorBar(1:length(nanmean(data)), nanmean(data), nansem(data), {'-', 'color', myIIColor, 'lineWidth', myLineWidth}, 1); hold on;
    title(['DD > BB (', subType ')']); hold on;
    xlim(myPupilPlotXlimVals); ylim([-.4 2.8]); grid on
    set(gca, 'xtick', [min(myPupilPlotXlimVals):1000:max(myPupilPlotXlimVals)], 'xticklabel', [(min(myPupilPlotXlimVals)- 10000)/1000     : (max(myPupilPlotXlimVals)-10000)/1000]);
    box on; axis square
    [H_T, P_T, CI_T, STATS_T] = ttest(data);
    try plot(find(P_T < myPThresh), 0.01, 'o', 'color', myTtestStarColor, 'markerfacecolor', myTtestStarColor); catch; end
    set(fh, 'name', [' Particpant-group: ', subType ])
    set(fh, 'position', [730 350 400 400])
    cd([myMotherPath , filesep, fullfile('Results_EST', 'EYE_ANAL', 'RUN1')]);
    set(gcf, 'paperunits', 'centimeters', 'paperposition', [0 0 15 15])
    print(['PupilDiff_DDhBB_Task_',num2str(task2Anal),'_subType_', subType ,'.tif'], '-dtiff', '-r300')
    fh = figure; set(gcf, 'color', [1 1 1])
    data = allSubsPupilStimLocked_IC - allSubsPupilStimLocked_CI;
    line([round(length(data)/2), round(length(data)/2)], [min(nanmean(data))-(max(nanmean(data))-min(nanmean(data)))*0.3, max(nanmean(data))+(max(nanmean(data))-min(nanmean(data)))*0.3], 'color', myEconStimOnColor ); hold on;
    shadedErrorBar(1:length(nanmean(data)), nanmean(data), nansem(data), {'-', 'color', myCCColor, 'lineWidth', myLineWidth}, 1); hold on;
    title(['D > B (', subType ')']); hold on;
    xlim(myPupilPlotXlimVals); ylim([-.4 2.8]); grid on
    set(gca, 'xtick', [min(myPupilPlotXlimVals):1000:max(myPupilPlotXlimVals)], 'xticklabel', [(min(myPupilPlotXlimVals)- 10000)/1000     : (max(myPupilPlotXlimVals)-10000)/1000]);
    box on; axis square
    [H_T, P_T, CI_T, STATS_T] = ttest(data);
    try plot(find(P_T < myPThresh), 0.01, 'o', 'color', myTtestStarColor, 'markerfacecolor', myTtestStarColor); catch; end
    set(fh, 'name', [' Particpant-group: ', subType ])
    set(fh, 'position', [1110 350 400 400])
    cd([myMotherPath , filesep, fullfile('Results_EST', 'EYE_ANAL', 'RUN1')]);
    set(gcf, 'paperunits', 'centimeters', 'paperposition', [0 0 15 15])
    print(['PupilDiff_DhB_Task_',num2str(task2Anal),'_subType_', subType ,'.tif'], '-dtiff', '-r300')
    cd(myMotherPath)
elseif task2Anal == 3
    fh = figure; set(gcf, 'color', [1 1 1])
    data = allSubsPupilStimLocked_CC - allSubsPupilStimLocked_CI;
    shadedErrorBar(1:length(nanmean(data)), nanmean(data), nansem(data), {'-', 'color', myCCColor, 'lineWidth', myLineWidth}, 1); hold on;
    title(['Target (15%) > Standard (', subType ')']); hold on;
    xlim(myPupilPlotXlimVals); ylim([-.2 .6]); grid on
    set(gca, 'xtick', [min(myPupilPlotXlimVals):1000:max(myPupilPlotXlimVals)], 'xticklabel', [(min(myPupilPlotXlimVals)- 10000)/1000     : (max(myPupilPlotXlimVals)-10000)/1000]);
    box on; axis square
    [H_T, P_T, CI_T, STATS_T] = ttest(data);
    try plot(find(P_T < myPThresh), 0.01, 'o', 'color', myTtestStarColor, 'markerfacecolor', myTtestStarColor); catch; end
    set(fh, 'name', [' Particpant-group: ', subType ])
    set(fh, 'position', [0 350 400 400])
    cd([myMotherPath , filesep, fullfile('Results_EST', 'EYE_ANAL', 'RUN1')]);
    set(gcf, 'paperunits', 'centimeters', 'paperposition', [0 0 15 15])
    print(['PupilDiff_ThS_Task_',num2str(task2Anal),'_subType_', subType ,'.tif'], '-dtiff', '-r300')
    fh = figure; set(gcf, 'color', [1 1 1])
    data = allSubsPupilStimLocked_IC - allSubsPupilStimLocked_CI;
    shadedErrorBar(1:length(nanmean(data)), nanmean(data), nansem(data), {'-', 'color', myICColor, 'lineWidth', myLineWidth}, 1); hold on;
    title(['Super Target (5%) > Standard (', subType ')']); hold on;
    xlim(myPupilPlotXlimVals); ylim([-.2 .6]); grid on
    set(gca, 'xtick', [min(myPupilPlotXlimVals):1000:max(myPupilPlotXlimVals)], 'xticklabel', [(min(myPupilPlotXlimVals)- 10000)/1000     : (max(myPupilPlotXlimVals)-10000)/1000]);
    box on; axis square
    [H_T, P_T, CI_T, STATS_T] = ttest(data);
    try plot(find(P_T < myPThresh), 0.01, 'o', 'color', myTtestStarColor, 'markerfacecolor', myTtestStarColor); catch; end
    set(fh, 'name', [' Particpant-group: ', subType ])
    set(fh, 'position', [360 350 400 400])
    cd([myMotherPath , filesep, fullfile('Results_EST', 'EYE_ANAL', 'RUN1')]);
    set(gcf, 'paperunits', 'centimeters', 'paperposition', [0 0 15 15])
    print(['PupilDiff_SThS_Task_',num2str(task2Anal),'_subType_', subType ,'.tif'], '-dtiff', '-r300')
    fh = figure; set(gcf, 'color', [1 1 1])
    data = allSubsPupilStimLocked_IC - allSubsPupilStimLocked_CC;
    shadedErrorBar(1:length(nanmean(data)), nanmean(data), nansem(data), {'-', 'color', 'g', 'lineWidth', myLineWidth}, 1); hold on;
    title(['Super Target (5%) > Target (15%)  (', subType ')']); hold on;
    xlim(myPupilPlotXlimVals); ylim([-.2 .6]); grid on
    set(gca, 'xtick', [min(myPupilPlotXlimVals):1000:max(myPupilPlotXlimVals)], 'xticklabel', [(min(myPupilPlotXlimVals)- 10000)/1000     : (max(myPupilPlotXlimVals)-10000)/1000]);
    box on; axis square
    [H_T, P_T, CI_T, STATS_T] = ttest(data);
    try plot(find(P_T < myPThresh), 0.01, 'o', 'color', myTtestStarColor, 'markerfacecolor', myTtestStarColor); catch; end
    set(fh, 'name', [' Particpant-group: ', subType ])
    set(fh, 'position', [730 350 400 400])
    cd([myMotherPath , filesep, fullfile('Results_EST', 'EYE_ANAL', 'RUN1')]);
    set(gcf, 'paperunits', 'centimeters', 'paperposition', [0 0 15 15])
    print(['PupilDiff_SThT_Task_',num2str(task2Anal),'_subType_', subType ,'.tif'], '-dtiff', '-r300')
    fh = figure; set(gcf, 'color', [1 1 1])
    data = allSubsPupilStimLocked_CI;
    shadedErrorBar(1:length(nanmean(data)), nanmean(data), nansem(data), {'-', 'color', myCIColor, 'lineWidth', myLineWidth}, 1); hold on;
    title(['Standard only (', subType ')']); hold on;
    xlim(myPupilPlotXlimVals); ylim([-.2 .6]); grid on
    set(gca, 'xtick', [min(myPupilPlotXlimVals):1000:max(myPupilPlotXlimVals)], 'xticklabel', [(min(myPupilPlotXlimVals)- 10000)/1000     : (max(myPupilPlotXlimVals)-10000)/1000]);
    box on; axis square
    [H_T, P_T, CI_T, STATS_T] = ttest(data);
    try plot(find(P_T < myPThresh), 0.01, 'o', 'color', myTtestStarColor, 'markerfacecolor', myTtestStarColor); catch; end
    set(fh, 'name', [' Particpant-group: ', subType ])
    set(fh, 'position', [1110 350 400 400])
    cd([myMotherPath , filesep, fullfile('Results_EST', 'EYE_ANAL', 'RUN1')]);
    set(gcf, 'paperunits', 'centimeters', 'paperposition', [0 0 15 15])
    print(['PupilDiff_SONly_Task_',num2str(task2Anal),'_subType_', subType ,'.tif'], '-dtiff', '-r300')
    cd(myMotherPath)
else
    disp('DANGER: mis-specified task')
end
