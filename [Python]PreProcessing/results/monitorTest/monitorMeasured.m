clear all;
close all;

addpath(genpath('../toolBox'));

Screen('Preference', 'SkipSyncTests', 1);

%%

SCREEN_YELLOW = 0;
SCREEN_BLUE = 1;

KbName('UnifyKeyNames');
key = [];
key.KEY_RETURN = KbName('a'); % for mac
key.KEY_ESCAPE = KbName('q'); % for mac

fileName = dir('../../../LEDcubeSimulation/tmp_test/*.csv');
for f = fileName'
    delete([f.folder '/' f.name])
end

cfg = [];
cfg.LUMINANCE_BACKGROUND = 0;
cfg.LUMINANCE_TEXT = ones(1,3)*80./255;

%% display setting
stereoMode = 10;
PsychDefaultSetup(2);

if ~IsWin
    % Assign left-eye view (the master window) to main display:
    scrnNum = 0;
    slaveScreen = 2;
else
    % Assign left-eye view (the master window) to main display:
    scrnNum = 1;
    slaveScreen = 2;
end

PsychImaging('PrepareConfiguration');
PsychImaging('AddTask', 'General', 'DualWindowStereo', slaveScreen);

bgColor=0;
[windowPtr, windowRect] = PsychImaging('OpenWindow', scrnNum, bgColor, [], [], [], stereoMode);

Screen('Flip', windowPtr);

while 1
    [ keyIsDown, seconds, keyCode ] = KbCheck;

    if keyCode(key.KEY_RETURN)
        break;
    end

    if keyCode(key.KEY_ESCAPE)
        Screen('CloseAll');
        Screen('ClearAll');
        ListenChar(0);
        return
    end
end

pause(1);

rect = [0 0 1920 1080];
iSide = 0;

iLight=1;

% light.light1 = [68 56 105 100 0 0];
% light.light2 = [68 43 82 100 0 82];
% light.light3 = [68 12 38 100 0 150];

% light.light3 = [95 0 1 21 0 75];

light.ipRGC_corrected = [89 33 0 30 0 70];
% light.light3 = [50 50 50 50 50 50];

lightName = fieldnames(light);

for mmName = lightName'
    for iProj = [SCREEN_YELLOW SCREEN_BLUE]


        [window_light.(mmName{1}).(['proj' num2str(iProj)]),screenRect] = Screen('OpenOffscreenWindow', scrnNum, cfg.LUMINANCE_BACKGROUND,[],[],32);

        if iProj == SCREEN_YELLOW
            Screen('FillRect', windowPtr, light.(mmName{1})(1:3)./255, windowRect);
        else
            Screen('FillRect', windowPtr, light.(mmName{1})(4:6)./255, windowRect);
        end

        Screen('CopyWindow',windowPtr, window_light.(mmName{1}).(['proj' num2str(iProj)]));

    end
end

for mmName = lightName'
    for iProj = [SCREEN_YELLOW SCREEN_BLUE]
        Screen('SelectStereoDrawBuffer', windowPtr, iProj);
        Screen('CopyWindow',window_light.(mmName{1}).(['proj' num2str(iProj)]), windowPtr);
    end

    Screen('DrawingFinished', windowPtr);
    Screen('Flip', windowPtr);

    while 1
        [ keyIsDown, seconds, keyCode ] = KbCheck;

        if keyCode(key.KEY_RETURN)
            break;
        end

        if keyCode(key.KEY_ESCAPE)
            Screen('CloseAll');
            Screen('ClearAl l');
            ListenChar(0);
            return
        end
    end
    pause(0.5)
end

sca;
ListenChar(0);
%%
% for iLight = 1:12
% 
%     disp(['Selected=' num2str(iLight)])
% 
%     load(['../../toolbox/light/lightIndData_Light' num2str(iLight) '.mat'])
%     lightName =  fieldnames(lightIndData);
% 
% 
%     disp(['Turned on'])
%     showLight();
% 
%     %     fileName = dir('../../../../LEDcubeSimulation/tmp_test/*.csv');
%     %
%     %     fName = ['../../../../LEDcubeSimulation/tmp_test/' num2str(iLight) '_' num2str(iPattern)];
%     %     mkdir(fName)
%     %     for f = fileName'
%     %         movefile([f.folder '/' f.name],fName)
%     %     end
%     %
%     %     for f = fileName'
%     %         delete([f.folder '/' f.name])
%     %     end
% end
% 
% sca;
% ListenChar(0);

%% audio stimulus settings

% data = [];
%
% fileName = '../../../../LEDcubeSimulation/data_LEDCube_20230711x40y25ipRGC109v2.json';
% str = fileread(fileName);
% data.control = jsondecode(str);
%
% fileName = '../../../../LEDcubeSimulation/data_LEDCube_20230711x40y25ipRGC186v2.json';
% str = fileread(fileName);
% data.ipRGC = jsondecode(str);
%
%
%
% for iSub = 1:length(dir('../../0_stimTest/results/s*'))
%
%     if iSub < 10
%         subName = ['s0' num2str(iSub)];
%     else
%         subName = ['s' num2str(iSub)];
%     end
%
%     fileName = ['../../0_stimTest/results/' subName '/stealthPupilBehavior*.mat'];
%     fileName = dir(fileName);
%
%     if ~isempty(fileName)
%
%         tmp_coeff = [];
%         for i = 1:length(fileName)
%             load([fileName(i).folder '/' fileName(i).name], 'cfg', 'param')
%
%             ind = find(param.ans == 1);
%
%             sz = [];
%             for j = 1:size(data.control.Yxy,1)
%                 sz = [sz; length(find(cfg.condition_frame.Color(ind) == j))];
%             end
%             tmp_coeff = [tmp_coeff sz];
%         end
%         [a I] = max(mean(tmp_coeff,2));
%
%         %     load([fileName(end).folder '/' fileName(end).name], 'param');
%         controlNum = I;
%     else
%         controlNum = 1;
%     end
%
%
%     fileName = ['../../1_Brightness/results/' subName '/stealthPupilBehavior*.mat'];
%     fileName = dir(fileName);
%
%     if ~isempty(fileName)
%         load([fileName(end).folder '/' fileName(end).name], 'cfg', 'param')
%     end
%
%     disp([subName ',selected=' num2str(controlNum) ',coeff=' num2str(param.coeff_adjusted)])
%
%
%     load(['../../toolbox/light/lightIndData_Light' num2str(controlNum) '.mat'])
%     light_name =  fieldnames(lightIndData);
%
% %     for lightName = fieldnames(lightIndData)'
%
%     showLight();
%     %     end
%     fileName = dir('../../../../LEDcubeSimulation/tmp_test/*.csv');
%
%     fNaame = ['../../../../LEDcubeSimulation/tmp_test/' subName];
%     mkdir(fName)
%     for f = fileName'
%         movefile([f.folder '/' f.name],fName)
%     end
% a
%     for f = fileName'
%         delete([f.folder '/' f.name])
%     end
%
%
% end
%
%
% sca;
% ListenChar(0);