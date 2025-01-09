clear all;
close all;

addpath(genpath('../toolBox'));

Screen('Preference', 'SkipSyncTests', 1);

%%

SCREEN_YELLOW = 1;
SCREEN_BLUE = 0;

KbName('UnifyKeyNames');
key = [];
key.KEY_RETURN = KbName('a'); % for mac
key.KEY_ESCAPE = KbName('q'); % for mac

fileName = dir('../../../../LEDcubeSimulation/tmp_test/*.csv');
for f = fileName'
    delete([f.folder '/' f.name])
end

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

cfg.VISUAL_DISTANCE = 100;
cfg.DOT_PITCH = 0.219; % 19inch (1920x1080)

lineWidth = round(pixel_size(cfg.DOT_PITCH, 1, cfg.VISUAL_DISTANCE));

%% audio stimulus settings

data = [];

fileName = '../../../../LEDcubeSimulation/data_LEDCube_20230711x40y25ipRGC109v2.json';
str = fileread(fileName);
data.control = jsondecode(str);

fileName = '../../../../LEDcubeSimulation/data_LEDCube_20230711x40y25ipRGC186v2.json';
str = fileread(fileName);
data.ipRGC = jsondecode(str);

for iLight = 1:12

    load(['../../toolbox/light/lightIndData_Light' num2str(iLight) '.mat'])
    light_name =  fieldnames(lightIndData);

    for iPattern = 1:1
        showLight();
        fileName = dir('../../../../LEDcubeSimulation/tmp_test/*.csv');

        fName = ['../../../../LEDcubeSimulation/tmp_test/' num2str(iLight) '_' num2str(iPattern)];
        mkdir(fName)
        for f = fileName'
            movefile([f.folder '/' f.name],fName)
        end

        for f = fileName'
            delete([f.folder '/' f.name])
        end


    end
end


sca;
ListenChar(0);