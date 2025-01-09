

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

% iLight=1;

%%
% for iSub = dir('../../1_Brightness/results/s*')'
% 
%     disp(['Loading ' iSub.name ])
% 
%     fileName = dir(['../../1_Brightness/results/' iSub.name '/stealthPupilBehavior*.mat']);
% 
%     load([fileName.folder '/' fileName.name], 'cfg', 'param')
%     load(['../../toolbox/light/lightIndData_Light' num2str(param.selected) '.mat'])
% 
% 
%     for iProj = [SCREEN_YELLOW SCREEN_BLUE]
%         Screen('SelectStereoDrawBuffer', windowPtr, iProj);
%         
%         tex = Screen('MakeTexture',windowPtr,lightIndData.control.(['proj' num2str(iProj)])(:,:,:,round(param.coeff_adjusted)) ./255);
%         Screen('DrawTexture', windowPtr, tex, [], windowRect);
% %         Screen('CopyWindow',windowPtr, window_light.(mmName{1}).(['proj' num2str(iProj)]));
% 
%     end
% 
%     Screen('DrawingFinished', windowPtr);
%     Screen('Flip', windowPtr);
% 
%     disp(['Disp: ' iSub.name ])
% 
% % end
% % 
% % 
% % for mmName = lightName'
% %     for iProj = [SCREEN_YELLOW SCREEN_BLUE]
% %         Screen('SelectStereoDrawBuffer', windowPtr, iProj);
% %         Screen('CopyWindow',window_light.(mmName{1}).(['proj' num2str(iProj)]), windowPtr);
% %     end
% % 
% %     Screen('DrawingFinished', windowPtr);
% %     Screen('Flip', windowPtr);
% 
%     while 1
%         [ keyIsDown, seconds, keyCode ] = KbCheck;
% 
%         if keyCode(key.KEY_RETURN)
%             break;
%         end
% 
%         if keyCode(key.KEY_ESCAPE)
%             Screen('CloseAll');
%             Screen('ClearAl l');
%             ListenChar(0);
%             return
%         end
%     end
%     pause(0.5)
% end
% 
% sca;
% ListenChar(0);

%%

for iLight = [3 9]
     disp(['Loading ' num2str(iLight) ])
    load(['../../toolbox/light/lightIndData_Light' num2str(iLight) '.mat'])


    for iProj = [SCREEN_YELLOW SCREEN_BLUE]
        Screen('SelectStereoDrawBuffer', windowPtr, iProj);
        
        tex = Screen('MakeTexture',windowPtr,lightIndData.control.(['proj' num2str(iProj)])(:,:,:,10) ./255 .*1);
        Screen('DrawTexture', windowPtr, tex, [], windowRect);
%         Screen('CopyWindow',windowPtr, window_light.(mmName{1}).(['proj' num2str(iProj)]));

    end

    Screen('DrawingFinished', windowPtr);
    Screen('Flip', windowPtr);

    disp(['Disp: ' num2str(iLight) ])

% end
% 
% 
% for mmName = lightName'
%     for iProj = [SCREEN_YELLOW SCREEN_BLUE]
%         Screen('SelectStereoDrawBuffer', windowPtr, iProj);
%         Screen('CopyWindow',window_light.(mmName{1}).(['proj' num2str(iProj)]), windowPtr);
%     end
% 
%     Screen('DrawingFinished', windowPtr);
%     Screen('Flip', windowPtr);

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
