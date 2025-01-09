

% if iSide == -1
%     tmp_rect = [0 0 rect(3)/2 rect(4)];
% %     tmp_rect = [0 0 rect(3)/3 rect(4)];
% elseif iSide == 0
%     tmp_rect = rect;
% %     tmp_rect = [rect(3)/3 0 rect(3)-rect(3)/3 rect(4)];
% elseif iSide == -2
%     tmp_rect = [0 0 rect(3) rect(4)/2];
% elseif iSide == 2
%     tmp_rect = [0 rect(4)/2 rect(3) rect(4)];
% else
% %     tmp_rect = [rect(3)-rect(3)/3 0 rect(3) rect(4)];
%     tmp_rect = [rect(3)/2 0 rect(3) rect(4)];
% end

      

for iProj = [SCREEN_YELLOW SCREEN_BLUE]

    Screen('SelectStereoDrawBuffer', windowPtr, iProj);

        tmp_image = lightIndData.(lightName{1}).(['proj' num2str(iProj)])(:,:,:,10)./255;

%         tmp_image = lightIndData.(lightName{2}).(['proj' num2str(iProj)])./255;
        

%         tmp_image = lightIndData.control.(['proj' num2str(iProj)])(:,:,:,round(param.coeff_adjusted))./255;

%         tmp_image = lightIndData.ipRGC.(['proj' num2str(iProj)])(:,:,:)./255;


%         tmp_iamage = lightIndData.control.(['proj' num2str(iProj)])(:,:,:,10)./255;

%     if iPattern == 1 %## draw ipRGC in the left
%         tmp_image =  lightIndData.ipRGC.(['proj' num2str(iProj)])./255;
%         width = size(tmp_image,1);
%         tmp_image(width/2:end,:,:) = lightIndData.control.(['proj' num2str(iProj)])(width/2:end,:,:,10) ./255;
%     else
%         tmp_image = lightIndData.control.(['proj' num2str(iProj)])(:,:,:,10) ./255;
%         tmp_image(width/2:end,:,:) =  lightIndData.ipRGC.(['proj' num2str(iProj)])(width/2:end,:,:)./255;
%     end

%     tmp_image(width/2-round(lineWidth):width/2+round(lineWidth),:) = 0;
    tmp_image = permute(tmp_image,[2 1 3]);


    gradtexture = Screen('MakeTexture',windowPtr, tmp_image);
    Screen('DrawTexture', windowPtr, gradtexture, [], windowRect);
end

%%
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

%%
% fileName = dir("../../../LEDcubeSimulation/tmp_test/*.csv");
% 
% M = readtable([fileName(end).folder '/' fileName(end).name]);
% 
% targetLight = M.Var2(M.Var1 >= 380 & M.Var1 <= 780);
% refLight = squeeze(sum(data.ipRGC.spectrum(iLight,:,:),2));
% refLight = refLight*iCoeff*0.5;

