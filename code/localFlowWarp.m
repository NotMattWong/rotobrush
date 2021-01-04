function [LocalWindows] = localFlowWarp(WarpedPrevFrame, CurrentFrame, LocalWindows, Mask, Width)
% LOCALFLOWWARP Calculate local window movement based on optical flow between frames.
    
    % create flow object and give it two frames
    flowObject = opticalFlowFarneback;
    flow = estimateFlow(flowObject,rgb2gray(WarpedPrevFrame));
    flow = estimateFlow(flowObject,rgb2gray(CurrentFrame));

    % now do each window
    [num, coords] = size(LocalWindows);
    NewLocalWindows = zeros(num, coords);
    NumWindows = length(LocalWindows);
    for i = 1:NumWindows
        % window for X
        centerX = LocalWindows(i,1);
        minX = centerX - (Width / 2);
        maxX = centerX + (Width / 2 - 1);
        X = minX:maxX;
        
        % window for Y
        centerY = LocalWindows(i,2);
        minY = centerY - (Width / 2);
        maxY = centerY + (Width / 2 - 1);
        Y = minY:maxY;
 
        Vx = flow.Vx(Y,X);
        Vy = flow.Vy(Y,X);
        [yf,xf] = find(Mask(Y,X));
        
        size_ = length(yf);
        foreX = [];
        foreY = [];

        for j = 1:size_
            foreX = [foreX Vx(yf(j),xf(j))];
            foreY = [foreY Vy(yf(j),xf(j))];
        end
        
        avgVx = sum(foreX)/length(foreX);
        avgVy = sum(foreY)/length(foreY);
        
        if isnan(avgVx) || isnan(avgVy)
            continue;
        end
        
        LocalWindows(i, 1) = round(LocalWindows(i, 1) + avgVx);
        LocalWindows(i, 2) = round(LocalWindows(i, 2) + avgVy);
    end

end

