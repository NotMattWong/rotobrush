function [mask, mask_outline, NewLocalWindows, ColorModels, ShapeConfidences] = ...
    updateModels(...
        NewLocalWindows, ...
        CurrentFrame, ...
        mask, ...
        mask_outline, ...
        WindowWidth, ...
        ColorModels, ...
        ShapeConfidences, ...
        ProbMaskThreshold, ...
        fcutoff, ...
        SigmaMin, ...
        R, ...
        A, ...
        BoundaryWidth ...
    )
% UPDATEMODELS: update shape and color models, and apply the result to generate a new mask.
% Feel free to redefine this as several different functions if you prefer.
    
    tmpShapeConf = initShapeConfidences(NewLocalWindows,ColorModels.Confidences,WindowWidth,SimgaMin,A,fcutoff,R,mask_outline);

    % Variables
    tmpNum = zeros(size(mask));
    tmpDenom = zeros(size(mask));
    pF = zeros(size(mask));
    
    % Looping through each window
    for i = 1:length(NewLocalWindows)
 
        % Center of window 
        win_x = NewLocalWindows(i,1);
        win_y = NewLocalWindows(i,2);
        
        % Range of window
        xRange = (win_x-(WindowWidth/2)):(win_x+(WindowWidth/2 - 1));   
        yRange = (win_y-(WindowWidth/2)):(win_y+(WindowWidth/2 - 1)); 

        % Distance for shape and color confidence formulas
        dist = bwdist(mask_outline(yRange,xRange));

%% Updating Color Model

        % Image and mask
        IMG = CurrentFrame(yRange,xRange,:);
        locMask = mask(yRange,xRange);
        fgData = [];
        bgData = [];
        
        for x = 1:length(xRange)
            for y = 1:length(yRange)
                if dist(y,x) < BoundaryWidth
                    continue
                end
                
                if locMask(y,x) == 1 && tmpShapeConf(i,y,x) >= 0.75
                    fgData(end+1,:) = IMG(y,x,:);
                elseif locMask(y,x) == 0 && tmpShapeConf(i,y,x) <= 0.25
                    bgData(end+1,:) = IMG(y,x,:);
                end
            end
        end
         
        if size(fgData,1) > 3
            foregroundGMM = fitgmdist(fgData, 3, 'RegularizationValue', 0.001, 'Options', statset('MaxIter',1500,'TolFun',1e-5));
        else
            foregroundGMM = ColorModels.foreGMM{i};
        end
 
        if size(bgData,1) > 3
            backgroundGMM = fitgmdist(bgData, 3, 'RegularizationValue', 0.001, 'Options', statset('MaxIter',1500,'TolFun',1e-5));   
        else
            backgroundGMM = ColorModels.backGMM{i};
        end
         
        % Reshaping IMG
        datafit = reshape(IMG,WindowWidth^2,3);

        % Calculating probability mask with old models
        pxFold = pdf(ColorModels.foreGMM{i},datafit);
        pxBold = pdf(ColorModels.backGMM{i},datafit);
        valOld = pxFold./(pxFold+pxBold);
        pCXold = reshape(valOld,WindowWidth,WindowWidth);

        % Calculating probability mask with new models
        pxFnew = pdf(foregroundGMM,datafit);
        pxBnew = pdf(backgroundGMM,datafit);
        valNew = pxFnew./(pxFnew+pxBnew);   
        pCXnew = reshape(valNew,WindowWidth,WindowWidth);

        [new,~] = find(pCXnew>ProbMaskThreshold);
        [old,~] = find(pCXold>ProbMaskThreshold);
        
        % If background data changes a lot, number of pixels classified as
        % foreground would be less so update the models and probabilities
        if (length(new)<length(old))
            weight = exp(-(dist.^2)/(WindowWidth*0.5)^2);
            denom = sum(sum(weight));
            numer = sum(sum(abs(locMask - pCXnew).*weight));
            ColorModels.Confidences{i} = 1 - numer/denom;
            ColorModels.foreGMM{i} = foregroundGMM;
            ColorModels.backGMM{i} = backgroundGMM;
            ColorModels.pcs(i,:,:) = pCXnew;
        else
            ColorModels.pcs(i,:,:) = pCXold;
        end
        
%% Updating Shape Model
        
        % Calculating the new shape confidence based on simga_s 
        if fcutoff < ColorModels.Confidences{i}
            SigmaS = SigmaMin + A*(ColorModels.Confidences{i} - fcutoff)^R;
        else
            SigmaS = SigmaMin;
        end

        ShapeConfidences(i,:,:) = 1 - exp(-(dist.^2)/((SigmaS)^2));

%% Merging color and shape confidences
        for x = 1:length(xRange)
            for y =1:length(yRange)
                pFx(i,y,x) = ShapeConfidences(i,y,x) * locMask(y,x) + (1-ShapeConfidences(i,y,x)) * ColorModels.pcs(i,y,x);
            end
        end
        
%% Merging windows
        % Calculating the final foreground probability for all pixels in the image
        for y = 1:length(yRange)
            for x = 1:length(xRange)
                dstFromCenter = 1/(sqrt((yRange(y)-win_y)^2 + (xRange(x)-win_x)^2)+0.1);
                tmpNum(yRange(y),xRange(x)) = tmpNum(yRange(y),xRange(x)) + pFx(i,y,x)*dstFromCenter;
                tmpDenom(yRange(y),xRange(x)) = tmpDenom(yRange(y),xRange(x))+ dstFromCenter;
            end
        end
        
    end 

%% Get final foreground probability mask

    % Setting up pF
    pF = (tmpNum)./(tmpDenom);
    pF(isnan(pF))=0;
    
    % Getting mask from pF if above threshold
    mask = (pF>ProbMaskThreshold);
    imshow(imoverlay(CurrentFrame, boundarymask(mask,8), 'red'));
    set(gcf,'visible','off')

    % Getting mask outline from the new mask generated above
    mask_outline = bwperim(mask,4);
    imshow(mask_outline);

end

