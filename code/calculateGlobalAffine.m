function [WarpedFrame, WarpedMask, WarpedMaskOutline, WarpedLocalWindows] = calculateGlobalAffine(IMG1,IMG2,Mask,Windows)
% CALCULATEGLOBALAFFINE: finds affine transform between two frames, and applies it to frame1, the mask, and local windows.
    % convert rgb to gray
    gray1 = (rgb2gray(IMG1));
    gray2 = (rgb2gray(IMG2));
    
    % match features
    [f1,valid1] = extractFeatures(gray1,detectSURFFeatures(gray1,'MetricThreshold',400));
    [f2,valid2] = extractFeatures(gray2,detectSURFFeatures(gray2,'MetricThreshold',400));
    indicies = matchFeatures(f1,f2);

    matchedPoints1 = valid1(indicies(:,1),:);
    matchedPoints2 = valid2(indicies(:,2),:);

    % transformation
    H = estimateGeometricTransform(matchedPoints2,matchedPoints1,'affine');

    figure; 
    showMatchedFeatures(gray1,gray2,matchedPoints1,matchedPoints2);

    % new size
    outputView = imref2d(size(IMG2));
    
    % apply transformation
    WarpedFrame = imwarp(IMG1, invert(H),'OutputView',outputView);
    WarpedMask = imwarp(Mask, invert(H),'OutputView',outputView);
    
    % move mask outline
    WarpedMaskOutline = bwperim(WarpedMask,4);
    % move windows
    WarpedLocalWindows = round(transformPointsForward(invert(H),Windows));
end

