function ShapeConfidences = initShapeConfidences(LocalWindows, ColorConfidences, WindowWidth, SigmaMin, A, fcutoff, R, MaskOutline)
% INITSHAPECONFIDENCES Initialize shape confidences.  ShapeConfidences is a struct you should define yourself.
    d = bwdist(MaskOutline);
    fs = [];
    
    for k=1:length(LocalWindows)
        fc = ColorConfidences{k};
        
        xcounter = 1;
        
        % Used -1 in loop to exclude 41st, not 100% sure if correctly done
        for x=LocalWindows(k,1)-(WindowWidth/2):LocalWindows(k,1)+(WindowWidth/2-1)
            ycounter = 1;
            for y=LocalWindows(k,2)-(WindowWidth/2):LocalWindows(k,2)+(WindowWidth/2-1)
                if fcutoff < fc
                    SigmaS = SigmaMin + A*(fc - fcutoff)^R;
                else
                    SigmaS = SigmaMin;
                end
                
                fs(k,ycounter,xcounter) = 1 - exp(-(d(y,x))^2/SigmaS^2);
                
                ycounter = ycounter + 1;
            end
            xcounter = xcounter + 1;
        end
    end
    ShapeConfidences = fs;
end
