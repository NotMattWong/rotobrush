function ColorModels = initColorModels(IMG, Mask, MaskOutline, LocalWindows, BoundaryWidth, WindowWidth)
    % INITIALIZAECOLORMODELS Initialize color models.  ColorModels is a struct you should define yourself.
    %
    % Must define a field ColorModels.Confidences: a cell array of the color confidence map for each local window.
    
    % DEFINE COLORMODELS STRUCT
    
    Confidences = {};
    foreGMM = {};
    backGMM = {};
    d = bwdist(MaskOutline);
    pcs = [];
    
    % Loop through windows
    for k=1:length(LocalWindows)
        F = [];
        B = [];
        
        
        % Loop through window pixels to initialize F and B
        % Used -1 in loop to exclude 41st, not 100% sure if correctly done

        for x=LocalWindows(k,1)-(WindowWidth/2):LocalWindows(k,1)+(WindowWidth/2-1)
            for y=LocalWindows(k,2)-(WindowWidth/2):LocalWindows(k,2)+(WindowWidth/2-1)
               
                if d(y,x) < BoundaryWidth
                    continue
                end

                if Mask(y,x) == 1
                   F(end+1,:) = IMG(y,x,:);
                else
                   B(end+1,:) = IMG(y,x,:);
                end
            end
        end
        
        Fgmm = fitgmdist(F, 3, 'RegularizationValue', .001, 'Options', statset('MaxIter',1500,'TolFun',1e-5));
        Bgmm = fitgmdist(B, 3, 'RegularizationValue', .001, 'Options', statset('MaxIter',1500,'TolFun',1e-5));
        
        sigma_s = WindowWidth/2;
    
        num_sum = 0;
        den_sum = 0;
        % Compute color confidence for each pixel in the window
        xcounter = 1;
        for x=LocalWindows(k,1)-(WindowWidth/2):LocalWindows(k,1)+(WindowWidth/2 - 1)
            ycounter = 1;
            for y=LocalWindows(k,2)-(WindowWidth/2):LocalWindows(k,2)+(WindowWidth/2 - 1)
                pf = pdf(Fgmm,[IMG(y,x,1), IMG(y,x,2), IMG(y,x,3)]);
                bf = pdf(Bgmm,[IMG(y,x,1), IMG(y,x,2), IMG(y,x,3)]);
                pc = pf/(pf+bf);
                pcs(k,ycounter,xcounter) = pc;
                num_sum = num_sum + abs(Mask(y,x)-pc) * exp((-d(y,x)^2)/sigma_s^2);
                den_sum = den_sum + exp((-d(y,x)^2)/sigma_s^2);
                ycounter = ycounter + 1;
            end
            xcounter = xcounter + 1;
        end
        fc = 1 - (num_sum/den_sum);
        Confidences{end+1} = fc;
        foreGMM{end+1} = Fgmm;
        backGMM{end+1} = Bgmm;
    end
    ColorModels.Confidences = Confidences;
    ColorModels.foreGMM = foreGMM;
    ColorModels.backGMM = backGMM;
    ColorModels.pcs = pcs;
end

