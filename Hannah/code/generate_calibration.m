function generate_calibration(path2smr)
    % If calibLeft and calibRight exist in path, do nothing
    % If they don't, then load the .smr file in do
    % for each channel in 5:8
    %   plot each channel vs channel 4
    %   compute the slope and the R^2 value
    % 
    % Keep the slope of the ch with the highest R^2 value (between 5 and 6 
    % on one side and between 7 and 8 on the other). Make 0 slopes for the
    % chanel with the lowest R^2.
    
    if (~ exist(fullfile(path2smr, 'calibLeft')))
        % generate calibLeft
    end
        
    if (~ exists(fullfile(path2smr, 'calibRight')))
        % generate calibLeft
    end
end

function MakeCalibrationFile(dat, ch1, ch2, fullpath)
    rsquare = ones(2,1);
    channels = [ch1, ch2];
    for i=1:2
        ch = channels(i);
        [fitobject, fog] = fit(dat(ch).data, dat(4).data, 'poly1');

        % Hanna's calibration file store results in a structure. I'm
        % copying her structure design (only the part that I need)
        if i==1
            scaleCh1 = fitobject.p1;
        else
            scaleCh2 = fitobject.p1;
        end
        rsquare(i) = fog.rsquare;
    end
    
    if rsquare(1) > rsquare(2)
        scaleCh2=0;
    else
        scaleCh1=0;
    end
    
    % save the calib file
    save(fullpath, 'scaleCh1','scaleCh2')
end
    