classdef dat<handle
    % DAT continuous-time AND discrete event data storage object
    %
    %   PROPERTIES
    %     'data' - matrix of numeric data, 1 row per sample, 1 column per
    %         channel. Can be measurement values (continuous) or event times
    %         (discrete)
    %     'chanlabel' - string labels for each channel
    %     'chanval' - numeric values associated with each
    %         channel (e.g. frequency at each bin for a spectrogram).
    %     'samplerate' - samples/second (continuous), 'event' for event channel
    %     'tstart','tend' - time of first/last sample; in seconds.
    %     'units' - descriptive string giving units of data.
    %     'nbadstart' and 'nbadend' - samples at start/end of data that are
    %         unreliable due to filtering edge effects)
    %
    %   METHODS
    %         obj = dat(data, chanlabel, chanval, samplerate, tstart, tend, units)
    %         chanind = datchanind(obj,chans)
    %         datout = datchan(obj, chans)
    %         data = datchandata(obj,chan)
    %         datout = datseg(obj,seg)
    %         data = datsegdata(obj, tstart, duration)
    %         datout = datsegmeans(obj, tstart, duration, shiftsegs,nonan)
    %         datout = dateventtocont(obj, npoints, tstart, tend)
    %         datrange = datrange(obj,chan)
    %         datout = datsmooth(obj,varargin)
    %         datout = datbin(obj, binwidth)
    %         datout = downsample(obj, varargin)         
    %         datout = upsample(obj, n)
    %         datout = resettime(obj);
    %         datout = dattimeind(obj,t) % returns data point at time t
    %         t = dattime(obj, seg)
    %         h = plot(obj, varargin)
    %               datfft
    %
    %   DEPENDENCIES
    %       tight_subplot.m
    %
    % Author: Hannah Payne 2012
    %
    properties
        data;
        chanlabel;
        chanval;
        samplerate;
        tstart;
        tend;
        units;
        nbadstart=0;
        nbadend=0;
        freq;
        PSD;
    end
    methods
        % Class constructor
        function obj = dat(varargin)
            %   dat(data, chanlabel, chanval, samplerate, tstart, tend,
            %   units)
            %     'data' - matrix of numeric data, 1 row per sample, 1 column per
            %         channel.
            %     'chanlabel' - string labels for each channel
            %     'chanval' - numeric values associated with each
            %         channel (e.g. frequency at each bin for a spectrogram).
            %     'samplerate' - derived, not claimed (samples/second, double-precision)
            %     'tstart','tend' - time of first/last sample; in seconds.
            %     'units' - descriptive string giving units of data.
            %
            %     'nbadstart' and 'nbadend' - (not set by user) samples at start/end of data that are
            %         unreliable due to filtering edge effects)
            %
            %    May be called with any or no input arguments
            %
            
            p = inputParser;
            addOptional(p,'data',[],@(x)isnumeric(x) || islogical(x));
            addOptional(p,'chanlabel',{},@(x) ischar(x) || iscell(x) ||isempty(x));
            addOptional(p,'chanval',[],@isnumeric);
            addOptional(p,'samplerate',1,@(x) isnumeric(x) || strcmp(x,'event'));
            addOptional(p,'tstart',0,@isnumeric);
            addOptional(p,'tend',[],@isnumeric);
            addOptional(p,'units',{},@(x) ischar(x) || iscell(x) ||isempty(x));
            
            parse(p,varargin{:});
            
            data = p.Results.data;
            chanlabels = p.Results.chanlabel;
            chanvals = p.Results.chanval;
            samplerates = p.Results.samplerate;
            tstarts = p.Results.tstart;
            tends = p.Results.tend;
            unitss = p.Results.units;
            
            if nargin>0
                if size(data,1)<size(data,2)
                    data = data';
                end
                if isempty(data)
                    nchans = 1;
                else
                    nchans = size(data,2);
                end
                
                % Preallocate object array
                obj(nchans) = dat;
                
                if ischar(chanlabels)
                    chanlabels = {chanlabels};
                end
                
                
                for i = 1:nchans
                    
                    if length(chanlabels)==nchans
                        ichanlabel = chanlabels{i};
                    else
                        ichanlabel = chanlabels;
                    end
                    if length(chanvals)==nchans
                        ichanval = chanvals(i);
                    else
                        ichanval = i;
                    end
                    if length(samplerates)==nchans
                        isamplerate = samplerates(i);
                    else
                        isamplerate = samplerates;
                    end
                    if length(tstarts)==nchans
                        itstart = tstarts(i);
                    else
                        itstart = tstarts;
                    end
                    if length(tends)==nchans
                        itend = tends(i);
                    else
                        if strcmp(isamplerate,'event')
                            itend = max(data);
                        else
                            itend = itstart + (length(data(:,i))-1)/isamplerate;
                        end
                    end
                    if length(unitss)==nchans
                        iunits = unitss(i);
                    else
                        iunits = unitss;
                    end
                    
                    
                    if ~isempty(data)
                        obj(i).data = data(:,i);
                    end
                    obj(i).chanlabel = ichanlabel;
                    obj(i).chanval = ichanval;
                    obj(i).samplerate = isamplerate;
                    obj(i).tstart = itstart;
                    obj(i).tend = itend;
                    obj(i).units = iunits;
                end
            end
            
        end
        
        % Methods
        function chanind = datchanind(obj,chans)
            % Get the index for the specified channel numbers or labels
            if isnumeric(chans)
                for i = 1:length(chans)
                    if find([obj.chanval]==chans(i))
                        chanind(i) = find([obj.chanval]==chans(i));
                    end
                end
            elseif ischar(chans) % Chanlabel specified (one)
                chanind = find(strcmp({obj.chanlabel},chans));     
            elseif iscell(chans) % Chanlabel specified (many)
                for i = 1:length(chans)
                    if find(strcmp({obj.chanlabel},chans{i}))
                        chanind(i) = find(strcmp({obj.chanlabel},chans{i}));
                    end
                end
            end
            chanind = chanind(chanind>0);
        end
        function datout = datchan(obj, chans)
            % Get a dat structure with just the specificied channels
            chaninds = datchanind(obj, chans);
            datout = obj(chaninds);
        end
        
        function datout = datdelete(obj,chans)
            chaninds = datchanind(obj,chans);
            datout = obj;
            datout(chaninds) = [];
        end
        
        function data = datchandata(obj,chan)
            % data = datchandata(obj,chan)
            % Get just the raw data from one channel
            
            if ~exist('chan','var') || length(obj)==1
                data = obj.data;
            else
                
                chanind = datchanind(obj,chan);
                datout = obj(chanind);
                if ~isempty(datout)
                    data = datout.data;
                else
                    error('No data or channel not found')
%                     data = [];
                end
            end
        end
        function datout = datseg(obj,seg, npoints)
            % datout = datseg(obj,seg)
            % Return dat object with only the given time segment of data
            
            if ~exist('npoints','var')
                npoints = [];
            end
            
            datout = obj;
            
            for i = 1:length(obj)
                
                if strcmp(obj(i).samplerate,'event')
                    mask = obj(i).data>seg(1) & obj(i).data<seg(2);
                    datout(i).data = obj(i).data(mask);
                    datout(i).tstart = seg(1);
                    datout(i).tend = seg(2);
                    
                else
                    
                    if isempty(npoints)
                        n = round((diff(seg)+2*eps)*obj(i).samplerate);
                    else
                        n= npoints;
                    end
                    
%                     startindex = round(obj(i).tstart + obj(i).samplerate*seg(1) + 1);
                    startindex = round((seg(1)-obj(i).tstart)* obj(i).samplerate + 1);
                    
%                     t = dattime(obj(i));
%                     startindex = find(t>=seg(1)-1e-4,1);
                    % Base sampling on total duration - so segments with the same duration
                    % and sample rate have the same length output
                    
                    stopindex = startindex + n-1;
                    
                    tstart = seg(1);
                    tend = seg(2);
                    if startindex < 1 || stopindex > length(obj(i).data)
                        %                         error('Time range exceeds data limits')
                        
                        
                        if length(obj)==1
                            if startindex < 1
                                startindex = 1;
                                tstart = obj.tstart;
                            end
                            if stopindex > length(obj(i).data)
                                stopindex = length(obj(i).data);
                                tend = obj.tend;
                            end
                        else
                            return
                        end
                    end
                    datout(i).data = obj(i).data(startindex:stopindex);
                    datout(i).tstart = tstart;%t(startindex);
                    datout(i).tend = tend;%t(stopindex); 11/1/13
                    
                end
            end
        end
        function data = datsegdata(obj, tstart, duration, alignsegs)
            % data = datsegdata(obj, tstart, duration, (alignsegs))
            %
            % Return the raw segmented data in a mtrix form with one
            % row per sample and one column per segment
            % since all segs must be the same length (for cont data), give tstart and
            % duration as input args
            %
            % For event data, returns a vector of all the event times,
            % aligned to each tstart
            %
            % DURATION can have two elements - [tbefore and tafter]
            
            if ~exist('alignsegs','var')
                alignsegs = 1;
            end
            if length(obj)>1
                error('Only input dat object with one channel - use datChan(obj, chan)')
            end
            
            if length(duration) == 2
                tbefore = duration(1);
                tafter = duration(2);
                duration = tafter - tbefore;
            elseif length(duration)==1
                tbefore = 0;
                tafter = duration;
            else
                tbefore  = 0;
                tafter = duration;
            end
            
            npoints = round( (tafter-tbefore)*obj(1).samplerate);
            
            data = [];
            for i = 1:length(tstart)
                                    
                if isnan(tstart(i))
                    data = [data NaN(npoints,1)];

                    continue
                end
                if obj.tend < tstart(i)+tafter || obj.tstart>tstart(i)+tbefore
                %*** HP 1/10/14
                    warning('Excluding segments that fall outside range')
%                     fprintf('trial = %g s, tbefore = %g, tafter = %g, obj.tstart = %g, obj.tend = %g\n\n',tstart(i), tbefore, tafter,obj.tstart, obj.tend)
                    continue
                end
                if length(tafter)>1
                datTemp = datseg(obj, [tstart(i)+tbefore tstart(i)+tafter(i)], npoints);
                else
                datTemp = datseg(obj, [tstart(i)+tbefore tstart(i)+tafter], npoints);
                end

                if strcmp(obj.samplerate,'event')
                    if alignsegs
                        data = [data; datTemp.data  - tstart(i)];
                    else
                        data = [data; datTemp.data];
                    end
                    else
                    data = [data datTemp.data(:)];
                end
                
            end
            
        end
        
        % inserted by PJ on 2/20/2015
        function ta = time_axis(obj)
            % return a time axis that can be used to plot(ta, obj.data)
            ta = obj.tstart:1/obj.samplerate:obj.tend-1/obj.samplerate;
        end
        
        function subdata = datsplit(obj, split_start_t, duration)
            % generate another dat object with data from obj starting at
            % "split_start_t" and lasting "duration"
            % split_start_t:   in seconds
            % duration: in seconds, if duration ==-1, goes to the end obj
            
            if duration==-1
                duration = obj.tend - split_start_t;
            end
            
            % find points in obj.data associated with tstart, tend
            startP = max(round((split_start_t - obj.tstart)*obj.samplerate), 1);
            endP = round(duration*obj.samplerate + startP - 1);

            % duplicate data before removing regions outsider ROI
            subdata = dat('data',obj.data(startP: endP), 'chanlabel', obj.chanlabel, ...
                'chanval', obj.chanval, 'samplerate', obj.samplerate, ...
                'tstart', split_start_t, 'tend', split_start_t + duration, ...
                'units', obj.units);
        end
        
        function [startT, duration] = detect_high_speed(obj, saccade_time, ...
                threshold, min_fixation_length)
            % identify regions with large velocities and return the startT
            % and duration of all those high velocity regions
            % This function returns two arrays, one with the start and
            % the other with the duration of each FEM segment detected.
            %
            % obj:          Hanna's object
            % saccade_time: How long do we expect a saccade or an artifact to last?
            %               in seconds
            % threshold:    in 0-1 units, actual threshold is computed by multiplying
            %               by the channel's range. If the difference in poisition
            %               at two frames separated by saccade_time/2, saccade_time/3
            %               is more than actual threshold then a high speed event
            %               is detected.
            %               Can be a single number or an array with
            %               individual thresholds (one per channel).
            % min_fixation_length:  minimum amount of quiet frames that I do accept
            %                       in between saccades or artifacts (in seconds)
            
            % just scale the "data range" by threshold
            thresh = range(obj.data) * threshold;
                
            % compute the diff between data separated by saccade_time (and
            % saccade_time/2, saccade_time/4).
            % saccade_diff is saccade_time*samplerate-1 shorter than data
            % A high number in saccade_diff means a saccade will happen shortly
            % or is under way
            saccade_diff = get_saccade_diff(obj.data, saccade_time, obj.samplerate);
                
            % get an array of frames where this channel had a strong velocity
            frames = 1:length(saccade_diff);
            frames = frames(saccade_diff > thresh);
            
            % change from frames to times, taking tstart into account and
            % appending two new 'frames', the beginning and end of the
            % experiment
            times = [obj.tstart obj.tstart + frames/obj.samplerate obj.tend];
            
            % discard regions of recording that are not quiet for at least
            % min_fixation_length. At this point, changing_times has more points
            % than the output.
            fem_length = diff(times);
            
            startT = zeros(length(fem_length),1);
            duration = zeros(length(fem_length),1);

            idx = 1;
            for i = 1:length(fem_length)
                
                % if interruption_distance is larger than min_fixation_length add
                % data from obj to startT, duration
                if fem_length(i) >= min_fixation_length
                    startT(idx) = times(i);
                    duration(idx) = times(i+1) - times(i);
                    idx = idx + 1;
                end
            end
            
            % remove preallocated elements from startT and duration that were not filled
            startT(idx:end) = [];
            duration(idx:end) = [];
        end
        
        function obj = calibrate_data(obj, path)
            % change channel from V to degrees (if units are in "V"olts.
            %
            % obj: the dat file to calibrate
            %
            % it tries to load files 'calibLeft.mat' and 'calibRight.mat'
            % (depending on the value of obj.chanval 
            % (5-6: calibLeft, 7-8: calibRight)
            % from path. If those files exists, it will use variable 
            % scaleCh1/2 to convert from V to deg.
            % calib.scaleCh1 -> channel 5 or 7
            % calib.scaleCh2 -> channel 6 or 8
            %
            % If calibRight and/or Left do not exist it will generate them
            

            if (strcmp(obj.units, 'V'))
                % get the calibration file name depending on ch
                if (obj.chanval > 6) % left eye are channels 5 and 6
                    % right eye are channels 7 and 8
                    calibfile = fullfile(path, 'calibRight.mat');
                else
                    calibfile = fullfile(path, 'calibLeft.mat');
                end
                
                if exist(calibfile,'file')
                    % load calibration data from path
                    calib = load(calibfile);
                else
                    % Calibrate channel based on the recording and
                    % store the calibration in calibfile
                    calib = generate_calibration(obj, calibfile);
                end
                
                if (obj.chanval==5 || obj.chanval==7)
                    % channels 5 and 7 are associated with scaleCh1
                    % of calibLeft and calibRight respectively
                    % and channels 6 and 8 with scaleCh2 of
                    % calibLeft and calibRight respectively
                    obj.data = obj.data * calib.scaleCh1;
                    obj.units = 'deg';
                else
                    obj.data = obj.data * calib.scaleCh2;
                    obj.units = 'deg';
                end
            end
        end
        function estimatePSD(obj)
            % compute power  spectral density estimation
            % I'm following the example in http://www.mathworks.com/help/signal/ug/psd-estimate-using-fft.html
            Fs = obj.samplerate;
            N = size(obj.data,1);
            xdft = fft(obj.data,[],1);
            xdft = xdft(1:N/2+1, :);
            psdx = (1/(Fs*N)) * abs(xdft).^2;
            psdx(2:end-1, :) = 2*psdx(2:end-1, :);
            obj.freq = 0:Fs/N:Fs/2;
            obj.PSD = psdx;
        end
        
        function plotPSD(obj)
            fig = figure();
            set(fig, 'Position', [0 0 1000 300]);
            plot(obj.freq, obj.PSD, 'k')
            xlabel('Freq (Hz)', 'Fontsize', 16);
            ylabel('Normalized power', 'Fontsize', 16);
            ylim([0 10]);
            xlim([0 50]);
            set(gca, 'Fontsize', 16);
            title('Nyquist = 500Hz', 'Fontsize', 16);
        end
        
        function [x, trajectory] = getEyeTrajectory(obj, cutOffFreq, varargin)
            % Filter data upto the cutOffFreq.
            % Returns both x and y such that you can plot(x, trajectory)
            % and overlay this with obj.plotSegment(segment, ch)

            %%%%%%%%%%% start and end are not implemented %%%%%%%%%%%%%%%%%
            %p = inputParser;
            %addOptional(p,'start',obj.tstart,@(x) x >= obj.tstart && x <= obj.tend);
            %addOptional(p,'end',obj.tend,@(x) x >= obj.tstart && x <= obj.tend);
            %parse(p,varargin{:});
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            Fs = obj.samplerate;
            N = size(obj.data,1);
            
            xdft = fft(obj.data,[],1);
            
            % matlab's fft start at freq = 0 and goes all the way to
            % Fs in N steps (each of length Fs/N). I want to convert 
            % cutOffFreq to a point in the freq
            % axis such that P * (Fs/N) = cutOffFreq
            P = cutOffFreq*N/Fs + 1;    % +1 because matlab index starts from 1.
            
            % Matlab's fft is redundant and includes both from freq=0 to
            % Nyqist freq (Fs/2) but also from Fs/2 to Fs (where
            % fft(i+Fs/2) = fft*(Fs/2-i) (second half is complex conjugate
            % of first half). When making zero points in the first half,
            % also make them in the 2nd half
            xdft(P:end-P+1, :) = 0;
            
            trajectory = real(ifft(xdft, [], 1));
            x = 0:size(trajectory,1)-1;
            x = x/Fs;
        end        
        % end of insertion by PJ on 2/20/2015
        
        function [datout, ncycles] = datsegmeans(obj, tstart, duration, shiftsegs,nonan, binwidth)
            % datout = datsegmeans(obj, tstart, duration, shiftsegs,nonan)
            % Returns mean values for all segments in a dat structure
            % Must specify tstart and duration
            % shiftsegs = 1 to shift all data so first point is 0
            % nonan = 1 ignores all segments containing any NaNs
            
            if ~exist('shiftsegs','var')
                shiftsegs = 0;
            end
            
            if ~exist('nonan','var')
                nonan = 0;
            end
            
            if ~exist('binwidth','var')
                binwidth = 1/obj(1).samplerate;
            end
            
            for i = 1:length(obj)
                datout(i) = obj(i);

                if isempty(obj(i).data)
                    continue;
                end
                
                if strcmp(obj(i).samplerate,'event')
                    segData = [];
                    nbins = round(duration/binwidth);
                    for j = 1:length(tstart)
%                         bins = tstart(j)-binwidth/2:binwidth:(tstart(j)+duration+binwidth/2);
                        
                        bins = (0:nbins-1)*binwidth + tstart(j);
                        counts = histc(obj(i).data,bins);
                        segData(:,j) = counts(1:end-1)/binwidth;
                        
                    end
                    datout(i).samplerate = 1/binwidth;
                    
                else
                    segData = datsegdata(obj(i), tstart, duration);

                    if shiftsegs
                        offsets = repmat(segData(1,:),[size(segData,1), 1]);
                        segData = segData-offsets;
                    end
                    if nonan
                        nansegs = find(sum(isnan(segData),1));
                        segData(:,nansegs)=[];
                    end
                end
                
                datout(i).data = nanmean(segData,2);
                ncycles = size(segData,2);
                
                if numel(duration)==1
                    datout(i).tstart = 0;
                    datout(i).tend = duration;
                else
                    datout(i).tstart = duration(1);
                    datout(i).tend = duration(2);
                end
            end
            
        end
        function datout = dateventtocont(obj, npoints, tstart, tend)
            % Make sure units of the event channel are in seconds!
            events = obj.data;
            nevents = length(events);
            
            % Max possible error in sample rate accuracy (samples/sec) to check
            maxerr = 1;
            
            if nargin<3
                currsamplerate = npoints;
                tstart = obj.tstart;
                tend = obj.tend;
            else
                currsamplerate = npoints/(tend-tstart);
            end
            t = (0:1/currsamplerate:(npoints-1)/currsamplerate)+tstart;
            
            inds = ones(1,nevents);
            
            for i = 1:nevents
                
                tindmin = floor((currsamplerate-maxerr)*(events(i) - tstart));
                tindmax = ceil((currsamplerate+maxerr)*(events(i) - tstart));
                
                if events(i)<tstart || events(i)>tend
                    fprintf('Truncating point %g (tstart = %g, tend = %g)', events(i),tstart,tend)
                    continue
                end
                
                if tindmin<1
                    tindmin = 1;
                elseif tindmax > npoints
                    tindmax = npoints;
                end
                
                tguess = t(tindmin:tindmax);
                
                diffs = tguess - events(i);
                
                diffs = abs(diffs);
                
                [~,I] = min(diffs);
                inds(i) = I + tindmin - 1;
            end
            
            currdata = false(npoints, 1);
            currdata(inds) = true;
            
            fprintf('Eventtocont\n Last input time:%.5f\n Last output time: %.5f\n',...
                obj.data(end),t(find(currdata,1,'last')))
            
            
            datout = dat(currdata, obj.chanlabel, obj.chanval, currsamplerate, tstart, tend, obj.units);
            
        end
        function datrange = datrange(obj,chan)
            % datRange Returns range (y-span) of data in specified channel
            %   chan is optional
            
            if exist('chan','var')
                datCurr = datchan(obj, chan);
                if numel(datCurr)~=1
                    error('Only pick one channel.')
                end
            elseif numel(obj)==1
                datCurr = obj;
            else
                error('Only pick one channel.')
            end
            datrange = [min(datCurr.data) max(datCurr.data)];
        end
        function datout = smooth(obj,varargin)
            % Smooths data in the dat structure DATIN
            % DATOUT = datsmooth(DATIN, [SPAN], [CHANS])
            %   DATIN - required input dat structure
            %   SPAN  - optional number of points to smooth
            %   CHANS - optional subset of channels to smooth (values
            %   or cell array of labels)
            
            
            p = inputParser;
            
            addRequired(p,'obj',@(x)isa(x,'dat'));
            addOptional(p,'span',5,@isnumeric);
            addOptional(p,'chans',[obj.chanval],@(x) isnumeric(x) || ischar(x) || iscell(x));
            
            parse(p,obj,varargin{:})
            
            obj = p.Results.obj;
            span = p.Results.span;
            chans = p.Results.chans;
            
            chanInd = datchanind(obj,chans);
            
            datout = obj;
            for i = 1:length(chanInd)
                if ~strcmp(obj(chanInd(i)).samplerate, 'event')
                    datout(chanInd(i)).data = smooth(obj(chanInd(i)).data, span,'moving');
                end
            end
        end
        function datout = datsmooth(obj,varargin)
            % Smooths data in the dat structure DATIN
            % DATOUT = datsmooth(DATIN, [SPAN], [CHANS])
            %   DATIN - required input dat structure
            %   SPAN  - optional number of points to smooth
            %   CHANS - optional subset of channels to smooth (values
            %   or cell array of labels)
            
           datout = smooth(obj,varargin{:});
        end
        function datout = datbin(obj, binwidth)
            % Returns a dat structure with binned data
            % Specify binwidth in seconds
            
            databinned = zeros(1,floor(length(obj.data)/(binwidth*obj.samplerate)));
            
            for i = 1:length(databinned)
                istart = round(1 + (i-1)*binwidth*obj.samplerate);
                istop = round(i*binwidth*obj.samplerate);
                databinned(i) = nanmean(obj.data(istart:istop));
                
            end
            
            datout = obj;
            datout.data = databinned;
            datout.samplerate = 1/binwidth;
        end
        function datout = downsample(obj, varargin)
            % datDownsample Downsample input signal
            %   Downsample by keeping every nth sample starting with the
            %   first
            
            p = inputParser;
            
            addRequired(p,'obj',@(x)isa(x,'dat'));
            addOptional(p,'n',10,@isnumeric);
            
            parse(p,obj,varargin{:})
            
            obj = p.Results.obj;
            n = p.Results.n;
            
            datout = obj;
            for i = 1:length(obj)
                if ~strcmp(obj(i).samplerate, 'event')
                    
                    datout(i).data = downsample(obj(i).data,n);
                    datout(i).data = datout(i).data(:);
                    datout(i).samplerate = datout(i).samplerate/n;
                end
            end
        end
        
         function datout = upsample(obj, n)
            % upsample Upsample input signal                           
            
            datout = obj;
            for i = 1:length(obj)
                if ~strcmp(obj(i).samplerate, 'event')
                    
                    datout(i).data = resample(obj(i).data,n,1);
                    datout(i).data = datout(i).data(:);
                    datout(i).samplerate = datout(i).samplerate*n;
                end
            end
        end
        
        function datout = resettime(obj,time0)
            if ~exist('time0','var')
                time0=0;
            end
            datout = obj;
            for i = 1:length(obj)
                datout(i).tend = obj(i).tend - obj(i).tstart + time0;      
                datout(i).tstart = time0;

            end
        end
        
        function t = dattime(obj, seg)
            % t =  DATTIME(obj)
            % t = DATTIME(obj,seg)
            %
            % obj must be single channel
            
            npoints = length(obj.data);
            t = (0:npoints-1)/obj.samplerate + obj.tstart;            
            if exist('seg','var')
                mask = t>=seg(1)-eps & t<seg(2)-eps;
                t(~mask) = [];
            end
        end
        
        
        function data = dattimeind(obj,t)            
            % data = dattimeind(obj,t)
            % Returns data point at specified time. obj must be 1 channel
            datt = dattime(obj);
            for i = 1:length(t)
            [~,temp] = min(abs(t(i)-datt));
            ind(i) = temp(1);
            end
            data = obj.data(ind);
        end
        
        function [freq, pwr] = datfft(obj)
            
            n = length(obj.data);
            nhalf = floor(n/2);
            
            x = fft(obj.data);
            pwr = abs(x(1:nhalf));
            freq = (0:nhalf-1)'*obj.samplerate/n;
        end
        
        function varargout = plot(obj, varargin)
            %  datPlot   Simple plotting function.
            %   datPlot(obj) plots all channels in the dat structure
            %
            %   datPlot(obj, [tstart tend], 'Range',[ymin ymax],'Overlaid',1)
            
            p = inputParser;
            
            addRequired(p,'obj',@(x)isa(x,'dat'));
            addOptional(p,'seg',[]);
            addParamValue(p,'Range',[]);
            addParamValue(p,'Overlaid',0);
            addParamValue(p,'Color',[]);
            
            parse(p,obj,varargin{:});
            
            obj = p.Results.obj;
            seg = p.Results.seg;
            plotrange = p.Results.Range;
            overlaid = p.Results.Overlaid;
            color = p.Results.Color;
            
            emptychans = cellfun(@isempty,{obj.data});
            
            if all(emptychans)
                return
            end
            
            obj = obj(~emptychans);
            
            if isempty(seg)
                datdisp = obj;
            else
                datdisp = datseg(obj, seg);
            end
            
            
            
            if numel(plotrange)==1
                plotrange = [-plotrange, plotrange];
            end
            
            nchans = length(datdisp);
            if isempty(color)
                colors = hsv(nchans);
            elseif length(color)==1
                colors = repmat(colorspec(color),nchans,1);
            elseif ~iscell(color)
                colors = color;
                if size(color,1)==1
                    colors = repmat(color,nchans,1);
                end
                
            else
                for i = 1:length(color)
                    colors(i,:) = colorspec(color{i});
                end
            end
            lightmask = sum(colors,2)>1.5;
            colors(lightmask,:) = max(0,colors(lightmask,:)-.5);   % Darken
            
            %% Main plotting loop
            for i = 1:nchans
                
                if isempty(datdisp(i))
                    continue
                end
                
                color = colors(i,:);
                
                % Invert the channel number for plotting
                inv = nchans-i+1;
                
                % Set up axis labels and limits
                if ~overlaid
                    haxis(inv) = subplot(nchans,1,i);
                    
                    s1 = get(gca,'position');
                    ylabel(datdisp(inv).chanlabel,'FontWeight','normal','FontSize',10);
                    %ylabel([datdisp(inv).chanlabel ' (' datdisp(inv).units ')'])
                    %                     set(get(gca,'YLabel'),'FontWeight','normal','FontSize',10)
                    set(gca,'FontWeight','normal','FontSize',10)
                    xlim([datdisp(inv).tstart datdisp(inv).tend])
                    if  strcmp(datdisp(inv).samplerate, 'event')
                        ylim([-1 1]);
                    elseif ~isempty(plotrange)
                        ylim(plotrange);
                    else
                        if ~strcmp(datdisp(inv).samplerate, 'event')
                            if islogical(datdisp(inv).data)
                                ylim([-.1 1.1])
                            elseif all(ismember(datdisp(inv).data,[-1 0 1]))
                                ylim([-1.1 1.1])
                            else
                                ylims = [prctile(datdisp(inv).data,.5) prctile(datdisp(inv).data,99.5)];
                                if ylims(2) > ylims(1)
                                ylim(ylims)
                                end
                            end
                        end
                    end
                    box off
                    set(gca,'position',s1);
                    if i<nchans
                        set(gca,'XColor','w')
                    else
                        xlabel('Time (s)')
                    end
                    
                else
                    xlabel('Time (s)')
                    haxis(1) = gca;
                    hold on;
                end
                
                
                % For normal time serieschannels
                if ~strcmp(datdisp(inv).samplerate, 'event')
                    t = dattime(datdisp(inv));
                    dispdata = datdisp(inv).data;
                    
                    h(inv) = line(t,dispdata,'LineWidth',1,'Color',color);
                    
                else % For event channels
                    y = zeros(1,length(datdisp(inv).data));
                    if ~isempty(y)
                        h(inv) = line(datdisp(inv).data,y,'Marker','+','LineStyle','none','Color',color);
                    end
                end
                
            end
            
            %% Allow interactive zooming on all subplots
%            dragzoom
            linkaxes(haxis,'x');
            
            %% Define output handle
            if nargout==1
                varargout = {h};
            end
            
        end % function datplot
        function datplot(obj,varargin)
            plot(obj,varargin);
        end
        
       
        % Standard arithmetic functions
        function num = mean(obj)
            for i = 1:length(obj)
                num(i) = mean(obj(i).data);
            end
        end
        
        function datout = datmean(obj)
            data = NaN(length(obj(1).data),length(obj));
            for i = 1:length(obj)
                if ~isempty(obj(i).data)
                data(:,i) = obj(i).data;
                end
            end
            datout = obj(1);
            datout.data = nanmean(data,2);
        end
        
        function datout = plus(obj,a)
            datout = obj;
            for i = 1:numel(obj)
                if isa(a,'dat')
                    if length(a) > 1
                        num = a(i).data;
                    else
                        num = a.data;
                    end
                else
                    num = a;
                end
                datout(i).data = obj(i).data + num;
            end
        end
        function datout = minus(obj,a)
            datout = plus(obj,-a);
        end
        function datout = uminus(obj)
            datout = obj;
            for i = 1:numel(obj)
                datout(i).data = -obj(i).data;
            end
        end
        function datout = times(obj,a)
            if ~isa(obj,'dat')
                temp = a;
                a = obj;
                obj = temp;
            end
            datout = obj;
            for i = 1:numel(obj)
                if isa(a,'dat')
                    datout(i).data = obj(i).data.*a(i).data;
                else
                    datout(i).data = obj(i).data .* a;
                end
            end
        end
        function datout = mtimes(obj, a)
            datout = times(obj, a);
        end
        function datout = mrdivide(obj,a)
            datout = obj;
            for i = 1:numel(obj)
                if isa(a,'dat')
                    datout(i).data = obj(i).data ./ a(i).data;
                else
                    datout(i).data = obj(i).data ./ a;
                end
            end
        end
        function datout = mldivide(obj,a)
            datout = mrdivide(obj,a);
        end
        function datout = ldivide(obj,a)
            datout = mrdivide(obj,a);
        end
        function datout = rdivide(obj,a)
            datout = mrdivide(obj,a);
        end
        
        function copy = low_pass_filter(obj, maxFreq)
            % low pass filter eye movements (up to 'maxFreq' frequency)
            N = length(obj(1).data);
            Fs = obj(1).samplerate;
            
            copy = obj;
            for ch=5:8
                xfft = fft(obj(ch).data);
                % x axis of xfft is in frequency. The difference between
                % consecutive points is Fs/N. I want a cutoff at maxFreq,
                % therefore Fs/N * cutoff = maxFreq or 
                % cutoff = maxFreq * N/Fs;
                cutoff = maxFreq * N / Fs;
                xfft(cutoff+1:end-cutoff) = 0;
                copy(ch).data = abs(ifft(xfft));
            end
        end
        
    end % METHODS
%{
    methods(Access = protected)
        % Override copyElement method:
        function cpObj = copyElement(obj)
            % Make a shallow copy of all four properties
            cpObj = copyElement@matlab.mixin.Copyable(obj);
        end
       
    end
%}
end % CLASSDEF DAT


function saccade_diff = get_saccade_diff(data, saccade_time, samplerate)
    % Compute the abs difference between frames separated by win_size which 
    % relate to "saccade_time" through win_size = saccade_time * samplerate
    % 
    % A high number in saccade_diff(i) means that a saccade is either
    % happening or will happen shortly
    win_size = saccade_time*samplerate;
    
    saccade_diff = zeros(length(data)-win_size+1, 3);
    
    s0 = length(data);
    s1 = round(win_size/2);
    s2 = s0-win_size+s1;
    s3 = round(win_size/4);
    s4 = s0-win_size+s3;
    saccade_diff(:,1) = abs(data(win_size:end) - data(1:end-win_size+1));
    saccade_diff(:,2) = abs(data(s1:s2) - data(1:end-win_size+1));
    saccade_diff(:,3) = abs(data(s3:s4) - data(1:end-win_size+1));
    saccade_diff = max(saccade_diff');
end

