classdef dat_segments<handle
    % an extension of Hanna's dat class
    %
    % a 'dat_segments' object is generated from Hanna's 'dat' object. 
    % startT and duration are arrays of the same length flagging into
    % regions of the original 'dat' recording
    % that we want to crop and save. Each sub 'dat' segment is stored as an
    % array in obj.segments. The rest of the fields just make it convenient
    % to work with dat_segment. 
    % 
    % samplerate:   since all channels in all segments have the same 
    %                'samplerate' I just generate a field obj.samplerate
    % segmentsN:    just the number of segments available (same as length(startT))
    % variance:     computed for each segment independently and then
    %               results are merged taking each segments length into
    %               consideration.
    % std:          idem variance
    % PSD:          (Power Spectrum Density) and freq can be used to run 
    %               'plot(freq, PSD) and compare different types of 
    %               recordings, for example with and without mice or with
    %               different stimulus protocols.
    
    properties
        segmentsN;
        segments = [dat dat];
        variance;
        std;
        PSD;
        freq;
        samplerate;
    end
    methods
        % Class constructor
        function obj = dat_segments(dat, startT, duration)            
            if length(startT) ~= length(duration)
                error, 'startT and endT should have the same number of points';
            end
            
            obj.segmentsN = length(startT);
            
            % Preallocate object array
            obj.segments(obj.segmentsN) = dat;
            
            obj.samplerate = dat.samplerate;
            
            for i = 1:length(startT)
                obj.segments(i) = dat.datsplit(startT(i), duration(i));
            end
            
        end
        
        function comparePlots(obj)
            % Make a plot of all selected segments in obj
                  
            for i = 1:obj.segmentsN
                segment = obj.segments(i);
                tdelta = 1/segment.samplerate;
                tstart = segment.tstart;
                tend = (length(segment.data) - 1)*tdelta + tstart;
                x = tstart-tdelta:tdelta:tend-tdelta;
                plot(x, segment.data, 'k');
            end
        end
        
        function output = cat_segments(obj, DCFlag)
            % concatenate all segments in obj
            % if DFflag is set then DC is kept
           
            pnts=0;
            for i=1:obj.segmentsN
                segment = obj.segments{i};
                pnts = pnts + length(segment(1).data);
            end
            
            % preallocate memory for channels 5-8
            output = zeros(pnts, 4);
            
            % fill cat with data from segments
            start_pnt = 1;
            for i=1:obj.segmentsN
                segment = obj.segments{i};
                pnts_being_added = length(segment(1).data); %channels 1-8 have the same number of points
                for ch=5:8
                    if DCFlag
                        output(start_pnt:start_pnt + pnts_being_added - 1, ch-4) = segment(ch).data;
                    else
                        output(start_pnt:start_pnt + pnts_being_added - 1, ch-4) = segment(ch).data - mean(segment(ch).data);
                    end
                end
                start_pnt = start_pnt + pnts_being_added;
            end
        end
        
        function get_variance(obj)
            % compute variance of signal across all segments and channels
            % obj: a list of dat objects.
            % Comput the variance for each channel and segment independently. Then
            % combine variances for different segments weighting by the segment
            % length.
            
            channelN = min(8, length(obj.segments{1}));
            segment_var = ones(channelN, obj.segmentsN);
            segment_length = zeros(channelN, obj.segmentsN);
            
            % Compute Variance for each independent segment
            for i = 1:obj.segmentsN
                segment = obj.segments{i};
                for ch = 1:channelN
                    segment_var(ch,i) = var(segment(ch).data);
                    segment_length(ch, i) = segment(ch).tend - segment(ch).tstart;
                end
            end
            
            % Merge all the variances for a given channel, taking the length of
            % each FEM into account
            obj.variance = sum(segment_var .* segment_length, 2)./sum(segment_length,2);
            obj.std = sqrt(obj.variance);
        end
        
        function scale_factor = get_calibration(obj)
            % concatenate all segments together. First find out how many
            % points arrays will ahve
            pnts = 0;
            for i=1:length(obj)
                pnts = pnts + length(obj(i).data);
            end
            
            channels = zeros(pnts, 5);   % chair is channel 1, Left is channels 2-3 and Right is channels 4-5
            pnts=1;
            
            for i=1:length(obj)
                segment = obj(i);
                new_pnts = length(segment(1).data);
                for ch=4:8
                    channels(pnts:new_pnts,ch-3) = segment(ch).data;
                end
            end
            
            scale_factor = ones(8,1);
            rsquare = ones(8,1);
            for ch = 4:8
                [fitobject, fog] = fit(channels(:,ch), channels(:,4), 'poly1');
                
                scale_factor(ch) = fitobject.p1;
                rsquare(ch) = fog.rsquare;
            end
            
            if rsquare(5) > rsquare(6)
                scale_factor(6)=0;
            else
                scale_factor(5)=0;
            end
            
            if rsquare(7) > rsquare(8)
                scale_factor(8)=0;
            else
                scale_factor(7)=0;
            end
        end
        
        function estimatePSD(obj)
            % compute power  spectral density estimation
            % I'm following the example in http://www.mathworks.com/help/signal/ug/psd-estimate-using-fft.html
            
            onesegment = obj.segments(1);
            Fs = onesegment(1).samplerate;
            x = obj.cat_segments(0);
            N = size(x,1);
            xdft = fft(x,[],1);
            xdft = xdft(1:N/2+1, :);
            psdx = (1/(Fs*N)) * abs(xdft).^2;
            psdx(2:end-1, :) = 2*psdx(2:end-1, :);
            obj.freq = 0:Fs/N:Fs/2;
            obj.PSD = psdx;
        end
        
        function plotPSD(obj, color, holdonflag)
            for ch=1:4
                subplot(4,1,ch)
                if holdonflag
                    hold on
                end
                plot(obj.freq, obj.PSD(:,ch), color)
                hold off
            end
            xlabel('Freq (Hz)');
            subplot(4,1,ch)
            title('Power Spectrum Density');
        end
        
        function [x, trajectory] = getEyeTrajectory(obj, cutOffFreq, segment, ch)
            % Filter data for the given segment and channel upto the cutOffFreq.
            % Returns both x and y such that you can plot(x, trajectory)
            % and overlay this with obj.plotSegment(segment, ch)
            
            one_segment = obj.segments{segment};
            data = one_segment(ch).data;
            Fs = obj.samplerate;
            N = size(data,1);
            
            xdft = fft(data,[],1);
            
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
        
        function plotSegment(obj, segment, ch, varargin)
            % just plot the desired segment and ch
            segment = obj.segments{segment};
            x = 0:length(segment(ch).data)-1;
            x = x/obj.samplerate;
            if isempty(varargin)
                plot(x, segment(ch).data);
            else
                plot(x, segment(ch).data, varargin{:});
            end
        end
        
        function testGetEyeTrajectory(obj, cutOffFreq, segment, ch)
           figure()
           obj.plotSegment(segment, ch, 'b');
           [x, trajectory] = obj.getEyeTrajectory(cutOffFreq, segment, ch);
           hold on
           plot(x, trajectory, 'r');
           hold off
        end
    end
end

