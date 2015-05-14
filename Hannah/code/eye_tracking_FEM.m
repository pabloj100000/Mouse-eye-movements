classdef eye_tracking_FEM
    properties
        startT;
        saccade_time;
        data;
        N;
        samplerate;
        units;
        threshold;
        min_fixation_length;
        var;
        std;
        channelN;
    end
    methods
        % Class constructor
        function FEM = eye_tracking_FEM(dat, threshold, min_fixation_length, ch_selection, varargin)
            % dat:  Hanna's object
            % threshold: changes grater than threshold * range(channel) are
            % trated as saccades
            % min_fixation_length: if fixation is less than
            % 'min_fixation_length' it is discarded
            % ch_selection: bit flag with as many values as channels. If a
            % channel is set it contributes to defining fixations. First 8
            % channels in Hanna's object are considered.
            
            p = inputParser;
            addOptional(p,'samplerate',dat(1).samplerate,@(x) isnumeric(x) || strcmp(x,'event'));
            addOptional(p,'units',dat(5).units,@(x) ischar(x) || iscell(x) ||isempty(x));
            addOptional(p, 'saccade_time', 0.03, @(x) isnumeric(x) && x>0);
            
            parse(p,varargin{:});

            FEM.samplerate = p.Results.samplerate;
            FEM.saccade_time = p.Results.saccade_time;
            FEM.units = p.Results.units;
            FEM.threshold = threshold;
            FEM.channelN = min(length(dat),8);      % channels 9 and 10 are different. DOn't look into those
            FEM.min_fixation_length = min_fixation_length;
            [FEM.data, FEM.startT, FEM.N] = detect_FEM_regions(dat, FEM, ...
                ch_selection);
            FEM.var = get_FEM_variance(FEM);
            FEM.std = sqrt(FEM.var);
        end
        
        
        function obj = plot(obj, index, varargin)
            % plot FEM for the given index.
            % if fig_ori is passed, the figure axis 5:8 (corresponding to
            % 1:4 in the plot I'm generating) get scaled accordingly. Data
            % should be the same.
            if index > obj.N
                error('Aborting: obj.N is less than requested index');
            end
                
            x = obj.startT{index}:1/obj.samplerate:obj.startT{index} + (length(obj.data{1}{index})-1)/obj.samplerate;
            for i = 1:obj.channelN
                subplot(obj.channelN,1,i)
                plot(x, obj.data{i}{index})
            end
            xlabel('Time (s)')
            
            %{
            if (varargin{1})
                targetX = xlim();
                yrange = {4,1};
                for i = 1:4
                    subplot(4,1,i)
                    yrange{i} = ylim();
                end
                
                figure(varargin{1});
                xlim(targetX)
                for i = 1:4
                    subplot(8,1,i)
                    yrange{i}
                    ylim(yrange{i})
                end
            end
            %}
        end
        
        function comparePlots(obj, FEM_ch, dat)
            % From Hanna's object "dat", plot the channel "dat_ch" in one
            % color and then 
            % add all the FEM regions for the corresponding channel from
            % FEM_ch
            % I'm assuming that dat_ch is the same as FEM_ch
            dat_ch = FEM_ch;
            dat(dat_ch).plot()
            hold on
            for i = 1:obj.N
                x = obj.startT{i}:1/obj.samplerate:obj.startT{i} + (length(obj.data{1}{i})-1)/obj.samplerate;
                % subtract 1/samplerate (1ms) to allign to Hanna's plot
                x = x - 1/obj.samplerate;
                plot(x, obj.data{FEM_ch}{i}, 'k');
            end
            hold off
        end
        
    end
end

function [FEM_data, startT, N] = detect_FEM_regions(dat, FEM, ...
    ch_selection)
    % for each eye and channel generate the still_starts/ends arrays. I
    % define still as anything below threshold
    % dat:          Hanna's object
    % threshold:    in 0-1 units, actual threshold is computed by multiplying
    %               by the channel's range. Can be a single number or an array with
    %               individual thresholds.
    % min_fixation_length:  minimum amount of quite frames that I do accept
    %                       as a fixation (in seconds)
    
    changing_frames = zeros(1000,1);
    idx = 1;
    for ch = 1:FEM.channelN
        if ~bitget(ch_selection, ch)
            continue
        end
        
        % compute the channel's threshold
        if (length(FEM.threshold) > 1)
            thresh = range(dat(ch).data(1:end)) * FEM.threshold(ch);
        else
            thresh = range(dat(ch).data(1:end)) * FEM.threshold;
        end

        % compute the diff between data separated by saccade_time.
        % saccade_diff is saccade_time*samplerate-1 shorter than data
        % A high number in saccade_diff means a saccade will happen shortly
        % or is under way
        saccade_diff = get_saccade_diff(dat(ch).data, FEM.saccade_time, dat(ch).samplerate);

        % get an array of frames where this channel had a strong velocity
        frames = 1:length(saccade_diff);
        frames = frames(saccade_diff > thresh);
        changing_frames(idx:idx+length(frames)-1) = frames;
    end

    % sort frames where large velocity was detected
    changing_frames = sort(changing_frames);

    % add frames corresponding to beggining and end of recording
    exp_length = dat(1).tend-dat(1).tstart;
    changing_frames = [1; changing_frames; exp_length*dat(1).samplerate];
    
    % discard regions of recording that are not quiet for at least
    % min_fixation_length
    interruption_distance = diff(changing_frames);

    FEM_data = cell(FEM.channelN, 1);
    startT = {};
    for i = length(interruption_distance):-1:1

        % if interruption_distance is larger than FEM.min_fixation_length add
        % data from dat to FEM_data
        if interruption_distance(i) >= FEM.min_fixation_length*dat(1).samplerate
%            offset = dat.tstart;
            FEM_starts = changing_frames(i)+1;
            FEM_ends = changing_frames(i+1);
            startT = [FEM_starts/dat(1).samplerate startT];
            for idx=1:FEM.channelN
                if isempty(FEM_data{idx})
                    FEM_data{idx} = {dat(idx).data(FEM_starts:FEM_ends)};
                else
                    FEM_data{idx} = [dat(idx).data(FEM_starts:FEM_ends) FEM_data{idx}];                
                end
            end
        end

    end

    N = length(FEM_data{1});
end

function FEM_var = get_FEM_variance(obj)
    % compute variance of signal (for each channel and FEM period)
    % in obj

    FEM_var = ones(obj.channelN,obj.N);
    FEM_length = zeros(obj.channelN, obj.N);

    % Compute Variance for each independent FEM segment
    for i = 1:obj.N
        for ch = 1:obj.channelN
            FEM_var(ch,i) = var(obj.data{ch}{i});
            FEM_length(ch, i) = length(obj.data{ch}{i});
        end
    end

    % Merge all the variances for a given channel, taking the length of
    % each FEM into account
    FEM_var = sum(FEM_var .* FEM_length, 2)./sum(FEM_length,2);
end

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

function dat_segments = exclude_artifacts(dat, saccade_time, ...
        threshold, min_fixation_length, ch_selection)
    % for each channel flagged by ch_selection from dat, identify regions
    % with large velocities and exclude those. 
    % This function returns a cell array of dat fragments once all high
    % velocity segments were excluded
    %
    % dat:          Hanna's object
    % saccade_time: How long do we expect a saccade or an artifact to last?
    %               in seconds
    % threshold:    in 0-1 units, actual threshold is computed by multiplying
    %               by the channel's range. Can be a single number or an array with
    %               individual thresholds (one per channel).
    % min_fixation_length:  minimum amount of quiet frames that I do accept
    %                       in between saccades or artifacts (in seconds)
    
    changing_frames = zeros(1000,1);
    idx = 1;
    dat_segments = cell(1000,1);   % this is where small dats will be appended to

    for ch = 1:FEM.channelN
        if ~bitget(ch_selection, ch)
            continue
        end
        
        % compute the channel's threshold
        if (length(FEM.threshold) > 1)
            thresh = range(dat(ch).data(1:end)) * threshold(ch);
        else
            thresh = range(dat(ch).data(1:end)) * threshold;
        end

        % compute the diff between data separated by saccade_time (and
        % saccade_time/2, saccade_time/4).
        % saccade_diff is saccade_time*samplerate-1 shorter than data
        % A high number in saccade_diff means a saccade will happen shortly
        % or is under way
        saccade_diff = get_saccade_diff(dat(ch).data, saccade_time, dat(ch).samplerate);

        % get an array of frames where this channel had a strong velocity
        frames = 1:length(saccade_diff);
        frames = frames(saccade_diff > thresh);
        changing_frames(idx:idx+length(frames)-1) = frames;
    end

    % sort frames where large velocity was detected
    changing_frames = sort(changing_frames);

    % change from frames to times, taking tstart into account
    changing_times = [dat(1).tstart; dat(1).tstart + changing_frames/dat(1).samplerate; dat(1).tend];
    
    % discard regions of recording that are not quiet for at least
    % min_fixation_length
    fem_length = diff(changing_times);

    for i = 1:length(fem_length)

        % if interruption_distance is larger than FEM.min_fixation_length add
        % data from dat to dat_segments
        if fem_length(i) >= min_fixation_length
            dat_segments{idx} = dat.datsplit(changing_times(i), changing_times(i+1));
        end
    end
end