classdef eye_tracking_FEM
    properties
        data;
        N;
        samplerate;
        units;
    end
    methods
        % Class constructor
        function obj = FEM(dat, threshold, min_fixation_length, varargin)
            
            p = inputParser;
            addOptional(p,'N',1,@(x) isnumeric(x)));
            addOptional(p,'samplerate',1,@(x) isnumeric(x) || strcmp(x,'event'));
            addOptional(p,'units',{},@(x) ischar(x) || iscell(x) ||isempty(x));
            addOptional(p,'threshold',1,@(x) isnumeric(x) && all(0<x) && all(x<1));
            addOptional(p,'min_fixation_length',1,@(x) isnumeric(x));
            
            parse(p,varargin{:});
            
            FEM.N = p.Results.chanlabel;
            FEM.samplerates = p.Results.samplerate;
            FEM.unitss = p.Results.units;
            FEM.data = detect_saccades(dat, threshold, min_fixation_length);
            FEM.threshold = threshold;
            FEM.min_fixation_length = min_fixation_length;
            
        end
        
        function FEM = detect_saccades(dat, threshold, min_fixation_length)
            % for each eye and channel generate the still_starts/ends arrays. I
            % define still as anything below threshold
            % dat: Hanna's object
            % threshold: in 0-1 units, actual threshold is computed by multiplying
            % by the channel's range. Can be a single number or an array with
            % individual thresholds.
            % min_fixation_length:  minimum amount of quite frames that I do accept
            % as a fixation (in seconds)
            
            changing_frames = [];
            
            for ch = 1:8
                % compute the channel's threshold
                if (length(threshold) > 1)
                    thresh = range(dat(ch).data) * threshold(ch);
                else
                    thresh = range(dat(ch).data) * threshold;
                end
                
                % compute the diff, diff is one frame shorter than data
                % ch_diff(1) is the velocity in between frames 2 and 1
                ch_diff = abs(diff(dat(ch).data));
                
                % get an array of frames where channel had a strong change
                frames = 1:length(ch_diff);
                changing_frames = [changing_frames, frames(ch_diff > thresh)];
            end
            
            % sort frames where large velocity was detected
            changing_frames = sort(changing_frames);
            
            % discard regions of recording that are not quiet for at least
            % min_fixation_length
            interruption_distance = diff(changing_frames);
            
            FEM = cell(4, 1);
            for i = length(interruption_distance):-1:1
                
                % if interruption_distance is larger than min_fixation_length add
                % the eye tracking position corresponding to left_ch1, left_ch2,
                % right_ch1, right_ch2 to FEM
                if interruption_distance(i) >= min_fixation_length*dat(1).samplerate
                    FEM_starts = changing_frames(i)+1;
                    FEM_ends = changing_frames(i+1);
                    for idx=1:4
                        if isempty(FEM{idx})
                            FEM{idx}{1} = dat(idx+4).data(FEM_starts:FEM_ends);
                        else
                            FEM{idx} = [FEM{idx} dat(idx+4).data(FEM_starts:FEM_ends)];
                        end
                    end
                end
                
            end
        end
        function plot_FEM(FEM, index)
            for i = 1:4
                subplot(4,1,i)
                plot(FEM{i}{index})
            end
        end

    end
end