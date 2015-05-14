classdef mouse<handle
    % a class to read and process data from Marcel De Jeu and 
    
    properties
        rate;
        tax;
        vert;
        hori;
    end
    methods
        % Class constructor
        function obj = mouse(file_path)
            % each file holds a struct that is named according to the mouse
            % I don't care about the particular mouse here and want to
            % access the data generically as just 'mouse'
            mouse = load(file_path);
            fname = fieldnames(mouse);  % this should have only one field
                                        % of the form "mouse#"
            
            obj.rate = 120;   % in Hz
            obj.vert = mouse.(fname{1}).calv;
            obj.hori = mouse.(fname{1}).calh;
            
            % mouse2 for example doesn't have the same number of vert and
            % horizontal points.
            fixLength(obj);

            obj.tax = (0:1/obj.rate:(size(obj.vert)-1)/obj.rate)';
            
            
        end
        
        function plot1(obj)
            % just plot both vert and horizontal eye movements against time
            % in the same plot
            s_plot_name = [inputname(1), '_HV_vs_T']
            
            % return the handle to the figure named s_plot_name
            h = findall(0, 'type', 'figure', 'name', s_plot_name);
            
            if isempty(h)
                h = figure('name', s_plot_name);
            end
            
            figure(h)
            title('Eye position');
            plot(obj.tax, obj.vert, 'r', obj.tax, obj.hori, 'b');
            
            xlabel('Time (s)');
            legend('vert', 'hori');
        end
        
        function plot2(obj, varargin)
            % plot vert vs hori in between the times given
            
            p  = inputParser;   % Create an instance of the inputParser class.
            p.addParameter('startT', 0, @(x)x>=0);
            p.addParameter('endT', obj.tax(end), @(x) x>=0 && x<=obj.tax(end));
            parse(p, varargin{:})
            startT = p.Results.startT;
            endT = p.Results.endT;

            % convert startT and endT to points in tax
            startP = find(obj.tax>=startT, 1);
            endP = find(obj.tax>=endT, 1);
            s_plot_name = [inputname(1), '_V_vs_H'];
            
            % get handle of corresponding figure
            h = findall(0, 'type', 'figure', 'name', s_plot_name);
            
            if isempty(h)
                h = figure('name', s_plot_name);
            end
            
            figure(h)
            plot(obj.hori(startP:endP), obj.vert(startP:endP), 'k');      
            xlabel('Horizontal Position (deg)');
            ylabel('Vertical Position (deg)');
        end
    end
end

function fixLength(obj)
    hori_pnts = length(obj.hori);
    vert_pnts = length(obj.vert);
    
    if (hori_pnts < vert_pnts)
        obj.vert(hori_pnts+1:end)=[];
    elseif (vert_pnts < hori_pnts)
        obj.hori(vert_pnts+1:end)=[];
    end
end