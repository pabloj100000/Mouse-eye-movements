function [mice_rec_segments, no_mice_data, mice_rec] = summary_20150319()
    chaninds = 1:12;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%% Mice Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Load Hanna's recording with mice.
    % Change channels 5-8 from V to deg.
    % Limit data to regions where Drum and Chair are still and light is On.
    % Exclude regions with high eye velocity.
    % Computer variance of eye position across selected segments.
    
    % {
    path1 = '/Users/jadz/Documents/Notebook/Matlab/Eye tracking Mice/data/20150108_H2_botheyes/20150108_H2_botheyes.smr';
    
    % Load mice data
    mice_rec = importSpike(path1, chaninds);

    % calibrate data, changing channels 5-8 from V to deg
    mice_rec = mice_rec.calibrate_data(5:8, fileparts(path1));

    % Chair and drum are still in between 500 and 680s
    % During that time, light is off between 491.4 and 557s.
    % I'm only using data with still chair and drum and light on.
    mice_rec = mice_rec.datsplit(557,123);

    % exclude regions with high speed eye movement
    [startT, duration] = detect_artifact_timing(mice_rec, .06, .05, 1, '0011110000');
    mice_rec_segments = dat_segments(mice_rec, startT, duration);

    % Compare recording with selected segments for further processing after
    % exclussion of high velocity segments
    figure(1), mice_rec_segments.comparePlots(mice_rec)
    
    % Plot power spectra of selected segments
    mice_rec_segments.estimatePSD();
    figure(2), mice_rec_segments.plotPSD('k', 0);

    % get the variance
    mice_rec_segments.get_variance();
    %}
    %%%%%%%%%%%%%%%%%%%%%%%%%% No Mice Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Load Hanna's recording with no mice.
    % Change channels 5-6 from V to deg (channels 7-8 are noise, no sensor
    % attached).
    % Limit data to regions where Drum and Chair are still
    % Computer variance of sensor position.
    % {
    chaninds=1:8;
    path2 = '/Users/jadz/Documents/Notebook/Matlab/Eye tracking Mice/data/20140820_HMC1512_075_200mm_3mmAway.smr';
    
    % Load data with no mice
    no_mice_data_ori = importSpike(path2 ,chaninds);
    
    % calibrate data
    no_mice_data_ori = no_mice_data_ori.calibrate_data(5:8, fileparts(path2));
    
    % For variance, limit experiment to where chair and drum are still
    no_mice_data = dat_segments(no_mice_data_ori, [0, 14.5], [4, -1]);
    
    % Compare recording with selected segments for further processing after
    % exclussion of high velocity segments
    figure(3), no_mice_data.comparePlots(no_mice_data_ori)

    % Plot power spectra of selected segments
    no_mice_data.estimatePSD();
    figure(2), no_mice_data.plotPSD('r', 1);
    for i=1:4
        subplot(4,1,i)
        xlim([0,150])
        ylim([0, .03])
    end
    
    % compute the variance
    no_mice_data.get_variance();
    
    
    %}
    %%%%%%%%%%%%%%%%%%%%%%%%%% Compare Variances %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (exist('no_mice_data', 'var') && exist('no_mice_data', 'var'))
        sprintf('\tVarNoMice\tVarMice\tsdNoMice\tsdMice')
        for ch=5:8
            [no_mice_data.variance(ch), mice_rec_segments.variance(ch), no_mice_data.std(ch), mice_rec_segments.std(ch)]
%            sprintf('\t%1.1e\t%1.1e\t%1.1e\t%1.1e', variance_no_mice(ch), ...
%                variance_mice(ch), sqrt(variance_no_mice(ch)), sqrt(variance_mice(ch)));
        end
    end
end
