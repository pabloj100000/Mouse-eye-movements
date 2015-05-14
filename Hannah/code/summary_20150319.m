function rec = summary_20150319()    
    %%%%%%%%%%%%%%%%%%%%%%%%%% Mice Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Load Hanna's recording with mice.
    % Change channels 5-8 from V to deg.
    % Limit data to regions where Drum and Chair are still and light is On.
    % Exclude regions with high eye velocity.
    % Computer variance of eye position across selected segments.
    
    % {
    path1 = '/Users/jadz/Documents/Notebook/Matlab/Eye tracking Mice/data/20150319_H2_stationary/20150319_H2_stationary.smr';
    

    % for these particular experiment, I only care about channel 5. 
    %   channels 1-4, have spourious data (chair and drum always stopped
    %   channel 5,  left eye that is somewhat linear
    %   channel 6,  left eye but less linear than channel 5
    %   channel 7-8, not recorded
    chaninds = 5;

    % Load mice data
    mice_rec = importSpike(path1, chaninds);

    % calibrate data, changing channels 5 and 6 from V to deg. Channels 7
    % and 8 have noise, where not recorded
    mice_rec = mice_rec.calibrate_data(fileparts(path1));

    % head fixed mouse with same gender mouse 1ft away
    rec(1) = mice_rec.datsplit(1, 660);
    
    % head fixed mouse with fixed stripes in front
    rec(2) = mice_rec.datsplit(830, 600);
    
    % head fixed mouse in the dark
    rec(3) = mice_rec.datsplit(1440, 600);
    
    % exclude regions with high speed eye movement from all recordings
    [startT, duration] = rec(1).detect_high_speed(.06, .05, 1);
    rec_segments(1) = dat_segments(rec(1), startT, duration);

    [startT, duration] = rec(2).detect_high_speed(.06, .05, 1);
    rec_segments(2) = dat_segments(rec(2), startT, duration);

    [startT, duration] = rec(3).detect_high_speed(.06, .05, 1);
    rec_segments(3) = dat_segments(rec(3), startT, duration);

    % Compare recording with selected segments for further processing after
    % exclussion of high velocity segments
    figure(1)
    subplot(3,1,1);
    plot(rec(1).time_axis, rec(1).data);
    hold on
    title('With same sex mouse at 1ft');
    rec_segments(1).comparePlots()
    hold off

    subplot(3,1,2);
    plot(rec(2).time_axis, rec(2).data);
    hold on
    title('With fixed background');
    rec_segments(2).comparePlots()
    hold off

    subplot(3,1,3);
    plot(rec(3).time_axis, rec(3).data);
    hold on
    title('In the dark');
    rec_segments(3).comparePlots()
    hold off
    
    xlabel('Time (s)')
    ylabel('deg')

end
