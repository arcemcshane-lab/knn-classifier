%% LOAD DATA
date = num2str(28);
cortical_areas = {'M1F'; 'M1U'; 'S1F'; 'S1U'};
load(strcat('201902', date, 'graphs.mat')) %hits for 6 markers for 2000 ms window trials
load(strcat('201902', date, '_CranialKinematics.mat')) %fully digitized Rocky trials

%% FIND MISALIGNED TRIALS

misaligned_trials = [];
misaligned_trials = 0;

for i = 1:length(graphme)
    
    if length(graphme{i}) ~= length(Kinematics.index{i})
        misaligned_trials = [misaligned_trials ; i];
    end
    
end

misaligned_trials(1) = [];

valid_trials = 1:length(graphme);
valid_trials(misaligned_trials) = [];

%% CALCULATE CONTACT ONSETS, OFFSETS
contact_onsets = [];
contact_offsets = [];

for i = valid_trials

    onset_offset_binary_1 = 0; % 0 = no contact during this timestep; 1 = yes contact during this timestep. init as 0
    contact_onsets{i} = [];
    contact_offsets{i} = [];

    for j = 1:size(graphme{i}, 1)
       if j > 1
           onset_offset_binary_1 = any(graphme{i}(j-1, :));
       end
       onset_offset_binary_2 = any(graphme{i}(j, :));
       if onset_offset_binary_1 - onset_offset_binary_2 == -1
           contact_onsets{i} = [contact_onsets{i} ; j];
       elseif onset_offset_binary_1 - onset_offset_binary_2 == 1
           contact_offsets{i} = [contact_offsets{i} ; j];
       end
    end
    
    if size(contact_onsets{i}, 1) ~= 0 && size(contact_offsets{i}, 1) ~= 0 % check for empty arrays

        if contact_onsets{i}(1) > contact_offsets{i}(1)
            contact_onsets{i} = [1 ; contact_onsets{i}];
        end

        if contact_onsets{i}(size(contact_onsets{i}, 1)) > contact_offsets{i}(size(contact_offsets{i}, 1))
            contact_offsets{i} = [contact_offsets{i} ; size(graphme{i}, 1) ];
        end
    
    end

    if size(contact_onsets{i}) ~= size(contact_offsets{i})
        error('Contact onsets and offsets not the same length. \nThis is caused by something I did not forsee and requires further testing.')
    end
    
    for j=1:size(contact_onsets{i}, 1) % Grab neuronal spiketimes
        contact_onsets{i}(j, 2) = Kinematics.index{i}(contact_onsets{i}(j, 1), 3);
        contact_offsets{i}(j, 2) = Kinematics.index{i}(contact_offsets{i}(j, 1), 3);
    end
end

clear onset_offset_binary_* i j

%% BEGIN FOR LOOP
% for area = 1:length(cortical_areas)
for area = 1:1
    %% LOAD DATA    
    NEV = load(strcat('201902', date, '_', cortical_areas{area}, '_sortedspikes.mat'));
    NEV_cell = struct2cell(NEV); %converts to cell for easier indexing
    
    
    %% MAKE SPIKETIMES
    
    spikenames = fields(NEV);
    for i = 1:length(NEV_cell)
        spiketimes.(string(spikenames(i))) = NEV_cell{i}.times;
    end
    spiketimes_cell = struct2cell(spiketimes);
   
    %% CALCULATE GRAPH STARTS AND ENDS
    spikes = [];
    
    for i = valid_trials

        starttime = Kinematics.index{i}(1, 3) / 30000;
        endtime = Kinematics.index{i}(size(Kinematics.index{i}, 1), 3) / 30000;

        for neuron=1:length(spiketimes_cell)
            spikes{neuron}{i} = spiketimes_cell{neuron}(find(spiketimes_cell{neuron} >= starttime & spiketimes_cell{neuron} <= endtime));
            if length(spikes{neuron}{i} > 0)
                spikes{neuron}{i} = spikes{neuron}{i} - spikes{neuron}{i}(1);
                spikes{neuron}{i} = spikes{neuron}{i} * 200; % convert to 200 Hz XRay data
            end
        end
    
    end
    
    %% GRAPH
    i = 1; %change to loop over all of valid_trials
    
    center_line = -3;
    
    hold on    
    yline(center_line);
    % pos = [x y w h]
    for j = 1:size(contact_onsets{i}, 1)
        pos = [contact_onsets{i}(j, 1) center_line-2 contact_offsets{i}(j, 1)-contact_onsets{i}(j, 1) 4];
        rectangle('Position',pos, 'FaceColor', 'w');
    end
    axis([0 2000 -6 length(spiketimes_cell) + 1]);
    
    for neuron = 1:length(spiketimes_cell)
        scatter(spikes{neuron}{i}, zeros(size(spikes{neuron}{i})) + neuron, '|');
    end
    
    legend({'Contact Events', 'Neuron Spiking'}, 'Location', 'northeastoutside')
    
    hold off
    
    %% POPULATE FR TABLES
    fr_before_after = [];
    fr_during = [];
    
    for i = valid_trials
    
        for neuron = 1:length(spiketimes_cell)

            fr_before_after{neuron}{i} = 0;
            fr_during{neuron}{i} = 0; 

            for j = 1:size(contact_onsets{i}, 1)

                % make starttime, endtime for ranges around each contact event. May
                % need starttime_1, endtime_1 and *_2 for each (before and after).
                % Calculate the FR, compare with fr_during

                if j == 1 % Special case for first contact event
                    if contact_onsets{i}(j, 1) ~= 1 % If first contact event does not start at beginning of trial
                        starttime_1 = Kinematics.index{i}(1, 3) / 30000;
                        endtime_1 = contact_onsets{i}(j, 2) / 30000;
                    else % If first contact event does start at beginning of trial
                        starttime_1 = 0;
                        endtime_1 = 0;
                    end
                    starttime_2 = contact_offsets{i}(j, 2) / 30000;
                    endtime_2 = contact_onsets{i}(j+1, 2) / 30000;
                elseif j == size(contact_onsets{i}, 1) % Special case for last contact event
                    if contact_offsets{i}(j, 1) ~= length(graphme{i}) % If last contact event does not go to the end of trial
                        starttime_2 = contact_offsets{i}(j, 2) / 30000;
                        endtime_2 = Kinematics.index{i}(length(graphme{i}), 3) / 30000;
                    else % If last contact event does go to the end of trail
                        starttime_2 = 0;
                        endtime_2 = 0;
                    end
                    starttime_1 = contact_offsets{i}(j-1, 2) / 30000;
                    endtime_1 = contact_onsets{i}(j, 2) / 30000;
                else % All other contact events
                    starttime_1 = contact_offsets{i}(j-1, 2) / 30000;
                    endtime_1 = contact_onsets{i}(j, 2) / 30000;

                    starttime_2 = contact_offsets{i}(j, 2) / 30000;
                    endtime_2 = contact_onsets{i}(j+1, 2) / 30000;
                end


                spikes_count_1 = length(find(spiketimes_cell{neuron} >= starttime_1 & spiketimes_cell{neuron} <= endtime_1));
                spikes_count_2 = length(find(spiketimes_cell{neuron} >= starttime_2 & spiketimes_cell{neuron} <= endtime_2));

                fr_before_after{neuron}{i} = [ fr_before_after{neuron}{i} ; (spikes_count_1 + spikes_count_2) / ((endtime_1 - starttime_1) + (endtime_2 - starttime_2)) ];

                starttime = contact_onsets{i}(j, 2) / 30000;
                endtime = contact_offsets{i}(j, 2) / 30000;        

                spikes_count = length(find(spiketimes_cell{neuron} >= starttime & spiketimes_cell{neuron} <= endtime));
                fr_during{neuron}{i} = [fr_during{neuron}{i} ; (spikes_count / (endtime - starttime))];
            end

            fr_before_after{neuron}{i}(1) = [];
            fr_during{neuron}{i}(1) = [];

        end
    
    end        
    
    %% CALCULATE GEOM INDEX FOR FR
    
    i = 1; % change to loop over all valid_trials
    
    geom_index = [];
    
    for i = valid_trials
    
        for neuron = 1:length(spiketimes_cell)

            geom_index{neuron}{i} = (fr_before_after{neuron}{i} - fr_during{neuron}{i}) ./ (fr_before_after{neuron}{i} + fr_during{neuron}{i});

        end
    
    end
    
    %% CALCULATE 95% CIs
    
    geom_stats = []; % (mean, SD)
    geom_stats.total = zeros(1, length(spiketimes_cell)); 
    
    for i = valid_trials
        
        if size(contact_onsets{i}, 1) > 0
        
            temp_mat = zeros(size(contact_onsets{i}, 1), 1);

            for neuron = 1:length(spiketimes_cell)

                temp_mat = horzcat(temp_mat, geom_index{neuron}{i});

            end

            temp_mat(:, 1) = [];

            geom_stats.total = [ geom_stats.total; temp_mat ];
            
        end
        
    end
    
    geom_stats.mean(1, :) = mean(geom_stats.total, 1, 'omitnan');
    geom_stats.std(1, :) = std(geom_stats.total, 1, 'omitnan');
    geom_stats.se(1, :) = geom_stats.mean ./ sqrt(geom_stats.std);
    geom_stats.ci(1, :) = geom_stats.mean - geom_stats.se*1.96;
    geom_stats.ci(2, :) = geom_stats.mean + geom_stats.se*1.96;
    geom_stats.ci(3, :) = geom_stats.se*1.96;
    geom_stats.total(1, :) = [];
    
    %% GRAPH GEOM INDICES
    
    i = 1; % change to iterate over valid_trials
    neuron = 1; % change to iterate over neurons
    
    x = 1:length(spiketimes_cell);
    
    fig = figure;
    fig = bar(x, geom_stats.mean);
    
    hold on
    
    er = errorbar(x, geom_stats.mean, -1 * geom_stats.ci(3, :), geom_stats.ci(3, :));
    er.Color = 'black';
    er.LineStyle = 'none';
    
    hold off
%% END FOR LOOP
end