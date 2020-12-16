%% LOAD DATA
date = 27;
date = num2str(date);
cortical_areas = {'M1F'; 'M1U'; 'S1F'; 'S1U'};
load(strcat('201902', date, 'graphs.mat')) %hits for 6 markers for 2000 ms window trials
load(strcat('201902', date, '_CranialKinematics.mat')) %fully digitized Rocky trials

max_window_size = 30; % frames; an int; set to 9999 for unbounded window sizing
offset_before_after_contact_events = 5; % frames; an int; how many frames ignored before/after each contact event

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
    disp(cortical_areas{area});
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
    
    title(sprintf('Spikes and Contact Events, %s', cortical_areas{area}));
    
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
                        starting_frame = contact_onsets{i}(j, 1) - min([max_window_size, contact_onsets{i}(j, 1)]) + 1;
                        starttime_1 = Kinematics.index{i}(starting_frame, 3) / 30000;
                        endtime_1 = contact_onsets{i}(j, 2) / 30000;
                    else % If first contact event does start at beginning of trial
                        starttime_1 = 0;
                        endtime_1 = 0;
                    end
                    starttime_2 = contact_offsets{i}(j, 2) / 30000;
                    ending_frame = contact_offsets{i}(j, 1) + min([max_window_size, contact_onsets{i}(j+1, 1) - contact_offsets{i}(j, 1)]);
                    endtime_2 = contact_onsets{i}(j+1, 2) / 30000;
                elseif j == size(contact_onsets{i}, 1) % Special case for last contact event
                    if contact_offsets{i}(j, 1) ~= length(graphme{i}) % If last contact event does not go to the end of trial
                        starttime_2 = contact_offsets{i}(j, 2) / 30000;
                        ending_frame = contact_offsets{i}(j, 1) + min([max_window_size, length(graphme{i}) - contact_offsets{i}(j, 1)]);
                        endtime_2 = Kinematics.index{i}(ending_frame, 3) / 30000;
                    else % If last contact event does go to the end of trail
                        starttime_2 = 0;
                        endtime_2 = 0;
                    end
                    starttime_1 = contact_offsets{i}(j-1, 2) / 30000;
                    endtime_1 = contact_onsets{i}(j, 2) / 30000;
                else % All other contact events
                    starting_frame = contact_onsets{i}(j, 1) - min([max_window_size, contact_onsets{i}(j, 1) - contact_offsets{i}(j-1, 1)]);
                    starttime_1 = Kinematics.index{i}(starting_frame, 3) / 30000;
                    endtime_1 = contact_onsets{i}(j, 2) / 30000;

                    starttime_2 = contact_offsets{i}(j, 2) / 30000;
                    ending_frame = contact_offsets{i}(j, 1) + min([max_window_size, contact_onsets{i}(j+1, 1) - contact_offsets{i}(j, 1)]);
                    endtime_2 = Kinematics.index{i}(ending_frame, 3) / 30000;
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
    
    geom_stats = []; % (mean, SD)
    
    geom_stats.index_per_trial = [];
    geom_stats.avg_fr_during = [];
    geom_stats.avg_fr_before_after = [];
    geom_stats.std_during = [];
    geom_stats.std_before_after = [];
    
    for i = valid_trials
    
        temp_1_avg = [];
        temp_1_std = [];
        temp_2_avg = [];
        temp_2_std = [];
        
        for neuron = 1:length(spiketimes_cell)

            geom_stats.index_per_trial{neuron}{i} = (fr_before_after{neuron}{i} - fr_during{neuron}{i}) ./ (fr_before_after{neuron}{i} + fr_during{neuron}{i});
            
            temp_1_avg = horzcat(temp_1_avg, mean(fr_during{neuron}{i}));
            temp_1_std = horzcat(temp_1_std, std(fr_during{neuron}{i}));
            temp_2_avg = horzcat(temp_2_avg, mean(fr_before_after{neuron}{i}));
            temp_2_std = horzcat(temp_2_std, std(fr_during{neuron}{i}));

        end
    
        geom_stats.avg_fr_during = [geom_stats.avg_fr_during; temp_1_avg];
        geom_stats.std_during = [geom_stats.std_during; temp_1_std];
        geom_stats.avg_fr_before_after = [geom_stats.avg_fr_before_after; temp_2_avg];
        geom_stats.std_before_after = [geom_stats.std_during; temp_2_std];
        
    end
    
    geom_stats.avg_fr_during = mean(geom_stats.avg_fr_during, 1, 'omitnan');
    geom_stats.avg_fr_before_after = mean(geom_stats.avg_fr_before_after, 1, 'omitnan');
    
    geom_stats.std_during = mean(geom_stats.std_during, 1, 'omitnan');
    geom_stats.std_before_after = mean(geom_stats.std_before_after, 1, 'omitnan');
    
    geom_stats.std_mean = mean([geom_stats.std_during ; geom_stats.std_before_after], 1, 'omitnan');
    
    geom_stats.geom_index(1, :) = (geom_stats.avg_fr_before_after - geom_stats.avg_fr_during) ./ (geom_stats.avg_fr_before_after + geom_stats.avg_fr_during);
    
    %% CALCULATE 95% CIs
    
%     geom_stats.total = zeros(1, length(spiketimes_cell)); 
%     
%     for i = valid_trials
%         
%         if size(contact_onsets{i}, 1) > 0
%         
%             temp_mat = zeros(size(contact_onsets{i}, 1), 1);
% 
%             for neuron = 1:length(spiketimes_cell)
% 
%                 temp_mat = horzcat(temp_mat, geom_index{neuron}{i});
% 
%             end
% 
%             temp_mat(:, 1) = [];
% 
%             geom_stats.total = [ geom_stats.total; temp_mat ];
%             
%         end
%         
%     end
    
%     geom_stats.mean(1, :) = mean(geom_stats.total, 1, 'omitnan');
%     geom_stats.std(1, :) = std(geom_stats.total, 1, 'omitnan');
    geom_stats.se(1, :) = geom_stats.geom_index ./ sqrt(geom_stats.std_mean);
    geom_stats.ci(1, :) = geom_stats.geom_index - geom_stats.se*1.96;
    geom_stats.ci(2, :) = geom_stats.geom_index + geom_stats.se*1.96;
    geom_stats.ci(3, :) = geom_stats.se*1.96;
%     geom_stats.total(1, :) = [];
    
    %% GRAPH GEOM INDICES
    
    i = 1; % change to iterate over valid_trials
    neuron = 1; % change to iterate over neurons
    
    x = 1:length(spiketimes_cell);
    
    fig = figure;
    fig = bar(x, geom_stats.geom_index);
    
    hold on
    
%     er = errorbar(x, geom_stats.geom_index, -1 * geom_stats.ci(3, :), geom_stats.ci(3, :));
    er = errorbar(x, geom_stats.geom_index, -1 * geom_stats.std_mean / sqrt(length(graphme)), geom_stats.std_mean / sqrt(length(graphme)));
    er.Color = 'black';
    er.LineStyle = 'none';
    title(sprintf('Geom Indices, %s', cortical_areas{area}));
    
    hold off
    
    %% GRAPH PER NEURON (as req'd by Fritzie)
    
    [geom_max, neuron_max] = max(geom_stats.geom_index);
    [geom_min, neuron_min] = min(geom_stats.geom_index);
    
    i = 1;
    
    fig1 = figure(1);
    hold on
    
    for i = valid_trials
    
        scatter(fr_during{neuron_max}{i}, fr_before_after{neuron_max}{i}, 'r');
        scatter(fr_during{neuron_min}{i}, fr_before_after{neuron_min}{i}, 'g');
        
    end
    
    xlabel('FR During Contact Events');
    ylabel('FR Before and After Contact Events');
    legend({sprintf('Max Geom: Neuron %i, Geom Index = %f', neuron_max, geom_max), sprintf('Min Geom: Neuron %i, Geom Index = %f', neuron_min, geom_min)})
    title(sprintf('Figure 1. All Iterations, %s', cortical_areas{area}));
    
    hold off
    
    fig3 = figure(3);    
    
    hold on
    
    scatter(geom_stats.avg_fr_during, geom_stats.avg_fr_before_after);
    xlabel('Average FR During Contact Events');
    ylabel('FR Before and After Contact Events');
    title(sprintf('Figure 3. Average FR During and After Contact Events, %s', cortical_areas{area}));
    x = max(geom_stats.avg_fr_during);
    y = max(geom_stats.avg_fr_before_after);
    plot(linspace(0, max(x, y)), linspace(0, max(x, y)));
    axis square
    
    hold off
    
    %% CLEARVARS
    
    clearvars -except area contact_* cortical_areas date graphme Kinematics max_window_size misaligned_trials offset_before_after_contact_events valid_trials
    
%% END FOR LOOP
end