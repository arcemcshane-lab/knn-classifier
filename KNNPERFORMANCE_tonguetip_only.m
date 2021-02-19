%% LOAD DATA
date = num2str(28);
cortical_areas = {'M1F'; 'M1U'; 'S1F'; 'S1U'};
load(strcat('201902', date, 'graphs.mat')) %hits for 6 markers for 2000 ms window trials
load(strcat('201902', date, '_CranialKinematics.mat')) %fully digitized Rocky trials
load(strcat('201902', date, 'contactbyregionsallmarkers.mat'))
% rng(42) % set seed for reproduction
sorting = false;
use_parallel = false;
select_neurons = 'all'; % either an int or 'all'

tic
%% BEGIN FOR LOOP
for area = 1:length(cortical_areas)
% for area = 1:1
    NEV = load(strcat('201902', date, '_', cortical_areas{area}, '_sortedspikes.mat'));
    NEV_cell = struct2cell(NEV); %converts to cell for easier indexing
    mult_cvmdlloss = [];
    
    %% FIND MISALIGNED TRIALS
    width = 0.05; % sliding window width (s)
    offset_contacts = 0; % sliding window offset for contactstable (s)
    offset_spikes = 0.05; % sliding window offset for spiketable (s)
    
    contacts_frames = [];
    spikes_frames = [];  
    total_frames = [];
    
    total_rows = 0;
    
    for i = 1:length(contactbyregionsallmarkers) % each contact event
      x = 0;
      y = 0;
%     for i = 3
        for j = (1+offset_contacts):size(contactbyregionsallmarkers{i}, 1)-(200)*width+offset_contacts
            x = x+1;
        end

%       for i = 3
        starttime = Kinematics.index{i}(1, 3) / 30000 + offset_spikes;
        endtime = Kinematics.index{i}(size(Kinematics.index{i}, 1), 3) / 30000 + offset_spikes;

        winstart = starttime:0.005:endtime - width;

        for j = 1:size(winstart, 2)
            y = y+1;
        end
        
     contacts_frames = [contacts_frames; x];
     spikes_frames = [spikes_frames; y];
     
     if x == y
         total_rows = total_rows + x;
         total_frames = [total_frames; x];
     end
     
    end
    
    misaligned_trials = contacts_frames - spikes_frames;
    misaligned_trials = find(misaligned_trials);
    
    valid_trials = 1:length(graphme);
    valid_trials(misaligned_trials) = [];
    
     %% REMOVE BAD SNR
    
    if sorting == true
        bad_snr = [];
        for i = 1:length(NEV_cell)
            if NEV_cell{i}.SNR < 2
                bad_snr = [bad_snr ; i];
            end
        end

        for i = 1:length(bad_snr)
            NEV_cell(bad_snr(i), :) = [];
        end   
    end
%% POPULATE CONTACTS
Contact.trialnames = Kinematics.trialnames;

for i = valid_trials
    Contact.waves{i} = [];
end

for marker = 1:6
% for marker = 6

    for i = valid_trials
%         Contact.waves{i} = [];
        trialnumber = i; % 1-78
        partial = 0.5; %0-1 CHANGE THIS TO ADJUST MOMENT OF CONTACT

        hits = sum(graphme{trialnumber}(:,marker),2);

        [~,locs] = findpeaks(hits, 'MinPeakHeight',100,'MinPeakDistance',50); %Peak coordinates

        waves = zeros(length(locs),5); %Frame of Start, Onset, Peak, Offset, End
        waves(:,3) = locs; %Assign Peak

        hits_inverse = hits * -1;

        for j = 1:length(hits_inverse)
            if hits_inverse(j) == 0
               hits_inverse(j) = NaN;
            end
        end 

        for j = 1:length(locs)
            % Look behind
            start_search = min([max([1, locs(j) - 20]), locs(j)]); % 0<locs(j-20)<locs(j)
            end_search = locs(j);
            [~, loc_min] = max(hits_inverse(start_search:end_search, 1));
            waves(j, 1) = loc_min + start_search;

            % Look ahead
            start_search = locs(j);
            end_search = min([max([locs(j), locs(j)+20]), size(hits_inverse, 1)]); % 0<locs(j)<locs(j+20);
            [~, loc_min] = max(hits_inverse(start_search:end_search, 1));
            waves(j, 5) = loc_min + start_search;
        end

        waves = [waves marker*ones(size(waves, 1), 1)]; % Add marker identity to waves

        Contact.waves{i} = [Contact.waves{i}; waves]; % matrix for 1 point convert to cells for multiple
    %     Contact.stats{i,1} = length(waves); % number of contact events
    %     Contact.stats{i,2} = mean(waves(:,5)-waves(:,1)); % mean contact frames
    %     Contact.stats{i,3} = median(waves(:,5)-waves(:,1)); % median contact frames
    end




for i = valid_trials
    for j = 1:size(Contact.waves{i}, 1) % however many waves in each trial
        clear('vector')
        vector = Kinematics.points{i}(Contact.waves{i}(j,1):Contact.waves{i}(j,5)-1,marker+27:marker+27+2); % 3 is trial 28:30 is x,y,z of 10
        Contact.vectors{i}{j} = vector;
    end
end

% Contact = rmfield(Contact,'vectors')


% Determine Location

for i = valid_trials % 78 trials
    for j = 1:length(Contact.vectors{i}) % however many trajectories
        if size(Contact.vectors{i},2) < 4
            Contact.vectors{i}{j} = horzcat(Contact.vectors{i}{j},zeros(length(Contact.vectors{i}{j}),1)); %add column of zeros
        end
        for k = 1:size(Contact.vectors{i}{j},1)
            Contact.vectors{i}{j}(k,4) = palate_locator(Contact.vectors{i}{j}(k,1),Contact.vectors{i}{j}(k,3)); %add contact region
%             Contact.vectors{i}{j}(k,5) = marker % TODO: Add marker identity
        end
    end
end

clear('i','j','k','locs','hits', 'loc_min', 'start_search', 'end_search', 'hits_inv','nframes','trialnumber','vector','waves');

%% FIND TIME FOR EACH CONTACT EVENT
knnindex = 1;
% for i = 1:40%length(nintypercent)%length(Contact.waves)-8 % 1 to 70
for i = valid_trials
    trial = i;
    for j = 1:size(Contact.waves{trial}, 1) % 1 to 6
        if any(Contact.waves{trial}(j,3)) % checks middle wave
            % look up in Kinematics
            knntable(knnindex,2) = Kinematics.index{trial}(find(Kinematics.index{trial}(:,2) == Contact.waves{trial}(j,3)),3);
            knnindex = knnindex+1;
        end
    end
end
%% COLLECT TRAINING DATASETS(CONTACT IDENTITY)
trainingtrials = [];

for i = valid_trials
    trial = i;
    for j = 1:length(Contact.vectors{trial}) % e.g. 1 to 6
        if any(Contact.vectors{trial}{j})
            trainingtrials = [trainingtrials, mode(Contact.vectors{1,trial}{1,j}(:,4))];
        end
    end
end
contactstable = trainingtrials';
%% MAKE SPIKETIMES
    
    spikenames = fields(NEV);
    for i = 1:length(NEV_cell)
        spiketimes.(string(spikenames(i))) = NEV_cell{i}.times;
    end
    spiketimes_cell = struct2cell(spiketimes);
%% CALCULATE SPIKES
spiketable = zeros(length(contactstable), length(spiketimes_cell));
for i = 1:length(contactstable) % each contact event
    timestamp = knntable(i,2)/30000;
    for j = 1:length(spiketimes_cell) % each channel
        spiketable(i,j) = length(find(spiketimes_cell{j} >= timestamp-width/2 & spiketimes_cell{j} <= timestamp+width/2));
    end
end
% clear('i', 'j', 'timestamp')
% Each row is the palate area, then timestamp for tongue tip, then 256
%% CALCULATE SPIKES 
% idx_spikes = 0;
% spiketable = zeros(total_rows, size(spiketimes_cell, 1));
% idx_spikes = 0;
%  
% for i = 1:length(valid_trials) % each contact event
%     starttime = Kinematics.index{valid_trials(i)}(1, 3) / 30000 + offset_spikes;
%     endtime = Kinematics.index{valid_trials(i)}(size(Kinematics.index{valid_trials(i)}, 1), 3) / 30000 + offset_spikes;
%     spikemat_temp = zeros(1, size(spiketimes_cell, 1));
% 
%     winstart = starttime:0.005:endtime - width;
% 
%     for j = 1:size(winstart, 2)
% 
%         idx_spikes = idx_spikes + 1;
% 
%         for k = 1:size(spiketimes_cell, 1)
% %                      spikemat_temp(j, k) = length(find(spiketimes_cell{k} >= winstart(j) & spiketimes_cell{k} <= winstart(j) + width));
% %             spikemat_temp(1, k) = length(find(spiketimes_cell{k} >= winstart(j) & spiketimes_cell{k} <= winstart(j) + width));
%         end
% 
%         x = x+1;
%         
% %         spiketable(idx_spikes, :) = spikemat_temp;
% 
%     end
% 
% %             spiketable = vertcat(spiketable, spikemat_temp);
% 
%         end
%% TRAIN ClassificationKNN Model

for run = 1:10
    if isempty(contactstable)
        cvmdlloss = 0;
    else
%     rng(42) % set seed for reproduction
    % responsetable = int2str(contactstable);
    % modeltable = array2table(knntable);
    % modeltable(:, 2) = []; % remove neuronal timesteps from model training data
    Mdl = fitcknn(spiketable, contactstable,  'Distance', 'euclidean', 'OptimizeHyperparameters', {'NumNeighbors'},...
        'HyperparameterOptimizationOptions',...
        struct('MaxObjectiveEvaluations',64, 'Verbose', 0, 'Repartition', true)); % train classifier model

    cvmodel=crossval(Mdl);
        cvmdlloss=kfoldLoss(cvmodel);
    end

mult_cvmdlloss = [mult_cvmdlloss cvmdlloss];

end

%     date = datestr(datetime(now, 'ConvertFrom', 'datenum'), 'mm_dd_yy_HHMM');
    filename = strcat('knnmodel_', date, '_', cortical_areas{area}, '_euclidean_distance');
    save(filename, 'area', 'mult_cvmdlloss'); % save labeled classifier model
    % TODO: Change to save Mdl, cvmodel
% saveLearnerForCoder(Mdl, filename); % save labeled classifier model
% disp(filename) % print filename for easy copy pasting, needed for next section
% clear('Mdl') % do not clear date to call same function again, or specify different date manually

end


%% GENERATE CLASSIFIER CODE
% % Must MANUALLY change name of saved knnmodel in knnpredictor.m
% codegen knnpredictor -args {coder.typeof(knntable, [Inf, 256], [1, 0])} 
%     % generate C/C++ code to reproduce generator
% %% USE CLASSIFIER FOR PREDICTIONS
% %rng(42);
% input_data = rand(10, 256) * 150; % change to something meaningful later
% predictions = knnpredictor_mex(input_data);
clearvars -except date cortical_areas graphme Kinematics contactbyregionsallmarkers sorting use_parallel area

end