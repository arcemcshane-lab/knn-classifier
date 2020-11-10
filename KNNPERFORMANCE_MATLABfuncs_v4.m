%% LOAD DATA
date = num2str(28);
cortical_areas = {'M1F'; 'M1U'; 'S1F'; 'S1U'};
load(strcat('201902', date, 'graphs.mat')) %hits for 6 markers for 2000 ms window trials
load(strcat('201902', date, '_CranialKinematics.mat')) %fully digitized Rocky trials
load(strcat('201902', date, 'contactbyregionsallmarkers.mat'))
% rng(42) % set seed for reproduction
sorting = false;
use_parallel = false;

tic
%% BEGIN FOR LOOP
for area = 1:length(cortical_areas)
% for area = 1:1
    NEV = load(strcat('201902', date, '_', cortical_areas{area}, '_sortedspikes.mat'));
    NEV_cell = struct2cell(NEV); %converts to cell for easier indexing
    
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
    %% MAKE SPIKETIMES
    
    spikenames = fields(NEV);
    for i = 1:length(NEV_cell)
        spiketimes.(string(spikenames(i))) = NEV_cell{i}.times;
    end
    spiketimes_cell = struct2cell(spiketimes);
   
    %% POPULATE CONTACTS AND SPIKES
%     contactstable = zeros(1, 6);
%     spiketable = zeros(1, size(spiketimes_cell, 1));
    contactstable = zeros(total_rows, 6);
    spiketable = zeros(total_rows, size(spiketimes_cell, 1));
    spikemat_temp = zeros(1, size(spiketimes_cell, 1));
    contacts_temp = zeros(1, 6);
    idx_contacts = 0;
    idx_spikes = 0;
    cum_frames = cumsum(total_frames);

    if use_parallel == true
        parfor i = 1:length(valid_trials) % each contact event
    %     for i = 4
            idx_contacts = cum_frames(i);
            idx_spikes = cum_frames(i);
            for j = (1+offset_contacts):size(contactbyregionsallmarkers{valid_trials(i)}, 1)-(200)*width+offset_contacts
%                 contactstable = vertcat(contactstable, any(contactbyregionsallmarkers{valid_trials(i)}(j:j+10, 13:18), 1));
%                 contacts_temp = vertcat(contactstable, any(contactbyregionsallmarkers{valid_trials(i)}(j:j+10, 13:18), 1));
                contactstable(idx_contacts, :) = any(contactbyregionsallmarkers{valid_trials(i)}(j:j+10, 13:18), 1);
            end

            starttime = Kinematics.index{valid_trials(i)}(1, 3) / 30000 + offset_spikes;
            endtime = Kinematics.index{valid_trials(i)}(size(Kinematics.index{valid_trials(i)}, 1), 3) / 30000 + offset_spikes;
            spikemat_temp = zeros(1, size(spiketimes_cell, 1));

            winstart = starttime:0.005:endtime - width;

            for j = 1:size(winstart, 2)

                for k = 1:size(spiketimes_cell, 1)
%                      spikemat_temp(j, k) = length(find(spiketimes_cell{k} >= winstart(j) & spiketimes_cell{k} <= winstart(j) + width));
                    spiketable(idx_spikes, :) = length(find(spiketimes_cell{k} >= winstart(j) & spiketimes_cell{k} <= winstart(j) + width));
                end
                
                idx_spikes = idx_spikes + 1;

            end

%             spiketable = vertcat(spiketable, spikemat_temp);
        end
    else
        for i = valid_trials % each contact event
%         for i = 67
            for j = (1+offset_contacts):size(contactbyregionsallmarkers{valid_trials(i)}, 1)-(200)*width+offset_contacts
%                 contactstable = vertcat(contactstable, any(contactbyregionsallmarkers{valid_trials(i)}(j:j+10, 13:18), 1));
                idx_contacts = idx_contacts + 1;
                contactstable(idx_contacts, :) = any(contactbyregionsallmarkers{valid_trials(i)}(j:j+200*width, 13:18), 1);
                
            end

            starttime = Kinematics.index{valid_trials(i)}(1, 3) / 30000 + offset_spikes;
            endtime = Kinematics.index{valid_trials(i)}(size(Kinematics.index{valid_trials(i)}, 1), 3) / 30000 + offset_spikes;
            spikemat_temp = zeros(1, size(spiketimes_cell, 1));

            winstart = starttime:0.005:endtime - width;

            for j = 1:size(winstart, 2)
                
                idx_spikes = idx_spikes + 1;

                for k = 1:size(spiketimes_cell, 1)
%                      spikemat_temp(j, k) = length(find(spiketimes_cell{k} >= winstart(j) & spiketimes_cell{k} <= winstart(j) + width));
                    spikemat_temp(1, k) = length(find(spiketimes_cell{k} >= winstart(j) & spiketimes_cell{k} <= winstart(j) + width));
                end
                
                spiketable(idx_spikes, :) = spikemat_temp;

            end

%             spiketable = vertcat(spiketable, spikemat_temp);

        end
    end

    contactstable(1, :) = [];
    spiketable(1, :) = [];
    
     %% REMOVE BAD SPIKETIMES
    % Note: Cannot calculate this before running spiketable because will
    % not account for spiking before and after trials
    
    if sorting == true
        suprameanspike = sum(spiketable((1:10:length(spiketable)),:))/0.050; % sum every tenth frame then divide by real world time
        suprameanspike = suprameanspike/length((1:10:length(spiketable))); % take mean

        bad_spikes = [];

        for i = 1:length(suprameanspike)
            if suprameanspike(i) < 2
                bad_spikes = [bad_spikes ; i];
            end
        end
        
        spiketable(:, bad_spikes) = [];
        
    end
    %% CLASSIFY TRAINING AND TEST DATA
%     rng(42);
%     training_indxs = datasample(1:length(contactstable), floor(length(contactstable)*0.9), 'Replace', false);
    
    predictor_train = spiketable;
    response_train = contactstable;
    
    %%
    
    % TRAIN ClassificationKNN Model
    response_train = int2str(response_train);    
    
    Mdl = fitcknn(predictor_train, response_train, 'Distance', 'euclidean', 'OptimizeHyperparameters', {'NumNeighbors'},...
        'HyperparameterOptimizationOptions',...
        struct('UseParallel',use_parallel, 'MaxObjectiveEvaluations',15, 'AcquisitionFunctionName','expected-improvement-plus','ShowPlots', false, 'Verbose', 0, 'Repartition', false)); % train classifier model
  elapsedtimeknn=toc;  
  
      cvmodel=crossval(Mdl);
    cvmdlloss=kfoldLoss(cvmodel);
   
%     date = datestr(datetime(now, 'ConvertFrom', 'datenum'), 'mm_dd_yy_HHMM');
    filename = strcat('knnmodel_', date, '_', cortical_areas{area}, '_euclidean_distance');
    save(filename); % save labeled classifier model
    %disp(filename) % print filename for easy copy pasting, needed for next section
%     clear('Mdl') % do not clear date to call same function again, or specify different date manually

% GENERATE CLASSIFIER CODE
    % Must MANUALLY change name of saved knnmodel in knnpredictor.m
    % codegen knnpredictor -args {coder.typeof(spiketable, [Inf, 256], [1, 0])} 
        % generate C/C++ code to reproduce generator
    % USE CLASSIFIER FOR PREDICTIONS
    % %rng(42);
    % input_data = rand(10, 256) * 10; % change to something meaningful later
    % predictions = knnpredictor_mex(input_data);
%     ypred = predict(Mdl, spiketable);
% END FOR LOOP

%%%%%%%%%%%
%rng(1); % For reproducibility
% CVKNNMdl = crossval(KNNMdl);
% classError = kfoldLoss(CVKNNMdl)


clearvars -except date cortical_areas graphme Kinematics contactbyregionsallmarkers sorting use_parallel area

end


% for i=1:length(inx)
%     tf(i)=strcmp(label(i,:),response_train(i,:));
% end
% 
% length(find(tf))
% Functions
% DETERMINE PALATE COORDINATES


