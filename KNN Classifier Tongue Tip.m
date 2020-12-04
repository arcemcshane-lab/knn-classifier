%% RUN EACH SECTION INDIVIDUALLY 
load('20190227graphs.mat') %hits for 6 markers for 2000 ms window trials
load('20190227_CranialKinematics.mat') %fully digitized Rocky trials
load('20190227contactbyregionsallmarkers.mat')
NEV = load('Ry201902270001.mat'); %electrode recordings of spikes by time
NEV_cell = struct2cell(NEV); %converts to cell for easier indexing
%load('RyNeuronSelectSorted.mat')
%knearestneighbors = 111;

%% POPULATE CONTACTS
Contact.trialnames = Kinematics.trialnames;

for i = 1:length(Kinematics.points) % 1 to 78
    trialnumber = i; % 1-78
    marker = 1; % 1-6
    partial = 0.5; %0-1 CHANGE THIS TO ADJUST MOMENT OF CONTACT

    hits = sum(graphme{trialnumber}(:,marker),2);

    [pks,locs] = findpeaks(hits, 'MinPeakHeight',100,'MinPeakDistance',50); %Peak coordinates

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
        [pk_min, loc_min] = max(hits_inverse(start_search:end_search, 1));
        waves(j, 1) = loc_min + start_search;
        
        % Look ahead
        start_search = locs(j);
        end_search = min([max([locs(j), locs(j)+20]), size(hits_inverse, 1)]); % 0<locs(j)<locs(j+20);
        [pk_min, loc_min] = max(hits_inverse(start_search:end_search, 1));
        waves(j, 5) = loc_min + start_search;
    end
    
    Contact.waves{i} = waves; % matrix for 1 point convert to cells for multiple
    Contact.stats{i,1} = length(waves); % number of contact events
    Contact.stats{i,2} = mean(waves(:,5)-waves(:,1)); % mean contact frames
    Contact.stats{i,3} = median(waves(:,5)-waves(:,1)); % median contact frames
end

for i = 1:length(Kinematics.points) %1:length(Kinematics.points) % 1 to 78
    for j = 1:size(Contact.waves{i}, 1) % however many waves in each trial
        clear('vector')
        vector = Kinematics.points{i}(Contact.waves{i}(j,1):Contact.waves{i}(j,5)-1,28:30); % 3 is trial 28:30 is x,y,z of 10
        Contact.vectors{i}{j} = vector;
    end
end

% Contact = rmfield(Contact,'vectors')

% Determine Location

for i = 1:length(Contact.vectors) % 78 trials
    for j = 1:length(Contact.vectors{i}) % however many trajectories
        if size(Contact.vectors{i},2) < 4
            Contact.vectors{i}{j} = horzcat(Contact.vectors{i}{j},zeros(length(Contact.vectors{i}{j}),1)); %add column of zeros
        end
        for k = 1:size(Contact.vectors{i}{j},1)
            Contact.vectors{i}{j}(k,4) = palate_locator(Contact.vectors{i}{j}(k,1),Contact.vectors{i}{j}(k,3)); %add contact region
        end
    end
end

clear('i','j','k','locs','pks','hits','pk_min', 'loc_min', 'start_search', 'end_search', 'hits_inv','nframes','trialnumber','vector','waves');

%% COLLECT TRAINING DATASETS(CONTACT IDENTITY)
trainingtrials = [];

for i = 1:length(graphme)%64:70%:length(Contact.vectors)-8 % 1 to 78, no 45 or 58, 64
    trial = i;
    for j = 1:length(Contact.vectors{trial}) % e.g. 1 to 6
        if any(Contact.vectors{trial}{j})
            trainingtrials = [trainingtrials, mode(Contact.vectors{1,trial}{1,j}(:,4))];
        end
    end
end

knntable = zeros(length(trainingtrials),length(NEV_cell)+2);
knntable(:,1) = trainingtrials';
%% FIND TIME FOR EACH CONTACT EVENT
knnindex = 1;
% for i = 1:40%length(nintypercent)%length(Contact.waves)-8 % 1 to 70
for i = 1:length(graphme)
    trial = i;
    for j = 1:size(Contact.waves{trial}, 1) % 1 to 6
        if any(Contact.waves{trial}(j,3)) % checks middle wave
            % look up in Kinematics
            knntable(knnindex,2) = Kinematics.index{trial}(find(Kinematics.index{trial}(:,2) == Contact.waves{trial}(j,3)),3);
            knnindex = knnindex+1;
        end
    end
end

%% PREPARE NEURAL DATA
spikenames = fields(NEV); %Extracts Spike Data from Struct and Sorts by Channel
for i = 1:length(spikenames)
    channel = NEV.(string(spikenames(i)));
    for j = 1:length(unique(channel(:,2)))
        unit_rows_inx = find(channel(:,2) == j); %takes slices with the correct unit number
        
        ch_ID = char(spikenames(i));
        
        spiketimes.(string(sprintf('%s_%03d_%02d', string(ch_ID(1:3)), i, j))) = channel(unit_rows_inx,3)';
    end
end

clear('ch_ID','channel','i','j','unit_rows_inx')
spiketimes_cell = struct2cell(spiketimes);
%% COUNT SPIKES FOR 51 MS WINDOW(NEURAL DATA FOR EACH CONTACT EVENT)
for i = 1:length(knntable) % each contact event
    timestamp = knntable(i,2)/30000;
    for j = 1:256 % each channel
        knntable(i,j+2) = length(find(spiketimes_cell{j} >= timestamp-25 & spiketimes_cell{j} <= timestamp+25));
    end
end
% clear('i', 'j', 'timestamp')
% Each row is the palate area, then timestamp for tongue tip, then 256
% channels of spikes
%% POPULATE CONTACTS
% TODO: PREALLOCATE

contactstable = zeros(1, 6);

for i = 1:length(contactbyregionsallmarkers) % each contact event
% for i = 1:2
    for j = 1:length(contactbyregionsallmarkers{i})-5*2
        contactstable = vertcat(contactstable, any(contactbyregionsallmarkers{i}(j:j+10, 13:18), 1));
    end
end

contactstable(1, :) = [];

%% COUNT SPIKES FOR EACH 50MS SLIDING WINDOW

spiketable = zeros(1, size(spiketimes_cell, 1));
spikemat_temp = zeros(1, size(spiketimes_cell, 1));

% x = 0;

width = 0.05; % sliding window width

for i = 1:length(graphme)
% for i = 1:2
    starttime = Kinematics.index{i}(1, 3) / 30000;
    endtime = Kinematics.index{i}(size(Kinematics.index{i}, 1), 3) / 30000;
    spikemat_temp = zeros(1, size(spiketimes_cell, 1));
    
    winstart = starttime:0.005:endtime - width;
    
    for j = 1:size(winstart, 2)
              
    
        for k = 1:size(spiketimes_cell, 1)
            spikemat_temp(j, k) = length(find(spiketimes_cell{k} >= winstart(j) & spiketimes_cell{k} <= winstart(j) + width));
        end
        
%         x = x+1;
        
    end
    
    spiketable = vertcat(spiketable, spikemat_temp);
end

spiketable(1, :) = [];
 
%% TRAIN ClassificationKNN Model
rng(42) % set seed for reproduction
responsetable = int2str(contactstable);
% modeltable = array2table(knntable);
% modeltable(:, 2) = []; % remove neuronal timesteps from model training data
% Mdl = fitcknn(modeltable, 'knntable1', 'OptimizeHyperparameters', 'auto',...
Mdl = fitcknn(spiketable, responsetable, 'OptimizeHyperparameters', 'auto',...
    'HyperparameterOptimizationOptions',...
    struct('MaxObjectiveEvaluations',256, 'Verbose', 0, 'Repartition', true)); % train classifier model
date = datestr(datetime(now, 'ConvertFrom', 'datenum'), 'mm_dd_yy_HHMM');
filename = strcat('trained_models/knnmodel_', date);
saveLearnerForCoder(Mdl, filename); % save labeled classifier model
disp(filename) % print filename for easy copy pasting, needed for next section
clear('Mdl') % do not clear date to call same function again, or specify different date manually

%% GENERATE CLASSIFIER CODE
% % Must MANUALLY change name of saved knnmodel in knnpredictor.m
% codegen knnpredictor -args {coder.typeof(knntable, [Inf, 256], [1, 0])} 
%     % generate C/C++ code to reproduce generator
% %% USE CLASSIFIER FOR PREDICTIONS
% %rng(42);
% input_data = rand(10, 256) * 150; % change to something meaningful later
% predictions = knnpredictor_mex(input_data);
%% Functions
%% DETERMINE PALATE COORDINATES
function [n] = palate_locator(x,z)
if x >= 3.6586
    if z >= -0.090
        n = 13;
    else
        n = 14;
    end
elseif x >= 2.5469
    if z >= -0.0711
        n = 15;
    else
        n = 16;
    end
elseif x >= 1.4064
    if z >=-0.0448
        n = 17;
    else
        n = 18;
    end
else
    n = NaN;
end
end
