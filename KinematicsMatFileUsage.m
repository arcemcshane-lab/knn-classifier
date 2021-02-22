%% Example Kinematics Mat File Usage

% Written By J.D. Laurence-Chasen 2/15/2021

%% Import single dataset

% I prefer keeping all of my datasets in one folder
datapath = 'C:\Users\jdlc7\Documents\PhD\Code\Datasets';
cd(datapath); % set working directory
cs = 'Cranium'; % use data in cranial coordinate system
load('20190228_Kinematics.mat'); % Load Rocky control dataset

% Access mandible pitch data from one trial (gape)
gape = Kinematics.TMJ{1}(:,6); % 6th column is mandible pitch (gape)

% Access anterior tongue point data from one trial
% ttptcols = column indices of tongue time (anterior midline marker)
ttptcols = find(contains(Kinematics.ColumnNames.points,'AnteriorM_'));
Kinematics.Cranium.points{1}(:,ttptcols)


% Access corresponding neural time
Kinematics.NeuralIndex{1}(:,3) % 3rd column is neural time in seconds


% Access tongue tip and gape data for all chews
% chews = row indices of all chew cycles in the GapeCycleInfo file
chews = find(contains(Kinematics.GapeCycleInfo.CycleType,'Chew'));

ttdata = {};
gapedata = {};
for c = 1:length(chews)
    
    % getPointData function will access the 3D point data from the kin mat
    % for the select rows of the given cycle. It will return all columns.
    %  
    % getGapeData will do the same, but for mandibular pitch
    %
    
    tmpdata = getPointData(chews(c),Kinematics);
    ttdata{c} = tmpdata(:,ttptcols);
    gapedata{c} = getGapeData(chews(c),Kinematics);
end

%% For loops to access multiple days of data

% Rocky and yosemite control and nerve block
dates = {'20190228' '20190227'; '20190509' '20190510'};% where rows are individuals and columns are conditions
individuals = {'Rocky' 'Yosemite'}; % storing names for figure titling purposes
conditions = {'Control' 'Nerve Block'}; % storing names for figure titling purposes

nind = length(individuals);
ncon = length(conditions);
data = {};
% Now each Kinematics mat will be stored in the variable 'data'
for i = 1:nind
    for c = 1:ncon
        
        load([dates{i,c} '_Kinematics.mat'])
        data{i,c} = Kinematics;
    end
end
clear Kinematics

%% EXAMPLE plot tongue tip mean trajectory + gape for all chews

for i = 1:nind % for each individual
    
    figure; % new figure
    
    for c = 1:ncon  % for each condition
      
      subplot(1,2,c)
      hold on
      % get chews IMPORTANT--modify this line to access anysubset of
      % gapecycleinfo file. I.e. CycleNumber < 5 --> only first 5 chews in
      % each sequence
      chews = find(contains(data{i,c}.GapeCycleInfo.CycleType,'Chew'));
      % get tongue tip columns for this dataset (sometimes they change!)
      ttptcols = find(contains(data{i,c}.ColumnNames.points,'AnteriorM_'));
        
      ttdata = [];
      gapedata = [];
        
        for cy = 1:length(chews)  % for each chew    
            tmppts = getPointData(chews(cy),data{i,c}); % get point data
            tmppts = tmppts(:,ttptcols); % just tongue tip data
            tmppts = resampleData(tmppts,100); % resample to 100 pts (percent of cycle)
            ttdata(:,:,cy) = tmppts; % store in array
            % repeat for gape
            tmpgape = getGapeData(chews(cy),data{i,c});
            tmpgape = resampleData(tmpgape,100);
            gapedata(:,:,cy) = tmpgape;
        end
        
        ttmean = nanmean(ttdata(:,1,:),3); %mean of X axis (hence the 1)
        tterr = nanstd(ttdata(:,1,:),0,3); %std
        gapemean = nanmean(gapedata(:,:,:),3);
        gapeerr = nanstd(gapedata(:,:,:),0,3);
        
        % Plot and add labels
        errorbar(ttmean,tterr)
        ylabel('Tongue Tip A-P Position (mm)')
        yyaxis right % plot gape on right axis
        errorbar(gapemean,gapeerr)
        ylabel('Jaw Pitch (deg)')
        title([individuals{i} ' ' conditions{c}])
    end
end

