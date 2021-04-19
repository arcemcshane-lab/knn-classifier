function ContactEvents = ContactEvents(Kinematics, locs, threshold, criteria)
%run these in command window
%datasets = 'C:\Users\Derrick Tang\Documents\MATLAB\Datasets\';
%load(strcat(datasets, '\rockylocsmm.mat'));
%load(strcat(datasets,'20190228_Kinematics.mat'));
%       TODO:
%       CE.Catalogue: Trial num, palate location, duration, trajectory,
%       angle, mean 3D position
%       CE.TypeA: Binary contact with criteria in Classifier ready format
%       CE.stats: number of CE's per marker, mean CE duration, mean CE per gape cycle
%       OOP verison
%       
arguments
    Kinematics struct
    locs double
    threshold double
    criteria double
end

superficialmarkers = {'AnteriorM_', 'AnteriorSuperficialR', 'AnteriorSuperficialL', 'IntermediateSuperficialR', 'IntermediateSuperficialM', 'IntermediateSuperficialL'};
% index marker x-y-z with find(contains(Kinematics.ColumnNames.points, superficialmarkers{1}));

ContactEvents.hits = cell(1,length(Kinematics.Cranium.points));
ContactEvents.Threshold = threshold;
ContactEvents.Criteria = criteria;
ContactEvents.TrialNames = Kinematics.TrialNames;
ContactEvents.Catalogue = cell(1,6);

for i = 1:1%length(Kinematics.Cranium.points)
    for j = 1:length(superficialmarkers)
        marker = Kinematics.Cranium.points{j}(:,find(contains(Kinematics.ColumnNames.points, superficialmarkers{j})));
        for k = 1:size(marker,1)
            if any(isnan(marker(k,:)))
                ContactEvents.hits{i}{k,j} = NaN;
            else
                ContactEvents.hits{i}{k,j} = DetermineHits(marker(k,:),locs,threshold);
            end
        end
        
    end
    ContactEvents.graphme{i} = ContactEvents.hits{i};
    % write a helper function to do this
    for j = 1:size(ContactEvents.graphme{i},1) % each frame
        for k = 1:size(ContactEvents.graphme{i},2) % each marker
            if isempty(ContactEvents.graphme{i}{j,k})
                ContactEvents.graphme{i}{j,k} = 0;
            elseif isnan(ContactEvents.graphme{i}{j,k})
                ContactEvents.graphme{i}{j,k} = 0;
            else
                ContactEvents.graphme{i}{j,k} = sum(ContactEvents.graphme{i}{j,k}(:,2));
            end
        end
    end
    ContactEvents.graphme{i} = cell2mat(ContactEvents.graphme{i});
    ContactEvents.graphme{i}(find(isnan(ContactEvents.graphme{i}))) = 0;
    
    for k = 1:size(ContactEvents.graphme{i},2) % each marker
        
        %[~,x] = max(ContactEvents.hits{1}{464,1}(:,2));
        hits = ContactEvents.graphme{i}(:,k);
        
        [pks,locs] = findpeaks(hits, 'MinPeakHeight',100,'MinPeakDistance',50); %Peak coordinates
        Peaky = pks;
        Peakx = locs;
        Starts = zeros(size(locs,2),1);
        Ends = zeros(size(locs,2),1);
        
        for l = 1:length(locs) % Assign Start and End
            for m = 1:locs(l)
                if hits(locs(l)-m+1) - hits(locs(l)-m) == 0 % Scan backwards
                    Starts(l,1) = locs(l)-m+1; % Start
                    break
                end
            end
            for m = 1:length(hits)-locs(l)-1
                if hits(locs(l)+m+1) - hits(locs(l)+m) == 0 % Scan forwards
                    Ends(l,1) = locs(l)+m; % End
                    break
                end
            end
        end
        Duration = Ends - Starts;
        Trial = repmat({Kinematics.TrialNames(i)},size(Duration,1),1);
        ContactEvents.Catalogue{k} = table(Trial,Starts,Peakx,Peaky,Ends,Duration);
    end
    
%     ContactEvents.binarycontact{i}(find(isnan(ContactEvents.binarycontact{i}))) = 0;
%     emptyIndex = cellfun('isempty', ContactEvents.binarycontact{i});
%     ContactEvents.binarycontact{i}(emptyIndex) = {0};
%     ContactEvents.contact(find(isnan(ContactEvents.contact))) = 0;
end
%%
% partial = 0.5
% for j = 1:length(locs) %Assign Onset
%     clear('x','y')
%     [x,y] = findpeaks(-1*(abs(hits(waves(j,1):waves(j,3))-(floor(hits(locs(j))/2)))),(waves(j,1):waves(j,3))); %Onset
%     if any(y)
%         waves(j,2) = y(1); %not robust to multiple peaks
%     else
%         waves(j,2) = NaN;
%     end
% end
% 
% for j = 1:length(locs) %Assign Offset
%         clear('x','y')
%     [x,y] = findpeaks(-1*(abs(hits(waves(j,3):waves(j,5))-(floor(hits(locs(j))/2)))),(waves(j,3):waves(j,5))); %Offset
%     if any(y)
%         waves(j,4) = y(1);
%     else
%         waves(j,4) = NaN;
%     end
% end
% Waves is complete
% 
% ContactEvents.waves{i} = waves; % matrix for 1 point convert to cells for multiple
% ContactEvents.stats{i,1} = length(waves); % number of contact events
% ContactEvents.stats{i,2} = mean(waves(:,3)-waves(:,1)); % mean contact frames
% ContactEvents.stats{i,3} = median(waves(:,3)-waves(:,1)); % median contact frames

end
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
% function tbl = DetermineHits(kin, ref, threshold)
% %Returns a tabulation of locators that fall within threshold
%     arguments
%         kin double % single xyz coordinate
%         ref double % locators
%         threshold double % distance in mm
%     end
%     tbl = tabulate(ref(find(vecnorm(ref(:,2:4) - kin(1,:),2,2) < threshold),1));
% end
%% https://www.mathworks.com/help/matlab/matlab_prog/techniques-for-improving-performance.html
%% LOAD ALL KINEMATIC DATA
% dates = {'20190228' '20190227'; '20190509' '20190510'};% where rows are individuals and columns are conditions
% individuals = {'Rocky' 'Yosemite'}; % storing names for figure titling purposes
% conditions = {'Control' 'NerveBlock'}; % storing names for figure titling purposes
%
% nind = length(individuals);
% ncon = length(conditions);
% data = {};
% datanames = {};
% % Now each Kinematics mat will be stored in the variable 'data'
% for i = 1:nind
%     for c = 1:ncon
%         load([dates{i,c} '_Kinematics.mat'])
%         data{i,c} = Kinematics;
%         datanames{i,c} = strcat(individuals{i},' ',conditions{c});
%     end
% end
% clear Kinematics nind ncon i c
% %%