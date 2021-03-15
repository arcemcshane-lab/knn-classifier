%% LOAD DATE AND TRIAL START AND END
datapath = 'C:\Users\Derrick Tang\Documents\MATLAB\20190227'; % wherever you keep your neural and kinematic data
cd(datapath);
dates = {'28','27'};
fr = cell(2,4); % first row is control, second is nb
cortical_areas = {'M1F'; 'M1U'; 'S1F'; 'S1U'};
regions = {'Caudal M1o', 'Rostral M1o', 'Areas 1/2', 'Areas 3a/3b'};

for k = 1:length(dates)
    date = dates{k};
    load(strcat('201902', date, '_Kinematics.mat')) % change this

    startend = zeros(length(Kinematics.NeuralIndex),2);
    for i = 1:length(Kinematics.NeuralIndex)
        startend(i,1) = min(Kinematics.NeuralIndex{i}(:,3)); % in seconds
        startend(i,2) = max(Kinematics.NeuralIndex{i}(:,3));
    end

    %%
    for area = 1:length(cortical_areas)
        NEV = load(strcat('201902', date, '_', cortical_areas{area}, '_sortedspikes.mat'));
        NEV = struct2cell(NEV); %converts to cell for easier indexing
        meanfr = zeros(length(NEV),length(Kinematics.NeuralIndex));
        for i = 1:length(NEV)
            for j = 1:length(Kinematics.NeuralIndex)
                meanfr(i,j) =  length(find(NEV{i}.times > startend(j,1) & NEV{i}.times < startend(j,2)));
                meanfr(i,j) = meanfr(i,j) / max(Kinematics.NeuralIndex{j}(:,3)); % spike count divided by frames
            end
        end
        fr{k,area} = meanfr;
    end
end

clearvars -except fr date Kinematics cortical_areas NEV startend regions

%% SCATTERPLOT
figure
for i = 1:size(fr,2)
    
    subplot(2,2,i)

    cy = ([1:size(fr{1,i},1)]);
    cx = zeros(size(fr{1,i},1),1);
    for j = 1:size(fr{1,i},1)
        cx(j) = mean(nonzeros(fr{1,i}(j,:)));
    end
    
    scatter(cx, cy, 'b','.')
    
    ny = ([1:size(fr{2,i},1)]);
    nx = zeros(size(fr{2,i},1),1);
    for j = 1:size(fr{2,i},1)
        nx(j) = mean(nonzeros(fr{2,i}(j,:)));
    end
    hold on
    scatter(nx, ny,'r' ,'.')
    title(regions{i});
    xlabel('Firing rate in spikes/second')
    ylabel('Unit number')
end
sgtitle('Overall Firing Rate by Cortical Region Irrespective of Kinematics')
legend({'Control','Nerve Block'})
%% HISTOGRAM
figure
for i = 1:size(fr,2)
    
    subplot(2,2,i)
    
    y = [0:ceil(max(cx)*1000)/20000:ceil(max(cx)*1000)/1000];    
    
    cx = zeros(size(fr{1,i},1),1);
    for j = 1:size(fr{1,i},1)
        cx(j) = mean(nonzeros(fr{1,i}(j,:)));
    end
    
    histogram(cx,y)
    
    nx = zeros(size(fr{2,i},1),1);
    for j = 1:size(fr{2,i},1)
        nx(j) = mean(nonzeros(fr{2,i}(j,:)));
    end
    
    hold on
    histogram(nx,y)
    
    title(regions{i});
    xlabel('Firing rate in spikes/second')
    ylabel('Unit number')
end
sgtitle('Overall Firing Rate by Cortical Region Irrespective of Kinematics')
legend({'Control','Nerve Block'})
%% BOXPLOT
figure
for i = 1:size(fr,2)
    
    subplot(2,2,i)

    cx = zeros(size(fr{1,i},1),1);
    for j = 1:size(fr{1,i},1)
        cx(j) = mean(nonzeros(fr{1,i}(j,:)));
    end
    g1 = repmat({'Control'}, size(fr{1,i},1),1);
    
    %boxplot(cx)
    
    ny = ([1:size(fr{2,i},1)]);
    nx = zeros(size(fr{2,i},1),1);
    for j = i:size(fr{2,i},1)
        nx(j) = mean(nonzeros(fr{2,i}(j,:)));
    end
    g2 = repmat({'Nerve Block'},size(fr{2,i},1),1);
    x = [cx; nx];
    
    g = [g1; g2];
    
    hold on
    boxplot(x,g)
    
    title(regions{i});
    ylabel('Firing rate in spikes/second')
end
sgtitle('Overall Firing Rate by Cortical Region Irrespective of Kinematics')