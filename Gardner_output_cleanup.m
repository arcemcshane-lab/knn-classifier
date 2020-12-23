% load and clean up output workspace from Gardner
folder = 'D:\SSD Backup\MatLab Code\Classifier\'; % folder containing classifier script workspaces
cortical_areas = {'M1F', 'M1U', 'S1F', 'S1U'};
clear('k_values','perf_values')
k_values = []; % initialize optimal k
perf_values = []; % initialize performace
%%
%for i = 1:1
    %matlabfile = strcat(folder,'KNN20190227_', cortical_areas{i}, '.mat');
    %disp(strcat('Loading',matlabfile,'. . .'));
    load(strcat(folder,'knnmodel_28_S1F_euclidean_distance_k1_tw100_lag0.mat'));
    disp('Succesfully loaded');
    k_values = [k_values ; Mdl.NumNeighbors];
    perf_values =  [perf_values ; 1-cvmdlloss];
%end
%clear('starttime','endtime','ans','area','cvmdlloss','contactstable','folder','graphme','i','j','k','Kinematics','NEV','NEV_cell','oldFolder')
%clear('spikemat_temp')
%%
load(strcat('D:\SSD Backup\MatLab Code\Classifier\KNN20190227_S1F.mat'))
k_values = [k_values ; Mdl.NumNeighbors];
perf_values =  [perf_values ; 1-cvmdlloss];
% add Mdl.Distance
% Control: eucl M1F all else corrleation
% eucl M1F cityblock M1U cosine S1U mahalanobis S1F

%% Make Performance Bar Graph
%subplot(1,2,1);
%perf = [perf_values(1),perf_values(5); perf_values(2),perf_values(6);perf_values(3),perf_values(7);perf_values(4),perf_values(8)];
perf = [0.69,0.74 ;0.74,0.77 ;0.61 ,0.73 ; 0.64,0.73 ];
b = bar(perf);
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(b(1).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
xtips2 = b(2).XEndPoints;
ytips2 = b(2).YEndPoints;
labels2 = string(b(2).YData);
text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
xlabel('Cortical Region')
ylabel('Performance')
xticks([1 2 3 4])
ylim([0 1]);
grid on
xticklabels({'M1F(≈36 units)'; 'M1U(≈190 units)'; 'S1F(≈27 units)'; 'S1U(≈33 units)'})
yline((1/2),'-r','LineWidth',2.0);
legend('Control','All Nerve Blocks','Chance = 1/2')
title('KNN Classifier Performance of Tongue Tip Contact by Cortical Region and Block Status')

%%
plot(nonerveblockprior,nerveblockpriors)

x = [];
for i = 1:length(nonerveblockprior)
    x(i) = nonerveblockclassnames(i) == nerveblockclassnames(i);
end
%subplot(1,2,2);
%title('Optimal K')
%uitable('Data', k_values, 'RowName', {'M1F(Control)'; 'M1U(Control)'; 'S1F(Control)'; 'S1U(Control)';'M1F(All Blocks)'; 'M1U(All Blocks)'; 'S1F(All Blocks)'; 'S1U(All Blocks)'},'ColumnName',{'Optimal K'}, 'Position', [100 20 500 300]);
