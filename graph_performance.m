%% LOAD CMVDLLOSS FROM EACH .MAT FILE
%areas = {'M1UOnly','M1FControl','S1UControl','S1FControl'};
areas = {'M1UOnly','M1FControl','S1UControl','S1FControl','M1UNerveBlock','M1FNerveBlock','S1UNerveBlock','S1FNerveBlock'};
perf = zeros(8,15);

for i=1:size(areas,2)
    files = dir (strcat('S:\',areas{i},'\*.mat')); % checks .mat files
    L = length(files);

    for j=1:L
        %fprintf(files(i).name)
        load(strcat('S:\',areas{i},'\',files(j).name),'cvmdlloss')
        perf(i,j) = 1-cvmdlloss;
    end
end

%% GRAPH WITH ERROR BARS
data = zeros(4,2);
SD = zeros(4,2);
x=1:8;
for i=1:4
   data(i,1) = mean(nonzeros(perf(i,:))); % control
   data(i,2) = mean(nonzeros(perf(i+4,:))); % nerveblock
   SD(i,1) = std(nonzeros(perf(i,:)));
   SD(i,2) = std(nonzeros(perf(i+4,:)));
end

bar(data)
hold on
% er = errorbar(x,data,SD,SD);
% er.Color = [0 0 0]; 
% er.LineStyle = 'none';
%%
%perf = [perf_values(1),perf_values(5); perf_values(2),perf_values(6);perf_values(3),perf_values(7);perf_values(4),perf_values(8)];
%perf = [0.69,0.74 ;0.74,0.77 ;0.61 ,0.73 ; 0.64,0.73 ];
b = bar(data,1);
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(round(b(1).YData,2));
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
xtips2 = b(2).XEndPoints;
ytips2 = b(2).YEndPoints;
labels2 = string(round(b(2).YData,2));
text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
ylabel('Performance')
xticks([1 2 3 4])
ylim([0 1]);
grid on
xticklabels({'     Rostral M1o\newline(34 units, n=11&8)'; '     Caudal M1o\newline(36 units, n=10&8)'; '       Areas 1/2\newline(33 units, n=10&8)'; '     Areas 3a/3b\newline(27 units, n=15&9)'})
% ≈
yline((1/2)^6,'-r','LineWidth',2.0);
%title('Mean KNN Classifier Performance of L-P Contact by Cortical Region and Block Status')

xBar=cell2mat(get(b,'XData')).' + [b.XOffset];  % compute bar centers
hold on
errorbar(xBar,data,SD,'k.');
legend('Control','All Nerve Blocks','Chance = (1/2)^6')
%set(get(get(errorbar,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% hold on
% er = errorbar(x,data,SD,SD);
% hold off
set(gcf,'color','w');
%% DETERMINE STATISTICAL SIGNIFICANCE
for i=1:4
    ztest(nonzeros(perf(i+4,:)),mean(nonzeros(perf(i,:))),std(nonzeros(perf(i,:))),'Alpha',0.001)
end