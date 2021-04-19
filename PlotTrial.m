function fig = PlotTrial(ContactEvents,trial)
% Takes CE struct and trial number and produces graph with
    arguments
        ContactEvents struct
        trial double
    end

x = (1:length(ContactEvents.graphme{trial}))/200;

fig = figure('NumberTitle','off','name',ContactEvents.TrialNames{trial},'units','normalized','outerposition',[0 0 1 1]);
xlabel('Time(ms)')
ylabel('Hits')

subplot(2,1,1); %Anterior
for i = 1:3
   y = ContactEvents.graphme{trial}(:,i);
   hold on
   plot(x,y)
end
title('Anterior Markers')
ylabel('Hits')
legend('AnteriorM','AnteriorSuperficialR','AnteriorSuperficialL')

subplot(2,1,2); %Intermediate
for i = 4:6
   y = ContactEvents.graphme{trial}(:,i);
   hold on
   plot(x,y)
end
title('Intermediate Markers')
ylabel('Hits')
legend('IntermediateSuperficialR','IntermediateSuperficialM','IntermediateSuperficialL')
end