% https://www.mathworks.com/help/matlab/matlab_prog/techniques-for-improving-performance.html
%% LOAD ALL KINEMATIC DATA
dates = {'20190228' '20190227'; '20190509' '20190510'};% where rows are individuals and columns are conditions
individuals = {'Rocky' 'Yosemite'}; % storing names for figure titling purposes
conditions = {'Control' 'Nerve Block'}; % storing names for figure titling purposes

nind = length(individuals);
ncon = length(conditions);
data = {};
datanames = {};
% Now each Kinematics mat will be stored in the variable 'data'
for i = 1:nind
    for c = 1:ncon
        
        load([dates{i,c} '_Kinematics.mat'])
        data{i,c} = Kinematics;
        datanames{i,c} = strcat(individuals{i},' ',conditions{c});
    end
end
clear Kinematics nind ncon i c
%%