function distancematrix = distancematrix(Kinematics,locs,markers)
% Accepts a Cranuim.points cell and an oral structure locators matrix, and
% matrix with specified tongue markers
% Returns a 3D euclidian distance matrix

npoints = length(markers);


for k = 1:size(Kinematics,1) % For each trial
    data = Kinematics{k}; % 
    nframes = length(data);
    
    distance = NaN(nframes,npoints,size(rockylocs,1)); % preallocate 3D distance matrix 2201 x 6 x Q
    
    for i = 1:npoints % for select markers
        for j = 1:nframes % for every frames
            if ~isnan(data(j,10))
                temppoint = data(j,((selectpoints(i)-1)*3+1):((selectpoints(i)-1)*3+3));
                temp = locs(:,2:4) - temppoint;
                tempdist = vecnorm(temp,2,2);
                distance(j,i,:) = tempdist;
            end
        end
    end
end

end