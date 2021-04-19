function tbl = DetermineHits(kin, ref, threshold)
% Takes a single xyz coordinate, a set of labeled locator coordinates and
% the desired distance in mm and returns a tabulation of hits.
    arguments
        kin double
        ref double
        threshold double
    end
    tbl = tabulate(ref(find(vecnorm(ref(:,2:4) - kin(1,:),2,2) < threshold),1));
end