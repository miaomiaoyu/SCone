function clusterpoints = d_doclustercorr(datamatrix, headortime, baseline, clusterthresh, npermutations)

% a generic function for doing nonparametric cluster correction of EEG data using t-tests
% based on Maris & Oostenveld (2007), J Neurosci Methods, 164: 177-190, doi: 10.1016/j.jneumeth.2007.03.024
% Inputs are:
%   datamatrix: a matrix of data with dimensions (participant,timepoint/sensor) for one-sample tests, or (participant, condition, timepoint/sensor) for paired tests
%   headortime: a flag that tells the script if you are cluster correcting across time (1) or across electrode position on the scalp (2)
%   baseline: the baseline value for t-tests (usually either 0 for ERPs or 0.5 for MVPA data)
%   clusterthresh: the alpha level for significance when comparing each cluster's t-value to the null distribution
%   npermutations: the number of permutations for creating the null distribution (at least 100, 10000 is fine)
% returns a cell array (clusterpoints) containing the indices of significant clusters
% DHB 5/10/16

s = size(datamatrix);           % work out how big the data matrix is
if length(s)==3    % this means we are doing paired t-tests, but since these are the same as one-sample tests we can just subtract the values
    datamatrix = squeeze(datamatrix(:,2,:) - datamatrix(:,1,:));
end

s = size(datamatrix);           % work out how big the data matrix is again

switch headortime               % create adjacency matrix for time points or electrodes
    case 1      % time
        adjacencymatrix = zeros(s(2));      % t*t matrix of zeros
        for n = 1:(s(2)-1)
            adjacencymatrix(n,n+1) = 1;      % adjacent time points get set to 1
            adjacencymatrix(n+1,n) = 1;      %
        end
    case 2      % head - note that this only works for the waveguard cap montage we use at York
        load Waveguard.mat;         % load in the schematic of the EEG cap
        threshdist = 0.18;          % arbitrary distance between electrodes that seems to work well in pairing up neighbouring electrodes
        
        for m = 1:64                % loop through all electrodes
            for n = 1:64            % calculate absolute distance from target electrode to each other electrode
                xy1 = lay.pos(channelmappings(m),:);
                xy2 = lay.pos(channelmappings(n),:);
                electrodedistances(m,n) = sqrt((xy1(1)-xy2(1)).^2 + (xy1(2)-xy2(2)).^2);
            end
        end
        
        adjacencymatrix = zeros(size(electrodedistances));
        adjacencymatrix(find(electrodedistances<threshdist)) = 1;   % electrode pairs that are less than the threshold distance away get set to 1
        adjacencymatrix(find(electrodedistances==0)) = 0;           % an electrode can't be paired with itself (set the major diagonal to 0)
end

for n = 1:s(2)          % now do t-tests for each individual sensor or time point
    temp = squeeze(datamatrix(:,n));
    [H,P,CI,STATS] = ttest(temp,baseline,'alpha',.05,'tail','right');
    alltvals(n) = STATS.tstat;
    allpvals(n) = P;
    allhvals(n) = H;
end

sigclustercount = 0;        % set up some empty variables and counters
clustercounter = 0;
allclusters = [];
clusterpoints = [];

for n = 1:s(2)                                  % cycle through all electrodes/time points
    if allhvals(n)                              % if electrode is significant look for others in the cluster
        clustercounter = clustercounter + 1;
        clusterlist = [];                       % clear the clusterlist variable
        clusterlist(1) = n;                     % add first electrode or time point to the cluster
        for m = 1:s(2)                          % cycle through all other electrodes or time points (including this one)
            if adjacencymatrix(n,m)             % see if the electrodes or time points are neighbours
                if allhvals(m)                  % see if the neighbouring value is significant
                    clusterlist(end+1) = m;     % add this electrode or time point to the cluster
                end
            end
        end
        allclusters{clustercounter} = clusterlist;  % put the cluster into a large cell array
    end
end

% now condense the clusters by seeing if any of them have overlapping electrodes or time points
ccount = 0;
if ~isempty(allclusters)                        % only do this if there are some clusters
    for n = 1:length(allclusters)               % do it for each cluster
        targetcluster = allclusters{n};         % allocate the target cluster to a variable
        if ~isempty(targetcluster)
            for m = (n+1):length(allclusters)   % only compare clusters we haven't got to yet, to save time
                compcluster = allclusters{m};   % allocate the comparison cluster to a variable
                if ~isempty(compcluster)
                    lia = ismember(targetcluster,compcluster);  % see if clusters have common elements
                    if sum(lia)
                        targetcluster(end+1:end+length(compcluster)) = compcluster;     % if they do, add comparison cluster to target cluster
                        allclusters{m} = [];                                            % then delete comparison cluster
                    end
                end
            end
            ccount = ccount + 1;
            condensedclusters{ccount} = unique(targetcluster);          % only unique values (so exclude repetitions of an electrode or time point)
            clustersizes(ccount) = length(unique(targetcluster));       % make a note of the cluster size
        end
    end
    
    [y,i] = max(clustersizes);                  % find the biggest cluster
    maxcluster = condensedclusters{i};          % save the cluster indices in a variable
    
    parfor n = 1:npermutations                     % repeat reshuffling npermutations times to generate a null distribution of summed t values
        allnullts = [];
        for elcounter = 1:length(maxcluster)    % do the following for each electrode or time point in the cluster
            tempdata = datamatrix(:,maxcluster(elcounter));                              % get the data that underlies the comparison for this observation
            permsigns = sign(randperm(length(tempdata))-(0.5+(length(tempdata)/2)))';    % generate a permuted list containing an equal number of +1 and -1s
            tempdiffs = (permsigns.*(tempdata - baseline)) + baseline;                   % flip half of the difference values about the baseline, then re-add the baseline
            allnullts(elcounter) = (mean(tempdiffs)-baseline)/(std(tempdiffs)/sqrt(length(tempdiffs)));  % manual calculation of t values is at least 10 times faster than calling ttest                                        % store the t value for this test
        end
        nulldist(n) = sum(allnullts);                                                    % sum the t values across all electrodes or time points in the cluster and add to the null distribution
    end
    
    for c = 1:ccount                                            % now let's see if each cluster has a significant pooled t value, compared to the null distribution
        sumtvals(c) = sum(alltvals(condensedclusters{c}));      % add up the real t values for all the elements in the cluster
        i = find(abs(sumtvals(c))<abs(nulldist));               % find all entries in the null distribution with a bigger t value
        
        if isempty(i)                                           % if there are none, the p value is 0
            clusterps(c) = 0;
        else                                                    % otherwise, the p value is the proportion of summed t values in the null distribution that are larger than the cluster summed t value
            clusterps(c) = length(i)/npermutations;
        end
        fprintf('\nClusterthresh = %.2d',clusterthresh);
        fprintf('\nClusterps = %.2d',clusterps(c));
        
        if clusterps(c)<clusterthresh                           % if the p value is smaller than the cluster threshold, store the cluster indices in a list for outputting
            sigclustercount = sigclustercount + 1;
            clusterpoints{sigclustercount} = condensedclusters{c};
        end
    end
end

end