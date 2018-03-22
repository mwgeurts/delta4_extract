function results = ComputeD4Stats(varargin)


%% Load input arguments
% Initialize return array
results = struct;

% Set default alpha for hypothesis testing
results.Alpha = 0.05;

% Set default number of datasets
results.Datasets = 1;

% Define variables to be statistically tested
results.Variables = {'GammaPassRateGlobal', 'GammaMeanGlobal', ...
    'GammaPassRateLocal', 'GammaMeanLocal',  'GammaPassRate', ...
    'MedianAbsDiff'};

% Define groups to be tested
results.Groups = {'Machine', 'Energy', 'PlanClass', 'Orientation'};

% If varargin{1} is a structure, append to results
if nargin > 0 && isstruct(varargin{1})
    names = fieldnames(varargin{1});
    for i = 1:length(names)
        results.(names{i}) = varargin{1}.(names{i});
    end
    
% Otherwise, if varargin{1} is a cell, assume it is the plan array
elseif nargin > 0 && iscell(varargin{1})
    results.Plans = varargin{1};
else
    error('ComputeD4Stats requires at least one input to execute');
end

% If additional Name/Value pair options exist, append them
if nargin > 2
    for i = 2:2:nargin
        results.(varargin{i}) = varargin{i+1};
    end
end

% Add dicom_tools submodule to search path
[path, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(path, 'libra'));

% Check if MATLAB can find mlochuber
if exist('mlochuber', 'file') ~= 2

    % If not, throw an error
    if exist('Event', 'file') == 2
        Event(['The libra submodule does not exist in the search path. ', ...
            'Use git clone --recursive or git submodule init followed by git ', ...
            'submodule update to fetch all submodules'], 'ERROR');
    else
        error(['The libra submodule does not exist in the search path. ', ...
            'Use git clone --recursive or git submodule init followed by git ', ...
            'submodule update to fetch all submodules']);
    end
end

%% Prepare data
% Verify results contains at least one plan
if size(results.Plans,1) < 1
    if exist('Event', 'file') == 2
        Event('No plans were found with matching dose data', 'WARN');
    else
        warning('No plans were found with matching dose data');
    end
    return;
end

% If there are more than one of each variable, append 'A' to the above 
% parameters (individual statistics are only reported for the first 
% reference dataset)
if results.Datasets > 1
    for i = 1:length(results.Variables)
        results.Variables{i} = [results.Variables{i}, char(1 + 64)];
    end
end

% Remove any variables that arent in the plans table
for i = 1:length(results.Variables)
    if ~ismember(results.Variables{i}, results.Plans.Properties.VariableNames)
        results.Variables{i} = '';
    end
end
results.Variables = results.Variables(~cellfun(@isempty, results.Variables));

% Define unique values for each group to be tested
groupVals = cell(1,4);
indices = cell(1,4);
for i = 1:length(results.Groups)
    
    % If the table column doesn't exist, ignore it
    if ~ismember(results.Groups{i}, results.Plans.Properties.VariableNames)
        groupVals{i} = cell(0);
        
    % Otherwise grab unique values    
    else
        groupVals{i} = unique(results.Plans.(results.Groups{i}));
    end
    
    % Also tore indices referring to each unique variable (this is used to
    % compute the combination array)
    indices{i} = 1:length(groupVals{i});
end

% Remove groups with only one unique value
results.Groups = results.Groups(cellfun(@length, groupVals) > 1);
groupVals = groupVals(cellfun(@length, groupVals) > 1);

% Compute all combinations of results.Groups
groupDims = cell(0);
for i = 1:length(results.Groups)
    a = nchoosek(1:length(results.Groups),i);
    for j = 1:size(a,1)
        groupDims{length(groupDims)+1} = a(j,:);
    end
end

% Compute permutations of each combination
e = cell(length(groupDims),length(results.Groups));
for i = 1:length(groupDims)
    for j = 1:length(groupDims{i})
        e{i,groupDims{i}(j)} = 1:length(groupVals{groupDims{i}(j)});
    end
end
e(cellfun(@isempty, e)) = {0};
groupCombs = zeros(1,length(results.Groups));
for i = 1:size(e,1)
    c = cell(1, length(e(i,:)));
    [c{:}] = ndgrid(e{i,:});
    c = cellfun(@(x) x(:), c, 'uniformoutput', false);
    groupCombs = [groupCombs; [c{:}]]; %#ok<AGROW>
end

%% Compute group statistics tables
% Initialize array
results.Stats = array2table(zeros(0,14), ...
    'VariableNames', {'Parameter', 'Group', 'N', 'Min', 'Median', 'Mean', ...
    'TrimMean', 'HuberM', 'Max', 'SD', 'MAD', 'Lilliefors', 'Skewness', ...
    'Kurtosis'});

% Temporarily suppress warnings (mainly for the Lilliefors test)
w = warning('off','all');

% Loop through each permutation
for i = 1:size(groupCombs,1)
    
    % Initialize index vector and group name
    match = true(size(results.Plans,1),1);
    g = cell(0);
    
    % Loop through each group
    for j = 1:size(groupCombs,2)
        
        % If non-zero, match only provided group
        if groupCombs(i,j) > 0
            if isnumeric(groupVals{j}(groupCombs(i,j)))
                match = match & (results.Plans.(results.Groups{j}) == ...
                    groupVals{j}(groupCombs(i,j)));
                g{length(g)+1} = [results.Groups{j}, '=', ...
                    num2str(groupVals{j}(groupCombs(i,j)))];
            else
                match = match & strcmp(results.Plans.(results.Groups{j}), ...
                    groupVals{j}(groupCombs(i,j)));
                g{length(g)+1} = [results.Groups{j}, '=', ...
                    groupVals{j}{groupCombs(i,j)}];
            end
        end
    end
    
    % If everything is matched, use 'all' for group name
    if isempty(g)
        g{1} = 'all';
    end
    
    % Loop through each parameter
    for j = 1:length(results.Variables)

        % If this parameter is a cell array
        if iscell(results.Plans.(results.Variables{j})(1))
            x = cell2mat(permute(results.Plans.(results.Variables{j})(match),[3 2 1]));
            if isempty(x)
                continue;
            end
            p = zeros(size(x,1),size(x,2));
            for k = 1:size(x,1)
                for l = 1:size(x,2)
                    try
                        [~, p(k,l)] = lillietest(squeeze(x(k,l,:)));
                    catch
                        p(k,l) = NaN;
                    end
                end
            end
            s = skewness(x,0,3);
            k = kurtosis(x,0,3);
            h = mlochuber(squeeze(x(k,l,:)), ...
                'k', 50, 'loc', 'median', 'sca', 'mad');
            results.Stats = [results.Stats; {results.Variables{j}, strjoin(g, ','), ...
                length(x), min(x,[],3), median(x,3), mean(x,3), ...
                trimmean(x,0.2,3), h, max(x,[],3), std(x,0,3), p, ...
                mad(squeeze(x(k,l,:))), s, k}];
        else
            x = results.Plans.(results.Variables{j})(match);
            if isempty(x)
                continue;
            end
            try
                [~, p] = lillietest(x);
            catch
                p = NaN;
            end
            s = skewness(x,0);
            k = kurtosis(x,0);
            h = mlochuber(x, 'k', 50, 'loc', 'median', 'sca', 'mad');
            results.Stats = [results.Stats; {results.Variables{j}, ...
                strjoin(g, ','), length(x), min(x), median(x), mean(x), ...
                trimmean(x,0.2),  h, max(x), std(x), mad(x), p, s, k}];
        end
    end
end

% Restore original warning state
warning(w);

%% Compute Levene's test for equal variance
results.Levene = array2table(zeros(length(results.Groups), ...
    length(results.Variables)), 'VariableNames', results.Variables, ...
    'RowNames', results.Groups);
for i = 1:length(results.Groups)
    if exist('Event', 'file') == 2
        Event(['Computing Kruskal-Wallis tests for group ', results.Groups{i}]);
    end
    for j = 1:length(results.Variables)
        try
            if iscell(results.Plans.(results.Variables{j})(1))
                x = cell2mat(permute(results.Plans.(results.Variables{j})(...
                    ismember(results.Plans.(results.Groups{i}), groupVals{i})),...
                    [3 2 1]));
                results.Levene{i,j} = vartestn(squeeze(x(1,1,:)), ...
                    results.Plans.(results.Groups{i})(ismember(results.Plans.(...
                    results.Groups{i}), groupVals{i})), 'Display', 'off', ...
                    'TestType', 'LeveneQuadratic');
            else
                results.Levene{i,j} = vartestn(results.Plans.(results.Variables{j})(...
                    ismember(results.Plans.(results.Groups{i}), groupVals{i})), ...
                    results.Plans.(results.Groups{i})(ismember(results.Plans.(...
                    results.Groups{i}), groupVals{i})), 'Display', 'off', ...
                    'TestType', 'LeveneQuadratic');
            end
        catch
            if exist('Event', 'file') == 2
                Event(['Error computing Levene statistic for ', ...
                   'parameter ', results.Variables{i}, ', group ', ...
                   results.Groups{i}], 'WARN');
            else
                warning(['Error computing Levene statistic for ', ...
                   'parameter ', results.Variables{i}, ', group ', ...
                   results.Groups{i}]);
            end
        end
    end
end

%% Compute plan Kruskal Wallis tests for each parameter, group
results.Kruskal = array2table(zeros(length(results.Groups), ...
    length(results.Variables)), 'VariableNames', results.Variables, ...
    'RowNames', results.Groups);
results.DunnRanks = array2table(zeros(0,11), 'VariableNames', {'Parameter', ...
    'GroupA', 'GroupB', 'RankA', 'SEA', 'RankB', 'SEB', ...
    'LowerCI', 'Diff', 'UpperCI', 'p'});
for i = 1:length(results.Groups)
    if exist('Event', 'file') == 2
        Event(['Computing Kruskal-Wallis tests for group ', results.Groups{i}]);
    end
    for j = 1:length(results.Variables)
        try
            if iscell(results.Plans.(results.Variables{j})(1))
                x = cell2mat(permute(results.Plans.(results.Variables{j})(...
                    ismember(results.Plans.(results.Groups{i}), groupVals{i})),...
                    [3 2 1]));
                [p,~,stats] = kruskalwallis(squeeze(x(1,1,:)), ...
                    results.Plans.(results.Groups{i})(ismember(results.Plans.(...
                    results.Groups{i}), groupVals{i})), 'off');
            else
                [p,~,stats] = kruskalwallis(results.Plans.(results.Variables{j})(...
                    ismember(results.Plans.(results.Groups{i}), groupVals{i})), ...
                    results.Plans.(results.Groups{i})(ismember(results.Plans.(...
                    results.Groups{i}), groupVals{i})), 'off');
            end
            results.Kruskal{i,j} = p;
            [d, m] = multcompare(stats, 'Alpha', results.Alpha, 'CType', ...
                'dunn-sidak', 'Display', 'off');
            for k = 1:size(d,1)
                results.DunnRanks = [results.DunnRanks; [results.Variables(j), ...
                    sprintf('%s=%s', results.Groups{i}, stats.gnames{d(k,1)}), ...
                    sprintf('%s=%s', results.Groups{i}, stats.gnames{d(k,2)}), ...
                    m(d(k,1),1), m(d(k,1),2), ...
                    m(d(k,2),1), m(d(k,2),2), d(k,3), d(k,4), d(k,5), ...
                    d(k,6)]];
            end
        catch
            if exist('Event', 'file') == 2
                Event(['Error computing Kruskal-Wallis statistic for ', ...
                   'parameter ', results.Variables{j}, ', group ', ...
                   results.Groups{i}], 'WARN');
            else
                warning(['Error computing Kruskal-Wallis statistic for ', ...
                   'parameter ', results.Variables{j}, ', group ', ...
                   results.Groups{i}]);
            end
        end
    end
end

%% Compute plan level multi-group ANOVA for each parameter
x = nan(size(results.Plans,1), length(results.Variables));
for i = 1:size(results.Plans,1)
    ingroup = true;
    for j = 1:length(results.Groups)
        if ~ismember(results.Plans.(results.Groups{j}), groupVals{j})
            ingroup = false;
            break;
        end
    end
    if ingroup    
        for j = 1:length(results.Variables)
            if iscell(results.Plans.(results.Variables{j})(i))
                x(i,j) = results.Plans.(results.Variables{j}){i}(1,1);
            else
                x(i,j) = results.Plans.(results.Variables{j})(i);
            end
        end
    end
end
g = cell(1, length(results.Groups));
c = [];
for i = 1:length(results.Groups) 
    g{i} = results.Plans.(results.Groups{i})(any(~isnan(x),2));
end
if results.Datasets == 1 && isfield(results.Plans, 'AbsVolume')
    g{length(g)+1} = results.Plans.AbsVolume(any(~isnan(x),2));
    terms = [results.Groups 'AbsVolume'];
    c(1) = length(g);
elseif results.Datasets > 1 && isfield(results.Plans, ['AbsVolume', char(65)])
    g{length(g)+1} = ...
        results.Plans.(['AbsVolume', char(65)])(any(~isnan(x),2));
    terms = [results.Groups 'AbsVolume'];
    c(1) = length(g);
else
    terms = results.Groups;
end
x = x(any(~isnan(x),2),:);
results.ANOVA = array2table(zeros(nchoosek(length(g),1) + ...
    nchoosek(length(g),2), length(results.Variables)), ...
    'VariableNames', results.Variables);
results.DunnMeans = array2table(zeros(0,11), 'VariableNames', {'parameter', ...
    'GroupA', 'GroupB', 'MeanA', 'SEA', 'MeanB', 'SEB', ...
    'LowerCI', 'Diff', 'UpperCI', 'p'});
for i = 1:length(results.Variables)
    if exist('Event', 'file') == 2
        Event(['Computing multi-group ANOVA for parameter ', ...
            results.Variables{i}]);
    end
    try
        [p,~,stats] = anovan(x(:,i), g, 'display', 'off', ...
            'model', 'interaction', 'continuous', c);
        results.ANOVA{:,i} = p;
        for j = 1:size(stats.terms,1)
            results.ANOVA.Properties.RowNames{j} = ...
                strjoin(terms(stats.terms(j,:) == 1), ',');
        end
        dims = groupDims;
        for j = 1:length(stats.grpnames)
            if length(stats.grpnames{j}) == 1
                for k = 1:length(dims)
                    if ~isempty(dims{k}) && ismember(j, dims{k})
                        dims{k} = {};
                    end
                end
            end
        end
        dims = dims(~cellfun(@isempty, dims));
        for j = 1:length(dims)
            [d, m, ~, n] = multcompare(stats, 'Alpha', results.Alpha, ...
                'CType', 'dunn-sidak', 'Display', 'off', ...
                'Dimension', dims{j});
            for k = 1:length(results.Groups)
                n = strrep(n, sprintf('X%i',k), results.Groups{k});
            end
            for k = 1:size(d,1)
                results.DunnMeans = [results.DunnMeans; ...
                    [results.Variables(i), n(d(k,1)), ...
                    n(d(k,2)), m(d(k,1),1), m(d(k,1),2), m(d(k,2),1), ...
                    m(d(k,2),1), d(k,3), d(k,4), d(k,5), d(k,6)]];
            end
        end
    catch
        if exist('Event', 'file') == 2
            Event(['Error computing ANOVA for parameter ', ...
                results.Variables{i}], ...
                'WARN');
        else
            warning(['Error computing ANOVA for parameter ', ...
                results.Variables{i}]);
        end
    end
end

%% Compute plan level MANOVA for each group
results.MANOVA = cell(length(results.Groups), 17);
for i = 1:length(results.Groups)
    if exist('Event', 'file') == 2
        Event(['Computing MANOVA for group ', results.Groups{i}]);
    end
    x = zeros(length(results.Plans.(results.Groups{i})(...
        ismember(results.Plans.(results.Groups{i}), groupVals{i}))), ...
        length(results.Variables));
    for j = 1:length(results.Variables)
        if iscell(results.Plans.(results.Variables{j})(1))
            y = cell2mat(permute(results.Plans.(results.Variables{j})(...
                ismember(results.Plans.(results.Groups{i}), ...
                groupVals{i})), [3 2 1]));
            x(:,j) = squeeze(y(1,1,:));
        else
            x(:,j) = results.Plans.(results.Variables{j})(...
                ismember(results.Plans.(results.Groups{i}), groupVals{i}));
        end
    end
    try
        [d,p,stats] = manova1(x, results.Plans.(results.Groups{i})(...
            ismember(results.Plans.(results.Groups{i}), groupVals{i})), ...
            results.Alpha);
        results.MANOVA{i,1} = d;
        results.MANOVA{i,2} = p;
        results.MANOVA{i,3} = stats.W;
        results.MANOVA{i,4} = stats.B;
        results.MANOVA{i,5} = stats.T;
        results.MANOVA{i,6} = stats.dfW;
        results.MANOVA{i,7} = stats.dfB;
        results.MANOVA{i,8} = stats.dfT;
        results.MANOVA{i,9} = stats.lambda;
        results.MANOVA{i,10} = stats.chisq;
        results.MANOVA{i,11} = stats.chisqdf;
        results.MANOVA{i,12} = stats.eigenval;
        results.MANOVA{i,13} = stats.eigenvec;
        results.MANOVA{i,14} = stats.canon;
        results.MANOVA{i,15} = stats.mdist;
        results.MANOVA{i,16} = stats.gmdist;
        results.MANOVA{i,17} = stats.gnames;
    catch
        if exist('Event', 'file') == 2
            Event(['Error computing MANOVA for group ', ...
                results.Groups{i}], 'WARN');
        else
            warning(['Error computing MANOVA for group ', results.Groups{i}]);
        end
    end
end
results.MANOVA = cell2table(results.MANOVA, 'RowNames', results.Groups, ...
    'VariableNames', {'Dimension', 'p', 'W', 'B', 'T', 'dfW', 'dfB', 'dfT', ...
    'lambda', 'chisq', 'chisqdf', 'eigenval', 'eigenvec', 'canon', 'mdist', ...
    'gmdist', 'gnames'});

%% Compute pairwise t-tests
for r = 2:results.Datasets
    
    % Initialize results table
    results.(['pairwise', char(r+64)]) = array2table(zeros(0,9), ...
        'VariableNames', {'parameter', 'group', 'N', 'se', ...
        'lowerCI', 'diff', 'upperCI', 'tstat', 'p'});
    
    % Loop through each permutation
    for i = 1:size(groupCombs,1)

        % Initialize index vector and group name
        match = true(size(results.Plans,1),1);
        g = cell(0);

        % Loop through each group
        for j = 1:size(groupCombs,2)

            % If non-zero, match only provided group
            if groupCombs(i,j) > 0
                if isnumeric(groupVals{j}(groupCombs(i,j)))
                    match = match & (results.Plans.(results.Groups{j}) == ...
                        groupVals{j}(groupCombs(i,j)));
                    g{length(g)+1} = [results.Groups{j}, '=', ...
                        num2str(groupVals{j}(groupCombs(i,j)))];
                else
                    match = match & strcmp(results.Plans.(results.Groups{j}), ...
                        groupVals{j}(groupCombs(i,j)));
                    g{length(g)+1} = [results.Groups{j}, '=', ...
                        groupVals{j}{groupCombs(i,j)}];
                end
            end
        end

        % If everything is matched, use 'all' for group name
        if isempty(g)
            g{1} = 'all';
        end

        % Loop through each parameter
        for j = 1:length(results.Variables)

            % If this parameter is a cell array
            if iscell(results.Plans.(results.Variables{j})(1))
                
                % Store matched refA and refB parameters
                x = cell2mat(permute(results.Plans...
                    .(results.Variables{j})(match),[3 2 1]));
                y = cell2mat(permute(results.Plans...
                    .([results.Variables{j}(1:end-1), ...
                    char(r+64)])(match),[3 2 1]));
                
                % Initialize arrays to store parameters 
                s = nan(size(x,1),size(x,2),6);
                for k = 1:size(x,1)
                    for l = 1:size(x,2)
                        try
                            [~, p, ci, stats] = ttest(x(k,l,:), y(k,l,:), ...
                                'Alpha', results.Alpha, 'Dim', 3, 'Tail', ...
                                'both');
                            s(k,l,1) = stats.sd/sqrt(length(x));
                            s(k,l,2) = ci(1);
                            s(k,l,3) = mean(x(k,l,:)-y(k,l,:),3);
                            s(k,l,4) = ci(2);
                            s(k,l,5) = stats.tstat;
                            s(k,l,6) = p;
                        catch
                        end
                    end
                end
                results.(['pairwise', char(r+64)]) = ...
                    [results.(['pairwise', char(r+64)]); ... 
                    {results.Variables{j}, strjoin(g, ','), length(x), ...
                    s(:,:,1), s(:,:,2), s(:,:,3), s(:,:,4), s(:,:,5), ...
                    s(:,:,6)}];
            else
                x = results.Plans.(results.Variables{j})(match);
                y = results.Plans.([results.Variables{j}(1:end-1),...
                    char(r+64)])(match);
                try
                    [~, p, ci, stats] = ttest(x, y, 'Alpha', results.Alpha, ...
                        'Tail', 'both');
                catch
                    ci = [NaN NaN];
                    p = NaN;
                    stats.tstat = NaN;
                    stats.sd = NaN;
                end
                results.(['pairwise', char(r+64)]) = ...
                    [results.(['pairwise', char(r+64)]); ...
                    {results.Variables{j}, strjoin(g, ','), length(x), ...
                    stats.sd/sqrt(length(x)), ci(1), mean(x-y), ci(2), ...
                    stats.tstat, p}];
            end
        end
    end
end

% Clear temporary variables
clear a c d e g i j k m n p r s w x dims indices ingroup match names stats  ...
    terms groupDims groupVals groupCombs;
