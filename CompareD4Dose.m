function results = CompareD4Dose(varargin)
% CompareD4Dose evaluates a single file or directory of Delta4 measurements 
% to RT dose volumes, including calculation of Gamma pass rates. The 
% function accepts either one or two sets of dose volume files; if two, the 
% function will conduct a pairwise comparison for each measurement and
% report the comparison statistics. In this manner, the tool can be used to
% perform a bulk comparison of how Delta4 results change with updates to a
% beam model.
%
% During excecution, this function will call the functions ScanDICOMPath
% and LoadDICOMDose from the dicom_tools submodule and function
% ParseD4Tables from this repository. It will attempt to match the patient
% last name, first name, ID, and plan name between each Delta4 measurement
% and dose volume. If a match is found, the results are stored in the
% return structure field results.plans. The matching is then repeated for
% each beam, using the beam name and storing the results to results.beams.
%
% As part of the analysis, this function accepts a Name/Pair configuration
% option 'planClasses' with an n x 2 cell array of plan class names and
% associated case insensitive regular expressions to match to the plan 
% name. If provided, the function will attempt to categorize each 
% measurement into a plan class, then compare each plan type. Note, each 
% plan is only matched to the first class that matches the regexp.
%
% The following required and optional inputs are available for this
% function:
%   varargin{1}: file/folder/list of files containing measurement files
%       (see ParseD4tables for information on compatible formats)
%   varargin{2}: file/folder/list of files containing reference DICOM RT
%       Plan and RT Dose files for each measurement file in varargin{1}
%   varargin{3} (optional): file/folder/list of files containing reference 
%       DICOM RT Plan and RT Dose files for each measurement for pairwise
%       comparison. If additional Name/Value parameters are passed to this
%       function, this input must be empty.
%   varargin{4:nargin}: Name/Value pairs of analysis sets, such as 
%       'Progress' => true or Gamma settings ('gammaAbs', 'gammaDta', 
%       'GammaRange'). See see below for a full list of options. 
%       Alternatively, a 'Report' option can be provided with the value set 
%       to a Delta4 report structure obtained from ParseD4Report. In this 
%       case, all relevant configuration options will be set to match what 
%       was set in the report. This is helpful when trying to replicate the 
%       same analysis as was performed in the Delta4 software. The 
%       name/value pair 'BeamAnalysis' => false can be provided to skip 
%       computation of individual beam statistics. Finally, a plan class
%       cell array can also be provided using 'planClasses', as described
%       above.
%
% Upon successfu completion, a structure is returned with the following
% fields:
%   gammaAbs: vector of absolute Gamma criteria evaluated, as a percent
%   gammaDta: vector of absolute Gamma criteria evaluated (together, with
%       gammaAbs, this produces a 2D table of all permutations), in mm
%   gammaRange: 2 element vector of the low and high percentage of the
%       reference dose to be included in Gamma statistics
%   gammaLimit: the number of DTAs to search when computing Gamma (for
%       example, is gammaLimit is 2 and gammaDta is 3, the code will search
%       in a 6 mm radius around each measurement point)
%   gammaRes: Gamma DTA search resolution (a value of 10 and gammaDta of 3
%       mm means that the code will evaluate every 0.3 mm).
%   refDose: a string indicating the type of reference dose. Can be
%       'measmax' (the reference dose will be the maximum measured dose) or
%       'planmax' (the reference dose will be the maximum reference dose)
%   absRange: 2 element vector of the low and high percentage of the
%       reference dose to be included in absolute dose difference
%   plans: table with one row for each plan matched to a dose volume
%   beams: table with one row for each beam matched to a dose volume
%   stats: table of the following parameter statistics (Gamma mean and pass 
%       rate, median dose difference) for the full dataset as well as each
%       subset: parameter, group, N (dataset size), min, median, mean,
%       trimmean (using 0.2 threshold), max, sd, lilliefors normality
%       p-value, skewness, and kurtosis. For this and other grouped
%       statistics, the "group" column contains a description of the
%       subset, where groups can be machine, beam energy, phantom
%       orientation, or plan class.
%   levene: table p-values for Levene's test comparing variance between 
%       groups, with each parameter as a column and each group category
%       (machine, energy, etc.) as a row.
%   kruskal: table p-values for Kruskal-Wallis rank test, with each 
%       parameter as a column and each group category (machine, energy, 
%       etc.) as a row.
%   dunnRanks: table of Dunn's multiple comparison test statistics based on 
%       Kruskal-Wallis groups. Columns include parameter, groupA, groupB, 
%       the standard error of each, median rank, rank difference, 
%       confidence intervals, and p-value.
%   anova: table of multi-group interaction ANOVA p-values, with each
%       parameter as a column and each interacting group as a row. In
%       addition to the groups used for basic stastical analysis, the
%       volume of the reference dose falling within the median absolute
%       range is included as a continuous variable.
%   dunnMeans: table of Dunn's multiple comparison test statistics based on 
%       the ANOVA groups. Columns include parameter, groupA, groupB, the 
%       standard error of each, mean, mean difference, confidence 
%       intervals, and p-value.
%   manova: table of p-values and MANOVA test statistics for each group.
%   pairwiseB: if a second reference dataset is provided, a table of two-
%       tailed Welch's t-test stastistics (p-value, diff, se, t-statistic)
%       between the two reference points.
%   
% The following example illustrates how to use this function:
%
% % Define folder containing Delta4 measurements
% meas = '../test_data/vmat_tests/';
% tps = '../test_data/tps_export/';
%
% % Execute evaluation using default criteria
% results = CompareD4Dose(meas, tps);
%
% % Execute evaluation, using folders above and 2%/2mm criteria
% results = CompareD4Dose(meas, tps, [], 'gammaAbs', 2, 'gammaDta', 2);
%
% % Extract a vector of pass rates for all plans using 'TrueBeam' and 6 MV
% pass = results.plans.gammaPassRateGlobal(strcmp(results.plans.machine, ...
%   'TrueBeam') & results.plans.energy == 6);
%
% % Find the mean local Gamma pass rate for all plans
% avg = results.stats.mean(strcmp(results.stats.parameter, ...
%       'gammaPassRateLocal') & strcmp(results.stats.group, 'all'));
%
% % Plot a scatter plot matric of the local pass rate, mean gamma, and 
% % median absolute difference grouped by machine
% gplotmatrix([results.plans.gammaPassRateLocal, ...
%   results.plans.gammaMeanLocal, results.plans.medianAbsDiff], [], 'bg', ...
%   '..', [],'on','', {'Pass Rate', 'Mean Gamma', 'Abs Diff'},  ...
%   {'Pass Rate', 'Mean Gamma', 'Abs Diff'});
%
% % Compare the pass rates between the original and a new TPS model
% pairwise = CompareD4Dose(meas, tps, '../test_data/new_model/');
%
% % Determine if the median absolute dose difference changed significantly
% % between the two models based on the t-test p-value
% pval = results.pairwiseB.p(strcmp(results.stats.parameter, ...
%       'medianAbsDiff') & strcmp(results.stats.group, 'all'));
%
% Author: Mark Geurts, mark.w.geurts@gmail.com
% Copyright (C) 2018 University of Wisconsin Board of Regents
%
% This program is free software: you can redistribute it and/or modify it 
% under the terms of the GNU General Public License as published by the  
% Free Software Foundation, either version 3 of the License, or (at your 
% option) any later version.
%
% This program is distributed in the hope that it will be useful, but 
% WITHOUT ANY WARRANTY; without even the implied warranty of 
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General 
% Public License for more details.
% 
% You should have received a copy of the GNU General Public License along 
% with this program. If not, see http://www.gnu.org/licenses/.

% Verify at least two inputs were provided
if nargin < 2
    error(['At least two inputs are required to execute CompareD4Dose: a ', ...
        'target measurement file/folder/list and a destination ', ...
        'file/folder/list']);
end
    
% Set default options
results.Progress = true;
results.alpha = 0.05;
results.absRange = [50 500];
results.gammaAbs = 3;
results.gammaDta = 3;
results.gammaRange = [20 500];
results.gammaLimit = 2;
results.gammaRes = 10;
results.BeamAnalysis = true;
results.refDose = 'measmax';
results.planClasses = cell(0,2);

% Update options
for i = 4:2:nargin
    results.(varargin{i}) = varargin{i+1};
end

% If 'Report' was set, copy parameters over
if isfield(results, 'Report')
    fields = fieldnames(results.Report);
    for i = 1:length(fields)
        if isfield(results, fields{i})
            results.(fields{i}) = results.Report.(fields{i});
        end
    end
    results = rmfield(results, 'Report');
    
    % Update refDose to 'planmax' as that is the default for ScandiDos
    results.refDose = 'planmax';
end

% Add dicom_tools submodule to search path
[path, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(path, 'dicom_tools'));

%% Scan reference dose folders for DICOM files
% Check if MATLAB can find ScanDICOMPath
if exist('ScanDICOMPath', 'file') ~= 2

    % If not, throw an error
    if exist('Event', 'file') == 2
        Event(['The dicom_tools submodule does not exist in the search path. ', ...
            'Use git clone --recursive or git submodule init followed by git ', ...
            'submodule update to fetch all submodules'], 'ERROR');
    else
        error(['The dicom_tools submodule does not exist in the search path. ', ...
            'Use git clone --recursive or git submodule init followed by git ', ...
            'submodule update to fetch all submodules']);
    end
end

% Display progress bar
if usejava('jvm') && feature('ShowFigureWindows') && results.Progress
    progress = waitbar(0, 'Scanning first reference dataset for RT Plans');
end

% Scan first reference folder
t = tic;
ref{1} = ScanDICOMPath(varargin{2}, 'Progress', false);

% Update waitbar
if exist('progress', 'var') && ishandle(progress) && results.Progress
    waitbar(0.1, progress, 'Scanning second reference dataset for RT Plans');
end

% Scan second reference folder
if nargin >= 3 && ~isempty(varargin{3})
    ref{2} = ScanDICOMPath(varargin{2}, 'Progress', false);
end

% Update waitbar
if exist('progress', 'var') && ishandle(progress) && results.Progress
    waitbar(0.2, progress, 'Scanning measurement files');
end

%% Scan directory of measurement files
% Log start of measurement folders
if exist('Event', 'file') == 2
    Event('Scanning for measurement files');
end

% Scan measurements folder directory contents
if iscell(varargin{1})
    meas = varargin{1};
elseif isfolder(varargin{1})
    meas = dir(fullfile(varargin{1}, '**'));
else
    meas = dir(varargin{1});
end

%% Initialize data structures
% Initialize plan list, beam list, and gamma arrays return fields
fixedNames = {'measurementFile', 'patient', 'ID', 'plan', 'machine', ...
    'energy', 'planClass', 'orientation'};
fixedUnits = {'', '', '', '', '', 'MV', '', ''};
varNames = {'referenceFile', 'referenceDose', 'gammaPassRateGlobal', ...
    'gammaMeanGlobal', 'gammaMaxGlobal', 'gammaPassRateLocal', ...
    'gammaMeanLocal', 'gammaMaxLocal', 'medianAbsDiff', 'absVolume', ...
    'gammaVolume'};
varUnits = {'', 'Gy', '%', '', '', '%', '', '', '%', 'cc', 'cc'};

% If only one reference exists
if length(ref) == 1
    v = varNames;
    u = varUnits;

% If multiple references exist
else
    
    % Initialize variable names
    v = cell(1,length(varNames)*length(ref));
    u = cell(1,length(varNames)*length(ref));
    for i = 1:length(ref)
        for j = 1:length(varNames)
            
            % Append a letter (A, B, C, ...) onto each variable name
            v{(i-1)*length(varNames)+j} = [varNames{j}, char(i+64)];
            u{(i-1)*length(varNames)+j} = varUnits{j};
        end
    end
end

% Create tables
results.plans = array2table(zeros(0,length(fixedNames)+length(v)), ...
    'VariableNames', horzcat(fixedNames, v));
results.plans.Properties.VariableUnits = horzcat(fixedUnits, u);
if results.BeamAnalysis
    results.beams = array2table(zeros(0,1+length(fixedNames)+length(v)), ...
        'VariableNames', horzcat(fixedNames, {'Beam'}, v));
    results.beams.Properties.VariableUnits = horzcat(fixedUnits, {''}, u);
end

% Clear and initialize GPU memory.  If CUDA is not enabled, or if the
% Parallel Computing Toolbox is not installed, this will error, and the
% function will automatically rever to CPU computation via the catch
% statement
try
    gpuDevice(1);
catch
    if exist('Event', 'file') == 2
        Event('GPU failed, will perform CPU interpolation', 'WARN');
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%% Process measurements %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop through measurement list
for i = 1:length(meas)

    %% Parse measurement
    % Update waitbar
    if exist('progress', 'var') && ishandle(progress) && results.Progress
        waitbar(0.2 + 0.7*i/length(meas), progress, ...
            sprintf('Comparing Delta4 measurement %i of %i', i, length(meas)));
    end
    
    % If the folder content is . or .. or a folder, skip to next file
    if strcmp(meas(i).name, '.') || strcmp(meas(i).name, '..') ...
            || meas(i).isdir == 1
        continue;
    end
    
    % Parse file extension
    [path, name, ext] = fileparts(fullfile(meas(i).folder, meas(i).name));
    
    % If extension matches potential report format, try to parse it
    if ismember(ext, {'.xlsx', '.xls', '.txt', '.csv'})
        try
            if exist('Event', 'file') == 2
                Event(['Parsing Delta4 measurement file ', name, ext]);
            end
            delta4 = ParseD4Tables(fullfile(path, [name, ext]), ...
                'CenterPhantom', true);
        
        % If parser fails, catch gracefully
        catch
            if exist('Event', 'file') == 2
                Event('Parsing failed, skipping measurement', 'WARN');
            else
                warning(['Parsing failed, skipping ', name, ext]);
            end
            continue;
        end
    else
        continue;
    end
    
    % Verify the plan dose exists
    if ~isfield(delta4, 'data') || isempty(delta4.data)
        continue;
    end

    % Split name
    nd = strsplit(delta4.name, ',');
    
    % Match plan class (will stop after the first match is found)
    class = 'Undefined';
    for j = 1:size(results.planClasses, 1)
        if ~isempty(regexpi(delta4.plan, results.planClasses{j,2}))
            class = results.planClasses{j,2};
            break;
        end
    end
    
    % Set orientation (if there are more than 400 diodes with same IEC Z 
    % value, this is a sagittal/coronal measurement)
    if sum(delta4.data(:,2) == mode(delta4.data(:,2))) > 400
        orient = 'Sag/Cor';
    else
        orient = 'Diagonal';
    end
    
    % Initialize plan IDs
    matched = zeros(1, length(ref));
    
    %% Match plan to reference dose files
    % Loop through each reference cell array
    for r = 1:length(ref)
        for j = 1:length(ref{r})
        
            % Split name
            nr = strsplit(ref{r}{j,5}, ',');
            
            % If reference is a plan and matches the ID, plan, and pt name
            if strcmp(ref{r}{j,10}, 'PLAN') && strcmpi(delta4.ID, ...
                    ref{r}{j,6}) && strcmpi(delta4.plan, ref{r}{j,9}) && ...
                    contains(nd{1}, nr{1}, 'IgnoreCase', true) && ...
                    (length(nd) < 2 || length(nr) < 2 || ...
                    contains(nd{2}, nr{2}, 'IgnoreCase', true)) && ...
                    contains(delta4.value, 'Absolute dose', 'IgnoreCase', ...
                    true)
                
                matched(r) = j;
                break;
            end
        end
    end
    
    % If all references were matched
    if all(matched)

        % Initialize results
        result = {[name, ext], delta4.name, delta4.ID, delta4.plan, ...
            ref{r}{matched(1),12}{1}, ...
            mean(cell2mat(ref{r}{matched(1),13})), class, orient};

        % Loop through each reference
        for j = 1:length(matched)
            
            % Load dose
            dose = LoadDICOMDose(ref{r}{matched(j),2}, ref{r}{matched(j),1});
        
            % Set reference value
            switch results.refDose
                case 'measmax'
                    refval = max(delta4.data(:,5));
                case 'planmax'
                    refval = max(max(max(dose.data)));
                otherwise
                    refval = 0;
            end
            
            % Compute dose volumes
            absvol = (sum(sum(sum(dose.data >= refval * ...
                results.absRange(1)/100))) - sum(sum(sum(dose.data > ...
                refval * results.absRange(2)/100)))) * prod(dose.width);
            gamvol = (sum(sum(sum(dose.data >= refval * ...
                results.gammaRange(1)/100))) - sum(sum(sum(dose.data > ...
                refval * results.gammaRange(2)/100)))) * prod(dose.width);
            
            % Log calculation
            if exist('Event', 'file') == 2
                Event(['Computing plan Gamma table between ', name, ext, ...
                    'and ', ref{r}{matched(j),1}]);
            end
            
            % Append reference file and gamma table
            result = [result, ref{r}{matched(j),1}, refval, ...
                gammaTable(delta4.data, dose, refval, results), absvol, ...
                gamvol]; %#ok<*AGROW>
        end
        
        % Append result onto results
        results.plans = [results.plans; result];
        
    % Otherwise, warn the user
    else
        if exist('Event', 'file') == 2
            Event(['No matching dose files were found for measurement ', ...
                name, ext], 'WARN');
        else
            warning(['No matching dose files were found for measurement ', ...
                name, ext])
        end
    end

    %% Match each beam to reference dose files
    % Next, if beams exist
    if results.BeamAnalysis && isfield(delta4, 'beams')
        
        % Loop through measured beam data
        for b = 1:length(delta4.beams)
            
            % Initialize beam IDs
            matched = zeros(1, length(ref));
            
            % If parse was successful, match plan to reference dose files
            for r = 1:length(ref)
                for j = 1:length(ref{r})
        
                    % Split name
                    nr = strsplit(ref{r}{j,5}, ',');

                    % If reference is a BEAM and matches the ID, plan, pt , 
                    % and beam name
                    if strcmp(ref{r}{j,10}, 'BEAM') && strcmpi(delta4.ID, ...
                            ref{r}{j,6}) && strcmpi(delta4.plan, ref{r}{j,9}) && ...
                            contains(nd{1}, nr{1}, 'IgnoreCase', true) && ...
                            (length(nd) < 2 || length(nr) < 2 || ...
                            contains(nd{2}, nr{2}, 'IgnoreCase', true)) && ...
                            strcmpi(delta4.beams{b}.name, ref{1}{j,11}) && ...
                            contains(delta4.beams{b}.value, ...
                            'Absolute dose', 'IgnoreCase', true)
                        
                        matched(r) = j;
                        break;
                    end
                end
            end
            
            % If all references were matched
            if all(matched)
                
                % Initialize results
                result = {[name, ext], delta4.name, delta4.ID, ...
                    delta4.plan, ref{r}{matched(1),12}, ...
                    ref{r}{matched(1),13}, class, orient, ...
                    delta4.beams{b}.name};

                % Loop through each reference
                for j = 1:length(matched)

                    % Load dose
                    dose = LoadDICOMDose(ref{r}{matched(j),2}, ...
                        ref{r}{matched(j),1});

                    % Set reference value
                    switch results.refDose
                        case 'measmax'
                            refval = max(delta4.beams{b}.data(:,5));
                        case 'planmax'
                            refval = max(max(max(dose.data)));
                        otherwise
                            refval = 0;
                    end
                    
                    % Compute dose volumes
                    absvol = (sum(sum(sum(dose.data >= refval * ...
                        results.absRange(1)/100))) - sum(sum(sum(dose.data > ...
                        refval * results.absRange(2)/100)))) * prod(dose.width);
                    gamvol = (sum(sum(sum(dose.data >= refval * ...
                        results.gammaRange(1)/100))) - sum(sum(sum(dose.data > ...
                        refval * results.gammaRange(2)/100)))) * prod(dose.width);

                    % Log calculation
                    if exist('Event', 'file') == 2
                        Event(['Computing beam Gamma table between ', ...
                            name, ext, 'and ', ref{r}{matched(j),1}]);
                    end
                    
                    % Append reference file and gamma table
                    result = [result, ref{r}{matched(j),1}, refval, ...
                        gammaTable(delta4.beams{b}.data, ...
                        dose, refval, results), absvol, gamvol]; %#ok<*AGROW>
                end

                % Append result onto results
                results.beams = [results.beams; result];
            end
        end
    end
end

% Clear temporary variables
clear b i j r t u v nd nr class delta4 dose path name ext fixedNames fixedUnits ... 
    matched meas refval result varNames varUnits orient;

%% %%%%%%%%%%%%%%%%%%%%%%% Compute statistics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if size(results.plans,1) < 1
    if exist('progress', 'var') == 1 && ishandle(progress)
        close(progress);
    end
    if exist('Event', 'file') == 2
        Event('No plans were found with matching dose data', 'WARN');
    else
        warning('No plans were found with matching dose data');
    end
    return;
end

if exist('progress', 'var') && ishandle(progress) && results.Progress
    waitbar(0.91, progress, 'Computing statistics');
end

% Define parameters to be statistically tested
params = {'gammaPassRateGlobal', 'gammaMeanGlobal', 'gammaPassRateLocal', ...
    'gammaMeanLocal', 'medianAbsDiff'};

% If there are more than one, append 'A' to the above parameters
% (individual statistics are only reported for the first reference dataset)
if length(ref) > 1
    for i = 1:length(params)
        params{i} = [params{i}, char(1 + 64)];
    end
end

% Define groups to be tested, and unique values in each
groups = {'machine', 'energy', 'planClass', 'orientation'};
indices = cell(1,4);
for i = 1:length(groups)
    
    % For energy, use beam energies (plans contain average energy)
    if strcmp(groups{i}, 'energy')
        if isfield(results, 'beams') && ~isempty(results.beams)
            groupVals{i} = unique(results.beams.(groups{i}));
        else
            groupVals{i} = unique(results.plans.(groups{i}));
        end
    else
        groupVals{i} = unique(results.plans.(groups{i}));
    end
    indices{i} = 1:length(groupVals{i});
end

% Remove groups with only one unique value
groups = groups(cellfun(@length, groupVals) > 1);
groupVals = groupVals(cellfun(@length, groupVals) > 1);

% Compute all combinations of groups
groupDims = cell(0);
for i = 1:length(groups)
    a = nchoosek(1:length(groups),i);
    for j = 1:size(a,1)
        groupDims{length(groupDims)+1} = a(j,:);
    end
end

% Compute permutations of each combination
e = cell(length(groupDims),length(groups));
for i = 1:length(groupDims)
    for j = 1:length(groupDims{i})
        e{i,groupDims{i}(j)} = 1:length(groupVals{groupDims{i}(j)});
    end
end
e(cellfun(@isempty, e)) = {0};
groupCombs = zeros(1,length(groups));
for i = 1:size(e,1)
    c = cell(1, length(e(i,:)));
    [c{:}] = ndgrid(e{i,:});
    c = cellfun(@(x) x(:), c, 'uniformoutput', false);
    groupCombs = [groupCombs; [c{:}]];
end

%% Compute group statistics tables
% Initialize array
results.stats = array2table(zeros(0,12), ...
    'VariableNames', {'parameter', 'group', 'N', 'min', 'median', 'mean', ...
    'trimmean', 'max', 'sd', 'lilliefors', 'skewness', 'kurtosis'});

% Temporarily suppress warnings (mainly for the Lilliefors test)
w = warning('off','all');

% Loop through each permutation
for i = 1:size(groupCombs,1)
    
    % Initialize index vector and group name
    match = true(size(results.plans,1),1);
    g = cell(0);
    
    % Loop through each group
    for j = 1:size(groupCombs,2)
        
        % If non-zero, match only provided group
        if groupCombs(i,j) > 0
            if isnumeric(groupVals{j}(groupCombs(i,j)))
                match = match & (results.plans.(groups{j}) == ...
                    groupVals{j}(groupCombs(i,j)));
                g{length(g)+1} = [groups{j}, '=', ...
                    num2str(groupVals{j}(groupCombs(i,j)))];
            else
                match = match & strcmp(results.plans.(groups{j}), ...
                    groupVals{j}(groupCombs(i,j)));
                g{length(g)+1} = [groups{j}, '=', ...
                    groupVals{j}{groupCombs(i,j)}];
            end
        end
    end
    
    % If everything is matched, use 'all' for group name
    if isempty(g)
        g{1} = 'all';
    end
    
    % Loop through each parameter
    for j = 1:length(params)

        % If this parameter is a cell array
        if iscell(results.plans.(params{j})(1))
            x = cell2mat(permute(results.plans.(params{j})(match),[3 2 1]));
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
            results.stats = [results.stats; {params{j}, strjoin(g, ','), ...
                length(x), min(x,[],3), median(x,3), mean(x,3), ...
                trimmean(x,0.2,3), max(x,[],3), std(x,0,3), p, s, k}];
        else
            x = results.plans.(params{j})(match);
            try
                [~, p] = lillietest(x);
            catch
                p = NaN;
            end
            s = skewness(x,0);
            k = kurtosis(x,0);
            results.stats = [results.stats; {params{j}, strjoin(g, ','), ...
                length(x), min(x), median(x), mean(x), trimmean(x,0.2), ...
                max(x), std(x), p, s, k}];
        end
    end
end

% Restore original warning state
warning(w);

% Move progress bar forward
if exist('progress', 'var') && ishandle(progress) && results.Progress
    waitbar(0.93, progress);
end

%% Compute Levene's test for equal variance
results.levene = array2table(zeros(length(groups), length(params)), ...
    'VariableNames', params, 'RowNames', groups);
for i = 1:length(groups)
    if exist('Event', 'file') == 2
        Event(['Computing Kruskal-Wallis tests for group ', groups{i}]);
    end
    for j = 1:length(params)
        try
            if iscell(results.plans.(params{j})(1))
                x = cell2mat(permute(results.plans.(params{j})(...
                    ismember(results.plans.(groups{i}), groupVals{i})),...
                    [3 2 1]));
                results.levene{i,j} = vartestn(squeeze(x(1,1,:)), ...
                    results.plans.(groups{i})(ismember(results.plans.(...
                    groups{i}), groupVals{i})), 'Display', 'off', ...
                    'TestType', 'LeveneQuadratic');
            else
                results.levene{i,j} = vartestn(results.plans.(params{j})(...
                    ismember(results.plans.(groups{i}), groupVals{i})), ...
                    results.plans.(groups{i})(ismember(results.plans.(...
                    groups{i}), groupVals{i})), 'Display', 'off', ...
                    'TestType', 'LeveneQuadratic');
            end
        catch
            if exist('Event', 'file') == 2
                Event(['Error computing Levene statistic for ', ...
                   'parameter ', params{i}, ', group ', groups{i}], 'WARN');
            else
                warning(['Error computing Levene statistic for ', ...
                   'parameter ', params{i}, ', group ', groups{i}]);
            end
        end
    end
end

if exist('progress', 'var') && ishandle(progress) && results.Progress
    waitbar(0.95, progress);
end

%% Compute plan Kruskal Wallis tests for each parameter, group
results.kruskal = array2table(zeros(length(groups), length(params)), ...
    'VariableNames', params, 'RowNames', groups);
results.dunnRanks = array2table(zeros(0,11), 'VariableNames', {'parameter', ...
    'groupA', 'groupB', 'rankA', 'seA', 'rankB', 'seB', ...
    'lowerCI', 'diff', 'upperCI', 'p'});
for i = 1:length(groups)
    if exist('Event', 'file') == 2
        Event(['Computing Kruskal-Wallis tests for group ', groups{i}]);
    end
    for j = 1:length(params)
        try
            if iscell(results.plans.(params{j})(1))
                x = cell2mat(permute(results.plans.(params{j})(...
                    ismember(results.plans.(groups{i}), groupVals{i})),...
                    [3 2 1]));
                [p,~,stats] = kruskalwallis(squeeze(x(1,1,:)), ...
                    results.plans.(groups{i})(ismember(results.plans.(...
                    groups{i}), groupVals{i})), 'off');
            else
                [p,~,stats] = kruskalwallis(results.plans.(params{j})(...
                    ismember(results.plans.(groups{i}), groupVals{i})), ...
                    results.plans.(groups{i})(ismember(results.plans.(...
                    groups{i}), groupVals{i})), 'off');
            end
            results.kruskal{i,j} = p;
            [d, m] = multcompare(stats, 'Alpha', results.alpha, 'CType', ...
                'dunn-sidak', 'Display', 'off');
            for k = 1:size(d,1)
                results.dunnRanks = [results.dunnRanks; [params(j), ...
                    sprintf('%s=%s', groups{i}, stats.gnames{d(k,1)}), ...
                    sprintf('%s=%s', groups{i}, stats.gnames{d(k,2)}), ...
                    m(d(k,1),1), m(d(k,1),2), ...
                    m(d(k,2),1), m(d(k,2),2), d(k,3), d(k,4), d(k,5), ...
                    d(k,6)]];
            end
        catch
            if exist('Event', 'file') == 2
                Event(['Error computing Kruskal-Wallis statistic for ', ...
                   'parameter ', params{j}, ', group ', groups{i}], 'WARN');
            else
                warning(['Error computing Kruskal-Wallis statistic for ', ...
                   'parameter ', params{j}, ', group ', groups{i}]);
            end
        end
    end
end

if exist('progress', 'var') && ishandle(progress) && results.Progress
    waitbar(0.97, progress);
end

%% Compute plan level multi-group ANOVA for each parameter
x = nan(size(results.plans,1), length(params));
for i = 1:size(results.plans,1)
    ingroup = true;
    for j = 1:length(groups)
        if ~ismember(results.plans.(groups{j}), groupVals{j})
            ingroup = false;
            break;
        end
    end
    if ingroup    
        for j = 1:length(params)
            if iscell(results.plans.(params{j})(i))
                x(i,j) = results.plans.(params{j}){i}(1,1);
            else
                x(i,j) = results.plans.(params{j})(i);
            end
        end
    end
end
g = cell(1, length(groups)+1);
for i = 1:length(groups) 
    g{i} = results.plans.(groups{i})(any(~isnan(x),2));
end
if length(ref) == 1
    g{length(groups)+1} = results.plans.absVolume(any(~isnan(x),2));
else
    g{length(groups)+1} = ...
        results.plans.(['absVolume', char(65)])(any(~isnan(x),2));
end
terms = [groups 'absVolume'];
x = x(any(~isnan(x),2),:);
results.anova = array2table(zeros(nchoosek(length(groups)+1,1) + ...
    nchoosek(length(groups)+1,2), length(params)), 'VariableNames', params);
results.dunnMeans = array2table(zeros(0,11), 'VariableNames', {'parameter', ...
    'groupA', 'groupB', 'meanA', 'seA', 'meanB', 'seB', ...
    'lowerCI', 'diff', 'upperCI', 'p'});
for i = 1:length(params)
    if exist('Event', 'file') == 2
        Event(['Computing multi-group ANOVA for parameter ', params{i}]);
    end
    try
        [p,~,stats] = anovan(x(:,i), g, 'display', 'off', ...
            'model', 'interaction', 'continuous', length(groups)+1);
        results.anova{:,i} = p;
        for j = 1:size(stats.terms,1)
            results.anova.Properties.RowNames{j} = ...
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
            [d, m, ~, n] = multcompare(stats, 'Alpha', results.alpha, 'CType', ...
                'dunn-sidak', 'Display', 'off', 'Dimension', dims{j});
            for k = 1:length(groups)
                n = strrep(n, sprintf('X%i',k), groups{k});
            end
            for k = 1:size(d,1)
                results.dunnMeans = [results.dunnMeans; [params(i), n(d(k,1)), ...
                    n(d(k,2)), m(d(k,1),1), m(d(k,1),2), m(d(k,2),1), ...
                    m(d(k,2),1), d(k,3), d(k,4), d(k,5), d(k,6)]];
            end
        end
    catch
        if exist('Event', 'file') == 2
            Event(['Error computing ANOVA for parameter ', params{i}], ...
                'WARN');
        else
            warning(['Error computing ANOVA for parameter ', params{i}]);
        end
    end
end

if exist('progress', 'var') && ishandle(progress) && results.Progress
    waitbar(0.99, progress);
end

%% Compute plan level MANOVA for each group
results.manova = cell(length(groups), 17);
for i = 1:length(groups)
    if exist('Event', 'file') == 2
        Event(['Computing MANOVA for group ', groups{i}]);
    end
    x = zeros(length(results.plans.measurementFile(...
        ismember(results.plans.(groups{i}), groupVals{i}))), length(params));
    for j = 1:length(params)
        if iscell(results.plans.(params{j})(1))
            y = cell2mat(permute(results.plans.(params{j})(...
                ismember(results.plans.(groups{i}), groupVals{i})), [3 2 1]));
            x(:,j) = squeeze(y(1,1,:));
        else
            x(:,j) = results.plans.(params{j})(...
                ismember(results.plans.(groups{i}), groupVals{i}));
        end
    end
    try
        [d,p,stats] = manova1(x, results.plans.(groups{i})(...
            ismember(results.plans.(groups{i}), groupVals{i})), ...
            results.alpha);
        results.manova{i,1} = d;
        results.manova{i,2} = p;
        results.manova{i,3} = stats.W;
        results.manova{i,4} = stats.B;
        results.manova{i,5} = stats.T;
        results.manova{i,6} = stats.dfW;
        results.manova{i,7} = stats.dfB;
        results.manova{i,8} = stats.dfT;
        results.manova{i,9} = stats.lambda;
        results.manova{i,10} = stats.chisq;
        results.manova{i,11} = stats.chisqdf;
        results.manova{i,12} = stats.eigenval;
        results.manova{i,13} = stats.eigenvec;
        results.manova{i,14} = stats.canon;
        results.manova{i,15} = stats.mdist;
        results.manova{i,16} = stats.gmdist;
        results.manova{i,17} = stats.gnames;
    catch
        if exist('Event', 'file') == 2
            Event(['Error computing MANOVA for group ', groups{i}], 'WARN');
        else
            warning(['Error computing MANOVA for group ', groups{i}]);
        end
    end
end
results.manova = cell2table(results.manova, 'RowNames', groups, ...
    'VariableNames', {'dimension', 'p', 'W', 'B', 'T', 'dfW', 'dfB', 'dfT', ...
    'lambda', 'chisq', 'chisqdf', 'eigenval', 'eigenvec', 'canon', 'mdist', ...
    'gmdist', 'gnames'});

%% Compute pairwise t-tests
for r = 2:length(ref)
    
    % Initialize results table
    results.(['pairwise', char(r+64)]) = array2table(zeros(0,9), ...
        'VariableNames', {'parameter', 'group', 'N', 'se', ...
        'lowerCI', 'diff', 'upperCI', 'tstat', 'p'});
    
    % Loop through each permutation
    for i = 1:size(groupCombs,1)

        % Initialize index vector and group name
        match = true(size(results.plans,1),1);
        g = cell(0);

        % Loop through each group
        for j = 1:size(groupCombs,2)

            % If non-zero, match only provided group
            if groupCombs(i,j) > 0
                if isnumeric(groupVals{j}(groupCombs(i,j)))
                    match = match & (results.plans.(groups{j}) == ...
                        groupVals{j}(groupCombs(i,j)));
                    g{length(g)+1} = [groups{j}, '=', ...
                        num2str(groupVals{j}(groupCombs(i,j)))];
                else
                    match = match & strcmp(results.plans.(groups{j}), ...
                        groupVals{j}(groupCombs(i,j)));
                    g{length(g)+1} = [groups{j}, '=', ...
                        groupVals{j}{groupCombs(i,j)}];
                end
            end
        end

        % If everything is matched, use 'all' for group name
        if isempty(g)
            g{1} = 'all';
        end

        % Loop through each parameter
        for j = 1:length(params)

            % If this parameter is a cell array
            if iscell(results.plans.(params{j})(1))
                
                % Store matched refA and refB parameters
                x = cell2mat(permute(results.plans.(params{j})(match),[3 2 1]));
                y = cell2mat(permute(results.plans...
                    .([params{j}(1:end-1), char(r+64)])(match),[3 2 1]));
                
                % Initialize arrays to store parameters 
                s = nan(size(x,1),size(x,2),6);
                for k = 1:size(x,1)
                    for l = 1:size(x,2)
                        try
                            [~, p, ci, stats] = ttest(x(k,l,:), y(k,l,:), ...
                                'Alpha', results.alpha, 'Dim', 3, 'Tail', ...
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
                    {params{j}, strjoin(g, ','), length(x), s(:,:,1), ...
                    s(:,:,2), s(:,:,3), s(:,:,4), s(:,:,5), s(:,:,6)}];
            else
                x = results.plans.(params{j})(match);
                y = results.plans.([params{j}(1:end-1), char(r+64)])(match);
                try
                    [~, p, ci, stats] = ttest(x, y, 'Alpha', results.alpha, ...
                        'Tail', 'both');
                catch
                    ci = [NaN NaN];
                    p = NaN;
                    stats.tstat = NaN;
                    stats.sd = NaN;
                end
                results.(['pairwise', char(r+64)]) = ...
                    [results.(['pairwise', char(r+64)]); ...
                    {params{j}, strjoin(g, ','), length(x), ...
                    stats.sd/sqrt(length(x)), ci(1), mean(x-y), ci(2), ...
                    stats.tstat, p}];
            end
        end
    end
end

%% Finish up
% Close waitbar
if exist('progress', 'var') == 1 && ishandle(progress)
    waitbar(1.0, progress);
    close(progress);
end

% Remove Progress and BeamAnalysis field from results
results = rmfield(results, 'Progress');
results = rmfield(results, 'BeamAnalysis');

% Log completion
if exist('Event', 'file') == 2
    Event(sprintf(['Comparison successfully completed in %0.3f seconds, ', ...
        'scanning %i plans and %i beams'], toc(t), size(results.plans,1), ...
        size(results.beams,1)));
end

% Clear temporary variables
clear i j k l x y t d p g params ingroup stats ref groups groupVals progress;

end

%% Compute Gamma Table Subfunction
function stats = gammaTable(meas, dose, refval, criteria)
% gammaTable is called by CompareD4Dose and computes a 2D gamma table based
% on a provided measurement array and reference dose volume, using the
% results gammaAbs, gammaDta, gammaRange, gammaRes, and gammaLimit fields
% from the .

% Compute mesh grid for reference dose
[meshA, meshB, meshC] = meshgrid(single(dose.start(2) + ...
    dose.width(2) * (size(dose.data,2) - 1):-dose.width(2):dose.start(2)), ...
    single(dose.start(1):dose.width(1):dose.start(1) + dose.width(1)...
    * (size(dose.data,1) - 1)), single(dose.start(3):dose.width(3):...
    dose.start(3) + dose.width(3) * (size(dose.data,3) - 1)));

% Interpolate dose at original position
try
    % Run GPU interp3 function to compute the dose values at the shifted 
    % measured coordinate points
    i = gather(interp3(gpuArray(meshA), gpuArray(meshB), ...
        gpuArray(meshC), gpuArray(single(dose.data)), ...
        gpuArray(meas(:,4)), gpuArray(meas(:,2)), gpuArray(-meas(:,3)), ...
        'linear', 0));
catch

    % Run CPU interp3 function to compute the dose values at the shifted 
    % measured coordinate points
    i = interp3(meshA, meshB, meshC, single(dose.data), meas(:,4), meas(:,2), ...
        -meas(:,3), '*linear', 0);
end

% Compute the median absolute dose difference
absDiff = (meas(:,5) - i)./refval * 100;
stats = cell(1,7);
stats{7} = median(absDiff(meas(:,5) > refval * criteria.absRange(1)/100 & ...
    meas(:,5) < refval * criteria.absRange(2)/100));

% Apply gamma range to measured values
meas = meas(meas(:,5) > refval * criteria.gammaRange(1)/100 & ...
    meas(:,5) < refval * criteria.gammaRange(2)/100,:);

% Compute shifts (in cm)
shifts = repmat(-criteria.gammaRes * criteria.gammaLimit:criteria.gammaRes * ...
    criteria.gammaLimit, size(meas,1), 1) / ...
    criteria.gammaRes * max(criteria.gammaDta)/10;

% Compute shifted 2D arrays for each position
measA = single(horzcat(repmat(meas(:,4), 1, size(shifts,2)) + shifts, ...
    repmat(meas(:,4), 1, 2*size(shifts,2))));
measB = single(horzcat(repmat(meas(:,2), 1, size(shifts,2)), ...
    repmat(meas(:,2),1, size(shifts,2)) + shifts, ...
    repmat(meas(:,2), 1, size(shifts,2))));
measC = single(horzcat(repmat(-meas(:,3), 1, 2*size(shifts,2)), ...
    repmat(-meas(:,3), 1, size(shifts,2)) + shifts));

% Interpolate dose at shifted position
try
    % Run GPU interp3 function to compute the dose values at the shifted 
    % measured coordinate points
    i = gather(interp3(gpuArray(meshA), gpuArray(meshB), ...
        gpuArray(meshC), gpuArray(single(dose.data)), ...
        gpuArray(measA), gpuArray(measB), gpuArray(measC), 'linear', 0));
catch

    % Run CPU interp3 function to compute the dose values at the shifted 
    % measured coordinate points
    i = interp3(meshA, meshB, meshC, single(dose.data), measA, measB, ...
        measC, '*linear', 0);
end

% Compute local gamma as min of each shifted position
loc = zeros(size(i,1), length(criteria.gammaAbs), length(criteria.gammaDta));
glob = zeros(size(i,1), length(criteria.gammaAbs), length(criteria.gammaDta));

% Loop through Abs, DTA criteria
for a = 1:length(criteria.gammaAbs)
    for d = 1:length(criteria.gammaDta)
        
        % Compute global gamma based on refval as min of each shifted position
        glob(:,a,d) = sqrt(min(((repmat(meas(:,5),1,size(i,2)) - i) ./ ...
            (refval * criteria.gammaAbs(a)/100)).^2 + ...
            (repmat(shifts,1,3)/(criteria.gammaDta(d)/10)).^2,[],2));
        
        % Compute local gamma as min of each shifted position
        loc(:,a,d) = sqrt(min(((repmat(meas(:,5),1,size(i,2)) - i) ./ ...
            (i * criteria.gammaAbs(a)/100)).^2 + ...
            (repmat(shifts,1,3)/(criteria.gammaDta(d)/10)).^2,[],2));
    end
end

% Compute pass rate, mean, and maximum Gamma statistics for global gamma
stats{1} = zeros(length(criteria.gammaAbs), length(criteria.gammaDta));
for a = 1:length(criteria.gammaAbs)
    for d = 1:length(criteria.gammaDta)
        stats{1}(a,d) = sum(glob(:,a,d) <= 1) / size(glob,1) * 100;
    end
end
stats{2} = squeeze(mean(glob,1));
stats{3} = squeeze(max(glob,[], 1));

% Compute pass rate, mean, and maximum Gamma statistics for local gamma
stats{4} = zeros(length(criteria.gammaAbs), length(criteria.gammaDta));
for a = 1:length(criteria.gammaAbs)
    for d = 1:length(criteria.gammaDta)
        stats{4}(a,d) = sum(loc(:,a,d) <= 1) / size(loc,1) * 100;
    end
end
stats{5} = squeeze(mean(loc,1));
stats{6} = squeeze(max(loc,[], 1));

% Clear temporary variables
clear measA measB measC a d i shifts glob loc meshA meshB meshC absDiff;

end
