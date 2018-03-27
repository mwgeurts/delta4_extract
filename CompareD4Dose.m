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
% ParseD4Tables from this repository. It will attempt to match the Patient
% last name, first name, ID, and Plan name between each Delta4 measurement
% and dose volume. If a match is found, the results are stored in the
% return structure field results.Plans. The matching is then repeated for
% each beam, using the beam name and storing the results to results.Beams.
%
% As part of the analysis, this function accepts a Name/Value configuration
% option 'PlanClasses' with an n x 2 cell array of Plan class names and
% associated case insensitive regular expressions to match to the Plan 
% name. If provided, the function will attempt to categorize each 
% measurement into a Plan class, then compare each Plan type. Note, each 
% Plan is only matched to the first class that matches the regexp.
%
% The following required and optional inputs are available:
%   varargin{1}: file/folder/list of files containing measurement files
%       (see ParseD4tables for information on compatible formats)
%   varargin{2}: file/folder/list of files containing reference DICOM RT
%       Plan and RT Dose files for each measurement file in varargin{1}
%   varargin{3} (optional): file/folder/list of files containing reference 
%       DICOM RT Plan and RT Dose files for each measurement for pairwise
%       comparison. If additional Name/Value parameters are passed to this
%       function, this input must be empty.
%   varargin{4:nargin}: Name/Value pairs of analysis sets, such as 
%       'Progress' => true or Gamma settings ('GammaAbs', 'GammaDTA', 
%       'GammaRange'). See see below for a full list of options. 
%       Alternatively, a 'Report' option can be provided with the value set 
%       to a Delta4 report structure obtained from ParseD4Report. In this 
%       case, all relevant configuration options will be set to match what 
%       was set in the report. This is helpful when trying to replicate the 
%       same analysis as was performed in the Delta4 software. The 
%       name/value pair 'BeamAnalysis' => false can be provided to skip 
%       computation of individual beam statistics. Finally, a Plan class
%       cell array can also be provided using 'PlanClasses', as described
%       above.
%
% Upon successful completion, a structure is returned with these fields:
%   Alpha: the significance level used for hypothesis testing and
%       confidence intervals
%   GammaAbs: vector of absolute Gamma criteria evaluated, as a percent
%   GammaDTA: vector of absolute Gamma criteria evaluated (together, with
%       GammaAbs, this produces a 2D table of all permutations), in mm
%   GammaRange: 2 element vector of the low and high percentage of the
%       reference dose to be included in Gamma statistics
%   GammaLimit: the number of DTAs to search when computing Gamma (for
%       example, is GammaLimit is 2 and GammaDTA is 3, the code will search
%       in a 6 mm radius around each measurement point)
%   GammaRes: Gamma DTA search resolution (a value of 10 and GammaDTA of 3
%       mm means that the code will evaluate every 0.3 mm).
%   RefDose: a string indicating the type of reference dose. Can be
%       'measmax' (the reference dose will be the maximum measured dose) or
%       'Planmax' (the reference dose will be the maximum reference dose)
%   AbsRange: 2 element vector of the low and high percentage of the
%       reference dose to be included in absolute dose difference
%   Plans: table with one row for each Plan matched to a dose volume
%   Beams: table with one row for each beam matched to a dose volume
%   Stats: MATLAB table of basic statistics (min, max, mean, etc.), along
%       with p-values from a Lilliefors normality test and a few robust
%       statistics (trimmed mean, Huber M location estimate, MAD) for each
%       variable being analyzed (Gamma Pass Rate, Median Absolute Dose
%       Difference, etc.). Results are also computed for each dataset group
%       (Machine, Energy, Plan Class, etc.) as well as each subgroup. The
%       Group table column describes the grouping; examples include 'all'
%       (all data) or 'Machine=TrueBeam,Energy=6FFF'.
%   Levene: table p-values for Levene's test comparing variance between 
%       groups, with each parameter as a column and each group category
%       (Machine, Energy, etc.) as a row.
%   Kruskal: table p-values for Kruskal-Wallis rank test, with each 
%       parameter as a column and each group category (Machine, Energy, 
%       etc.) as a row.
%   DunnRanks: table of Dunn's multiple comparison test statistics based on 
%       Kruskal-Wallis groups. Columns include parameter, groupA, groupB, 
%       the standard error of each, median rank, rank difference, 
%       confidence intervals, and p-value.
%   ANOVA: table of multi-group interaction ANOVA p-values, with each
%       parameter as a column and each interacting group as a row. In
%       addition to the groups used for basic stastical analysis, the
%       volume of the reference dose falling within the median absolute
%       range is included as a continuous variable.
%   DunnMeans: table of Dunn's multiple comparison test statistics based on 
%       the ANOVA groups. Columns include parameter, groupA, groupB, the 
%       standard error of each, mean, mean difference, confidence 
%       intervals, and p-value.
%   MANOVA: table of p-values and MANOVA test statistics for each group.
%   PairwiseB: if a second reference dataset is provided, a table of two-
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
% results = CompareD4Dose(meas, tps, [], 'GammaAbs', 2, 'GammaDTA', 2);
%
% % Extract a vector of pass rates for all Plans using 'TrueBeam' and 6 MV
% pass = results.Plans.GammaPassRateGlobal(strcmp(results.Plans.Machine, ...
%   'TrueBeam') & results.Plans.Energy == 6);
%
% % Find the mean local Gamma pass rate for all Plans
% avg = results.Stats.mean(strcmp(results.Stats.Variable, ...
%       'GammaPassRateLocal') & strcmp(results.Stats.Group, 'all'));
%
% % Plot a scatter plot matric of the local pass rate, mean gamma, and 
% % median absolute difference grouped by Machine
% gplotmatrix([results.Plans.GammaPassRateLocal, ...
%   results.Plans.GammaMeanLocal, results.Plans.MedianAbsDiff], [], 'bg', ...
%   '..', [],'on','', {'Pass Rate', 'Mean Gamma', 'Abs Diff'},  ...
%   {'Pass Rate', 'Mean Gamma', 'Abs Diff'});
%
% % Compare the pass rates between the original and a new TPS model
% pairwise = CompareD4Dose(meas, tps, '../test_data/new_model/');
%
% % Determine if the median absolute dose difference changed significantly
% % between the two models based on the t-test p-value
% pval = results.pairwiseB.p(strcmp(results.stats.parameter, ...
%       'MedianAbsDiff') & strcmp(results.stats.group, 'all'));
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
results.ComputeStats = true;
results.BeamAnalysis = true;
results.IgnoreName = false;
results.IgnoreID = false;
results.AbsRange = [50 500];
results.GammaAbs = 3;
results.GammaDTA = 3;
results.GammaRange = [20 500];
results.GammaLimit = 2;
results.GammaRes = 10;
results.RefDose = 'measmax';
results.PlanClasses = cell(0,2);

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
    
    % Update RefDose to 'Planmax' as that is the default for ScandiDos
    results.RefDose = 'Planmax';
end

% Add dicom_tools submodule to search path
[path, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(path, 'dicom_tools'));

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

%% Scan reference dose folders for DICOM files
% Display progress bar
if usejava('jvm') && feature('ShowFigureWindows') && results.Progress
    progress = waitbar(0, 'Scanning first reference dataset for RT Plans');
end

% Scan first reference folder
if exist('Event', 'file') == 2
    Event('Scanning first reference dataset for RT Plans');
    t = tic;
end

% Store cell array of DICOM files as ref{1}
ref{1} = ScanDICOMPath(varargin{2}, 'Progress', false);

% Scan second reference folder
if nargin >= 3 && ~isempty(varargin{3})
    
    % Log action
    if exist('Event', 'file') == 2
        Event('Scanning second reference dataset for RT Plans');
        t = tic;
    end
    
    % Update waitbar
    if exist('progress', 'var') && ishandle(progress) && results.Progress
        waitbar(0.1, progress, 'Scanning second reference dataset for RT Plans');
    end
    
    % Store second cell array of DICOM files as ref{2}
    ref{2} = ScanDICOMPath(varargin{2}, 'Progress', false);
end

%% Scan directory of measurement files
% Update waitbar
if exist('progress', 'var') && ishandle(progress) && results.Progress
    waitbar(0.2, progress, 'Scanning measurement files');
end

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
% Initialize Plan list, beam list, and gamma arrays return fields
fixedNames = {'MeasurementFile', 'Patient', 'ID', 'Plan', 'Machine', ...
    'Energy', 'PlanClass', 'Orientation'};
fixedUnits = {'', '', '', '', '', 'MV', '', ''};
varNames = {'ReferenceFile', 'ReferenceDose', 'GammaPassRateGlobal', ...
    'GammaMeanGlobal', 'GammaMaxGlobal', 'GammaPassRateLocal', ...
    'GammaMeanLocal', 'GammaMaxLocal', 'MedianAbsDiff', 'AbsVolume', ...
    'GammaVolume'};
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
results.Datasets = length(ref);

% Create tables
results.Plans = array2table(zeros(0,length(fixedNames)+length(v)), ...
    'VariableNames', horzcat(fixedNames, v));
results.Plans.Properties.VariableUnits = horzcat(fixedUnits, u);
if results.BeamAnalysis
    results.Beams = array2table(zeros(0,1+length(fixedNames)+length(v)), ...
        'VariableNames', horzcat(fixedNames, {'Beam'}, v));
    results.Beams.Properties.VariableUnits = horzcat(fixedUnits, {''}, u);
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
    if ismember(lower(ext), {'.xlsx', '.xls', '.txt', '.csv'})
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
    if ~isfield(delta4, 'Data') || isempty(delta4.Data)
        continue;
    end

    % Split name
    nd = strsplit(delta4.Name, ',');
    
    % Match plan class (will stop after the first match is found)
    class = 'Undefined';
    for j = 1:size(results.PlanClasses, 1)
        if ~isempty(regexpi(delta4.Plan, results.PlanClasses{j,2}))
            class = results.PlanClasses{j,2};
            break;
        end
    end
    
    % Set Orientation (if there are more than 400 diodes with same IEC Z 
    % value, this is a sagittal/coronal measurement)
    if sum(delta4.Data(:,2) == mode(delta4.Data(:,2))) > 400
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
            if strcmp(ref{r}{j,10}, 'PLAN') && (strcmpi(delta4.ID, ...
                    ref{r}{j,6}) || results.IgnoreID) && ...
                    strcmpi(delta4.Plan, ref{r}{j,9}) && ...
                    ((contains(nd{1}, nr{1}, 'IgnoreCase', true) && ...
                    (length(nd) < 2 || length(nr) < 2 || ...
                    contains(nd{2}, nr{2}, 'IgnoreCase', true))) || ...
                    results.IgnoreName) && contains(delta4.Value, ...
                    'Absolute dose', 'IgnoreCase', true)
                
                matched(r) = j;
                break;
            end
        end
    end
    
    % If all references were matched
    if all(matched)

        % If all beams have the same energy, use that, otherwise report Mixed
        if length(unique(ref{r}{matched(1),13})) == 1
            energy = ref{r}{matched(1),13}{1};
        else
            energy = 'Mixed';
        end
        
        % Initialize results
        result = {[name, ext], delta4.Name, delta4.ID, delta4.Plan, ...
            ref{r}{matched(1),12}{1}, energy, class, orient};

        % Loop through each reference
        for j = 1:length(matched)
            
            % Load dose
            dose = LoadDICOMDose(ref{r}{matched(j),2}, ref{r}{matched(j),1});
        
            % Set reference value
            switch results.RefDose
                case 'measmax'
                    refval = max(delta4.Data(:,5));
                case 'Planmax'
                    refval = max(max(max(dose.data)));
                otherwise
                    refval = 0;
            end
            
            % Compute dose volumes
            absvol = (sum(sum(sum(dose.data >= refval * ...
                results.AbsRange(1)/100))) - sum(sum(sum(dose.data > ...
                refval * results.AbsRange(2)/100)))) * prod(dose.width);
            gamvol = (sum(sum(sum(dose.data >= refval * ...
                results.GammaRange(1)/100))) - sum(sum(sum(dose.data > ...
                refval * results.GammaRange(2)/100)))) * prod(dose.width);
            
            % Log calculation
            if exist('Event', 'file') == 2
                Event(['Computing plan Gamma table between ', name, ext, ...
                    'and ', ref{r}{matched(j),1}]);
            end
            
            % Append reference file and gamma table
            result = [result, ref{r}{matched(j),1}, refval, ...
                gammaTable(delta4.Data, dose, refval, results), absvol, ...
                gamvol]; %#ok<*AGROW>
        end
        
        % Append result onto results
        results.Plans = [results.Plans; result];
        
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
        for b = 1:length(delta4.Beams)
            
            % Initialize beam IDs
            matched = zeros(1, length(ref));
            
            % If parse was successful, match Plan to reference dose files
            for r = 1:length(ref)
                for j = 1:length(ref{r})
        
                    % Split name
                    nr = strsplit(ref{r}{j,5}, ',');

                    % If reference is a BEAM and matches the ID, Plan, pt , 
                    % and beam name
                    if strcmp(ref{r}{j,10}, 'BEAM') && (strcmpi(delta4.ID, ...
                            ref{r}{j,6}) || results.IgnoreID) && ...
                            strcmpi(delta4.Plan, ref{r}{j,9}) && ...
                            ((contains(nd{1}, nr{1}, 'IgnoreCase', true) && ...
                            (length(nd) < 2 || length(nr) < 2 || ...
                            contains(nd{2}, nr{2}, 'IgnoreCase', true))) || ...
                            results.IgnoreName) && ...
                            strcmpi(delta4.Beams{b}.Name, ref{1}{j,11}) && ...
                            contains(delta4.Beams{b}.Value, ...
                            'Absolute dose', 'IgnoreCase', true)
                        
                        matched(r) = j;
                        break;
                    end
                end
            end
            
            % If all references were matched
            if all(matched)
                
                % Initialize results
                result = {[name, ext], delta4.Name, delta4.ID, ...
                    delta4.Plan, ref{r}{matched(1),12}, ...
                    ref{r}{matched(1),13}, class, orient, ...
                    delta4.Beams{b}.Name};

                % Loop through each reference
                for j = 1:length(matched)

                    % Load dose
                    dose = LoadDICOMDose(ref{r}{matched(j),2}, ...
                        ref{r}{matched(j),1});

                    % Set reference value
                    switch results.RefDose
                        case 'measmax'
                            refval = max(delta4.Beams{b}.Data(:,5));
                        case 'Planmax'
                            refval = max(max(max(dose.data)));
                        otherwise
                            refval = 0;
                    end
                    
                    % Compute dose volumes
                    absvol = (sum(sum(sum(dose.data >= refval * ...
                        results.AbsRange(1)/100))) - sum(sum(sum(dose.data > ...
                        refval * results.AbsRange(2)/100)))) * prod(dose.width);
                    gamvol = (sum(sum(sum(dose.data >= refval * ...
                        results.GammaRange(1)/100))) - sum(sum(sum(dose.data > ...
                        refval * results.GammaRange(2)/100)))) * prod(dose.width);

                    % Log calculation
                    if exist('Event', 'file') == 2
                        Event(['Computing beam Gamma table between ', ...
                            name, ext, 'and ', ref{r}{matched(j),1}]);
                    end
                    
                    % Append reference file and gamma table
                    result = [result, ref{r}{matched(j),1}, refval, ...
                        gammaTable(delta4.Beams{b}.Data, ...
                        dose, refval, results), absvol, gamvol]; %#ok<*AGROW>
                end

                % Append result onto results
                results.Beams = [results.Beams; result];
            end
        end
    end
end

% Clear temporary variables
clear b i j r u v nd nr class delta4 dose path name ext fixedNames ... 
    fixedUnits matched meas refval result varNames varUnits orient absvol ...
    gamvol ref energy;

%% Compute Statistics
if results.ComputeStats
    if exist('progress', 'var') == 1 && ishandle(progress)
        waitbar(0.95, progress, 'Computing statistics');
    end
    results = ComputeD4Stats(results);
end

%% Finish up
% Close waitbar
if exist('progress', 'var') == 1 && ishandle(progress)
    waitbar(1.0, progress);
    close(progress);
end

% Remove execution flags from results
results = rmfield(results, 'Progress');
results = rmfield(results, 'BeamAnalysis');
results = rmfield(results, 'ComputeStats');
results = rmfield(results, 'IgnoreName');
results = rmfield(results, 'IgnoreID');


% Log completion
if exist('Event', 'file') == 2
    Event(sprintf(['Comparison successfully completed in %0.3f seconds, ', ...
        'scanning %i Plans and %i beams'], toc(t), size(results.Plans,1), ...
        size(results.Beams,1)));
end

% Clear temporary variables
clear t progress;

end

%% Compute Gamma Table Subfunction
function stats = gammaTable(meas, dose, refval, criteria)
% gammaTable is called by CompareD4Dose and computes a 2D gamma table based
% on a provided measurement array and reference dose volume, using the
% results GammaAbs, GammaDTA, GammaRange, GammaRes, and GammaLimit fields
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
stats{7} = median(absDiff(meas(:,5) > refval * criteria.AbsRange(1)/100 & ...
    meas(:,5) < refval * criteria.AbsRange(2)/100));

% Apply gamma range to measured values
meas = meas(meas(:,5) > refval * criteria.GammaRange(1)/100 & ...
    meas(:,5) < refval * criteria.GammaRange(2)/100,:);

% Compute shifts (in cm)
shifts = repmat(-criteria.GammaRes * criteria.GammaLimit:criteria.GammaRes * ...
    criteria.GammaLimit, size(meas,1), 1) / ...
    criteria.GammaRes * max(criteria.GammaDTA)/10;

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
loc = zeros(size(i,1), length(criteria.GammaAbs), length(criteria.GammaDTA));
glob = zeros(size(i,1), length(criteria.GammaAbs), length(criteria.GammaDTA));

% Loop through Abs, DTA criteria
for a = 1:length(criteria.GammaAbs)
    for d = 1:length(criteria.GammaDTA)
        
        % Compute global gamma based on refval as min of each shifted position
        glob(:,a,d) = sqrt(min(((repmat(meas(:,5),1,size(i,2)) - i) ./ ...
            (refval * criteria.GammaAbs(a)/100)).^2 + ...
            (repmat(shifts,1,3)/(criteria.GammaDTA(d)/10)).^2,[],2));
        
        % Compute local gamma as min of each shifted position
        loc(:,a,d) = sqrt(min(((repmat(meas(:,5),1,size(i,2)) - i) ./ ...
            (i * criteria.GammaAbs(a)/100)).^2 + ...
            (repmat(shifts,1,3)/(criteria.GammaDTA(d)/10)).^2,[],2));
    end
end

% Compute pass rate, mean, and maximum Gamma statistics for global gamma
stats{1} = zeros(length(criteria.GammaAbs), length(criteria.GammaDTA));
for a = 1:length(criteria.GammaAbs)
    for d = 1:length(criteria.GammaDTA)
        stats{1}(a,d) = sum(glob(:,a,d) <= 1) / size(glob,1) * 100;
    end
end
stats{2} = squeeze(mean(glob,1));
stats{3} = squeeze(max(glob,[], 1));

% Compute pass rate, mean, and maximum Gamma statistics for local gamma
stats{4} = zeros(length(criteria.GammaAbs), length(criteria.GammaDTA));
for a = 1:length(criteria.GammaAbs)
    for d = 1:length(criteria.GammaDTA)
        stats{4}(a,d) = sum(loc(:,a,d) <= 1) / size(loc,1) * 100;
    end
end
stats{5} = squeeze(mean(loc,1));
stats{6} = squeeze(max(loc,[], 1));

% Clear temporary variables
clear measA measB measC a d i shifts glob loc meshA meshB meshC absDiff;

end
