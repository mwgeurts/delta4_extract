function results = CompareD4Dose(varargin)
% CompareD4Dose evaluated a file, list, or directory of Delta4 measurements 
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
%       same analysis as was performed in the Delta4 software. Finally, the 
%       name/pair 'BeamStats' => false can be provided to skip computation 
%       of individual beam statistics.
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
%
% The following example illustrates how to use this function:
%
% % Define folder containing Delta4 measurements
% meas = '../test_data/VMAT Tests/';
% tps = '../test_data/RayStation/';
%
% % Execute evaluation, using folders above and 2%/2mm criteria
% results = CompareD4Dose(meas, tps, 'gammaAbs', 2, 'gammaDta', 2);
%
% % Extract pass rates for all plans using 'TrueBeam' and 6 MV
% pass = results.plans.gammaPassRateGlobal(strcmp(results.plans.machine, ...
%   'TrueBeam') & results.plans.energy == 6);
%
% % Compute percentage of plans with 95% pass rate or greater
% sum(pass > 1)/length(pass)*100
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
results.gammaAbs = 3;
results.gammaDta = 3;
results.gammaRange = [20 500];
results.gammaLimit = 2;
results.gammaRes = 10;
results.BeamStats = true;
results.refDose = 'measmax';
results.absRange = [50 500];

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
fixedNames = {'measurementFile', 'patient', 'ID', 'plan', 'machine', 'energy'};
varNames = {'referenceFile', 'referenceDose', 'gammaPassRateGlobal', ...
    'gammaMeanGlobal', 'gammaMaxGlobal', 'gammaPassRateLocal', ...
    'gammaMeanLocal', 'gammaMaxLocal', 'medianAbsDiff'};

% If only one reference exists
if length(ref) == 1
    v = varNames;

% If multiple references exist
else
    
    % Initialize variable names
    v = cell(1,length(varNames)*length(ref));
    for i = 1:length(ref)
        for j = 1:length(varNames)
            
            % Append a letter (A, B, C, ...) onto each variable name
            v{(i-1)*length(varNames)+j} = [varNames{j}, char(i + 64)];
        end
    end
end

% Create tables
results.plans = array2table(zeros(0,length(fixedNames)+length(v)), ...
    'VariableNames', horzcat(fixedNames, v));
if results.BeamStats
    results.beams = array2table(zeros(0,1+length(fixedNames)+length(v)), ...
        'VariableNames', horzcat(fixedNames, {'Beam'}, v));
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

%% Process measurements
% Loop through measurement list
for i = 1:length(meas)

    % Update waitbar
    if exist('progress', 'var') && ishandle(progress) && results.Progress
        waitbar(0.2 + 0.8*i/length(meas), progress, ...
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

    % Split name
    nd = strsplit(delta4.name, ',');
    
    % Initialize plan IDs
    matched = zeros(1, length(ref));
    
    % If parse was successful, match plan to reference dose files
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
            ref{r}{matched(1),12}{1}, mean(cell2mat(ref{r}{matched(1),13}))};

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
            
            % Log calculation
            if exist('Event', 'file') == 2
                Event(['Computing plan Gamma table between ', name, ext, ...
                    'and ', ref{r}{matched(j),1}]);
            end
            
            % Append reference file and gamma table
            result = [result, ref{r}{matched(j),1}, refval, gammaTable(delta4.data, ...
                dose, refval, results)]; %#ok<*AGROW>
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

    % Next, match each beam to reference dose files
    if results.BeamStats && isfield(delta4, 'beams')
        
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
                    ref{r}{matched(1),13}, delta4.beams{b}.name};

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

                    % Log calculation
                    if exist('Event', 'file') == 2
                        Event(['Computing beam Gamma table between ', ...
                            name, ext, 'and ', ref{r}{matched(j),1}]);
                    end
                    
                    % Append reference file and gamma table
                    result = [result, ref{r}{matched(j),1}, refval, ...
                        gammaTable(delta4.beams{b}.data, ...
                        dose, refval, results)]; %#ok<*AGROW>
                end

                % Append result onto results
                results.beams = [results.beams; result];
            end
        end
    end
end

%% Finish up
% Close waitbar
if exist('progress', 'var') == 1 && ishandle(progress)
    close(progress);
end

% Remove Progress and BeamStats field from results
results = rmfield(results, 'Progress');
results = rmfield(results, 'BeamStats');

% Log completion
if exist('Event', 'file') == 2
    Event(sprintf(['Comparison successfully completed in %0.3f seconds, ', ...
        'scanning %i plans and %i beams'], toc(t), size(results.plans,1), ...
        size(results.beams,1)));
end

% Clear temporary variables
clear i t meas ref{1} ref{2} refval progress;

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
