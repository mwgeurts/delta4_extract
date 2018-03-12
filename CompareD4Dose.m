function results = CompareD4Dose(varargin)
% CompareD4Dose evaluated a file, list, or directory of Delta4 measurements 
% to RT dose volumes, including calculation of Gamma pass rates. The 
% function accepts either one or two sets of dose volume files; if two, the 
% function will conduct a pairwise comparison for each measurement and
% report the comparison statistics. In this manner, the tool can be used to
% perform a bulk comparison of how Delta4 results change with updates to a
% beam model.
%
% During excecution, this function will call the function ScanDICOMPath
% from the dicom_tools submodule.
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
%       'GammaRange'). See see code documentation for a full list of 
%       options. Alternatively, a 'Report' option can be provided with the
%       value set to a Delta4 report structure obtained from ParseD4Report.
%       In this case, all relevant configuration options will be set to
%       match what was set in the report. This is helpful when trying to
%       replicate the same analysis as was performed in the Delta4
%       software.

% Verify at least two inputs were provided
if nargin < 2
    error(['At least two inputs are required to execute CompareD4Dose: a ', ...
        'target measurement file/folder/list and a destination ', ...
        'file/folder/list']);
end
    
% Set default options
results.Progress = true;
results.gammaAbs = 1:5;
results.gammaDta = 1:5;
results.gammaRange = [20 500];
results.gammaRes = 0.1;
results.absRange = [50 500];

% Update options
for i = 4:nargin
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

% Display progress bar
if usejava('jvm') && feature('ShowFigureWindows') && results.Progress
    progress = waitbar(0, 'Scanning first reference dataset for RT Plans');
end

% Scan first reference folder
if iscell(varargin{2})
    refA = ScanDICOMPath(varargin{2}, 'Progress', false);
end

% Update waitbar
if exist('progress', 'var') && ishandle(progress) && results.Progress
    waitbar(0.1, progress, 'Scanning first reference dataset for RT Plans');
end

% Scan second reference folder
if nargin >= 3 && ~isempty(varargin{3})
    refB = ScanDICOMPath(varargin{2}, 'Progress', false);
else
    refB = [];
end

% Update waitbar
if exist('progress', 'var') && ishandle(progress) && results.Progress
    waitbar(0.2, progress);
end

% Scan the directory for DICOM files
if exist('Event', 'file') == 2
    Event('Scanning for measurement files');
    t = tic;
end

% Scan measurements folder directory contents
if iscell(varargin{1})
    meas = varargin{1};
elseif isfolder(varargin{1})
    meas = dir(fullfile(varargin{1}, '**'));
else
    meas = {varargin{1}};
end

% Initialize plan list, beam list, and gamma arrays return fields
if isempty(refB)
    results.plans = array2table(zeros(0,5), 'VariableNames', ...
        {'Patient', 'ID', 'Plan', 'MeasurementFile', 'ReferenceFile'});
    results.beams = array2table(zeros(0,6), 'VariableNames', ...
        {'Patient', 'ID', 'Plan', 'Beam', 'MeasurementFile', ...
        'ReferenceFile'});

    results.planGamma = zeros(0, length(results.gammaAbs), ...
        length(results.gammaDta));
    results.beamGamma = zeros(0, length(results.gammaAbs), ...
        length(results.gammaDta));
else
    results.plans = array2table(zeros(0,6), 'VariableNames', ...
        {'Patient', 'ID', 'Plan', 'MeasurementFile', 'ReferenceFileA', ...
        'ReferenceFileB'});
    results.beams = array2table(zeros(0,7), 'VariableNames', ...
        {'Patient', 'ID', 'Plan', 'Beam', 'MeasurementFile', ...
        'ReferenceFileA', 'ReferenceFileB'});

    results.planGamma = zeros(0, length(results.gammaAbs), ...
        length(results.gammaDta), 2);
    results.beamGamma = zeros(0, length(results.gammaAbs), ...
        length(results.gammaDta), 2);
end

% Loop through measurement list
for i = 1:length(meas)

    % Update waitbar
    if exist('progress', 'var') && ishandle(progress) && results.Progress
        waitbar(0.2 + 0.8*i/length(meas), progress, ...
            'Comparing measurement files to reference dose volumes');
    end
    
    % If the folder content is . or .. or a folder, skip to next file
    if strcmp(list(i).name, '.') || strcmp(list(i).name, '..') ...
            || list(i).isdir == 1
        continue;
    end
    
    % Parse file extension
    [path, name, ext] = fileparts(fullfile(varargin{1}, list(i).name));
    
    % If extension matches potential report format, try to parse it
    if ismember(ext, {'.xlsx', '.xls', '.txt', '.csv'})
        try
            if exist('Event', 'file') == 2
                Event(['Parsing Delta4 measurement file ', name, ext]);
            end
            delta4 = ParseD4Tables([path, name, ext]);
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
    na = strsplit(delta4.name, ',');
    
    % If parse was successful, match to reference dose file
    for j = 1:size(refA, 1)

        % Split name
        nb = strsplit(refA{j,5}, ',');
        
        % If first, last, ID, and plan matches
        if strcmp(refA{j,3}, 'RTDOSE') && strcmpi(delta4.ID, refA{j,6}) && ...
                strcmpi(delta4.plan, refA{j,9}) && ...
                contains(na{1}, nb{1}, 'IgnoreCase', true) && ...
                (length(na) < 2 || length(nb) < 2 || ...
                contains(na{2}, nb{2}, 'IgnoreCase', true))
                
            % If this is the plan dose and measured plan dose exists
            if isempty(refB) && strcmp(refA{j,10}, 'PLAN') && ...
                    isfield(delta4, 'data') && isfield(delta4, 'value') && ...
                    contains(delta.value, 'Absolute dose', 'IgnoreCase', ...
                    true)
                
                % Load reference dose
                doseA = LoadDICOMDose(fullfile(refA{j,2}, refA{j,1}));
                
                % Add plan to plans list
                results.plans = [results.plans; {delta4.name, delta4.ID, ... 
                    delta4.plan, [name, ext], refA{j,1}}];
                
                % Compute gamma
                results.planGamma(size(results.planGamma,1)+1,:,:) = ...
                    gammaTable(delta4.data, doseA);
                
            % Otherwise, if this is a beam dose
            elseif isempty(refB) && strcmp(refA{j,10}, 'BEAM') && ...
                    isfield(delta4, 'beams')
           
                % Loop through measured beam data
                for k = 1:length(delta4.beams)
                   
                    % If beam names match
                    if strcmpi(delta4.beams{k}.name, refA{j,11}) && ...
                            contains(delta.value, 'Absolute dose', ...
                            'IgnoreCase', true)
                    
                        % Load reference dose
                        doseA = LoadDICOMDose(fullfile(refA{j,2}, ...
                            refA{j,1}));
                        
                        % Add beam to beams list
                        results.beams = [results.plans; {delta4.name, ...
                            delta4.ID,  delta4.plan, delta4.beams{k}.name, ...
                            [name, ext], refA{j,1}}];

                        % Compute gamma
                        results.beamGamma(size(results.beamGamma,1)+1,:,:) = ...
                            gammaTable(delta4.beams{k}.data, doseA);
                    end
                end
                
            % Otherwise, if refB is not empty, also match a refB file
            
            
            
            
            
            
            
            end
        end
    end
end

% Close waitbar
if exist('progress', 'var') == 1 && ishandle(progress)
    close(progress);
end

% Remove Progress field from results
results = rmfield(results, 'Progress');

% Log completion
if exist('Event', 'file') == 2
    Event(sprintf(['Comparison succefullly completed in %0.3f seconds'], ...
        toc(t)));
end

% Clear temporary variables
clear i t meas refA refB progress;

end

%% Compute Gamma Subfunction
function gamma = gammaTable(meas, dose)
% gammaTable is called by CompareD4Dose and computes a 2D gamma table based
% on a provided measurement array and reference dose volume, using the
% results gammaAbs, gammaDta, gammaRange, and gammaRange fields as
% criteria.

% Compute mesh grid for reference dose
[meshX, meshY, meshZ] = meshgrid(single(dose.start(2) + ...
    dose.width(2) * (size(dose.data,2) - 1):-dose.width(2):dose.start(2)), ...
    single(dose.start(1):dose.width(1):dose.start(1) + dose.width(1)...
    * (size(dose.data,1) - 1)), single(dose.start(3):dose.width(3):...
    dose.start(3) + dose.width(3) * (size(dose.data,3) - 1)));


end