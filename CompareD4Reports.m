function results = CompareD4Reports(varargin)


% Verify at least one input was provided
if nargin < 1
    error(['At least one input is required to execute CompareD4Reports, ', ...
        'either a string containing the target folder path or a cell array', ...
        ' of files']);
end

% Set default options
results.Progress = true;
results.ComputeStats = true;
results.BeamAnalysis = true;
results.PlanClasses = cell(0,2);

% Update options
for i = 2:2:nargin
    results.(varargin{i}) = varargin{i+1};
end

% Display progress bar
if usejava('jvm') && feature('ShowFigureWindows') && results.Progress
    progress = waitbar(0, 'Scanning folder for Delta4 reports');
end

% Scan first reference folder
if exist('Event', 'file') == 2
    Event('Scanning folder for Delta4 reports');
    t = tic;
end

% Store reports cell array
if iscell(varargin{1})
    reports = varargin{1};
elseif isfolder(varargin{1})
    reports = dir(fullfile(varargin{1}, '**'));
else
    reports = dir(varargin{1});
end

% Define variable names from ParseD4Report that should be stored
results.Plans = array2table(zeros(0, 33), ...
    'VariableNames', {'MeasurementFile', 'Title', 'Patient', 'ID', ...
    'Clinic', 'Plan', 'PlanDate', 'PlanUser', 'PlanClass', 'MeasDate', ...
    'MeasUser', 'ReviewStatus', 'ReviewDate', 'ReviewUser', 'Comments', ...
    'Machine', 'Temperature', 'Reference', 'NormDose', 'AbsPassRate', ...
    'DTAPassRate', 'GammaPassRate', 'MedianAbsDiff', 'Energy', 'AbsRange', ...
    'AbsPassLimit', 'DTARange', 'DTAPassLimit', 'GammaRange', 'GammaAbs', ...
    'GammaDTA', 'GammaPassLimit', 'GammaTable'});
results.Plans.Properties.VariableUnits = {'', '', '', '', '', '', '', '', ...
    '', '', '', '', '', '', '', '', 'C', '', 'Gy', '%', '%', '%', '%', 'MV', ...
    '%,%', '%,%', '%/mm', '%,mm', '%,%', '%', 'mm', '%,-', '%'};

% If beam data should also be stored
if results.BeamAnalysis
    results.Beams = array2table(zeros(0, 13), 'VariableNames', ...
        {'MeasurementFile', 'Patient', 'ID', 'Plan', 'PlanClass', 'Name', ...
        'Energy', 'DailyCF', 'NormDose', 'AbsPassRate', 'DTAPassRate', ...
        'GammaPassRate', 'MedianAbsDiff'});
    results.Beams.Properties.VariableUnits = {'', '', '', '', '', '', 'MV', ...
        '', 'Gy', '%', '%', '%', '%'};
end

%% Parse reports
% Loop through measurement list
for i = 1:length(reports)

    % Update waitbar
    if exist('progress', 'var') && ishandle(progress) && results.Progress
        waitbar(0.9*i/length(reports), progress, ...
            sprintf('Parsing Delta4 report %i of %i', i, length(reports)));
    end
    
    % If the folder content is . or .. or a folder, skip to next file
    if strcmp(reports(i).name, '.') || strcmp(reports(i).name, '..') ...
            || reports(i).isdir == 1
        continue;
    end
    
    % Parse file extension
    [path, name, ext] = fileparts(fullfile(reports(i).folder, reports(i).name));
    
    % If extension matches report format, try to parse it
    if strcmpi(ext, '.pdf')
        try
            if exist('Event', 'file') == 2
                Event(['Parsing Delta4 report file ', name, ext]);
            end
            delta4 = ParseD4Report(fullfile(path, [name, ext]));
        
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
        
    % Verify the patient and plan names at least exist
    if ~isfield(delta4, 'Plan') || isempty(delta4.Plan) || ...
            ~isfield(delta4, 'Name') || isempty(delta4.Name)
        continue;
    end

    % Initialize plan row
    r = cell(1, length(results.Plans.Properties.VariableNames));
    r{1} = [name, ext];
    
    % Loop through each variable, storing the result if it exists
    for j = 2:length(results.Plans.Properties.VariableNames)
        
        % If report contains that field, store it
        if isfield(delta4, results.Plans.Properties.VariableNames{j})
            
            % If input is a datetime, store as a string
            if isdatetime(delta4.(...
                    results.Plans.Properties.VariableNames{j}))
                r{j} = datestr(delta4.(...
                    results.Plans.Properties.VariableNames{j}));
            
            % Otherwise, if a table, store as an array
            elseif istable(delta4.(...
                    results.Plans.Properties.VariableNames{j}))
                r{j} = table2array(delta4.(...
                    results.Plans.Properties.VariableNames{j}));
            
            else
                r{j} = delta4.(results.Plans.Properties.VariableNames{j});
            end
        end
        
        % If value is a datetime, convert it to a string
        if isdatetime(r{j})
            r{j} = datestr(r{j});
        end
    end
    
    % Store patient name
    r{strcmp(results.Plans.Properties.VariableNames, ...
        'Patient')} = delta4.Name;
    
    % Match plan class (will stop after the first match is found)
    class = 'Undefined';
    for j = 1:size(results.PlanClasses, 1)
        if ~isempty(regexpi(delta4.Plan, results.PlanClasses{j,2}))
            class = results.PlanClasses{j,2};
            break;
        end
    end
    r{strcmp(results.Plans.Properties.VariableNames, 'PlanClass')} = class;
    
    % If all beams have the same energy, use that, otherwise report Mixed
    if isfield(delta4, 'Beams') && istable(delta4.Beams) && ...
            size(delta4.Beams,1) > 0
        if length(unique(delta4.Beams.Energy)) == 1
            r{strcmp(results.Plans.Properties.VariableNames, 'Energy')} = ...
            	delta4.Beams.Energy{1};
        else
            r{strcmp(results.Plans.Properties.VariableNames, 'Energy')} = ...
                'Mixed';
        end
    end
    
    % Append to plan array
    results.Plans = [results.Plans; r];
end

% Clear temporary variables
clear i j r class delta4 path name ext reports;

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

% Log completion
if exist('Event', 'file') == 2
    Event(sprintf(['Comparison successfully completed in %0.3f seconds, ', ...
        'scanning %i Plans and %i beams'], toc(t), size(results.Plans,1), ...
        size(results.Beams,1)));
end

% Clear temporary variables
clear t progress;


