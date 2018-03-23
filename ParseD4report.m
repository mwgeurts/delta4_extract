function delta4 = ParseD4Report(content)
% ParseD4Report reads a ScandiDos Delta4 report into a MATLAB structure.
% The function input argument can either be a string containing a file path
% and/or Name corresponding to a report PDF file, or a cell array of text 
% data from that PDF. If a PDF file, this function will call XpdfText to
% extract the file contents (via the xpdf_tools submodule).
%
% Upon successful completion, this function will return a structure
% containing the following fields:
%   Title: string containing report title
%   Name: string containing patient name
%   ID: string containing patient ID
%   Clinic: cell array containing clinic name and address
%   Plan: string containing plan name
%   PlanDate: planned datetime
%   PlanUser: string containing planned user (if present)
%   MeasDate: measured datetime
%   MeasUser: string containing planned user (if present)
%   Comments: cell array of comments
%   Temperature: double containing temperature
%   Reference: string containing reference (i.e. 'Planned Dose')
%   NormDose: double containing normalization dose, in Gy
%   AbsPassRate: double containing Absolute pass rate, as a percentage
%   DTAPassRate: double containing DTA pass rate, as a percentage
%   GammaPassRate: double containing Gamma pass rate, as a percentage
%   MedianAbsDiff: double containing median dose deviation, as a percentage
%   Beams: table of each beams containing the following 
%       columns: Name, DailyCF, NormDose, AbsPassRate, DTAPassRate, 
%       GammaPassRate, and MedianAbsDiff.
%   AbsRange: 2 element vector of dose deviation range, as percentages
%   AbsPassLimit: 2 element vector of dose deviation acceptance criteria,
%       as percentages (i.e. [90 3] means 90% within 3%)
%   DTARange: 2 element vector of DTA range, in %/mm
%   DTAPassLimit: 2 element vector of DTA acceptance criteria, in % and mm
%   GammaRange: 2 element vector of Gamma range
%   GammaAbs: double containing Gamma Absolute criterion as a percentage
%   GammaDTA: double containing Gamma DTA criterion in mm
%   GammaPassLimit: 2 element vector of Gamma acceptance criteria (i.e.
%       [95 1] means 95% less than 1)
%   GammaTable: structure containing Gamma Index Evaluation table (if
%       present) with the following fields: DTA, Abs, PassRate
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

% If content is a file Name (with PDF extension), read file contents in
if ~iscell(content) && endsWith(content, '.pdf', 'IgnoreCase',true)
    
    % Add xpdf_tools submodule to search path
    [path, ~, ~] = fileparts(mfilename('fullpath'));
    addpath(fullfile(path, 'xpdf_tools'));

    % Check if MATLAB can find XpdfText
    if exist('XpdfText', 'file') ~= 2

        % If not, throw an error
        if exist('Event', 'file') == 2
            Event(['The xpdf_tools submodule does not exist in the search path. ', ...
                'Use git clone --recursive or git submodule init followed by git ', ...
                'submodule update to fetch all submodules'], 'ERROR');
        else
            error(['The xpdf_tools submodule does not exist in the search path. ', ...
                'Use git clone --recursive or git submodule init followed by git ', ...
                'submodule update to fetch all submodules']);
        end
    end

    % Read PDF text from pages 1&2 into contents
    pages = XpdfText(content);
    content = horzcat(pages{1}, pages{2});
    
    % Clear temporary variables
    clear path pages;
end

% Log start
if exist('Event', 'file') == 2
    Event('Parsing data from Delta4 report');
    tic;
end

% Initialize empty return variable
delta4 = struct;

% If Plan report is from version April 2016 or later
if length(content{7}) >= 7 && strcmp(content{7}(1:7), 'Clinic:')
    
    % Store Title, patient Name, and ID
    delta4.Title = strtrim(content{1});
    delta4.Name = strtrim(content{3});
    delta4.ID = strtrim(content{5});

    % Initialize row counter
    r = 6;
    
else 
    % Store Title and patient Name
    fields = strsplit(content{1}, '   ');
    if ~isempty(fields)
        delta4.Title = strtrim(fields{1});
    else
        delta4.Title = '';
    end
    if length(fields) > 1
        delta4.Name = strtrim(fields{2});
    else
        delta4.Name = '';
    end
    for i = 3:length(fields)
        delta4.Name = [delta4.Name, ' ', strtrim(fields{i})];
    end

    % Store patient ID
    delta4.ID = strtrim(content{3});

    % Initialize row counter
    r = 4;
end

% Loop through rows until Clinic info is found
while r < length(content)
    
    % If row starts with 'Clinic:'
    if length(content{r}) >= 7 && strcmp(content{r}(1:7), 'Clinic:')
        content{r} = content{r}(8:end);
        delta4.Clinic = cell(0);
        break;
    else
        r = r + 1;
    end
end

% Store Clinic contact info, followed by Plan Name
while r < length(content)
    
    % If row starts with 'Plan:'
    if length(content{r}) >= 5 && strcmp(content{r}(1:5), 'Plan:')
        delta4.Plan = strtrim(content{r}(6:end));
        break;
    else
        if ~isempty(content{r})
            delta4.Clinic = vertcat(delta4.Clinic, strtrim(content{r}));
        end
        r = r + 1;
    end
end

% Loop through rows until Planned date info is found
while r < length(content)
    
    % If row starts with 'Planned:'
    if length(content{r}) >= 8 && strcmp(content{r}(1:8), 'Planned:')
        
        % Store Planned date
        fields = strsplit(content{r});
        delta4.PlanDate = datetime([fields{2}, ' ', fields{3}, ' ', ...
            fields{4}], 'InputFormat', 'M/d/yyyy h:m a');
        
        % Store user, if present
        if length(fields) > 4
            delta4.PlanUser = fields{5};
        end
        
        break;
    else
        r = r + 1;
    end
end

% Loop through rows until measured date info is found
while r < length(content)
    
    % If row starts with 'Measured:'
    if length(content{r}) >= 9 && strcmp(content{r}(1:9), 'Measured:')
        
        % Store measured date
        fields = strsplit(content{r});
        delta4.MeasDate = datetime([fields{2}, ' ', fields{3}, ' ', ...
            fields{4}], 'InputFormat', 'M/d/yyyy h:m a');
        
        % Store user, if present
        if length(fields) > 4
            delta4.MeasUser = fields{5};
        end
        
        break;
    else
        r = r + 1;
    end
end

% Loop through rows until reviewed status info is found
while r < length(content)
    
    % If row starts with 'Accepted:' or 'Rejected:' or 'Failed:'
    if length(content{r}) >= 9 && (strcmp(content{r}(1:9), 'Accepted:') || ...
                strcmp(content{r}(1:9), 'Rejected:') || ...
                strcmp(content{r}(1:7), 'Failed:'))
        
        % Store measured date
        fields = strsplit(content{r});
        delta4.ReviewStatus = fields{1}(1:end-1);
        delta4.ReviewDate = datetime([fields{2}, ' ', fields{3}, ' ', ...
            fields{4}], 'InputFormat', 'M/d/yyyy h:m a');
        
        % Store user, if present
        if length(fields) > 4
            delta4.ReviewUser = fields{5};
        end
        
        % Otherwise, move to next row
        r = r + 1;
    
    % Otherwise, stop if row starts with 'Comments:'
    elseif length(content{r}) >= 9 && strcmp(content{r}(1:9), 'Comments:')
        
        content{r} = content{r}(10:end);
        break;
        
    % Otherwise, move to next row
    else
        r = r + 1;
    end
end

% Store Comments and look for treatment summary
delta4.Comments = cell(0);
while r < length(content)
    
    % If row is Treatment Summary
    if ~isempty(regexp(content{r}, 'Treatment Summary', 'ONCE'))
        break;
    else
        if ~isempty(content{r})
            delta4.Comments = vertcat(delta4.Comments, ...
                strtrim(content{r}));
        end
        r = r + 1;
    end
end

% Look for and store radiation device
while r < length(content)
    
    % If row starts with 'Radiation Device:'
    if length(content{r}) >= 17 && ...
            strcmp(content{r}(1:17), 'Radiation Device:')
        delta4.Machine = strtrim(content{r}(18:end));
        break;
    else
        r = r + 1;
    end
end

% Look for and store Temperature
while r < length(content)
    
    % If row starts with 'Temperature:'
    if length(content{r}) >= 12 && strcmp(content{r}(1:12), 'Temperature:')
        fields = regexp(content{r}(13:end), '([0-9\.]+)', 'tokens');
        
        if ~isempty(fields)
            delta4.Temperature = str2double(fields{1}(1));
        end
        break;
    else
        r = r + 1;
    end
end

% Look for and store dose Reference
while r < length(content)
    
    % If row starts with 'Reference:'
    if length(content{r}) >= 10 && ...
            strcmp(content{r}(1:10), 'Reference:')
        delta4.Reference = strtrim(content{r}(11:end));
        break;
    else
        r = r + 1;
    end
end

% Look for and store fraction statistics
while r < length(content)
    
    % If row starts with 'Fraction' or 'Composite'
    if length(content{r}) >= 9 && ...
            (strcmp(content{r}(1:8), 'Fraction') || ...
            strcmp(content{r}(1:9), 'Composite'))
        fields = regexp(content{r}(9:end), ['([0-9\.]+) +(c?Gy) +([0-9\.]', ...
            '+)% +([0-9\.]+)% +([0-9\.]+)% +(-?[0-9\.]+)%'], 'tokens');
        if strcmp(fields{1}(2), 'cGy')
            delta4.NormDose = str2double(fields{1}(1))/100;
        else
            delta4.NormDose = str2double(fields{1}(1));
        end
        delta4.AbsPassRate = str2double(fields{1}(3));
        delta4.DTAPassRate = str2double(fields{1}(4));
        delta4.GammaPassRate = str2double(fields{1}(5));
        delta4.MedianAbsDiff = str2double(fields{1}(6));
        r = r + 1;
        break;
    else
        r = r + 1;
    end
end

% Initialize Beams table
delta4.Beams = array2table(zeros(0,8), 'VariableNames', {'Name', 'Energy', ...
    'DailyCF', 'NormDose', 'AbsPassRate', 'DTAPassRate', 'GammaPassRate', ...
    'MedianAbsDiff'});
delta4.Beams.Properties.VariableUnits = {'', 'MV', '', 'Gy', '%', '%', ...
    '%', '%'};

% Look for and store beam statistics
while r < length(content)
    
    % If row is 'Histograms'
    if ~isempty(regexp(content{r}, 'Histograms', 'ONCE'))
        break
    else
        
        % If beam data exists
        if ~isempty(regexp(content{r}, ['([0-9]+) MV([0-9 MV,F\.]+)([0-9\.]+) ', ...
                '+([0-9\.]+) +(c?Gy) +([0-9\.]+)% +([0-9\.]+)% ', ...
                '+([0-9\.]+)% +(-?[0-9\.]+)%'], 'ONCE'))
            
            % Parse beam data
            fields = regexp(content{r}, ['([0-9]+) MV([0-9 MV,F\.]+)([0-9\.]+) ', ...
                '+([0-9\.]+) +(c?Gy) +([0-9\.]+)% +([0-9\.]+)% ', ...
                '+([0-9\.]+)% +(-?[0-9\.]+)%'], 'tokens');
            
            % Store number and FFF flag
            e = [fields{1}{1}, regexprep([fields{1}{2}, fields{1}{3}], ...
                '[0-9\., ]+', '')];
            
            % Store CF (correcting for bad parsing)
            cf = regexprep([fields{1}{2}, fields{1}{3}], ...
                '[^0-9\.]+', '');
            cf = str2double(cf);
            
            % Convert norm dose to Gy
            if strcmp(fields{1}{5}, 'cGy')
                n = str2double(fields{1}{4})/100;
            else
                n = str2double(fields{1}{4});
            end
            
            % Append beam row
            delta4.Beams = [delta4.Beams; [{regexp(strtrim(content{r}), ...
                '\S+', 'match', 'once')}, e, cf, n, ...
                str2double(fields{1}{6}), str2double(fields{1}{7}), ...
                str2double(fields{1}{8}), str2double(fields{1}{9})]];
        end
        r = r + 1;
    end
end

% Look for and store Gamma table
while r < length(content)
    
    % If row starts with 'Dose Deviation', skip ahead
    if startsWith(content{r}, 'Dose Deviation')
        break;
    end
    
    % If row contains 'Gamma Index Evaluations'
    if ~isempty(regexp(content{r}, 'Gamma\s+Index\s+Evaluations', 'ONCE'))
        
        % Initialize gamma table
        dta = [];
        data = [];
        r = r + 1;
        
        % Loop through Gamma table
        while r < length(content)
            
            % If table row exists
            if ~isempty(regexp(content{r}, ...
                    '([0-9\.]+ mm)\s+([0-9\.]+\s+)+', 'ONCE'))
                fields = regexp(content{r}, '([0-9\.]+)', 'tokens');
                dta(length(dta)+1) = str2double(fields{1}{1}); %#ok<*AGROW>
                data(size(data,1)+1,:) = zeros(1, length(fields)-1);
                for i = 1:size(data,2)
                    data(size(data,1),i) = ...
                        str2double(fields{1+i}{1});
                end
                
            elseif ~isempty(regexp(content{r}, '([0-9\.]+ %)+', 'ONCE'))
                fields = regexp(content{r}, '([0-9\.]+) %', 'tokens');
                delta4.GammaTable = array2table(data);
                delta4.GammaTable.Properties.UserData.GammaDTA = dta;
                delta4.GammaTable.Properties.UserData.GammaAbs = ...
                    zeros(1, length(fields));
                for i = 1:length(fields)
                    delta4.GammaTable.Properties.UserData.GammaAbs(i) = ...
                        str2double(fields{i}{1});
                end
                break;
                
            elseif ~isempty(regexp(content{r}, 'Dose\s+Deviation', 'ONCE'))
                break;
            end
            r = r + 1;
        end
        
        r = r + 1;
        break;
    else
        r = r + 1;
    end
    
    % Clear temporary variables
    clear dta data;
end

% Look for and store dose deviation parameters
while r < length(content)
    
    % If row starts with 'Dose Deviation'
    if length(content{r}) >= 14 && ...
            strcmp(content{r}(1:14), 'Dose Deviation')
        
        % Parse dose data
        fields = regexp(content{r}(15:end), ['([0-9\.]+)%[^0-9]+([0-9\.]', ...
            '+)%[^0-9]+([0-9\.]+)%[^0-9]+([0-9\.]+)%'], 'tokens');
        
        % Store dose data
        delta4.AbsRange(1) = str2double(fields{1}(1));
        delta4.AbsRange(2) = str2double(fields{1}(2));
        delta4.AbsPassLimit(1) = str2double(fields{1}(3));
        delta4.AbsPassLimit(2) = str2double(fields{1}(4));
        r = r + 1;
        break;
    else
        r = r + 1;
    end
end

% Look for and store DTA parameters
while r < length(content)
    
    % If row starts with 'Dist to Agreement'
    if length(content{r}) >= 17 && ...
            strcmp(content{r}(1:17), 'Dist to Agreement')
        
        % Parse DTA data
        fields = regexp(content{r}(18:end), ['([0-9\.]+)%[^0-9]+([0-9\.]', ...
            '+)%[^0-9]+([0-9\.]+)'], 'tokens');
        
        % Store DTA data
        delta4.DTARange(1) = str2double(fields{1}(1));
        delta4.DTARange(2) = inf;
        delta4.DTAPassLimit(1) = str2double(fields{1}(2));
        delta4.DTAPassLimit(2) = str2double(fields{1}(3));
        r = r + 1;
        break;
    else
        r = r + 1;
    end
end

% Look for and store Gamma Index parameters
while r < length(content)
    
    % If row starts with 'Gamma Index'
    if length(content{r}) >= 11 && ...
            strcmp(content{r}(1:11), 'Gamma Index')
        
        % Parse Gamma data
        fields = regexp(content{r}(12:end), ['([0-9\.]+)%[^0-9]+([0-9\.]', ...
            '+)%[^0-9]+([0-9\.]+)%[^0-9]+([0-9\.]+)[^0-9]+([0-9\.]+)', ...
            '%[^0-9]+([0-9\.]+)'], 'tokens');
        
        % Store Gamma data
        delta4.GammaRange(1) = str2double(fields{1}(1));
        delta4.GammaRange(2) = str2double(fields{1}(2));
        delta4.GammaAbs = str2double(fields{1}(3));
        delta4.GammaDTA = str2double(fields{1}(4));
        delta4.GammaPassLimit(1) = str2double(fields{1}(5));
        delta4.GammaPassLimit(2) = str2double(fields{1}(6));
        break;
    else
        r = r + 1;
    end
end


% Clear temporary variables
clear fields;

% Log finish
if exist('Event', 'file') == 2
    Event(sprintf('Delta4 report parsed successfully in %0.3f seconds', toc));
end