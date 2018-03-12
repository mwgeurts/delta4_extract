function delta4 = ParseD4Report(content)
% ParseD4Report reads a ScandiDos Delta4 report into a MATLAB structure.
% The function input argument can either be a string containing a file path
% and/or name corresponding to a report PDF file, or a cell array of text 
% data from that PDF. If a PDF file, this function will call XpdfText to
% extract the file contents (via the xpdf_tools submodule).
%
% This function also searches the comments field for two possible inputs:
% phantom name and delivered vs. expected MU (for TomoTherapy plans). The
% function looks for a comment line containing the word "delta4", and if
% present, will store the line under the "phantom" structure field. Second,
% if a line contains two numbers separated by / (example, 3453/3456), the
% first number will be stored as the cumulativeMU field, with the second
% as the expectedMU field.
%
% Upon successful completion, this function will return a structure
% containing the following fields:
%   title: string containing report title
%   name: string containing patient name
%   ID: string containing patient ID
%   clinic: cell array containing clinic name and address
%   plan: string containing plan name
%   planDate: planned datetime
%   planUser: string containing planned user (if present)
%   measDate: measured datetime
%   measUser: string containing planned user (if present)
%   comments: cell array of comments
%   phantom: string containing phantom name (if stored in comments),
%       otherwise 'Unknown'
%   cumulativeMU: double containing cumulative MU (if stored in comments)
%   expectedMU: double containing expected MU (if stored in comments)
%   machine: string containing radiation device
%   temperature: double containing temperature
%   reference: string containing reference (i.e. 'Planned Dose')
%   normDose: double containing normalization dose, in Gy
%   absPassRate: double containing absolute pass rate, as a percentage
%   dtaPassRate: double containing DTA pass rate, as a percentage
%   gammaPassRate: double containing Gamma pass rate, as a percentage
%   doseDev: double containing median dose deviation, as a percentage
%   beams: cell array of structures for each beam containing the following 
%       fields: name, dailyCF, normDose, absPassRate, dtaPassRate, 
%       gammaPassRate, and doseDev.
%   absRange: 2 element vector of dose deviation range, as percentages
%   absPassLimit: 2 element vector of dose deviation acceptance criteria,
%       as percentages (i.e. [90 3] means 90% within 3%)
%   dtaRange: 2 element vector of DTA range, in %/mm
%   dtaPassLimit: 2 element vector of DTA acceptance criteria, in % and mm
%   gammaRange: 2 element vector of Gamma range
%   gammaAbs: double containing Gamma absolute criterion as a percentage
%   gammaDta: double containing Gamma DTA criterion in mm
%   gammaPassLimit: 2 element vector of Gamma acceptance criteria (i.e.
%       [95 1] means 95% less than 1)
%   gammaTable: structure containing Gamma Index Evaluation table (if
%       present) with the following fields: dta, abs, passRate
%
% Author: Mark Geurts, mark.w.geurts@gmail.com
% Copyright (C) 2017 University of Wisconsin Board of Regents
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

% If content is a file name (with PDF extension), read file contents in
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

% If plan report is from version April 2016 or later
if length(content{7}) >= 7 && strcmp(content{7}(1:7), 'Clinic:')
    
    % Store title, patient name, and ID
    delta4.title = strtrim(content{1});
    delta4.name = strtrim(content{3});
    delta4.ID = strtrim(content{5});

    % Initialize row counter
    r = 6;
    
else 
    % Store title and patient name
    fields = strsplit(content{1}, '   ');
    delta4.title = strtrim(fields{1});
    delta4.name = strtrim(fields{2});
    for i = 3:length(fields)
        delta4.name = [delta4.name, ' ', strtrim(fields{i})];
    end

    % Store patient ID
    delta4.ID = strtrim(content{3});

    % Initialize row counter
    r = 4;
end

% Loop through rows until clinic info is found
while r < length(content)
    
    % If row starts with 'Clinic:'
    if length(content{r}) >= 7 && strcmp(content{r}(1:7), 'Clinic:')
        content{r} = content{r}(8:end);
        delta4.clinic = cell(0);
        break;
    else
        r = r + 1;
    end
end

% Store clinic contact info, followed by plan name
while r < length(content)
    
    % If row starts with 'Plan:'
    if length(content{r}) >= 5 && strcmp(content{r}(1:5), 'Plan:')
        delta4.plan = strtrim(content{r}(6:end));
        break;
    else
        if ~isempty(content{r})
            delta4.clinic = vertcat(delta4.clinic, strtrim(content{r}));
        end
        r = r + 1;
    end
end

% Loop through rows until planned date info is found
while r < length(content)
    
    % If row starts with 'Planned:'
    if length(content{r}) >= 8 && strcmp(content{r}(1:8), 'Planned:')
        
        % Store planned date
        fields = strsplit(content{r});
        delta4.planDate = datetime([fields{2}, ' ', fields{3}, ' ', ...
            fields{4}], 'InputFormat', 'M/d/yyyy h:m a');
        
        % Store user, if present
        if length(fields) > 4
            delta4.planUser = fields{5};
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
        delta4.measDate = datetime([fields{2}, ' ', fields{3}, ' ', ...
            fields{4}], 'InputFormat', 'M/d/yyyy h:m a');
        
        % Store user, if present
        if length(fields) >= 4
            delta4.measUser = fields{5};
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
        delta4.reviewStatus = fields{1}(1:end-1);
        delta4.reviewDate = datetime([fields{2}, ' ', fields{3}, ' ', ...
            fields{4}], 'InputFormat', 'M/d/yyyy h:m a');
        
        % Store user, if present
        if length(fields) > 4
            delta4.reviewUser = fields{5};
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

% Store comments and look for treatment summary
delta4.comments = cell(0);
while r < length(content)
    
    % If row is Treatment Summary
    if ~isempty(regexp(content{r}, 'Treatment Summary', 'ONCE'))
        break;
    else
        if ~isempty(content{r})
            delta4.comments = vertcat(delta4.comments, ...
                strtrim(content{r}));
        end
        r = r + 1;
    end
end

% Initialize unknown phantom
delta4.phantom = 'Unknown';

% Search for specific tags in comments
for i = 1:length(delta4.comments)
    
    % If phantom name is in the comments
    if contains(delta4.comments{i}, 'delta4', 'IgnoreCase',true)
        delta4.phantom = delta4.comments{i};
   
    % If cumulative/expected MU are in comments
    elseif regexp(delta4.comments{i}, '([0-9]+)[ ]?/[ ]?([0-9]+)') > 0
        fields = regexp(delta4.comments{i}, ...
            '([0-9]+)[ ]?/[ ]?([0-9]+)', 'tokens');
        delta4.cumulativeMU = str2double(fields{1}(1));
        delta4.expectedMU = str2double(fields{1}(2));
    end
end

% Look for and store radiation device
while r < length(content)
    
    % If row starts with 'Radiation Device:'
    if length(content{r}) >= 17 && ...
            strcmp(content{r}(1:17), 'Radiation Device:')
        delta4.machine = strtrim(content{r}(18:end));
        break;
    else
        r = r + 1;
    end
end

% Look for and store temperature
while r < length(content)
    
    % If row starts with 'Temperature:'
    if length(content{r}) >= 12 && strcmp(content{r}(1:12), 'Temperature:')
        fields = regexp(content{r}(13:end), '([0-9\.]+)', 'tokens');
        
        if ~isempty(fields)
            delta4.temperature = str2double(fields{1}(1));
        end
        break;
    else
        r = r + 1;
    end
end

% Look for and store dose reference
while r < length(content)
    
    % If row starts with 'Reference:'
    if length(content{r}) >= 10 && ...
            strcmp(content{r}(1:10), 'Reference:')
        delta4.reference = strtrim(content{r}(11:end));
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
            delta4.normDose = str2double(fields{1}(1))/100;
        else
            delta4.normDose = str2double(fields{1}(1));
        end
        delta4.absPassRate = str2double(fields{1}(3));
        delta4.dtaPassRate = str2double(fields{1}(4));
        delta4.gammaPassRate = str2double(fields{1}(5));
        delta4.doseDev = str2double(fields{1}(6));
        r = r + 1;
        break;
    else
        r = r + 1;
    end
end

% Initialize beams counter
b = 0;

% Look for and store beam statistics
while r < length(content)
    
    % If row is 'Histograms'
    if ~isempty(regexp(content{r}, 'Histograms', 'ONCE'))
        break
    else
        
        % If beam data exists
        if ~isempty(regexp(content{r}, ['([0-9\.]+) +([0-9\.]+) +(c?Gy) ', ...
                '+([0-9\.]+)% +([0-9\.]+)% +([0-9\.]+)% +(-?[0-9\.]+)%'], ...
                'ONCE'))
            
            % Increment beam counter
            b = b + 1;
            
            % Parse beam name
            delta4.beams{b,1}.name = regexp(strtrim(content{r}), '\S+', ...
                'match', 'once');
            
            % Parse beam data
            fields = regexp(content{r}, ['([0-9\.]+) +([0-9\.]+) +(c?Gy) ', ...
                '+([0-9\.]+)% +([0-9\.]+)% +([0-9\.]+)% +(-?[0-9\.]+)%'], ...
                'tokens');
            
            % Store beam data
            delta4.beams{b,1}.dailyCF = str2double(fields{1}(1));
            if strcmp(fields{1}(3), 'cGy')
                delta4.beams{b,1}.normDose = str2double(fields{1}(1))/100;
            else
                delta4.beams{b,1}.normDose = str2double(fields{1}(1));
            end
            delta4.beams{b,1}.absPassRate = str2double(fields{1}(4));
            delta4.beams{b,1}.dtaPassRate = str2double(fields{1}(5));
            delta4.beams{b,1}.gammaPassRate = str2double(fields{1}(6));
            delta4.beams{b,1}.doseDev = str2double(fields{1}(7));
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
        delta4.gammaTable.dta = [];
        delta4.gammaTable.passRate = [];
        r = r + 1;
        
        % Loop through Gamma table
        while r < length(content)
            
            % If table row exists
            if ~isempty(regexp(content{r}, ...
                    '([0-9\.]+ mm)\s+([0-9\.]+\s+)+', 'ONCE'))
                fields = regexp(content{r}, '([0-9\.]+)', 'tokens');
                delta4.gammaTable.dta(length(delta4.gammaTable.dta)+1) = ...
                    str2double(fields{1}{1}); %#ok<*AGROW>
                delta4.gammaTable.passRate(size(delta4.gammaTable.passRate,1)+1,:) = ...
                    zeros(1, length(fields)-1);
                for i = 1:size(delta4.gammaTable.passRate,2)
                    delta4.gammaTable.passRate(size(...
                        delta4.gammaTable.passRate,1),i) = ...
                        str2double(fields{1+i}{1});
                end
                
            elseif ~isempty(regexp(content{r}, '([0-9\.]+ %)+', 'ONCE'))
                fields = regexp(content{r}, '([0-9\.]+) %', 'tokens');
                delta4.gammaTable.abs = zeros(1, length(fields));
                for i = 1:length(fields)
                    delta4.gammaTable.abs(i) = str2double(fields{i}{1});
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
        delta4.absRange(1) = str2double(fields{1}(1));
        delta4.absRange(2) = str2double(fields{1}(2));
        delta4.absPassLimit(1) = str2double(fields{1}(3));
        delta4.absPassLimit(2) = str2double(fields{1}(4));
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
        delta4.dtaRange(1) = str2double(fields{1}(1));
        delta4.dtaRange(2) = inf;
        delta4.dtaPassLimit(1) = str2double(fields{1}(2));
        delta4.dtaPassLimit(2) = str2double(fields{1}(3));
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
        delta4.gammaRange(1) = str2double(fields{1}(1));
        delta4.gammaRange(2) = str2double(fields{1}(2));
        delta4.gammaAbs = str2double(fields{1}(3));
        delta4.gammaDta = str2double(fields{1}(4));
        delta4.gammaPassLimit(1) = str2double(fields{1}(5));
        delta4.gammaPassLimit(2) = str2double(fields{1}(6));
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