function delta4 = ParseD4Tables(content, varargin)
% ParseD4Tables reads a ScandiDos Delta4 tabular export into a structure.
% Either a file name or a cell array of exported data can be passed to this
% function. If a file is provided, it can be a Microsoft Excel spreadsheet
% (where each beam is saved to a different sheet), a CSV file, or a
% white-space delimited text file (where each beam is appended after
% another).
%
% Note, if the input data contains both plan and beam-level data, they will
% both be stored. If only beam level data is stored, however, the tool will 
% sum the individual beam values and store as the plan value only if
% absolute data is stored.
%
% Upon successful completion, the function will return a structure
% containing the following fields:
%   Name: string containing patient name
%   ID: string containing patient ID
%   Plan: string containing plan name
%   Data: n x 5 matrix of plan level data (or sum of beam data if absolute 
%       dose), where column 1 is distance from phantom center, column 2 is
%       the IEC X position, column 3 is IEC Y, column 4 is IEC Z, and
%       column 5 is the diode data. If the data is stored in cGy, it is
%       converted to Gy.
%   Value: string describing the stored data (i.e. Absolute dose in [Gy])
%   Beams: cell array of structures containing the following fields: name,
%       value, data
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

% Set default options
opt.CenterPhantom = false;

% Parse provided options
for i = 2:2:nargin
    opt.(varargin{i-1}) = varargin{i};
end

% If content is a file name with Excel extension
if ~iscell(content) && (endsWith(content, '.xlsx', 'IgnoreCase', true) || ...
        endsWith(content, '.xls', 'IgnoreCase', true))
    
    % Get file information
    [status, sheets] = xlsfinfo(content);
    
    % If file could not be read, throw an error
    if strcmp(status, '')
        if exist('Event', 'file') == 2
            Event(sheets, 'ERROR');
        else
            error(sheets);
        end
    end
    
    % Otherwise, loop through each sheet
    for i = 1:length(sheets)
        
        % Read in sheet contents
        [num, txt] = xlsread(content, sheets{i}, '', 'basic');
        
        % Store header data
        for j = 1:size(txt,1)
            
            % Store patient name and ID
            if startsWith(txt{j,1}, 'Patient:')
                fields = regexp(strtrim(txt{j,1}(9:end)), ...
                    '^(\S+)\s+(.+)$', 'tokens');
                delta4.ID = fields{1}{1};
                delta4.Name = fields{1}{2};
                
            % Store plan name
            elseif startsWith(txt{j,1}, 'Plan:')
                delta4.Plan = strtrim(txt{j,1}(6:end));
                
            % Store beam name
            elseif startsWith(txt{j,1}, 'Beam:')
                delta4.Beams{i}.Name = strtrim(txt{j,1}(6:end));
            
            % Store value specifier
            elseif startsWith(txt{j,1}, 'Intensity values:')
                delta4.Beams{i}.Value = strtrim(txt{j,1}(18:end));
                break;
            end
        end
        
        % If a beam name was not present, assume this is the plan data
        if ~isfield(delta4.Beams{i}, 'Name')
            delta4.Data = [];
            for j = 5:size(num,1)
                for k = 2:size(num,2)
                    if ~isnan(num(j,k))
                        delta4.Data(size(delta4.Data,1)+1,:) = ...
                            [num(1,k) num(2,k) num(1,j) num(3,k) num(j,k)];
                    end
                end
            end
            
        % Otherwise, store as beam data    
        else
            delta4.Beams{i}.Data = [];
            for j = 5:size(num,1)
                for k = 2:size(num,2)
                    if ~isnan(num(j,k))
                        delta4.Beams{i}.Data(size(delta4.Beams{i}.Data,1)+1,:) = ...
                            [num(1,k) num(2,k) num(j,1) num(3,k) num(j,k)];
                    end
                end
            end
        end
        
        % Scale data to Gy
        if contains(delta4.Beams{i}.Value, '[cGy]', 'IgnoreCase', true)
            delta4.Beams{i}.Data(:,5) = delta4.Beams{i}.Data(:,5)/100;
        end
    end
    
    % Clear temporary variables
    clear i j k status sheets num txt;
   
% Otherwise, if content is a file name of a text file
elseif ~iscell(content) && (endsWith(content, '.txt', 'IgnoreCase',true) ...
        || endsWith(content, '.csv', 'IgnoreCase',true))
    
    % Open file handle and read in contents
    fid = fopen(content, 'r');
    
    % Store file contents as cell array
    content = cell(0);
    while ~feof(fid)
        content{length(content)+1} = fgetl(fid);
    end
    fclose(fid);
    
    % Clear temporary variables
    clear fid;
    
% Otherwise, if content is some unknown file extension
else
   if exist('Event', 'file') == 2
        Event(['Unknown format: ', content], 'ERROR');
   else
        error(['Unknown format: ', content]);
   end
end

% If text contents exist
if iscell(content)
    
    % Initialize return variable
    delta4.Beams = cell(0);

    % Loop through content cell array, parsing
    i = 1;
    while i <= length(content)

        % Store patient name and ID
        if startsWith(content{i}, 'Patient:')
            fields = regexp(strtrim(content{i}(9:end)), ...
                '^(\S+)\s+(.+)$', 'tokens');
            delta4.ID = fields{1}{1};
            delta4.Name = fields{1}{2};
            b = false;

        % Store plan name
        elseif startsWith(content{i}, 'Plan:')
            delta4.Plan = strtrim(content{i}(6:end));

        % Store beam name
        elseif startsWith(content{i}, 'Beam:')
            delta4.Beams{length(delta4.Beams)+1}.Name = ...
                strtrim(content{i}(6:end));
            b = true;

        % Store value specifier
        elseif startsWith(content{i}, 'Intensity values:')
            delta4.Beams{length(delta4.Beams)}.Value = ...
                strtrim(content{i}(18:end));

        % Store distances
        elseif startsWith(content{i}, 'Distance')
            d = cellfun(@str2double, strsplit(strtrim(content{i}(9:end)), ...
                {'\s', ','}, 'CollapseDelimiters', false, 'DelimiterType', ...
                'RegularExpression'));

        % Store IEC X
        elseif startsWith(content{i}, 'X (iec-left)')
            x = cellfun(@str2double, strsplit(strtrim(content{i}(13:end)), ...
                {'\s', ','}, 'CollapseDelimiters', false, 'DelimiterType', ...
                'RegularExpression'));

        % Store IEC Z
        elseif startsWith(content{i}, 'Z (iec-up)')
            z = cellfun(@str2double, strsplit(strtrim(content{i}(11:end)), ...
                {'\s', ','}, 'CollapseDelimiters', false, 'DelimiterType', ...
                'RegularExpression'));

        % Store data
        elseif startsWith(content{i}, 'Y (iec-head)')
            i = i + 1;

            % Loop through data lines
            data = [];
            while i <= length(content)
                if regexp(content{i}, '^[0-9]')
                    row = cellfun(@str2double, strsplit(strtrim(content{i}), ...
                        {'\s', ','}, 'CollapseDelimiters', false, ...
                        'DelimiterType', 'RegularExpression'));
                    for j = 1:length(row)-1
                        if ~isnan(row(j+1))
                            data(size(data,1)+1,:) = ...
                                [d(j) x(j) row(1) z(j) row(j+1)]; %#ok<*AGROW>
                        end
                    end
                elseif startsWith(content{i}, 'Patient:')
                    break;
                end
                i = i + 1;
            end

            % Scale data to Gy
            if contains(delta4.Beams{length(delta4.Beams)}.Value, '[cGy]', ...
                    'IgnoreCase', true)
                data(:,5) = data(:,5)/100;
            end

            % Store beam or plan data 
            if b
                delta4.Beams{length(delta4.Beams)}.Data = data;
            else
                delta4.Data = data;
            end
            i = i - 1;
        end

        % Increment counter
        i = i + 1;
    end
end

% Remove empty beam names
for i = 1:length(delta4.Beams)
    if ~isfield(delta4.Beams{i}, 'Name')
        if isfield(delta4.Beams{i}, 'Value')
            delta4.Value = delta4.Beams{i}.Value;
        end
        delta4.Beams{i} = []; 
    end
end
delta4.Beams = delta4.Beams(~cellfun(@isempty, delta4.Beams));

% If a plan data is not present, compute it
if ~isfield(delta4, 'Data') && length(delta4.Beams) >= 1 ...
        && isfield(delta4.Beams{1}, 'Data') && ...
        isfield(delta4.Beams{1}, 'Value') && contains(delta4.Beams{1}.Value, ...
        'Absolute dose', 'IgnoreCase', true)
    delta4.Data = delta4.Beams{1}.Data;
    delta4.Value = delta4.Beams{1}.Value;
    for i = 2:length(delta4.Beams)
        if isfield(delta4.Beams{i}, 'Data')
            delta4.Data(:,5) = delta4.Data(:,5) + delta4.Beams{i}.Data(:,5);
        end
    end
end

% If center flag is enabled, center the phantom axially
if opt.CenterPhantom == 1
    
    % Center plan data
    if isfield(delta4, 'data')
        delta4.Data(:,2) = delta4.Data(:,2) - ...
            (max(delta4.Data(:,2)) + min(delta4.Data(:,2)))/2;
        delta4.Data(:,4) = delta4.Data(:,4) - ...
            (max(delta4.Data(:,4)) + min(delta4.Data(:,4)))/2;
    end
    
    % Center beams
    for i = 1:length(delta4.Beams)
        if isfield(delta4.Beams{i}, 'data')
            delta4.Beams{i}.Data(:,2) = delta4.Beams{i}.Data(:,2) - ...
                (max(delta4.Beams{i}.Data(:,2)) + ...
                min(delta4.Beams{i}.Data(:,2)))/2;
            delta4.Beams{i}.Data(:,4) = delta4.Beams{i}.Data(:,4) - ...
                (max(delta4.Beams{i}.Data(:,4)) + ...
                min(delta4.Beams{i}.Data(:,4)))/2;
        end
    end
end