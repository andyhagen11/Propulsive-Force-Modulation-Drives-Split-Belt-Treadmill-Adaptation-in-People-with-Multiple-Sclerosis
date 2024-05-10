%% Split-belt Joint Adaptation Analysis
% Originating authors: Isaac Bast and Andrew Hagen
% Last Updated: 12/24/2023

% This script uses data from a structure to calculate ROM (or range) and Peak adaptation time series (by gait cycle) for all joints and forces in all phases of a split-belt treadmill adaptation paradigm
%   see Vicon_Data_Grabber in repository for structure organization and creation using Vicon Nexus outputs
% Creates table of timeseries data for all participants
% Makes plots of each joint/force variable over all phases for SelectedParticipants

%% Script Initialization

% Make sure you have GaitCycle_Data_Interpolated.mat loaded. This will take a bit so its best to do this manually rather than have it be loaded up each time the script is run

% Don't accidentally clear the mat structure 
    clearvars -except GaitCycle_Data_Interpolated

Set Path 
    RootPath = 'R:\SBMS Joint Coordination\Joint Data\';
    cd(RootPath)

% Define the list of specific participants from the structure to use in the analysis.
SelectedParticipants = {'SBMS_004','SBMS_010','SBMS_011','SBMS_014','SBMS_026','SBMS_027','SBMS_040','SBMS_012','SBMS_018','SBMS_021','SBMS_022','SBMS_028','SBMS_030','SBMS_033', 'SBMS_006','SBMS_009','SBMS_013','SBMS_016','SBMS_019','SBMS_023','SBMS_024','SBMS_031','SBMS_037','SBMS_002','SBMS_003','SBMS_007','SBMS_017','SBMS_029','SBMS_034','SBMS_035','SBMS_036','SBMS_041'};

ParticipantNames = cell(1, numel(SelectedParticipants));

% Initialize a cell array to store the results tables for each participant
ROMCellArray = cell(12,numel(SelectedParticipants));
PeaksCellArray = cell(12,numel(SelectedParticipants));
% Only 8 rows if using Forces
% ROMCellArray = cell(4,numel(SelectedParticipants));
% PeaksCellArray = cell(8,numel(SelectedParticipants));


% Iterate through the selected participants
for i = 1:numel(SelectedParticipants)
    participant = SelectedParticipants{i};

%% Create List of Expressions to Iterate, Execute
Subjects = fieldnames(GaitCycle_Data_Interpolated);
Periods = fieldnames(GaitCycle_Data_Interpolated.SBMS_001); 
Limbs = fieldnames(GaitCycle_Data_Interpolated.SBMS_001.Baseline);
Angles = fieldnames(GaitCycle_Data_Interpolated.SBMS_001.Baseline.More_Affected);
Phases = fieldnames(GaitCycle_Data_Interpolated.SBMS_001.Baseline.More_Affected.Hip_Angles);
Planes = fieldnames(GaitCycle_Data_Interpolated.SBMS_001.Baseline.More_Affected.Hip_Angles.Stance);

    % Create expressions for the current participant
    expr = ['GaitCycle_Data_Interpolated.', participant, '.'];


    % Nested loops to iterate through periods, limbs, angles, phases, and planes
    % Look at Vicon_Data_Grabber in repository for how this structure is organized
    h = 1;
    for j = 1:(numel(Periods) - 1) % Excludes Adaptation 2
        for k = 1:numel(Limbs) % 2 Limbs
            for l = 7 % 6 Angles, excludes Force data - ADUJUST THIS BASED ON IF YOU WANT JOINT ANKLES OR FORCE 
                for m = 1 % Phases of gait: stance, swing
                    for n = 1 % Only looks at X-axis (sagittal plane angles)

                        tempName = strcat(expr, Periods{j}, '.', Limbs{k}, '.', Angles{l}, '.', Phases{m}, '.', Planes{n});
                        Expression(h, 1) = string(tempName);
                        h = h + 1;
                    end
                end
            end
        end
    end

%% Create AllCycles Structure to store all gait cycles of variables found in Expression 
    % Also create an interpolated version so we can compare time series across subjects with different numbers of gait cycles
% Initialize AllCycles and AllCyclesInterp structures
AllCycles = struct();
AllCyclesInterp = struct();
FieldNames = arrayfun(@(s) replace(s, '.', '_'), Expression, 'UniformOutput', false);

for q = 1:numel(FieldNames)
    TSFieldNames{q} = FieldNames{q}{1}(38:end);
end

CalcIDX = numel(Expression);

% Initialize a cell array to store the interpolated data for each field
InterpolatedData = cell(CalcIDX, 1);

for p = 1:CalcIDX
    AllCycles.(TSFieldNames{p}) = (eval(Expression(p)));
    
    % Define X as a common set of data points
    X = 1:(size(eval(Expression(p), 1))-5); % Crop first five gait cycles to allow for treadmill to speed up

    if contains(Expression(p), 'Baseline') 
        InterpolatedDataSize = 100;
    elseif contains(Expression(p), 'Adaptation1') 
        InterpolatedDataSize = 500;
    elseif contains(Expression(p), 'Catch')
        InterpolatedDataSize = 50;
    end
    
    % Loop through each column and interpolate each one
    for Col = 1:size(eval(Expression{p}), 2)
        ColData = eval(Expression{p});
        ColData = ColData(:, Col); %Interpolate each row iteratively
        ColData = ColData(6:end,:); % Crop first five gait cycles to allow for treadmill to speed up
        InterpolatedData{p}(:, Col) = interp1(X, ColData, linspace(1, numel(X), InterpolatedDataSize), 'linear', 'extrap');
    end
end

% Store the interpolated data in the AllCyclesInterp structure
for p = 1:CalcIDX
    AllCyclesInterp.(TSFieldNames{p}) = InterpolatedData{p};
end


%% Horizontally concatenate stance and swing phases of all variables

FieldNames = fieldnames(AllCyclesInterp);

% Initialize a structure to store concatenated data
ConcatenatedData = struct();

% Iterate over the field names
for q = 1:length(FieldNames)
    FieldName = FieldNames{q};
    
    % Split the field name into parts using underscores
    Parts = strsplit(FieldName, '_');
    
    % Construct the common part of the field name (excluding 'Stance' or 'Swing' and the last part)
    CommonPart = strjoin(Parts(1:end-3), '_');
    
    % Check if the common part exists in the concatenatedData structure
    if ~isfield(ConcatenatedData, CommonPart)
        ConcatenatedData.(CommonPart) = [];
    end
    
    % Check if it's 'Stance_X' or 'Swing_X' and concatenate accordingly
    if strcmp(Parts{end-1}, 'Stance') || strcmp(Parts{end-1}, 'Swing')
        ConcatenatedData.(CommonPart) = [ConcatenatedData.(CommonPart), AllCyclesInterp.(FieldName)];
    end
end

% Reinterpolate across each row
   FieldNames = fieldnames (ConcatenatedData) ;
    for r = 1:size(FieldNames)
        % Find length of row
        X = 1:(size(ConcatenatedData.(FieldNames{r}), 2));

        % Loop through each row and interpolate each one
        for Row = 1: size(ConcatenatedData.(FieldNames{r}), 1)
            RowData = ConcatenatedData.(FieldNames{r});
            RowData = RowData(Row, :); % Interpolate each row iteratively
            ConcatInterpolatedData{r}(Row, :) = interp1(X, RowData, linspace(1, numel(X), 200), 'linear', 'extrap');
        end
    end

% Store the interpolated data by rewriting into the ConcatenatedData structure
    for r = 1:size(FieldNames)
        ConcatenatedData.(FieldNames{r}) = ConcatInterpolatedData{r};
    end

%% Then vertically concatenate all phases together for each variable 

    % Get the field names from concatenatedData
    FieldNames = fieldnames(ConcatenatedData);
    
    % Initialize a new structure to store vertically concatenated data
    VerticalConcatenatedData = struct();
    
    % Iterate over the field names in concatenatedData
    for r = 1:numel(FieldNames)
        FieldName = FieldNames{r};
        
        % Split the field name into parts using underscores
        Parts = strsplit(FieldName, '_');
        
        % Extract the suffix part (e.g., "More_Affected_Hip")
        suffix = strjoin(Parts(2:end), '_');
        
        % Check if the suffix exists in verticalConcatenatedData
        if isfield(VerticalConcatenatedData, suffix)
            % Vertically concatenate the matrix to the existing field
            VerticalConcatenatedData.(suffix) = [VerticalConcatenatedData.(suffix); ConcatenatedData.(FieldName)];
%                 VerticalConcatenatedData.(suffix) = [VerticalConcatenatedData.(suffix); AllCyclesInterp.(FieldName)]; % - use this line if not horizontally concatenating above
        else
            % Create a new field with the matrix
            VerticalConcatenatedData.(suffix) = ConcatenatedData.(FieldName);
%                 VerticalConcatenatedData.(suffix) = AllCyclesInterp.(FieldName); %- use this line if not horizontally concatenating above
        end
    end

%% Create ROM and Peak adaptation time series cell array
% Can only do joint angles, forces (or impluse) one at a time in this version. Will need to rerun and adjust variable 'l' above to switch output.
  
    FieldNames = fieldnames(VerticalConcatenatedData);
    
    for s = 1:numel(FieldNames)
        FieldName = FieldNames{s};

%        % IF USING JOINT ANGLES UNCOMMENT
        ROMCellArray{s,i}  = range(VerticalConcatenatedData.(FieldName), 2); % Calculate ROM for each gait cycle
        if contains((FieldName), 'Knee') % Calculate max only for knee but min for everything else
            PeaksCellArray{s,i} = max(VerticalConcatenatedData.(FieldName),[],2);
        elseif contains((FieldName), 'Ankle') % Calculate max and min for dorsiflexion and plantar flexion
            if contains((FieldName), 'More_Affected')
            PeaksCellArray{(13),i} = max(VerticalConcatenatedData.(FieldName)(:,20:100),[],2); % Dorsiflexion
            PeaksCellArray{s,i} = min(VerticalConcatenatedData.(FieldName)(:,100:120),[],2); % Plantarflexion
            elseif contains((FieldName), 'Less_Affected')
            PeaksCellArray{(14),i} = max(VerticalConcatenatedData.(FieldName)(:,20:100),[],2); % Dorsiflexion
            PeaksCellArray{s,i} = min(VerticalConcatenatedData.(FieldName)(:,100:120),[],2); % Plantarflexion
            end
        else
            PeaksCellArray{s,i} = min(VerticalConcatenatedData.(FieldName),[],2);
        end
% 
         % IF USING FORCES UNCOMMENT
%         ROMCellArray{s,i}  = range(VerticalConcatenatedData.(FieldName), 2); % Calculate amplitude for each force profile
%         if contains((FieldName), 'X') % Calculate max and min for anterior GRF
%             if contains((FieldName), 'More_Affected')
%             PeaksCellArray{5,i} = max(VerticalConcatenatedData.(FieldName),[],2); % Braking
%             PeaksCellArray{(s),i} = min(VerticalConcatenatedData.(FieldName),[],2); % Propulsion
%             elseif contains((FieldName), 'Less_Affected')
%             PeaksCellArray{6,i} = max(VerticalConcatenatedData.(FieldName),[],2); % Braking
%             PeaksCellArray{(s),i} = min(VerticalConcatenatedData.(FieldName),[],2); % Propulsion
%             end
%         elseif contains((FieldName), 'Z')
%             if contains((FieldName), 'More_Affected')
%             PeaksCellArray{(s),i} = min((VerticalConcatenatedData.(FieldName)(:,1:(end/2))),[],2); % First Peak
%             PeaksCellArray{(7),i} = min((VerticalConcatenatedData.(FieldName)(:,(end/2):end)),[],2); % Second Peak 
%             elseif contains((FieldName), 'Less_Affected')
%             PeaksCellArray{(s),i} = min((VerticalConcatenatedData.(FieldName)(:,1:(end/2))),[],2); % First Peak
%             PeaksCellArray{(8),i} = min((VerticalConcatenatedData.(FieldName)(:, (end/2):end)),[],2); % Second Peak 
%             end
%         end

%        UNCOMMENT IF WANT TO CALCULATE FORCE IMPULSE
%        Make cell array with both braking and prop - used ROM for easy fit into script
%          PropArea = [];
%          BrakeArea = [];
%            x = VerticalConcatenatedData.(FieldName)';
%             for a = x  % loop through each column (gait cycle)
%                 mx_x = max(a(1:50));
%                 mx_xloc = find(a == mx_x);
%                 
%                 mn_x = min(a(mx_xloc:100));
%                 mn_xloc = find(a == mn_x);
% 
%                 if sum(a == 0) > 20 || ~any(conv(double(a(mx_xloc:end) < 0), ones(1, 5), 'valid') == 5)  % Check if 'a' contains more than 20 zeros or or has less than 5 consecutive negative values after max(meaning missing/invalid gait cycle)
%                 PropArea = [PropArea, NaN];
%                 BrakeArea = [BrakeArea, NaN];
%                 else
%                 
%                     for b = mx_xloc:100
%                         if all(a(b:b+5) < 0)
%                             ZeroSpot = b;
%                             break
%                         end
%                     end                 
%                     tp = trapz(a(ZeroSpot:100));
%                     PropArea = [PropArea,tp];
%                     tb = trapz(a(1:ZeroSpot));
%                     BrakeArea = [BrakeArea,tb];
%                 end
%             end
% 
%         if contains((FieldName), 'More_Affected')
%         ROMCellArray{1,i} = PropArea';
%         ROMCellArray{2,i} = BrakeArea';
%         elseif contains((FieldName), 'Less_Affected')
%         ROMCellArray{3,i} = PropArea';
%         ROMCellArray{4,i} = BrakeArea';
%         end
    end
end

%% Export combined individual data for each participant

%  Rename FieldNames because we added Ankle Dorsiflexion to peaks cells
   AdditionalFields = {'More_Affected_Ankle_Dorsi'; 'Less_Affected_Ankle_Dorsi'};
% AdditionalFields = {'More_Affected_Forces_X_Braking'; 'Less_Affected_Forces_X_Braking'; 'More_Affected_Forces_Z_Late'; 'Less_Affected_Forces_Z_Late'};
   FieldNamesPeaks = [FieldNames; AdditionalFields];
% FieldNames = {'MA_Propulsion_Impulse', 'MA_Braking_Impulse', 'LA_Propulsion_Impulse', 'LA_Braking_Impulse'};

ExportIndivData = inputdlg('Save and Export to Individual Data xlsx? Y/N: ');
if strcmp(ExportIndivData,'Y')
    ExportDir = uigetdir('Select a directory to save the XLSX files');
        
    TransposedROMData = cell2table(ROMCellArray', 'VariableNames', FieldNames);
    FullFilePath1 = fullfile(ExportDir, 'Individual_ROM_Table.xlsx');
    writetable(TransposedROMData, FullFilePath1, 'Sheet', 'Combined_Data');
    for t = 1:numel(FieldNames) % Create a sheet for each variable
        TempROMColData = TransposedROMData(:,t);  % Extract the i-th column
        SheetName = [char(FieldNames(t))];  % Name for the sheet
        writetable(TempROMColData, FullFilePath1, 'Sheet', SheetName,'WriteMode', 'append');
    end

    TransposedPeaksData = cell2table(PeaksCellArray', 'VariableNames', FieldNamesPeaks);
    FullFilePath2 = fullfile(ExportDir, 'Individual_Peaks_Table.xlsx');
    writetable(TransposedPeaksData, FullFilePath2, 'Sheet', 'Combined_Data');
    for t = 1:numel(FieldNamesPeaks) % Create a sheet for each variable
        TempPeaksColData = TransposedPeaksData(:,t); 
        SheetName = [char(FieldNamesPeaks(t))];  % Name for the sheet
        writetable(TempPeaksColData, FullFilePath2, 'Sheet', SheetName,'WriteMode', 'append');
    end

end 

%% Make Group Data

% Initialize cell arrays for average ROM and Peaks then calculate groupwise means for each gait cycle number
AverageROMCellArray = cell(numel(FieldNames),1);
%AveragePeaksCellArray = cell(numel(FieldNamesPeaks),1);
SE_ROMCellArray = cell(numel(FieldNames), 1);
%SE_PeaksCellArray = cell(numel(FieldNamesPeaks), 1);

for u = 1:numel(FieldNames)
    ROMData = ROMCellArray(u, :);
    AverageROM = mean(cat(2, ROMData{:}), 2,'omitnan');
    StdErrorROM = std(cat(2, ROMData{:}), 0, 2, 'omitnan') / sqrt(numel(ROMData));
    AverageROMCellArray{u} = AverageROM;
    SE_ROMCellArray{u} = StdErrorROM;
end


for v = 1:numel(FieldNamesPeaks)
    PeaksData = PeaksCellArray(v, :);
    AveragePeaks = mean(cat(2, PeaksData{:}), 2);
    StdErrorPeaks = std(cat(2, PeaksData{:}), 0, 2) / sqrt(numel(PeaksData));
    AveragePeaksCellArray{v} = AveragePeaks;
    SE_PeaksCellArray{v} = StdErrorPeaks;
end

% Convert the cell arrays to tables
    ROMTable = array2table(cell2mat(AverageROMCellArray'), 'VariableNames', FieldNames);
    PeaksTable = array2table(cell2mat(AveragePeaksCellArray'), 'VariableNames', FieldNamesPeaks);
    SE_ROMTable = array2table(cell2mat(SE_ROMCellArray'), 'VariableNames', FieldNames);
    SE_PeaksTable = array2table(cell2mat(SE_PeaksCellArray'), 'VariableNames', FieldNamesPeaks);

% Create rolling mean of tables
    RM_ROMTable = array2table(movmean(ROMTable{:,:},10), 'VariableNames', FieldNames);
    RM_PeaksTable = array2table(movmean(PeaksTable{:,:},10), 'VariableNames', FieldNamesPeaks);


% Smooth SE curves
    SE_ROM_SmoothedUpper = cell(size(SE_ROMCellArray)); % Smooth ROM SE Upper Bound
    SE_ROM_SmoothedLower = cell(size(SE_ROMCellArray)); % Smooth ROM SE Lower Bound

for w = 1:numel(SE_ROMCellArray) 
    x = AverageROMCellArray{w} + SE_ROMCellArray{(w)}; % Upper
    y = AverageROMCellArray{w} - SE_ROMCellArray{(w)}; % Lower
    % Smooth Upper Bound
    SmoothedDataUpper = fit((1:numel(x))', x, 'smoothingspline', 'SmoothingParam', 0.5);
    SmoothedValuesUpper = feval(SmoothedDataUpper, (1:numel(x))'); % Evaluate the fit over the range of values
    SE_ROM_SmoothedUpper{w} = SmoothedValuesUpper;
    % Smooth Lower Bound
    SmoothedDataLower = fit((1:numel(y))', y, 'smoothingspline', 'SmoothingParam', 0.5);
    SmoothedValuesLower = feval(SmoothedDataLower, (1:numel(y))'); % Evaluate the fit over the range of values
    SE_ROM_SmoothedLower{w} = SmoothedValuesLower;
end

SE_Peaks_SmoothedUpper = cell(size(SE_PeaksCellArray)); % Smooth Peaks SE Upper Bound
SE_Peaks_SmoothedLower = cell(size(SE_PeaksCellArray)); % Smooth Peaks SE Lower Bound

for w = 1:numel(SE_PeaksCellArray) 
    x = AveragePeaksCellArray{w} + SE_PeaksCellArray{w}; % Upper
    y = AveragePeaksCellArray{w} - SE_PeaksCellArray{w}; % Lower
    % Smooth Upper Bound
    SmoothedDataUpper = fit((1:numel(x))', x, 'smoothingspline', 'SmoothingParam', 0.5);
    SmoothedValuesUpper = feval(SmoothedDataUpper, (1:numel(x))'); % Evaluate the fit over the range of values
    SE_Peaks_SmoothedUpper{w} = SmoothedValuesUpper;
    % Smooth Lower Bound
    SmoothedDataLower = fit((1:numel(y))', y, 'smoothingspline', 'SmoothingParam', 0.5);
    SmoothedValuesLower = feval(SmoothedDataLower, (1:numel(y))'); % Evaluate the fit over the range of values
    SE_Peaks_SmoothedLower{w} = SmoothedValuesLower;
end

%% Export and Plot

ExportData = inputdlg('Save and Export to Group Data xlsx? Y/N: ');
 if strcmp(ExportData,'Y')
    ExportDir = uigetdir('Select a directory to save the XLSX files');
        % Create the full file paths for each table with XLSX file extension
        FullFilePath3 = fullfile(ExportDir, 'ROMTable.xlsx');
        FullFilePath4 = fullfile(ExportDir, 'RM_ROMTable.xlsx');
        FullFilePath5 = fullfile(ExportDir, 'PeaksTable.xlsx');
        FullFilePath6 = fullfile(ExportDir, 'RM_PeaksTable.xlsx');
    
        % Use writetable to export each table to the selected location with sheet names
        writetable(ROMTable, FullFilePath3, 'Sheet', 'ROMTable');
        writetable(RM_ROMTable, FullFilePath4, 'Sheet', 'RM_ROMTable');
        writetable(PeaksTable, FullFilePath5, 'Sheet', 'PeaksTable');
        writetable(RM_PeaksTable, FullFilePath6, 'Sheet', 'RM_PeaksTable');
    
        disp('Tables exported to the selected location as XLSX files.');
 end

%% Make plots of each joint data over all phases for selected participants
ExportPlot = inputdlg('Create plots? Y/N: ');
 if strcmp(ExportPlot,'Y')

    ROMPlot = figure(1);
    set(gcf, 'position', get(0, 'ScreenSize'))
    sgtitle('Adaptation Time Series for Impulse')
    set(sgtitle, 'FontWeight', 'bold');
    hold on

    Titles = { 'MA Hip ROM','LA Hip ROM', 'MA Knee ROM','LA Knee ROM','MA Ankle ROM','LA Ankle ROM'};
    Colors = {[0.259, 0.259, 0.0259],[0.376, 0.490, 0.545]};  % Colors for plotting
    AffectedType = [1, 2, 1, 2, 1, 2];  % 1 for More Affected, 2 for Less Affected
    Indices = [1, 7, 2, 8, 3, 9];  % Indices corresponding to ROMTable and RM_ROMTable

    GapSize = 30;  % Set the size of the gap
            legend('Baseline','Early Adapt','Late Adapt','Early Catch','Late Catch')
        set(legend,'Location', 'Best')
    
    for i = 1:numel(Indices)
        subplot(3, 2, i);
        hold on
    
        title(Titles{i},'FontWeight','bold','FontSize',12);
        subtitle('           Baseline                                                 Adaptation                                             Post-Adaptation')
        ylabel('Degrees');
        xlabel('Gait Cycle Number');
        xlim([1, 510]);  % Adjust the x-axis limit
        xticks([0,50,100,130,180,230,280,330,380,430,460,510]);  % Adjust the x-ticks limit
        xticklabels({'0', '50', '100', '0', '50', '100', '150', '', '450', '500', '0', '50'});
        ColorIndex = AffectedType(i);

        % Plot the first segment which is only Baseline
        plot(ROMTable{1:100, Indices(i)}, 'o', 'MarkerSize', 4, 'Color', Colors{ColorIndex});
        plot(RM_ROMTable{1:97, Indices(i)}, 'Color', Colors{ColorIndex}, 'LineWidth', 3);
        plot(1:100, (SE_ROM_SmoothedUpper{Indices(i)}(1:100)), 'LineWidth', 1, 'Color', [0.7, 0.7, 0.7]);
        plot(1:100, (SE_ROM_SmoothedLower{Indices(i)}(1:100)), 'LineWidth', 1, 'Color', [0.7, 0.7, 0.7]);
    
        % Add a gap by plotting NaN between Baseline and Adapt
        plot([100, 100] + [0, GapSize], ylim, 'Color', 'w');

        % Plot the second segment starting at Adapt and crop 300:500
        plot((101 + GapSize):(400 + GapSize), ROMTable{([101:300,501:600]), Indices(i)}, 'o', 'MarkerSize', 4, 'Color', Colors{ColorIndex});
        plot((103 + GapSize):(397 + GapSize), RM_ROMTable{[103:300,501:597], Indices(i)}, 'Color', Colors{ColorIndex}, 'LineWidth', 3);
        plot((101 + GapSize):(400 + GapSize), (SE_ROM_SmoothedUpper{Indices(i)}([101:300,501:600])), 'LineWidth', 1, 'Color', [0.7, 0.7, 0.7]);
        plot((101 + GapSize):(400 + GapSize), (SE_ROM_SmoothedLower{Indices(i)}([101:300,501:600])), 'LineWidth', 1, 'Color', [0.7, 0.7, 0.7]);
    
        % Add a gap by plotting NaN between Adapt and Catch
        plot([(400 + GapSize), (400 + GapSize)] + [0, GapSize], ylim, 'Color', 'w');
        
        % Plot the third segment which is Catch
        plot((401 + 2 * GapSize):(450+ 2 * GapSize), ROMTable{601:end, Indices(i)}, 'o', 'MarkerSize', 4, 'Color', Colors{ColorIndex});
        plot((405 + 2 * GapSize):(450+ 2 * GapSize), RM_ROMTable{605:end, Indices(i)}, 'Color', Colors{ColorIndex}, 'LineWidth', 3);
        plot((401 + 2 * GapSize):(450+ 2 * GapSize), (SE_ROM_SmoothedUpper{Indices(i)}(601:650)), 'LineWidth', 1, 'Color', [0.7, 0.7, 0.7]);
        plot((401 + 2 * GapSize):(450+ 2 * GapSize), (SE_ROM_SmoothedLower{Indices(i)}(601:650)), 'LineWidth', 1, 'Color', [0.7, 0.7, 0.7]);

        % Add vertical dotted line between Baseline and Adapt
        line([100+GapSize/2, 100+GapSize/2], ylim, 'Color', 'k', 'LineStyle', '--');
        % Add vertical dotted line between Adapt and Catch
        line([400+1.5*GapSize, 400+1.5*GapSize], ylim , 'Color', 'k', 'LineStyle', '--');

        % Extract y-limits for the subplot
        ylimits = ylim;
        
          % Add colored shaded bars with opacity based on subplot y-limits
        fill([86, 100, 100, 86], [ylimits(1), ylimits(1), ylimits(2), ylimits(2)], [0, 0.278, 0.671], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
        fill([131, 140, 140, 131], [ylimits(1), ylimits(1), ylimits(2), ylimits(2)], [0.875, 0.455, 0.318], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
        fill([415, 430, 430, 415], [ylimits(1), ylimits(1), ylimits(2), ylimits(2)], [0.431, 0.149, 0.055], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
        fill([461, 470, 470, 461], [ylimits(1), ylimits(1), ylimits(2), ylimits(2)], [0.010, 0.545, 0.341], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
        fill([495, 510, 510, 495], [ylimits(1), ylimits(1), ylimits(2), ylimits(2)], [0, 0.278, 0.165], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    end

    PeaksPlot = figure(2);
        set(gcf, 'position', get(0, 'ScreenSize'))
        sgtitle('Adaptation Time Series for Peak Angles of Each Joint')
        set(sgtitle, 'FontWeight', 'bold');
        hold on
        
        Titles = {'MA Peak Hip Extension','LA Peak Hip Extension','MA Peak Knee Flexion','LA Peak Knee Flexion','MA Peak Ankle Plantarflexion','LA Peak Ankle Plantarflexion','MA Peak Ankle Dorsiflexion','LA Peak Ankle Dorsiflexion'};    
        Colors = {[0.259, 0.259, 0.0259],[0.376, 0.490, 0.545]};  % Colors for plotting
        AffectedType = [1, 2, 1, 2, 1, 2, 1, 2];  % 1 for More Affected, 2 for Less Affected
        Indices = [1, 7, 2, 8, 3, 9, 13,14];  % Indices corresponding to PeaksTable and RM_PeaksTable
        GapSize = 30;  % Set the size of the gap
    
    for i = 1:numel(Indices)
        subplot(4, 2, i);
        hold on
    
        title(Titles{i},'FontWeight','bold','FontSize',12);
        subtitle('           Baseline                                                 Adaptation                                             Post-Adaptation')
        ylabel('Degrees');
        xlabel('Gait Cycle Number');
        xlim([1, 510]);  % Adjust the x-axis limit
        xticks([0,50,100,130,180,230,280,330,380,430,460,510]);  % Adjust the x-ticks limit
        xticklabels({'0', '50', '100', '0', '50', '100', '150', '', '450', '500', '0', '50'});
    
        % Choose color based on affected type
        ColorIndex = AffectedType(i);

        % Plot the first segment which is only Baseline
        plot(PeaksTable{1:100, Indices(i)}, 'o', 'MarkerSize', 4, 'Color', Colors{ColorIndex});
        plot(RM_PeaksTable{1:97, Indices(i)}, 'Color', Colors{ColorIndex}, 'LineWidth', 3);
        plot(1:100, (SE_Peaks_SmoothedUpper{Indices(i)}(1:100)), 'LineWidth', 1, 'Color', [0.7, 0.7, 0.7]);
        plot(1:100, (SE_Peaks_SmoothedLower{Indices(i)}(1:100)), 'LineWidth', 1, 'Color', [0.7, 0.7, 0.7]);
    
        % Add a gap by plotting NaN between Baseline and Adapt
        plot([100, 100] + [0, GapSize], ylim, 'Color', 'w');

        % Plot the second segment starting at Adapt and crop 300:500
        plot((101 + GapSize):(400 + GapSize), PeaksTable{([101:300,501:600]), Indices(i)}, 'o', 'MarkerSize', 4, 'Color', Colors{ColorIndex});
        plot((103 + GapSize):(397 + GapSize), RM_PeaksTable{[103:300,501:597], Indices(i)}, 'Color', Colors{ColorIndex}, 'LineWidth', 3);
        plot((101 + GapSize):(400 + GapSize), (SE_Peaks_SmoothedUpper{Indices(i)}([101:300,501:600])), 'LineWidth', 1, 'Color', [0.7, 0.7, 0.7]);
        plot((101 + GapSize):(400 + GapSize), (SE_Peaks_SmoothedLower{Indices(i)}([101:300,501:600])), 'LineWidth', 1, 'Color', [0.7, 0.7, 0.7]);
    
        % Add a gap by plotting NaN between Adapt and Catch
        plot([(400 + GapSize), (400 + GapSize)] + [0, GapSize], ylim, 'Color', 'w');
        
        % Plot the third segment which is Catch
        plot((401 + 2 * GapSize):(450+ 2 * GapSize), PeaksTable{601:end, Indices(i)}, 'o', 'MarkerSize', 4, 'Color', Colors{ColorIndex});
        plot((405 + 2 * GapSize):(447+ 2 * GapSize), RM_PeaksTable{605:end-3, Indices(i)}, 'Color', Colors{ColorIndex}, 'LineWidth', 3);
        plot((401 + 2 * GapSize):(450+ 2 * GapSize), (SE_Peaks_SmoothedUpper{Indices(i)}(601:650)), 'LineWidth', 1, 'Color', [0.7, 0.7, 0.7]);
        plot((401 + 2 * GapSize):(450+ 2 * GapSize), (SE_Peaks_SmoothedLower{Indices(i)}(601:650)), 'LineWidth', 1, 'Color', [0.7, 0.7, 0.7]);

        % Add vertical dotted line between Baseline and Adapt
        line([100+GapSize/2, 100+GapSize/2], ylim, 'Color', 'k', 'LineStyle', '--');
        % Add vertical dotted line between Adapt and Catch
        line([400+1.5*GapSize, 400+1.5*GapSize], ylim , 'Color', 'k', 'LineStyle', '--');

        % Extract y-limits for the subplot
        ylimits = ylim;
        
          % Add colored shaded bars with opacity based on subplot y-limits
        fill([86, 100, 100, 86], [ylimits(1), ylimits(1), ylimits(2), ylimits(2)], [0, 0.278, 0.671], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
        fill([131, 140, 140, 131], [ylimits(1), ylimits(1), ylimits(2), ylimits(2)], [0.875, 0.455, 0.318], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
        fill([415, 430, 430, 415], [ylimits(1), ylimits(1), ylimits(2), ylimits(2)], [0.431, 0.149, 0.055], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
        fill([461, 470, 470, 461], [ylimits(1), ylimits(1), ylimits(2), ylimits(2)], [0.010, 0.545, 0.341], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
        fill([495, 510, 510, 495], [ylimits(1), ylimits(1), ylimits(2), ylimits(2)], [0, 0.278, 0.165], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    end 
 end

Forces Plots
    ROMPlot = figure(1);
        set(gcf, 'position', get(0, 'ScreenSize'))
        sgtitle('Adaptation Time Series for Force Amplitude')
        set(sgtitle, 'FontWeight', 'bold');
        hold on
        
        titles = {
            'More Affected Anterior-Posterior Force Amplitude Adaptation',
            'Less Affected Anterior-Posterior Force Amplitude Adaptation',
            'More Affected Vertical Force Amplitude Adaptation',
            'Less Affected Vertical Force Amplitude Adaptation'
        };
        
        colors = {'r', 'b'};  % Colors for plotting
        
        indices = [1, 3, 2, 4];  % Indices corresponding to ROMTable and RM_ROMTable
        
        for i = 1:4
            subplot(2, 2, i);
            hold on
            plot(ROMTable{:, indices(i)}, [colors{1} 'o'], 'MarkerSize', 4);
            plot(RM_ROMTable{:, indices(i)}, colors{1}, 'LineWidth', 3);
            title(titles{i});
            ylabel('Raw Force (N)');
            xlabel('Gait Cycle Number');
            xlim([1, 650]);
            xticks([92, 105, 592, 605, 642]);
            xticklabels({'Baseline', 'Early Adaptation', 'Late Adaptation', 'Early Post-Adaptation', 'Late Post-Adaptation'});
            
            % Add error bands
            plot(1:numel(SE_ROMCellArray{indices(i)}), (ROMTable{indices(i)} + SE_ROMCellArray{indices(i)}), 'LineWidth', 1, 'Color', [0.7, 0.7, 0.7]);
            plot(1:numel(SE_ROMCellArray{indices(i)}), (ROMTable{indices(i)} - SE_ROMCellArray{indices(i)}), 'LineWidth', 1, 'Color', [0.7, 0.7, 0.7]);
        end
   PeaksPlot = figure(2);
        set(gcf, 'position', get(0, 'ScreenSize'))
        sgtitle('Adaptation Time Series for Peak Forces')
        set(sgtitle, 'FontWeight', 'bold');
        hold on
        
        titles = {
            'More Affected Propulsion Force Peak Adaptation',
            'Less Affected Propulsion Force Peak Adaptation',
            'More Affected Braking Force Peak Adaptation',
            'Less Affected Braking Force Peak Adaptation',
            'More Affected Vertical Force Early Peak Adaptation',
            'Less Affected Vertical Force Early Peak Adaptation',
            'More Affected Vertical Force Late Peak Adaptation',
            'Less Affected Vertical Force Late Peak Adaptation'
        };
        
        colors = {'r', 'b'};  % Colors for plotting
        
        indices = [1, 3, 5, 6, 2, 4, 7, 8];  % Indices corresponding to PeaksTable and RM_PeaksTable
        
        for i = 1:8
            subplot(4, 2, i);
            hold on
            plot(PeaksTable{:, indices(i)}, [colors{1} 'o'], 'MarkerSize', 4);
            plot(RM_PeaksTable{:, indices(i)}, colors{1}, 'LineWidth', 3);
            title(titles{i});
            ylabel('Raw Force (N)');
            xlabel('Gait Cycle Number');
            
            % Add error bands
            plot(1:numel(SE_PeaksCellArray{indices(i)}), (PeaksTable{:, indices(i)} + SE_PeaksCellArray{indices(i)}), 'LineWidth', 1, 'Color', [0.7, 0.7, 0.7]);
            plot(1:numel(SE_PeaksCellArray{indices(i)}), (PeaksTable{:, indices(i)} - SE_PeaksCellArray{indices(i)}), 'LineWidth', 1, 'Color', [0.7, 0.7, 0.7]);
        end
