%% Split-belt Treadmill Joint Profile Analysis
% Originating authors: Isaac Bast and Andrew Hagen
% Last Updated: 12/17/2023

% This script uses data from a structure to calculate the maximum value, minimum value, ROM, max location, and min location for all joints and forces in all phases of a split-belt treadmill adaptation paradigm
%   see Vicon_Data_Grabber in repository for structure organization and creation using Vicon Nexus outputs
% Outputs a table for the results of all participants in SelectedParticipants
% Plots average gait cycle profiles for all participants

%% Script Initialization

% Make sure you have GaitCycle_Data_Interpolated.mat loaded. This will take a bit so its best to do this manually rather than have it be loaded up each time the script is run

% Don't accidentally clear the mat structure 
    clearvars -except GaitCycle_Data_Interpolated

% Set Path 
%     RootPath = 'R:\SBMS Joint Coordination\Joint Data\';
%     cd(RootPath)

% Define the list of specific participants you want to be included in analysis (first layer of fields in the structure)
SelectedParticipants = {'SBMS_004','SBMS_010','SBMS_011','SBMS_014','SBMS_026','SBMS_027','SBMS_040','SBMS_012','SBMS_018','SBMS_021','SBMS_022','SBMS_028','SBMS_030','SBMS_033', 'SBMS_006','SBMS_009','SBMS_013','SBMS_016','SBMS_019','SBMS_023','SBMS_024','SBMS_031','SBMS_037','SBMS_002','SBMS_003','SBMS_007','SBMS_017','SBMS_029','SBMS_034','SBMS_035','SBMS_036','SBMS_041'};
ParticipantNames = cell(1, numel(SelectedParticipants));

% Initialize a cell array to store the results tables for each participant
ResultsTables = cell(1, numel(SelectedParticipants));

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

    % Initialize a cell array to store results for the current participant
    ResultsCell = cell(1, 270);  

    % Nested loops to iterate through periods, limbs, angles, phases, and planes
    % Look at Vicon_Data_Grabber in repository for how this structure is organized
    h = 1;
    for j = 1:(numel(Periods) - 1) % Excludes Adaptation 2
        for k = 1:numel(Limbs) % 2 Limbs
            for l = 1:(numel(Angles) - 1) % 6 Angles, excludes Force data
                for m = 1:numel(Phases) % Phases of gait: stance, swing
                    for n = 1 % Only looks at X-axis (sagittal plane angles)

                        tempName = strcat(expr, Periods{j}, '.', Limbs{k}, '.', Angles{l}, '.', Phases{m}, '.', Planes{n});
                        Expression(h, 1) = string(tempName);
                        h = h + 1;
                    end
                end
            end
        end
    end
   
 % Create variable names to assign to average output tables
    baselineVarNamesIDX = find(contains(Expression,'Baseline'));
    baselineVarNames_Long = Expression(baselineVarNamesIDX);
    baselineVarNames = erase(baselineVarNames_Long,'GaitCycle_Data_Interpolated.');
    
    adaptVarNamesIDX = find(contains(Expression,'Adaptation1'));
    adaptVarNames_Long = Expression(adaptVarNamesIDX);
    adaptVarNames = erase(adaptVarNames_Long,'GaitCycle_Data_Interpolated.');
    
    catchVarNamesIDX = find(contains(Expression,'Catch'));
    catchVarNames_Long = Expression(catchVarNamesIDX);
    catchVarNames = erase(catchVarNames_Long,'GaitCycle_Data_Interpolated.');
    
    clear baselineVarNamesIDX baselineVarNames_Long adaptVarNamesIDX adaptVarNames_Long catchVarNamesIDX catchVarNames_Long
    
    % Define iteration variables for average calculation loop
    avgMixedOutput = zeros(120,100);
    CalcIDX = numel(Expression);
    
    for p = 1: CalcIDX
        tempRaw = eval(Expression(p)); % Pull data from Struct into a temporary variable
        [tempRows,tempCols] = size(tempRaw);
    
        % Crop to window of data based on Period, calculate average for output,
        % write to mixed output matrix with known location
    
        if contains(Expression(p),'Baseline') == 1
            tempWindow = tempRaw(tempRows-14:tempRows,:);
            for bIDX = 1:100
                tempBaselineAvg = mean(nonzeros(tempWindow(:,bIDX)));
                avgMixedOutput(p,bIDX) = tempBaselineAvg;
            end
        elseif contains(Expression(p),'Adaptation1') == 1
            tempWindowEarly = tempRaw(6:15,:);
            tempWindowLate = tempRaw(tempRows-14:tempRows,:);
    
            for aIDX = 1:100
                tempEAdaptAvg = mean(nonzeros(tempWindowEarly(:,aIDX)));
                avgMixedOutput(p,aIDX) = tempEAdaptAvg;
                
                tempLAdaptAvg = mean(nonzeros(tempWindowLate(:,aIDX)));
                avgMixedOutput(p+48,aIDX) = tempLAdaptAvg;
            end
    
        elseif contains(Expression(p),'Catch')
            tempWindowEarly = tempRaw(6:15,:);
            tempWindowLate = tempRaw(tempRows-14:tempRows,:);
    
            for cIDX = 1:100
                tempECatchAvg = mean(nonzeros(tempWindowEarly(:,cIDX)));
                avgMixedOutput(p,cIDX) = tempECatchAvg;
    
                tempLCatchAvg = mean(nonzeros(tempWindowLate(:,cIDX)));
                avgMixedOutput(p+48,cIDX) = tempLCatchAvg;
            end
        end
    end
    clear tempName tempRows tempCols tempBaselineAvg  tempEAdaptAvg tempECatchAvg tempLAdaptAvg tempLCatchAvg tempRaw tempWindow tempWindowEarly tempWindowLate
    clear aIDX bIDX cIDX
    
    % Write Avg Output data to tables divided by Period, with rows named to
    % match all descriptors
    
    BaselineAvgTable = array2table(avgMixedOutput(1:24,:), 'RowNames', baselineVarNames);
    EAdaptAvgTable = array2table(avgMixedOutput(25:48,:), 'RowNames', adaptVarNames);
    ECatchAvgTable = array2table(avgMixedOutput(49:72,:), 'RowNames', catchVarNames);
    LAdaptAvgTable = array2table(avgMixedOutput(73:96,:), 'RowNames', adaptVarNames);
    LCatchAvgTable = array2table(avgMixedOutput(97:120,:), 'RowNames', catchVarNames);
    
    
    %% Make calculations and assign them to a table
    
    % Define the row locations for stance and swing for each joint
    StanceRowLocations = {1, 13, 3, 15, 5, 17};
    SwingRowLocations = {2, 14, 4, 16, 6, 18};
    JointNames = {'MA_Hip', 'LA_Hip', 'MA_Knee', 'LA_Knee', 'MA_Ankle', 'LA_Ankle'};
    phases = {'Baseline', 'EAdapt', 'LAdapt', 'ECatch', 'LCatch'};
    
    % Initialize a cell array to store the results for each joint and phase (9 variables per joint per phase)
    ResultsCell = cell(1, 9 * length(JointNames) * length(phases));
    
    % Loop over your iterations
    for q = 1:length(JointNames)
        % Loop over each phase (Baseline, EAdapt, LAdapt, ECatch, LCatch)
        for r = 1:length(phases)
            % Extract the row of interest from each table for each joint
            StanceData = eval([phases{r}, 'AvgTable{StanceRowLocations{q},:}']);
            SwingData = eval([phases{r}, 'AvgTable{SwingRowLocations{q},:}']);
            AllData = horzcat(StanceData, SwingData);

            CycleData{q,r} = AllData;
    
            % Calculate the maximum value, minimum value, ROM, max location, and min location for the current joint and phase
            eval([phases{r}, '_', JointNames{q}, '_ROM = range(AllData);']);
            eval([phases{r}, '_', JointNames{q}, '_Stance_Max = max(StanceData, [], 2);']); 
            eval([phases{r}, '_', JointNames{q}, '_Stance_Min = min(StanceData, [], 2);']);
            eval([phases{r}, '_', JointNames{q}, '_Swing_Max = max(SwingData, [], 2);']);
            eval([phases{r}, '_', JointNames{q}, '_Swing_Min = min(SwingData, [], 2);']);
            eval(['[~, loc] = max(StanceData(:, 20:100), [], 2); ', phases{r}, '_', JointNames{q}, '_Stance_Max_Loc = loc;']);
            eval(['[~, loc] = min(StanceData(:, 80:100), [], 2); ', phases{r}, '_', JointNames{q}, '_Stance_Min_Loc = loc;']);
            eval(['[~, loc] = max(SwingData(:, 20:100), [], 2); ', phases{r}, '_', JointNames{q}, '_Swing_Max_Loc = loc;']);
            eval(['[~, loc] = min(SwingData(:, 1:20), [], 2); ', phases{r}, '_', JointNames{q}, '_Swing_Min_Loc = loc;']);
 
            % Only use the specific windows (e.g. (:, 20:100)) for ankle - it wont work for hip or knee
%             eval([phases{r}, '_', JointNames{q}, '_ROM = range(AllData);']);
%             eval([phases{r}, '_', JointNames{q}, '_Stance_Max = max(StanceData(:, 20:100), [], 2);']); 
%             eval([phases{r}, '_', JointNames{q}, '_Stance_Min = min(StanceData(:, 80:100), [], 2);']);
%             eval([phases{r}, '_', JointNames{q}, '_Swing_Max = max(SwingData(:, 20:100), [], 2);']);
%             eval([phases{r}, '_', JointNames{q}, '_Swing_Min = min(SwingData(:, 1:20), [], 2);']);
%             eval(['[~, loc] = max(StanceData(:, 20:100), [], 2); ', phases{r}, '_', JointNames{q}, '_Stance_Max_Loc = loc + 19;']);
%             eval(['[~, loc] = min(StanceData(:, 80:100), [], 2); ', phases{r}, '_', JointNames{q}, '_Stance_Min_Loc = loc + 79;']);
%             eval(['[~, loc] = max(SwingData(:, 20:100), [], 2); ', phases{r}, '_', JointNames{q}, '_Swing_Max_Loc = loc + 19;']);
%             eval(['[~, loc] = min(SwingData(:, 1:20), [], 2); ', phases{r}, '_', JointNames{q}, '_Swing_Min_Loc = loc;']);

    
            % Store the results in the cell array
            % In the index calculations, use 9 to represent the number of variables per joint and phase
            index = (q - 1) * 9 * length(phases) + (r - 1) * 9 + 1;
            ResultsCell{1, index} = eval([phases{r}, '_', JointNames{q}, '_ROM']);
            ResultsCell{1, index + 1} = eval([phases{r}, '_', JointNames{q}, '_Stance_Max']);
            ResultsCell{1, index + 2} = eval([phases{r}, '_', JointNames{q}, '_Stance_Min']);
            ResultsCell{1, index + 3} = eval([phases{r}, '_', JointNames{q}, '_Swing_Max']);
            ResultsCell{1, index + 4} = eval([phases{r}, '_', JointNames{q}, '_Swing_Min']);
            ResultsCell{1, index + 5} = eval([phases{r}, '_', JointNames{q}, '_Stance_Max_Loc']);
            ResultsCell{1, index + 6} = eval([phases{r}, '_', JointNames{q}, '_Stance_Min_Loc']);
            ResultsCell{1, index + 7} = eval([phases{r}, '_', JointNames{q}, '_Swing_Max_Loc']);
            ResultsCell{1, index + 8} = eval([phases{r}, '_', JointNames{q}, '_Swing_Min_Loc']);
        end
    end
    
    % Convert the cell array to a table
    VariableNames = cell(1, 9 * length(JointNames) * length(phases));
    for q = 1:length(JointNames)
        for r = 1:length(phases)
            index = (q - 1) * 9 * length(phases) + (r - 1) * 9 + 1;
            VariableNames{index} = [phases{r}, '_', JointNames{q}, '_ROM'];
            VariableNames{index + 1} = [phases{r}, '_', JointNames{q},'_Stance_Max_Ang'];
            VariableNames{index + 2} = [ phases{r}, '_', JointNames{q},'_Stance_Min_Ang'];
            VariableNames{index + 3} = [phases{r}, '_', JointNames{q},'_Swing_Max_Ang'];
            VariableNames{index + 4} = [phases{r}, '_', JointNames{q},'_Swing_Min_Ang'];
            VariableNames{index + 5} = [phases{r}, '_', JointNames{q},'_Stance_Max_Ang_Loc'];
            VariableNames{index + 6} = [phases{r}, '_', JointNames{q},'_Stance_Min_Ang_Loc'];
            VariableNames{index + 7} = [phases{r}, '_', JointNames{q},'_Swing_Max_Ang_Loc'];
            VariableNames{index + 8} = [phases{r}, '_', JointNames{q},'_Swing_Min_Ang_Loc'];
        end
    end
  
    % Store the participant name
    ParticipantNames{i} = participant;

    % Create a results table for the current participant
    ResultsTable = array2table(cell2mat(ResultsCell), 'VariableNames', VariableNames);

    % Store the table and cycle data in the cell arrays
    ResultsTables{i} = ResultsTable;
    ResultsCycleData{i} = CycleData;

    % Clear variables for the next iteration
    clear ResultsCell Expression VariableNames

end

% Combine all participants' results into one table
CombinedResultsTable = vertcat(ResultsTables{:});
CombinedResultsTable.Properties.RowNames = ParticipantNames;

% Calculate the means of each column 
Means = mean(CombinedResultsTable{:, :}); 

% Create a "Mean" row with NaN for any non-numeric columns
MeanRow = array2table(Means, 'VariableNames', CombinedResultsTable.Properties.VariableNames, 'RowNames', {'Means'});
VariableNames = CombinedResultsTable.Properties.VariableNames;
MeanRow.Properties.VariableNames = VariableNames;

% Append the "Mean" row to the bottom of the table
CombinedResultsTable = [CombinedResultsTable; MeanRow];

ExportData = inputdlg('Save and Export to Data xlsx? Y/N: ');
 if strcmp(ExportData,'Y')
        % Prompt the user to select a directory and specify a file name
        [fileName, filePath] = uiputfile('*.xlsx', 'Save Table as Excel', 'JointProfileAnalysis.xlsx');

        % Check if the user canceled the operation
        if isequal(fileName, 0) || isequal(filePath, 0)
         disp('File save canceled.');
        else
        % Construct the full file path
        fullFilePath = fullfile(filePath, fileName);

        % Save the table to the selected location and sheet
        writetable(CombinedResultsTable, fullFilePath, 'WriteRowNames', true);  % Include row names

        disp(['Table saved to ' fullFilePath]);
        end
 end

% Create mean gait cycle profiles for all selected participants
AllCycleData = cell(6,5);
MeanCycleData = cell (6,5);
for s = 1:numel(SelectedParticipants)
    for row = 1:6
        for col = 1:5
            AllCycleData{row, col} = vertcat(AllCycleData{row, col}, ResultsCycleData{s}{row, col});
            MeanCycleData{row, col} = mean(AllCycleData{row, col},'omitnan');
        end
    end
end

CycleData_Smoothed = cell(size(MeanCycleData)); % Smooth ROM SE Lower Bound

for t = 1:numel(CycleData_Smoothed) 
    x = MeanCycleData{t}'; 
    SmoothedData = fit((1:numel(x))', x, 'smoothingspline', 'SmoothingParam', 0.1);
    SmoothedValues= feval(SmoothedData, (1:numel(x))'); % Evaluate the fit over the range of values
    CycleData_Smoothed{t} = SmoothedValues;
end

% Make plots of each joint data over all phases for selected participants
ExportPlot = inputdlg('Create plots? Y/N: ');
 if strcmp(ExportPlot,'Y')
        CyclePlot = figure(1);
            TimeVector = 1:200;
            set(gcf,'position',get(0, 'ScreenSize'))
            sgtitle('Adaptation Time Series for Gait Cycle Joint Profile')
            set(sgtitle, 'FontWeight', 'bold');    
            hold on

        subplot (3,2,1) % More affected hip
        hold on
            plot(TimeVector, CycleData_Smoothed{1,1}, LineWidth=1.75, Color=[0, 0.278, 0.671])
            plot(TimeVector, CycleData_Smoothed{1,2}, LineWidth=1.75, Color=[0.875, 0.455, 0.318])
            plot(TimeVector, CycleData_Smoothed{1,3}, LineWidth=1.75, Color=[0.431, 0.149, 0.055])
            plot(TimeVector, CycleData_Smoothed{1,4}, LineWidth=1.75, Color=[0.010, 0.545, 0.341])
            plot(TimeVector, CycleData_Smoothed{1,5}, LineWidth=1.75, Color=[0, 0.278, 0.165])
            title('MA Hip Joint Profiles')
            legend('Baseline','Early Adapt','Late Adapt','Early Post-Adapt','Late Post-Adapt')
            set(legend,'Location', 'best', fontsize=12)
            ylabel('Degrees');
            xlabel('Gait Cycle Percentage');
            xticks([0,20,40,60,80,100,120,140,160,180,200]);  % Adjust the x-ticks limit
            xticklabels({'0', '10', '20', '30', '40', '50', '60', '70', '80', '90', '100'});
        subplot(3,2,2) % Less affected hip
        hold on
            plot(TimeVector, CycleData_Smoothed{2,1}, LineWidth=1.75, Color=[0, 0.278, 0.671])
            plot(TimeVector, CycleData_Smoothed{2,2}, LineWidth=1.75, Color=[0.875, 0.455, 0.318])
            plot(TimeVector, CycleData_Smoothed{2,3}, LineWidth=1.75, Color=[0.431, 0.149, 0.055])
            plot(TimeVector, CycleData_Smoothed{2,4}, LineWidth=1.75, Color=[0.010, 0.545, 0.341])
            plot(TimeVector, CycleData_Smoothed{2,5}, LineWidth=1.75, Color=[0, 0.278, 0.165])
            title('LA Hip Joint Profiles')
            ylabel('Degrees');
            xlabel('Gait Cycle Percentage');
            xticks([0,20,40,60,80,100,120,140,160,180,200]);  % Adjust the x-ticks limit
            xticklabels({'0', '10', '20', '30', '40', '50', '60', '70', '80', '90', '100'});
        subplot(3,2,3) % More affected knee
        hold on
            plot(TimeVector, CycleData_Smoothed{3,1}, LineWidth=1.75, Color=[0, 0.278, 0.671])
            plot(TimeVector, CycleData_Smoothed{3,2}, LineWidth=1.75, Color=[0.875, 0.455, 0.318])
            plot(TimeVector, CycleData_Smoothed{3,3}, LineWidth=1.75, Color=[0.431, 0.149, 0.055])
            plot(TimeVector, CycleData_Smoothed{3,4}, LineWidth=1.75, Color=[0.010, 0.545, 0.341])
            plot(TimeVector, CycleData_Smoothed{3,5}, LineWidth=1.75, Color=[0, 0.278, 0.165])
            title('MA Knee Joint Profiles')
            ylabel('Degrees');
            xlabel('Gait Cycle Percentage');
            xticks([0,20,40,60,80,100,120,140,160,180,200]);  % Adjust the x-ticks limit
            xticklabels({'0', '10', '20', '30', '40', '50', '60', '70', '80', '90', '100'});
        subplot(3,2,4) % Less affected knee
        hold on
            plot(TimeVector, CycleData_Smoothed{4,1}, LineWidth=1.75, Color=[0, 0.278, 0.671])
            plot(TimeVector, CycleData_Smoothed{4,2}, LineWidth=1.75, Color=[0.875, 0.455, 0.318])
            plot(TimeVector, CycleData_Smoothed{4,3}, LineWidth=1.75, Color=[0.431, 0.149, 0.055])
            plot(TimeVector, CycleData_Smoothed{4,4}, LineWidth=1.75, Color=[0.010, 0.545, 0.341])
            plot(TimeVector, CycleData_Smoothed{4,5}, LineWidth=1.75, Color=[0, 0.278, 0.165])
            title('LA Knee Joint Profiles')
            ylabel('Degrees');
            xlabel('Gait Cycle Percentage');
            xticks([0,20,40,60,80,100,120,140,160,180,200]);  % Adjust the x-ticks limit
            xticklabels({'0', '10', '20', '30', '40', '50', '60', '70', '80', '90', '100'});
        subplot(3,2,5) % More affected ankle
        hold on
            plot(TimeVector, CycleData_Smoothed{5,1}, LineWidth=1.75, Color=[0, 0.278, 0.671])
            plot(TimeVector, CycleData_Smoothed{5,2}, LineWidth=1.75, Color=[0.875, 0.455, 0.318])
            plot(TimeVector, CycleData_Smoothed{5,3}, LineWidth=1.75, Color=[0.431, 0.149, 0.055])
            plot(TimeVector, CycleData_Smoothed{5,4}, LineWidth=1.75, Color=[0.010, 0.545, 0.341])
            plot(TimeVector, CycleData_Smoothed{5,5}, LineWidth=1.75, Color=[0, 0.278, 0.165])
            title('MA Ankle Joint Profiles')
            ylabel('Degrees');
            xlabel('Gait Cycle Percentage');
            xticks([0,20,40,60,80,100,120,140,160,180,200]);  % Adjust the x-ticks limit
            xticklabels({'0', '10', '20', '30', '40', '50', '60', '70', '80', '90', '100'});
        subplot(3,2,6) % less affected ankle
        hold on
            plot(TimeVector, CycleData_Smoothed{6,1}, LineWidth=1.75, Color=[0, 0.278, 0.671])
            plot(TimeVector, CycleData_Smoothed{6,2}, LineWidth=1.75, Color=[0.875, 0.455, 0.318])
            plot(TimeVector, CycleData_Smoothed{6,3}, LineWidth=1.75, Color=[0.431, 0.149, 0.055])
            plot(TimeVector, CycleData_Smoothed{6,4}, LineWidth=1.75, Color=[0.010, 0.545, 0.341])
            plot(TimeVector, CycleData_Smoothed{6,5}, LineWidth=1.75, Color=[0, 0.278, 0.165])
            title('LA Ankle Joint Profiles')
            ylabel('Degrees');
            xlabel('Gait Cycle Percentage');
            xticks([0,20,40,60,80,100,120,140,160,180,200]);  % Adjust the x-ticks limit
            xticklabels({'0', '10', '20', '30', '40', '50', '60', '70', '80', '90', '100'});
 end
%  