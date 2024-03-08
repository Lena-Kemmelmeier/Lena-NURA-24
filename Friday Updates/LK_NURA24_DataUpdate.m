% Lena NURA '24 - data update, output into Excel

clear all;
clearvars;
clc;

path = '/Users/lena/Desktop/Lena-NURA-24/';
cd(path);

% update this!
latestSubID = 18;

% subs = 1:latestSubID; % when latestSubID was 8
% subs = subs(subs ~= 6); % remove 6, did not pass Ishihara's

subs = 12:latestSubID; % 12 was first participant with changed times
subs = subs(subs ~= 13); % remove 13, did not pass Ishihara's


% ------------------

qualtrics_excel = "Lena NURA - Spring '24.xlsx";
excel_file = "Lena NURA '24 Data.xlsx";
all_data_sheet = "All Data";
demo_sheet = "Demographics";
group_sheet = "mTBI vs Control";
time_sheet = "Performance Over Time";

% update this!
output_name = sprintf('%s/FridayUpdateData/3.7.24.mat',path);

%import data
fullData = readtable(qualtrics_excel);

allN = height(fullData); %total number of responses

% only keep rows of data for valid subject IDs (must be members of subs)
subs_cell = num2cell(subs);
subs_cell_str = cellfun(@num2str, subs_cell, 'UniformOutput', false);
fullData = fullData(ismember(fullData.SubID, subs_cell_str), :);


% Check demographics requirements
keep = (fullData.NDisease ==2 & fullData.Age >= 14 & fullData.Vision ==2);
data = fullData(keep,:);


% store subIDs for those who did not meet demographics requirements
removedSubsDemo = fullData.SubID(~keep,:);
removedSubsTask = {};
keptSubsDemo = fullData.SubID(keep,:);
mTBI_subs = {};
control_subs = {};
% removedSubs{end + 1} = '7'; % test

% Check if any of these participants that met demographics requirements did
% not meet the performance minimum (at-chance)

path_to_ana = strcat(path,'Analyzed/');
addpath(path_to_ana);

for i = 1:length(keptSubsDemo)
    
    ID = cell2mat(keptSubsDemo(i));

    path_to_file = strcat(path,'Analyzed/',ID,'_RecallRecog_ana.mat');

    if exist(path_to_file,'file') == 2 % check if this file exists
        load(strcat(ID,'_RecallRecog_ana.mat'));

        if (all_recogBlockAcc <= 0.5 || all_recogIntermixAcc <= 0.5 || aveError_recallIntermixed >= 90 || aveError_recallBlocked >= 90)
            disp(strcat('Participant #',ID,' dropped for at-chance performance'));
            
            removedSubsTask{end + 1,1} = cell2mat(keptSubsDemo(i)); % keep track of those removed this for this reason
        end

    else
        noData = strcat('Error! No behavioral data for participant #',ID);
        error(noData);
    end

end

% get rid of the demographics data for those who were at chance
removedSubsTask = cellstr(removedSubsTask);
rows_to_remove = ismember(data.SubID, removedSubsTask);
data = data(~rows_to_remove, :);

% get rid of these subs from keptSubsDemo too
keptSubsDemo = keptSubsDemo(~rows_to_remove,:);

N = height(data); % this is our new N


% add path to folder w/ excel sheet
addpath('/Users/lena/Desktop/Friday Data Updates')

% Load in data for each participant, write into Excel (WM data & questionnaire score/group category)
for i = 1:length(keptSubsDemo)

    ID = cell2mat(keptSubsDemo(i));
    load(strcat(ID,'_RecallRecog_ana.mat'));

    BIS_score = data.SC1(i);
    age = fullData.Age(i);

    % get their gender
    switch fullData.Gender(i)
        case 1
            gender = "Female";
        case 2
            gender = "Male";
        case 3
            gender = "Other";
    end

  
    % write in data into All Data sheet
    cell_range1 = strcat('B',num2str(i + 1),':AC',num2str(i + 1));
    
    writematrix(str2num(ID), excel_file, 'Sheet',all_data_sheet, 'Range',strcat('A',num2str(i + 1)));
    writematrix(str2double(excel_data(end,1:end)), excel_file, 'Sheet',all_data_sheet, 'Range',cell_range1);
    writematrix(BIS_score, excel_file, 'Sheet',all_data_sheet, 'Range',strcat('AD',num2str(i + 1)));


    if (data.PreviousTBI(i) == 1) % they are mTBI
        status = 'mTBI';
        mTBI_subs{end + 1} = ID;
    else % they are control
        status = 'Control';
        control_subs{end + 1} = ID;
    end

    writematrix(status, excel_file, 'Sheet',all_data_sheet, 'Range',strcat('AE',num2str(i + 1)));


    % write in data into Demographics sheet
    writematrix(str2num(ID), excel_file, 'Sheet',demo_sheet, 'Range',strcat('A',num2str(i + 1)));
    writematrix(age, excel_file, 'Sheet',demo_sheet, 'Range',strcat('B',num2str(i + 1)));
    writematrix(gender, excel_file, 'Sheet',demo_sheet, 'Range',strcat('C',num2str(i + 1)));
    writematrix(status, excel_file, 'Sheet',demo_sheet, 'Range',strcat('D',num2str(i + 1)));

    % write in performance over time sheet into Performance over Time
    writematrix(str2num(ID), excel_file, 'Sheet',time_sheet, 'Range',strcat('A',num2str(i + 1)));
    writematrix((performance_time.blockRecall(:,2))',excel_file,'Sheet',time_sheet,'Range',strcat('B',num2str(i + 1),':E',num2str(i + 1)));
    writematrix((performance_time.mixRecall(:,2))',excel_file,'Sheet',time_sheet,'Range',strcat('F',num2str(i + 1),':I',num2str(i + 1)));
    writematrix((performance_time.blockRecog(:,2))',excel_file,'Sheet',time_sheet,'Range',strcat('J',num2str(i + 1),':M',num2str(i + 1)));
    writematrix((performance_time.mixRecog(:,2))',excel_file,'Sheet',time_sheet,'Range',strcat('N',num2str(i + 1),':Q',num2str(i + 1)));
    writematrix((performance_time.blockRecog(:,3))',excel_file,'Sheet',time_sheet,'Range',strcat('R',num2str(i + 1),':U',num2str(i + 1)));
    writematrix((performance_time.mixRecog(:,3))',excel_file,'Sheet',time_sheet,'Range',strcat('V',num2str(i + 1),':Y',num2str(i + 1)));

end


% mTBI vs. control split file

% write mTBI info
for i = 1:length(mTBI_subs)

    ID = cell2mat(mTBI_subs(i));
    ID_dbl = str2double(ID);
    load(strcat(ID,'_RecallRecog_ana.mat'));

    % Assuming 'SubID' is the column name
    subid_numeric = cellfun(@str2double, data.SubID); % Convert cell array to numeric array
    
    subData = data(subid_numeric == ID_dbl, :);
    subData = subData(1,:);
    


    BIS_score = subData.SC1;
    age = subData.Age;

    % get their gender
    switch subData.Gender
        case 1
            gender = "Female";
        case 2
            gender = "Male";
        case 3
            gender = "Other";
    end

  
    % write in data into mTBI vs. control sheet
    cell_range1 = strcat('B',num2str(i + 2),':AC',num2str(i + 2));

    
    writematrix(str2num(ID), excel_file, 'Sheet',group_sheet,'Range',strcat('A',num2str(i + 2)));
    writematrix(str2double(excel_data(end,1:end)), excel_file,'Sheet', group_sheet,'Range',cell_range1);
    writematrix(BIS_score, excel_file, 'Sheet',group_sheet, 'Range',strcat('AD',num2str(i + 2)));

end



% write control info
for i = 1:length(control_subs)

    ID = cell2mat(control_subs(i));
    ID_dbl = str2double(ID);
    load(strcat(ID,'_RecallRecog_ana.mat'));

    % Assuming 'SubID' is the column name
    subid_numeric = cellfun(@str2double, data.SubID); % Convert cell array to numeric array
    
    subData = data(subid_numeric == ID_dbl, :);
    subData = subData(1,:);
    


    BIS_score = subData.SC1;
    age = subData.Age;

    % get their gender
    switch subData.Gender
        case 1
            gender = "Female";
        case 2
            gender = "Male";
        case 3
            gender = "Other";
    end

  
    % write in data into mTBI vs. control sheet
    cell_range1 = strcat('B',num2str(i + 23),':AC',num2str(i + 23));
    
    writematrix(str2num(ID), excel_file, 'Sheet',group_sheet,'Range',strcat('A',num2str(i + 23)));
    writematrix(str2double(excel_data(end,1:end)), excel_file,'Sheet', group_sheet,'Range',cell_range1);
    writematrix(BIS_score, excel_file, 'Sheet',group_sheet, 'Range',strcat('AD',num2str(i + 23)));



end


disp('Done!')
