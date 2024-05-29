% Lena NURA '24 - data update, output into Excel

clear all;
clearvars;
clc;

path = '/Users/lena/Desktop/Lena-NURA-24/';
cd(path);

% update this!
latestSubID = 80; % this was the last participant we had according to the researcher spreadsheet 
startSubID = 24; % this was the first participant we had the 45 degrees for


subs = startSubID:latestSubID;


% removing outliers - 2 SD away from the mean in overall accuracy, RT, or error
subs = subs(subs ~= 48);
subs = subs(subs ~= 42);
subs = subs(subs ~= 73);


% remoce any participants who did not pass Ishihara's or the practice after 3 attempts
subs = subs(subs ~= 37); % remove 37, did not pass Ishihara's
subs = subs(subs ~= 38); % remove 38, did not pass practice after 3 attempts
subs = subs(subs ~= 40); % remove 40, did not pass Ishihara's
subs = subs(subs ~= 49); % remove 49, did not pass Ishihara's
subs = subs(subs ~= 50); % remove 50, did not pass practice after 3 attempts
subs = subs(subs ~= 55); % remove 55, stopped responding for the intermixed trials for recog
subs = subs(subs ~= 59); % remove 59, did not pass Ishihara's
subs = subs(subs ~= 60); % remove 60, did not pass Ishihara's
subs = subs(subs ~= 62); % remove 62, stopped responding for the intermixed trials for recog
subs = subs(subs ~= 70); % remove 70, stopped responding for the intermixed trials for recog
subs = subs(subs ~= 70); % remove 74, did not pass Ishihara's
subs = subs(subs ~= 77); % remove 74, did not pass Ishihara's

% 32 dropped automatically later in the script due to chance performance
% 33 dropped later for chance performance
% 72 dropped for chance performance


% ------------------

qualtrics_excel = "Lena NURA - Spring '24.xlsx";
excel_file = "Lena NURA '24 Data.xlsx";
all_data_sheet = "All Data";
demo_sheet = "Demographics";
group_sheet = "mTBI vs Control";
time_sheet = "Performance Over Time";

% update this!
output_name = sprintf('%s/FridayUpdateData/5.28.24.mat',path);

%import data
fullData = readtable(qualtrics_excel);

allN = height(fullData); %total number of responses

% only keep rows of data for valid subject IDs (must be members of subs)
subs_cell = num2cell(subs);
subs_cell_str = cellfun(@num2str, subs_cell, 'UniformOutput', false);
fullData = fullData(ismember(fullData.SubID, subs_cell_str), :);


% Check demographics requirements

% participant 28 does actually qualify - there was a typo in the
% demographics sheet that made them say they did meet a piece of exclusion
% criteria even when they did not (confirmed over email)

% Convert SubID values to strings for comparison
% subs_cell_str = string(subs);

% Filter out rows based on conditions
keep = (fullData.NDisease == 2 & fullData.Age >= 14 & fullData.Vision == 2) | strcmp(fullData.SubID, '28');

% keep = ((fullData.NDisease == 2 & fullData.Age >= 14 & fullData.Vision == 2) | fullData.SubID == 28);
surveyData = fullData(keep,:);


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
    path_to_raw = strcat(path,'Raw/',ID,'_RecallRecog_raw.mat');

    if exist(path_to_file,'file') == 2 % check if this file exists
        load(strcat(ID,'_RecallRecog_ana.mat'));

        if (all_recogBlockAcc <= 0.5 || all_recogIntermixAcc <= 0.5 || aveError_recallIntermixed >= 90 || aveError_recallBlocked >= 90)
            disp(strcat('Participant #',ID,' dropped for at-chance performance'));
            
            removedSubsTask{end + 1,1} = cell2mat(keptSubsDemo(i)); % keep track of those removed this for this reason
        end


        addpath("Raw/")
        load(strcat(ID,'_RecallRecog_raw.mat'));
        
        blockRecall = block_data.recall;
        blockRecog = block_data.recog;

        intermixRecall = intermix_data.recall;
        intermixRecog = intermix_data.recog;

        allRecall = vertcat(blockRecall,blockRecall);
        allRecog = vertcat(blockRecog, intermixRecog);

        % calcualate average accuracy for all recog
        
            % remove recog trials with no response
            allRecogData = allRecog(cat(1,allRecog(:,3)) > -1,:);

            % calculate
            all_RecogAcc = mean(allRecogData(:,3)); % does not discirminate between match and mismatch trials

            
        % calculate median rt for all recog
            all_RecogRT = median(allRecogData(:,4));
        
       
        % calculate average degree error for all recog
        
            % remove recall trials with invalid response
            alllRecallData = allRecall(cat(1,allRecall(:,4)) < 500,:);

            % calculate
            all_AvgRecallErr = mean(abs(alllRecallData(:,4)));


       % add these to excel_data
       excel_data(1,29) = "Average Recall Error";
       excel_data(2,29) = all_AvgRecallErr;
       excel_data(1,30) = "All Recog Accuracy";
       excel_data(2,30) = all_RecogAcc;
       excel_data(1,31) = "All Recog RT";
       excel_data(2,31) = all_RecogRT;
    
       save(path_to_file,'excel_data',"-append");

        
    else
        noData = strcat('Error! No behavioral data for participant #',ID);
        error(noData);
    end

end

% get rid of the demographics data for those who were at chance
removedSubsTask = cellstr(removedSubsTask);
rows_to_remove = ismember(surveyData.SubID, removedSubsTask);
surveyData = surveyData(~rows_to_remove, :);

% get rid of these subs from keptSubsDemo too
keptSubsDemo = keptSubsDemo(~rows_to_remove,:);

N = height(surveyData); % this is our new N


% add path to folder w/ excel sheet
addpath('/Users/lena/Desktop/Friday Data Updates')

% Load in data for each participant, write into Excel (WM data & questionnaire score/group category)
for i = 1:length(keptSubsDemo)

    ID = cell2mat(keptSubsDemo(i));
    load(strcat(ID,'_RecallRecog_ana.mat'));

    BIS_score = surveyData.SC1(i);
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
    cell_range1 = strcat('B',num2str(i + 1),':AF',num2str(i + 1));

    writematrix(str2num(ID), excel_file, 'Sheet',all_data_sheet, 'Range',strcat('A',num2str(i + 1)));
    writematrix(str2double(excel_data(end,1:end)), excel_file, 'Sheet',all_data_sheet, 'Range',cell_range1);
    writematrix(BIS_score, excel_file, 'Sheet',all_data_sheet, 'Range',strcat('AG',num2str(i + 1)));


    if (surveyData.PreviousTBI(i) == 1) % they are mTBI
        status = 'mTBI';
        mTBI_subs{end + 1} = ID;
    else % they are control
        status = 'Control';
        control_subs{end + 1} = ID;
    end

    writematrix(status, excel_file, 'Sheet',all_data_sheet, 'Range',strcat('AH',num2str(i + 1)));


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
    subid_numeric = cellfun(@str2double, surveyData.SubID); % Convert cell array to numeric array

    subData = surveyData(subid_numeric == ID_dbl, :);
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
    cell_range1 = strcat('B',num2str(i + 2),':AF',num2str(i + 2));


    writematrix(str2num(ID), excel_file, 'Sheet',group_sheet,'Range',strcat('A',num2str(i + 2)));
    writematrix(str2double(excel_data(end,1:end)), excel_file,'Sheet', group_sheet,'Range',cell_range1);
    writematrix(BIS_score, excel_file, 'Sheet',group_sheet, 'Range',strcat('AG',num2str(i + 2)));

end



% write control info
for i = 1:length(control_subs)

    ID = cell2mat(control_subs(i));
    ID_dbl = str2double(ID);
    load(strcat(ID,'_RecallRecog_ana.mat'));

    % Assuming 'SubID' is the column name
    subid_numeric = cellfun(@str2double, surveyData.SubID); % Convert cell array to numeric array

    subData = surveyData(subid_numeric == ID_dbl, :);
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
    cell_range1 = strcat('B',num2str(i + 23),':AF',num2str(i + 23));

    writematrix(str2num(ID), excel_file, 'Sheet',group_sheet,'Range',strcat('A',num2str(i + 23)));
    writematrix(str2double(excel_data(end,1:end)), excel_file,'Sheet', group_sheet,'Range',cell_range1);
    writematrix(BIS_score, excel_file, 'Sheet',group_sheet, 'Range',strcat('AG',num2str(i + 23)));



end


disp('Done!')
