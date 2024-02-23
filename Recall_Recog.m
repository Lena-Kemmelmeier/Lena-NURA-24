% L Kemmelmeier - Combined color recall (delayed estimation) and recognition (change detection) task
% Based off the code by A and B (add these), task design: Zhang and Luck 2008

clear;
root = '/Users/mblab/Desktop/Lena-NURA-24/';
addpath(root); % just to be sure...
addpath(root, char(root + "/Raw"), char(root + "/Data"), char(root + "/Analyzed"))

% GUI
noResponseGUI = 'Program aborted, missing participant info';
incorrectGUI = 'Program aborted, incorrect mouse or order settings entered';
prompt = {'Enter participant number:','Enter mouse settings:','Enter order settings:'};
dlg_title = 'New Participant';
num_lines = 1;
answer = inputdlg(prompt,dlg_title,num_lines); % presents box to enter data into

if isempty(answer{1,1} | answer{2,1} | answer{3,1})
    error(noResponseGUI); % throw an error
elseif ((answer{2,1} ~= '0' && answer{2,1} ~= '1') || (answer{3,1} ~= '0' && answer{3,1} ~= '1'))
    error(incorrectGUI); % throw an error
end 

% if everything is good, initialize variables
ID = answer{1,1}; % subID
mouseCondition = str2num(answer{2,1}); % 0 = right key is match, left key is mismatch; 1 = left key is match, right key is mismatch
a = str2num(answer{3,1}); % 1 = intermixed trial set; 0 = blocked trial set
b = ~a; % if a is 1, then b (second condition) is 0, and vice versa (1 = intermixed trial set; 0 = blocked trial set)
conditionSeq = [a b]; % make condition sequence


% Check for duplicated file/pre-existing participant ID

duplicateSubID = 'Program aborted, files with the same subID already exist';
test_path = char(root + "/Raw" + '/' + strcat(ID,'_RecallRecog_raw.mat'));

if exist(test_path,'file') == 2 % 2 if file exists, 0 if it doesn't
    error(duplicateSubID);
end

output_name = sprintf('%s/Raw/%s_RecallRecog_raw',root,ID);
data_name = sprintf('%s/Data/%s_RecallRecog_data',root,ID);
ana_name = sprintf('%s/Analyzed/%s_RecallRecog_ana',root,ID);




  try % prepare stim, run task
    prepareEnvironment;

    [window, rect] = openWindow(); % added by LK
    
    % Code below is from Conjuctions_Exp4.m
    % Turn on alpha blending. this makes drawing the stims much easier.
    % Type Screen('BlendingFunction?') for more info.
    Screen('BlendFunction',window.onScreen,GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
    
    prefs = getPreferences(); % prefs.nTrialsRecall must always be even!


    nTrialsTotal = prefs.nTrialsRecall * 2; % this is the total number of trials for each intermixed and blocked conditions

    % Base trial lists
    recallTrials = zeros(nTrialsTotal/2, 1) + 1; % 1 denotes recall trials
    recogMatchTrials = zeros(nTrialsTotal/4, 1) + 2; % 2 denotes recognition trials where the probe DOES match the sample array
    recogMismatchTrials = zeros(nTrialsTotal/4, 1) + 3; % 3 denotes recognition trials where the probe DOES NOT match the sample aaray
    recogTrials = vertcat(recogMismatchTrials, recogMatchTrials); % combine the match and mismatch recog trials

    % Create the intermixed sequence
    intermixedTrialList = vertcat(recallTrials,recogTrials); % combine the recall and recog trials
    intermixedTrialList = intermixedTrialList(randperm(length(intermixedTrialList))); % scramble the order of the recall and recog trials
    
    % Settings for block
    numBlocks = 8;
    trialsPerBlock = 25; % change to 25 once we're done!

    % Determine the order of the recognition and recall blocks
    recallBlocks = zeros(numBlocks/2) + 5; % 5 denotes recall block
    recogBlocks = zeros(numBlocks/2) + 6; % 6 denotes recognition block
    blockList = vertcat(recallBlocks, recogBlocks); % combine the recognition and recall blocks
    blockList = blockList(randperm(length(blockList))); % scramble the order of the different blocks

    scrambledRecogTrials = recogTrials(randperm(length(recogTrials))); % the trials for recognition blocks should be pseudo-randomized in whether they are match/mismatch

    
    % Create the blocked sequence
    blockedTrialList = [];
    recogTrialCtr = 0;

    for i = 1:numBlocks % iterate through the 8 blocks, build the trialList
    
        % check the type for the block
        if blockList(i) == 5 % this is a recall block
            
            % add a recall trial to the trial list for every number of trials there are in a block
            for j = 1:trialsPerBlock
               blockedTrialList = [blockedTrialList; 1]; % 1 denotes a recall trial
               
            end

        else % this is recog block
            
            % add a recog trial to the trial list for every number of trials there are in a block
             for j = 1:trialsPerBlock
                
                % we are adding both match/mismatch trials (pseudorandomized order)
                % 2 is a match trial, 3 is a mismatch trial
                recogTrialCtr = recogTrialCtr + 1;
                blockedTrialList = [blockedTrialList; scrambledRecogTrials(recogTrialCtr)];

             end
        end
    end
    


    %--------------------------------------------------------
    % Preallocate some vectors to control experimental params
    %--------------------------------------------------------
    
    % Stimulus parameters:
    stim.tStructure = ones(2,prefs.nTrialsRecall)*-99;
    stim.col = nan(prefs.nTrialsRecall,3); % RGB vals for the to be remembered color

    
    % Response params
    stim.colorOffset = nan(1,prefs.nTrialsRecall);
    stim.recallRT = nan(1,prefs.nTrialsRecall); % good, check how it's being calculated though
    stim.reportedColorDeg = nan(1,prefs.nTrialsRecall); % good
    stim.reportedColorRad = nan(1,prefs.nTrialsRecall); % good

        
    % scramble trial order
    stim.rndInd = randperm(prefs.nTrialsRecall);
    % check this! what is prefs.nTtialsRecall? shou
     
    % create stimulus sequence
    cond.oriChan = prefs.nChans;
    cond.colChan = prefs.nChans;
    [tStructure] = buildTrialStructure(cond,prefs);
    
    stim.tStructure = tStructure(:,stim.rndInd);
    
    stim.oriChan = tStructure(1,stim.rndInd);
    stim.colChan = tStructure(2,stim.rndInd);

    
    % save separate stim structs for each - save the data we have that is
    % the same for both condition
    intermix_stim = stim;
    block_stim = stim;

    data.mouseOnWheel = []; % inialize these
    data.isValidResponse = []; % inialize these

    % create the data structs - these will save the experiment data for
    % each condition separately
    intermix_data = data;
    block_data = data;


    % initalize recall stim
    intermix_stim = initRecallStim(prefs, intermix_stim, window);
    block_stim = initRecallStim(prefs, block_stim, window);


    % Determine how many items there are for each trial, and how long they
    % are shown for (change this! - do not need prefs/fullfactorial design mess)
    retentionInterval = prefs.retentionIntervals(prefs.fullFactorialDesign(prefs.order(1), 2));



    % Put up instructions and wait for keypress.
    instruct(window, mouseCondition);
    returnToFixation(window, prefs.fixationSize);



    % Get rects for each item.
    rects = cell(1, max(prefs.setSizes));
    for i = 1:max(prefs.setSizes)
      rects{i} = circularArrayRects([0, 0, prefs.squareSize, prefs.squareSize], ...
        i, prefs.radius, window.centerX, window.centerY)';
    end



    
    for s = 1:length(conditionSeq) % start condition loop
        currCondition = conditionSeq(s);

        if s == length(conditionSeq)
            % next condition screen
            conditionScreen(window);
        end

        % Keep track of which recall/recog trial we are on - reset this at
        % the beginning of both the intermixed and condition block
        recallCtr = 1; 
        recogCtr = 1;

        if currCondition == 1 % start of intermixed condition
            % disp('Test = this is the intermixed condition!') % check
            blockCtr = 0;
 
            for t = 1:nTrialsTotal % start trial loop
              theCloser;

              if (mod(t-1,trialsPerBlock) == 0)
                  intermixScreen(window);
              end
        
              % get the type for this trial (recall = 1, recog match = 2, recog mismatch = 3)
              trialType = intermixedTrialList(t, 1);
        
              % Pick which one of the three squares to test (this is randomized, not pseudorandomized); for both recall and recog trials
              itemToTest_mix(t) = randsample(prefs.nItems,1); % edited by LK
        
        
        
              % 400 - 500 ms SOA fixation
                jitter = prefs.possTimes(randi(randi([1, numel(prefs.possTimes)])));
                returnToFixation(window, prefs.fixationSize);
                drawFixation(window, prefs.fixationSize); % added, for some reason this was needed too...
            
                Screen('Flip', window.onScreen);
            
                WaitSecs(jitter);
                trial_jitter_intermix{t,1} = jitter; % keep track of jitter
        
        
        
              % Select which colors will be shown - and which one will be tested
              if trialType == 1 % this is a recall trial
        
                    % Pick two random colors (values are degrees/indices), these are the colors that WILL NOT be tested
                  
                    % num_numbers = 3; % Number of random numbers
                    min_value = 1; % Minimum value
                    max_value = 360; % Maximum value
                    % min_distance = 90; % Minimum distance between numbers (22.5 degrees rounded up)
                    
                    % The first number is the number that will be tested (generated earlier)
                    random_numbers = block_stim.color(recallCtr);
                    % disp(block_stim.color(recallCtr))
                    
                    % Generate the remaining random numbers
                    for i = 2:prefs.nItems
                        while true
                            % Generate a new random number
                            new_number = randi([min_value, max_value]);
                            
                            % Check if the new number satisfies the minimum distance condition
                            min_distance = min(abs(random_numbers(1:i-1) - new_number), 360 - abs(random_numbers(1:i-1) - new_number));
                            if min_distance >= 50
                                % If satisfied, add the new number to the array
                                random_numbers(i) = new_number;
                                break;

                            end
                        end
                    end
            
            
                    % Display the array of random numbers
                    % disp(random_numbers);
        
            
                  if itemToTest_mix(t) == 2
                      % the color to be tested should swap places with the color originially in the second index
                      random_numbers([2 1]) = random_numbers([1 2]);
        
                  elseif itemToTest_mix(t) == 3
                      % the color to be tested should swap places with the color originially in the first index
                      random_numbers([3 1]) = random_numbers([1 3]);
              
                  end
        
                  intermix_colorsInDegrees{t} = random_numbers;
                
                  colorsToDisplay = zeros(3,3);
                  colorsToDisplay = prefs.originalcolorwheel(intermix_colorsInDegrees{t}, :)'; % get the colors from the color wheel using the 'angle' (index w/i color wheel)
        
                    % disp(colorsToDisplay(:,:)); % check
              else % this is a recognition trial (either match or mismatxh)
                
                  % Pick three random colors (values are degrees/indices), this INCLUDES the square that WILL be tested          
                  % each color must be at least 15 degrees apart from one another (don't want colors to look exactly the same!)
        
                    % num_numbers = 3; % Number of random numbers
                    min_value = 1; % Minimum value
                    max_value = 360; % Maximum value
                    % min_distance = 90; % Minimum distance between numbers (22.5 degrees rounded up)

                    % Generate the first random number
                    random_numbers = randi([min_value, max_value]);

                    % Generate the remaining random numbers
                    for i = 2:prefs.nItems
                        while true
                            % Generate a new random number
                            new_number = randi([min_value, max_value]);
                            
                            % Check if the new number satisfies the minimum distance condition
                            min_distance = min(abs(random_numbers(1:i-1) - new_number), 360 - abs(random_numbers(1:i-1) - new_number));
                            if min_distance >= 50
                                % If satisfied, add the new number to the array
                                random_numbers(i) = new_number;
                                break;

                            end
                        end
                    end

            
                    intermix_colorsInDegrees{t} = random_numbers;
            
                    % Display the array of random numbers
                    % disp(random_numbers);
        
        
                  colorsToDisplay = zeros(3,3);
        
                  colorsToDisplay = prefs.originalcolorwheel(intermix_colorsInDegrees{t}, :)'; % get the colors from the color wheel using the 'angle' (index w/i color wheel)
        
        
              end
        
        
        
        
        
              % Show the three colored squares
        
              % Draw/post the stimulus and wait.
              drawFixation(window, prefs.fixationSize);
              Screen('FillRect', window.onScreen, colorsToDisplay, rects{prefs.nItems});
              Screen('Flip', window.onScreen);
              WaitSecs(prefs.stimulusDuration);
        
        
              % Remove stimulus, return to blank, wait for delay period to pass.
              returnToFixation(window,prefs.fixationSize);
              WaitSecs(retentionInterval);
        
        
        
        
              % Display the response screen.
        
              if trialType == 1 % probe this trial with recall
                  
                  colorsOfTest = repmat([120 120 120], prefs.nItems, 1); % rgb [120 120 120] is dark gray
                  colorsOfTest(itemToTest_mix(t), :) = [145 145 145]; % rgb [145 145 145] is light gray -- this is for the color to be tested
                  drawFixation(window, prefs.fixationSize);
                  Screen('FillRect', window.onScreen, colorsOfTest', rects{prefs.nItems});
        
        
                  [intermix_data,intermix_stim,prefs] = recallProbe(intermix_data,prefs,intermix_stim,recallCtr, t, window, colorsOfTest, itemToTest_mix, rects, data_name);
        
                  [x,y,buttons] = GetMouse(window.onScreen);
        
                  while any(buttons)
                    [x,y,buttons] = GetMouse(window.onScreen);
                  end

                  trial_colors = intermix_colorsInDegrees{t};
                  intermix_data.recall(recallCtr,6:8) =  trial_colors; % saves the three sqaures colors (in degrees)

                  % increment the recallCtr
                  if recallCtr < length(recallTrials)
                    recallCtr = recallCtr + 1;
                  end

        
              else % probe this trial with recognition (trialType is 2 or 3)
        
                  intermix_data = recogProbe(intermix_data, prefs, trialType, t, window, intermix_colorsInDegrees, colorsToDisplay, itemToTest_mix, rects, mouseCondition, recogCtr, data_name);
        
                  trial_colors = intermix_colorsInDegrees{t};
                  intermix_data.recog(recogCtr,10:12) =  trial_colors; % saves the three squares colors (in degrees)

                  % increment the recogCtr
                  if recogCtr < length(recogTrials)
                    recogCtr = recogCtr + 1;
                  end

        
              end % end of trial condition

              % if (mod(t,trialsPerBlock) == 0 &&  t < nTrialsTotal) % we are at the end of a block, but it's not the last block
              % 
              %   breakBetween(window);
              % 
              % end
        
            end % end of trial loop

            % end of intermixed condition

        else % start of blocked condition
            blockCtr = 0;

            % disp('Test = this is the blocked condition!') % check

            for t = 1:nTrialsTotal % start trial loop
                theCloser;


              if(mod(t-1,trialsPerBlock) == 0)
                blockCtr = blockCtr + 1;

                blockType = blockList(blockCtr);

                if blockType == 5 % recall
                    recallScreen(window);
                else % show recog screen
                    recogScreen(window);
                end
              end
        
              % get the type for this trial (recall = 1, recog match = 2, recog mismatch = 3)
              trialType = blockedTrialList(t, 1);
        
              % Pick which one of the three squares to test (this is randomized, not pseudorandomized); for both recall and recog trials
              itemToTest_block(t) = randsample(prefs.nItems,1); % edited by LK
        
        
        
              % 400 - 500 ms SOA fixation
                jitter = prefs.possTimes(randi(randi([1, numel(prefs.possTimes)])));
                returnToFixation(window, prefs.fixationSize);
                drawFixation(window, prefs.fixationSize); % added, for some reason this was needed too...
            
                Screen('Flip', window.onScreen);
            
                WaitSecs(jitter);
                trial_jitter_block{t,1} = jitter; % keep track of jitter
        
        
        
              % Select which colors will be shown - and which one will be tested
              if trialType == 1 % this is a recall trial
        
                    % Pick two random colors (values are degrees/indices), these are the colors that WILL NOT be tested
                  
                    % num_numbers = 3; % Number of random numbers
                    min_value = 1; % Minimum value
                    max_value = 360; % Maximum value
                    % min_distance = 90; % Minimum distance between numbers (22.5 degrees rounded up)
                    
                    % The first number is the number that will be tested (generated earlier)
                    random_numbers = block_stim.color(recallCtr);
                    % disp(block_stim.color(recallCtr))
                    
                    % Generate the remaining random numbers
                    for i = 2:prefs.nItems
                        while true
                            % Generate a new random number
                            new_number = randi([min_value, max_value]);
                            
                            % Check if the new number satisfies the minimum distance condition
                            min_distance = min(abs(random_numbers(1:i-1) - new_number), 360 - abs(random_numbers(1:i-1) - new_number));
                            if min_distance >= 50
                                % If satisfied, add the new number to the array
                                random_numbers(i) = new_number;
                                break;

                            end
                        end
                    end
            
            
                    % Display the array of random numbers
                    % disp(random_numbers);
        
            
                  if itemToTest_block(t) == 2
                      % the color to be tested should swap places with the color originially in the second index
                      random_numbers([2 1]) = random_numbers([1 2]);
        
                  elseif itemToTest_block(t) == 3
                      % the color to be tested should swap places with the color originially in the first index
                      random_numbers([3 1]) = random_numbers([1 3]);
              
                  end
        
                  block_colorsInDegrees{t} = random_numbers;
                
                  colorsToDisplay = zeros(3,3);
                  colorsToDisplay = prefs.originalcolorwheel(block_colorsInDegrees{t}, :)'; % get the colors from the color wheel using the 'angle' (index w/i color wheel)
        
                    % disp(colorsToDisplay(:,:)); % check


              else % this is a recognition trial (either match or mismatxh)
                
                  % Pick three random colors (values are degrees/indices), this INCLUDES the square that WILL be tested          
                  % each color must be at least 15 degrees apart from one another (don't want colors to look exactly the same!)
        
                    % num_numbers = 3; % Number of random numbers
                    min_value = 1; % Minimum value
                    max_value = 360; % Maximum value
                    % min_distance = 90; % Minimum distance between numbers (22.5 degrees rounded up)

                    % Generate the first random number
                    random_numbers = randi([min_value, max_value]);

                    % Generate the remaining random numbers
                    for i = 2:prefs.nItems
                        while true
                            % Generate a new random number
                            new_number = randi([min_value, max_value]);
                            
                            % Check if the new number satisfies the minimum distance condition
                            min_distance = min(abs(random_numbers(1:i-1) - new_number), 360 - abs(random_numbers(1:i-1) - new_number));
                            if min_distance >= 50
                                % If satisfied, add the new number to the array
                                random_numbers(i) = new_number;
                                break;

                            end
                        end
                    end

            
                    block_colorsInDegrees{t} = random_numbers;
            
                    % Display the array of random numbers
                    % disp(random_numbers);
        
        
                  colorsToDisplay = zeros(3,3);
        
                  colorsToDisplay = prefs.originalcolorwheel(block_colorsInDegrees{t}, :)'; % get the colors from the color wheel using the 'angle' (index w/i color wheel)
        
        
              end
        
        
        
        
        
              % Show the three colored squares
        
              % Draw/post the stimulus and wait.
              drawFixation(window, prefs.fixationSize);
              Screen('FillRect', window.onScreen, colorsToDisplay, rects{prefs.nItems});
              Screen('Flip', window.onScreen);
              WaitSecs(prefs.stimulusDuration);
        
        
              % Remove stimulus, return to blank, wait for delay period to pass.
              returnToFixation(window, prefs.fixationSize);
              WaitSecs(retentionInterval);
        
        
        
        
              % Display the response screen.
        
              if trialType == 1 % probe this trial with recall
                  
                  colorsOfTest = repmat([120 120 120], prefs.nItems, 1); % rgb [120 120 120] is dark gray
                  colorsOfTest(itemToTest_block(t), :) = [145 145 145]; % rgb [145 145 145] is light gray -- this is for the color to be tested
                  drawFixation(window, prefs.fixationSize);
                  Screen('FillRect', window.onScreen, colorsOfTest', rects{prefs.nItems});
        
        
                  [block_data,block_stim,prefs] = recallProbe(block_data,prefs,block_stim,recallCtr, t, window, colorsOfTest, itemToTest_block, rects, data_name);
        
                  [x,y,buttons] = GetMouse(window.onScreen);
        
                  while any(buttons)
                    [x,y,buttons] = GetMouse(window.onScreen);
                  end
                  
                  trial_colors = block_colorsInDegrees{t};
                  block_data.recall(recallCtr,6:8) =  trial_colors; % saves the three sqaures colors (in degrees)

                  % increment the recallCtr
                  if recallCtr < length(recallTrials)
                    recallCtr = recallCtr + 1;
                  end

        
              else % probe this trial with recognition (trialType is 2 or 3)
        
                  block_data = recogProbe(block_data, prefs, trialType, t, window, block_colorsInDegrees, colorsToDisplay, itemToTest_block, rects, mouseCondition, recogCtr, data_name);
                  
                  trial_colors = block_colorsInDegrees{t};
                  block_data.recog(recogCtr,10:12) =  trial_colors; % saves the three squares colors (in degrees)

                  % increment the recogCtr
                  if recogCtr < length(recogTrials)
                    recogCtr = recogCtr + 1;
                  end
        
              end

              % if (mod(t,trialsPerBlock) == 0 &&  t < nTrialsTotal) % we are at the end of a block, but it's not the last block
              % 
              %   breakBetween(window);
              % 
              % end
        
            end % end of trial loop
            
   
        end % end of blocked condition

    end % end of condition loop

    endScreen(window);

    % save raw variables
    save(output_name);

    postpareEnvironment;

    
    % analyze data


    % separate data out
    recogIntermixData = intermix_data.recog;
    recogBlockData = block_data.recog;

    recallIntermixData = intermix_data.recall;
    recallBlockData = block_data.recall;

    % recog - remove trials with no response
    recogIntermixData = recogIntermixData(cat(1,recogIntermixData(:,3)) > -1,:);
    recogBlockData = recogBlockData(cat(1,recogBlockData(:,3)) > -1,:);

    % recall - remove trials with invalid response
    recallIntermixData = recallIntermixData(cat(1,recallIntermixData(:,4)) < 500,:);
    recallBlockData = recallBlockData(cat(1,recallBlockData(:,4)) < 500,:);

    % recognition - calculate overall accuracies for intermixed vs. blocked

    % intermixed
    match_recogIntermixAcc = mean(recogIntermixData(recogIntermixData(:,2) ==2,3));
    mismatch_recogIntermixAcc = mean(recogIntermixData(recogIntermixData(:,2) ==3,3));
    all_recogIntermixAcc = mean(recogIntermixData(:,3));

    % blocked
    match_recogBlockAcc = mean(recogBlockData(recogBlockData(:,2) ==2,3));
    mismatch_recogBlockAcc = mean(recogBlockData(recogBlockData(:,2) ==3,3));
    all_recogBlockAcc = mean(recogBlockData(:,3));

    % recognition - reaction times for intermixed vs. blocked

    % intermixed
    match_recogIntermixRT = median(recogIntermixData(recogIntermixData(:,2) ==2,4));
    mismatch_recogIntermixRT = median(recogIntermixData(recogIntermixData(:,2) ==3,4));
    all_recogIntermixRT = median(recogIntermixData(:,4));

    corMatch_recogIntermixRT = median(recogIntermixData(recogIntermixData(:,2)==2 & recogIntermixData(:,3)==1,4));
    corMismatch_recogIntermixRT = median(recogIntermixData(recogIntermixData(:,2)==3 & recogIntermixData(:,3)==1,4));
    corAll__recogIntermixRT = median(recogIntermixData(recogIntermixData(:,3)==1,4));

    % blocked
    match_recogBlockRT = median(recogBlockData(recogBlockData(:,2) ==2,4));
    mismatch_recogBlockRT = median(recogBlockData(recogBlockData(:,2) ==3,4));
    all_recogBlockRT = median(recogBlockData(:,4));
    
    corMatch_recogBlockRT = median(recogBlockData(recogBlockData(:,2)==2 & recogBlockData(:,3)==1,4));
    corMismatch_recogBlockRT = median(recogBlockData(recogBlockData(:,2)==3 & recogBlockData(:,3)==1,4));
    corAll__recogBlockRT = median(recogBlockData(recogBlockData(:,3)==1,4));


    % recall trials - get average and standard deviation of error

    % intermixed
    aveError_recallIntermixed = mean(abs(recallIntermixData(:,4)));
    aveError_recallBlocked = mean(abs(recallBlockData(:,4)));

    stdError_recallIntermixed = std(abs(recallIntermixData(:,4)));
    stdError_recallBlocked = std(abs(recallBlockData(:,4)));


    
    % get number of hits for intermixed
    intermixed_hits = height(recogIntermixData(recogIntermixData(:,2) == 2 & recogIntermixData(:,3) == 1));
    intermixed_hit_rate = intermixed_hits/height(recogIntermixData(recogIntermixData(:,2) == 2));

    % get number of FAs for intermixed
    intermixed_fas = height(recogIntermixData(recogIntermixData(:,2) == 3 & recogIntermixData(:,3) == 0));
    intermixed_fa_rate = intermixed_fas/height(recogIntermixData(recogIntermixData(:,2) == 3));

    % get the K for intermixed
    K_intermixed = 3*(intermixed_hit_rate - intermixed_fa_rate);


    % get number of hits for blocked
    blocked_hits = height(recogBlockData(recogBlockData(:,2) == 2 & recogBlockData(:,3) == 1));
    blocked_hit_rate = blocked_hits/height(recogBlockData(recogBlockData(:,2) == 2));

    % get number of FAs for blocked
    blocked_fas = height(recogBlockData(recogBlockData(:,2) == 3 & recogBlockData(:,3) == 0));
    blocked_fa_rate = blocked_fas/height(recogBlockData(recogBlockData(:,2) == 3));

    % get the K for blocked
    K_blocked = 3*(blocked_hit_rate - blocked_fa_rate);



    % store this info in final data
    %block - recog
    final_data(1,1) = "Recog Acc (match - block)"; final_data(1,2) = match_recogBlockAcc;
    final_data(2,1) = "Recog Acc (mismatch - block)"; final_data(2,2) = mismatch_recogBlockAcc;
    final_data(3,1) = "Recog Acc (all - block)"; final_data(3,2) = all_recogBlockAcc;
    final_data(4,1) = "Recog Correct RT (match - block)"; final_data(4,2) = corMatch_recogBlockRT;
    final_data(5,1) = "Recog Correct RT (mismatch - block)"; final_data(5,2) = corMismatch_recogBlockRT;
    final_data(6,1) = "Recog Correct RT (all - block)"; final_data(6,2) = corAll__recogBlockRT;
    final_data(7,1) = "Recog RT (match - block"; final_data(7,2) = match_recogBlockRT;
    final_data(8,1) = "Recog RT (mismatch - block)"; final_data(8,2) = mismatch_recogBlockRT;
    final_data(9,1) = "Recog RT (all - block)"; final_data(9,2) = all_recogBlockRT; 
    %block - recall
    final_data(10,1) = "Recall Error - Average (block)"; final_data(10,2) = aveError_recallBlocked;
    final_data(11,1) = "Recall Error - Standard Deviation (block)"; final_data(11,2) = stdError_recallBlocked;

    %intermix - recog
    final_data(12,1) = "Recog Acc (match - intermixed)"; final_data(12,2) = match_recogIntermixAcc;
    final_data(13,1) = "Recog Acc (mismatch - intermixed)"; final_data(13,2) = mismatch_recogIntermixAcc;
    final_data(14,1) = "Recog Acc (all - intermixed)"; final_data(14,2) = all_recogIntermixAcc;
    final_data(15,1) = "Recog Correct RT (match - intermixed)"; final_data(15,2) = corMatch_recogIntermixRT;
    final_data(16,1) = "Recog Correct RT (mismatch - intermixed)"; final_data(16,2) = corMismatch_recogIntermixRT;
    final_data(17,1) = "Recog Correct RT (all - intermixed)"; final_data(17,2) = corAll__recogIntermixRT;
    final_data(18,1) = "Recog RT (match - intermixed"; final_data(18,2) = match_recogIntermixRT;
    final_data(19,1) = "Recog RT (mismatch - intermixed)"; final_data(19,2) = mismatch_recogIntermixRT;
    final_data(20,1) = "Recog RT (all - intermixed)"; final_data(20,2) = all_recogIntermixRT; 
    %intermix - recall
    final_data(21,1) = "Recall Error - Average (intermixed) "; final_data(21,2) = aveError_recallIntermixed;
    final_data(22,1) = "Recall Error - Standard Deviation (intermixed)"; final_data(22,2) = stdError_recallIntermixed;

    %block - K calculations
    final_data(23,1) = "Recog hit rate (block)"; final_data(23,2) = blocked_hit_rate;
    final_data(24,1) = "Recog false alarm rate (block)"; final_data(24,2) = blocked_fa_rate;
    final_data(25,1) = "Recog K (block)"; final_data(25,2) = K_blocked;

    %intermix - K calculations
    final_data(26,1) = "Recog hit rate (intermix)"; final_data(26,2) = intermixed_hit_rate;
    final_data(27,1) = "Recog false alarm rate (intermix)"; final_data(27,2) = intermixed_fa_rate;
    final_data(28,1) = "Recog K (intermix)"; final_data(28,2) = K_intermixed;

    excel_data = final_data'; %flip easier to copy into excel


    save(ana_name,'excel_data','final_data','corAll__recogBlockRT','match_recogIntermixAcc','mismatch_recogIntermixAcc','all_recogIntermixAcc','match_recogBlockAcc','mismatch_recogBlockAcc','all_recogBlockAcc','match_recogIntermixRT','mismatch_recogIntermixRT','all_recogIntermixRT','corMatch_recogIntermixRT','corMismatch_recogIntermixRT','corAll__recogIntermixRT','match_recogBlockRT','mismatch_recogBlockRT','all_recogBlockRT','corMatch_recogBlockRT','corMismatch_recogBlockRT','corAll__recogIntermixRT','aveError_recallIntermixed','aveError_recallBlocked','stdError_recallIntermixed','stdError_recallBlocked');
    save(output_name)


  catch
    postpareEnvironment;
    psychrethrow(psychlasterror);

  end % end try/catch

 


% function definitions

function offsets = circularArrayOffsets(n, radius, rotation)
  degreeStep = 360/n;
  offsets = [sind(0:degreeStep:(360-degreeStep) + rotation)'.* radius, ...
             cosd(0:degreeStep:(360-degreeStep) + rotation)'.* radius];
end

function rects = circularArrayRects(rect, nItems, radius, centerX, centerY)
  coor = circularArrayOffsets(nItems, radius, 0) + repmat([centerX centerY], nItems, 1);
  rects = [coor(:, 1)-rect(3)/2 , coor(:, 2)-rect(3)/2, coor(:, 1)+rect(3)/2, coor(:, 2)+rect(3)/2];
end

function drawColorWheel(window, prefs, stim, t)
  rind = squeeze(stim.rimLocs(t,:));

  for i = 1:360 % draw the color wheel probe
    color = prefs.colorwheel;
    Screen('FrameArc',window.onScreen,color(i,:),rind,(i-1)*1+90,prefs.penArc,prefs.rimThick,prefs.rimThick);
  end

end

function data = recogProbe(data, prefs, trialType, t, window, colorsInDegrees, colorsToDisplay, itemToTest, rects, mouseCondition, recogCtr, data_name)

    if trialType == 2 % match trial

        colorsOfTest = colorsToDisplay'; % show the same colors again
        % each row represents a color in colorsOfTest
        % each column represents a color in colorsToDisplay 
        % ' operator accounts for the difference

        trialColorInDegrees = colorsInDegrees{1,t}; % select the color index/degree for the square that will be tested
        targetColorInDegrees = trialColorInDegrees(itemToTest(t));

        % disp(colorsOfTest);

        % make the other two squares blend in w/ the background
        allIndicesArr = 1:3;
        nonTestedColors = allIndicesArr(allIndicesArr ~= itemToTest(t)); % get row #s that aren't tested
        colorsOfTest(nonTestedColors, :) = [127.5 127.5 127.5; 127.5 127.5 127.5]; % % matches the 127.5 intensity of background?  
        
    else % mismatch trial - target square will differ in color by 180 degrees

        colorsOfTest = colorsToDisplay';
        % disp(colorsOfTest)
        
        trialColorInDegrees = colorsInDegrees{1,t}; % select the color index/degree for the square that will be tested
        targetColorInDegrees = trialColorInDegrees(itemToTest(t));
        % disp(targetColorInDegrees)

        % offset the color by 180 degrees
        diffColorInDegrees = targetColorInDegrees - 180;
        % disp(diffColorInDegrees)
        if diffColorInDegrees < 0 % color indices/degree calues cannot be negative, correct by one full rotation if it is
            diffColorInDegrees = diffColorInDegrees + 360;
        elseif diffColorInDegrees == 0 % we cannot have 180 degrees - 180 degrees because color wheel is 1:360
            diffColorInDegrees = 360;
        end

        diffColor = prefs.originalcolorwheel(diffColorInDegrees, :); % 


        % store this new color
        colorsOfTest(itemToTest(t),:) = diffColor;

        % make the other two squares blend in w/ the background
        allIndicesArr = 1:3;
        nonTestedColors = allIndicesArr(allIndicesArr ~= itemToTest(t)); % get row #s that aren't tested
        colorsOfTest(nonTestedColors, :) = [window.bcolor; window.bcolor]; % match the non-tested squares to background
        % the reason we make them blend in and not delete them is because it would mess w/ the location of the shown probe square

    end
    
    % show the test/probe array
    WaitSecs(0.9);
    drawFixation(window, prefs.fixationSize);
    Screen('FillRect', window.onScreen, colorsOfTest', rects{prefs.nItems});
    Screen('Flip', window.onScreen);


    % collect a response - limited to 3000 ms
    startTime = GetSecs();
    MousePress = 0;

    while ((GetSecs() - startTime < 3) && ~MousePress)
    
        [~, ~, buttons] = GetMouse();
        MousePress = any(buttons);

        if(buttons(1) || buttons(2)) % left or right mouse was clicked
            RT = GetSecs() - startTime;
            % disp(RT);
            break;
        else
            RT = -1; % junk RT
        end
    end

    % what was pressed?
    if buttons(1) % clicked left mouse

        if mouseCondition == 1
            click = 2; % said it was match
        else % mouseCondition = 0
            click = 3; % said it was mismatxh
        end

    elseif buttons(2) % clicked right mouse

        if mouseCondition == 1
            click = 3; % said it was mismatch
        else
            click = 2; % said it was match
        end

    else
        click = 0; % no response
    end
    
    % was this response correct?
    if click == trialType % correct response!
        cor = 1;
    elseif click > 0 % responded, but incorrect
        cor = 0;
        % disp("incorrect!")
    else
        cor = -4; % otherwise, no response
        % disp("no response!")
    end

    % store the data
    data.recog(recogCtr, 1) = recogCtr; % trial number w/i condition
    data.recog(recogCtr, 2) = trialType; % 2 = match, 3 = mismatch trial
    data.recog(recogCtr, 3) = cor; % correct = 1 , incorrect = 0; no response = -4
    data.recog(recogCtr, 4) = RT; % RT in seconds
    data.recog(recogCtr, 5) = mouseCondition; % mouse condition (0 or 1)
    data.recog(recogCtr, 6) = click; % participant's judgement (2 = match, 3 = mismatch, 0 = no response)
    data.recog(recogCtr, 7) = buttons(1); % left mouse clicked?
    data.recog(recogCtr, 8) = buttons(2); % right mouse clicked?
    data.recog(recogCtr, 9) = itemToTest(t); % was square 1, 2, or 3 probed? (randomized)
    data.recog(recogCtr, 13) = targetColorInDegrees; % probe color

    if trialType == 3
        data.recog(recogCtr, 14) = diffColorInDegrees; % if it is a mismatch trial, this is actually the probe color
    end

    % save the data, preferences, stim info
    save(data_name, 'data','prefs','colorsOfTest');

    % end of trial

end


function  [data, stim, prefs] = recallProbe(data, prefs, stim, recallCtr, t, window, colorsOfTest, itemToTest, rects, data_name) % presents a rotated color wheel, gets response, calculates offset
    
    HideCursor;

    % allow the color wheel to rotate
    prefs.colorWheelLocations = [cosd([1:360]).*prefs.colorWheelRadius; ...
        sind([1:360]).*prefs.colorWheelRadius];
    prefs.ind(recallCtr,1) = ceil(randsample(360,1));  % setting prefs.ind == 0 would be NO change in color wheel position.
    ind = wshift('1d',1:360,prefs.ind(recallCtr,1));
    prefs.colorwheel = prefs.originalcolorwheel(ind,:);

    % draw the color wheel
    drawColorWheel(window, prefs, stim, recallCtr);
    
    % set the mouse to the center of the screen
    SetMouse(window.centerX*2,window.centerY*2,window.onScreen); % added by LK, better...  
    % [badX,badY,buttons] = GetMouse(window.onScreen); % this was here originally, I commented out - LK
    ShowCursor('Arrow');
    rtStart = GetSecs;

    [~,~,buttons] = GetMouse(window.onScreen); % added by LK

    % Collect a response for the recall
    data.mouseOnWheel(recallCtr) = false; % is the mouse w/i bounds (touching the color wheel)?
    data.isValidResponse(recallCtr) = false; % did the participant CLICK on a part of the color wheel? - LK
    while (~any(buttons))
        
        % draw the color wheel
        drawColorWheel(window, prefs, stim, recallCtr)

        [x,y,buttons] = GetMouse(window.onScreen); % get the location of the mouse

        % get the index for the color we are closest to
        % the wheel is not continuous, it is made up of 360 colors (we can't be between colors)
        [~, minDistanceIndex] = min(sqrt((prefs.colorWheelLocations(1, :)+window.centerX - x).^2 + (prefs.colorWheelLocations(2, :)+ window.centerY - y).^2));
        mouseDistanceFromCenter = sqrt((abs(window.centerX - x))^2 + ((abs(window.centerY - y))^2)); % how far away is the mouse from the center?

        
        % is the mouse within bounds (hovering over the color wheel)?
        if (mouseDistanceFromCenter > (prefs.colorWheelRadius - 1) && mouseDistanceFromCenter < (prefs.colorWheelRadius + prefs.rimThick))
          data.mouseOnWheel(recallCtr) = true;
        else
          data.mouseOnWheel(recallCtr) = false; % added this - LK (if mouse is not on wheel, we should not see a colored square)
        end

        % if the mouse is within bounds, make sure the target square shows which color they are hovering over
        if(data.mouseOnWheel(recallCtr)) % make sure mouse is ON the color wheel for the target square to change color
          colorsOfTest(itemToTest(t), :) = prefs.colorwheel(minDistanceIndex,:); % this will change the color of the original dark gray square
        else
          colorsOfTest(itemToTest(t), :) = [145 145 145]; % else, light gray - the mouse is out of bounds
        end
        
        % draw the fixation and three squares (two light gray, target square is either dark gray (default) or the color that the mouse is hovering over)
        drawFixation(window, prefs.fixationSize);
        Screen('FillRect', window.onScreen, colorsOfTest', rects{prefs.nItems});
        
        % draw the color wheel - this ensures it stays up
        drawColorWheel(window, prefs, stim, recallCtr)
        Screen('Flip', window.onScreen);

    end


  while any(buttons) % wait for release
    [~,~,buttons] = GetMouse(window.onScreen);
    if data.mouseOnWheel(recallCtr) % if the final click is on the color wheel, then the participant's answer for this trial is valid
        data.isValidResponse(recallCtr) = true;
    end
  end

  rtEnd = GetSecs;
  HideCursor;

    
    
    % get color offset
    if(data.isValidResponse(recallCtr)) % if the response is valid, calculate the offset

        if minDistanceIndex+prefs.ind(recallCtr)>360
            stim.reportedColor(recallCtr) = minDistanceIndex+prefs.ind(recallCtr)-360;
        else
            stim.reportedColor(recallCtr) = minDistanceIndex+prefs.ind(recallCtr);
        end
        
        stim.recallRT(recallCtr) = rtEnd-rtStart;
        
        %Gives actual reported color as number of degrees (already been
        %corrected for color wheel shift).
        stim.reportedColorDeg(recallCtr) = stim.reportedColor(recallCtr);
        stim.reportedColorRad(recallCtr) = deg2rad(stim.reportedColor(recallCtr));
        
        colOffset = stim.reportedColorDeg(recallCtr)- stim.color(recallCtr);
        
        if colOffset > 180
            colOffset = colOffset-360;
        elseif colOffset < -180
            colOffset = colOffset+360;
        end
        
        % Calculate the color offset!
        stim.colorOffset(recallCtr) = colOffset;

        % save the offset
        data.recall(recallCtr, 4) = stim.colorOffset(recallCtr); % color offset in degrees
        
        % save the RT
        data.recall(recallCtr, 5) = stim.recallRT(recallCtr);
        
        
        % this is a check, what is our color offset/degree of error? should always work, especially w/ rotating the wheel - LK
            % text1 = num2str(stim.colorOffset(recallCtr));
            % tCenter1 = [window.centerX-RectWidth(Screen('TextBounds', window.onScreen, text1))/2 window.centerY-120];
            % 
            % Screen('FillRect',window.onScreen,window.bcolor); % added by LK, make background gray
            % 
            % Screen('DrawText', window.onScreen, text1, tCenter1(1), tCenter1(2),stim.col(recallCtr,:), [0 0 0]);
            % Screen('DrawingFinished', window.onScreen);
            % Screen('Flip', window.onScreen);
            % 
            % WaitSecs(.5);

    else % invalid trial
        data.recall(recallCtr, 4) = 500; % color offset - junk value if invalid trial
        data.recall(recallCtr, 5) = - 1; % RT - junk value if invalid trial
    end

    % store the data
    data.recall(recallCtr, 1) = recallCtr; % trial number within condition
    data.recall(recallCtr, 2) = data.isValidResponse(recallCtr); % was this participant's response valid?
    data.recall(recallCtr, 3) = data.mouseOnWheel(recallCtr); % did the participant end on the wheel?
    % color offset & RT saved earlier because only accurately calculated if the response is valid
    data.recall(recallCtr,9) = stim.color(recallCtr); % store the probe color in degrees
    data.recall(recallCtr,10) = stim.reportedColorDeg(recallCtr); % what color did they click on in degrees?

    if(data.recall(recallCtr, 4) == 500)
        data.recall(recallCtr,10) = 500; % what color did they click on in degrees? - junk value if invalid trial
    end

    % save the data, the stim, and the preferences
    save(data_name, 'data','stim','prefs');

    % end of trial
end

function drawFixation(window, fixationSize)
  Screen('DrawDots', window.onScreen, [window.centerX, window.centerY], fixationSize, 255);
end

function prefs = getPreferences(prefs)
  prefs.nTrialsPerCondition = 200;
  prefs.setSizes = 3; % modified by LK: we are only using a set size of 3 - get rid of this....
  prefs.retentionIntervals = 0.9; % modified by LK: 900 ms delay period
  prefs.stimulusDuration = .1; % modified by LK: 100 ms encoding period
  prefs.possTimes = 0.4:0.001:0.5; % possible times for SOA
  prefs.squareSize = 90; % size of each stimulus object, in pixels (originally 75 - LK)
  prefs.radius = 180; % what is this? - LK
  prefs.fixationSize = 10; % originally 3 - LK
  prefs.nItems = 3;

  % Colorwheel details.
  prefs.colorWheelRadius = 400; % determines size of color wheel - LK


   % Dimensions for Stims
    prefs.rimSize = 1000; % how big is the color wheel?? - LK
    prefs.rimThick = 100; % how thick is the rim of the color wheel? - LK
    prefs.handX = prefs.rimThick;
    prefs.handY = prefs.rimSize-2;
    prefs.faceSize = prefs.rimThick*2;
    prefs.penArc = 2;
    
    % FEM
    prefs.nChans = 4;
    
    
    % % Stimulus Orientations for FEM -  - finite element method? (LK)
    
    % Stimulus positions for FEM
    prefs.arcChannels = linspace(0,360-360/prefs.nChans,prefs.nChans);
    prefs.arcChanVar = (prefs.arcChannels(2)-prefs.arcChannels(1))/2;
    
    % Stimulus Colors for FEM
    prefs.colChannels = linspace(1,361-360/prefs.nChans,prefs.nChans);
    prefs.colChanVar = (prefs.colChannels(2)-prefs.colChannels(1))/2;
    prefs.colChanVar = prefs.colChanVar-0.5; % correction because you can't be halfway between colors
    

  % Loads in 360 different colors for color wheel - old
  % prefs.originalcolorwheel = load('colorwheel360.mat', 'fullcolormatrix'); % columns: R, G, B values, respectively
  % prefs.originalcolorwheel = prefs.originalcolorwheel.fullcolormatrix;

  % Loads in 360 different colors for the color wheel (these were converted from CIE L*a*b* colour space)
  % Circle was centered in the colour space at (L = 70, a = 20, b = 38) with a radius of 60
  % All colors have equal luminance
  prefs.originalcolorwheel = load( 'originalcolorwheelCIE.mat', 'fullcolormatrix'); % columns: R, G, B values, respectively
  prefs.originalcolorwheel = prefs.originalcolorwheel.fullcolormatrix;

  % Randomize trial order of full factorial design order

  % get this of this stuff...
  prefs.fullFactorialDesign = fullfact([length(prefs.setSizes), ...
    length(prefs.retentionIntervals), ...
    prefs.nTrialsPerCondition]);

  prefs.order = Shuffle(1:length(prefs.fullFactorialDesign));
  prefs.nTrials = length(prefs.fullFactorialDesign);
  prefs.nTrialsRecall = 100; 

end


% This function is from Conjunctions_Exp4.m
function stim = initRecallStim(prefs, stim, window)
    
    for t = 1:prefs.nTrialsRecall % start trial loop

        % DEFINE LOCATIONS OF STIMULI
        [discMat,faceMat,handMat,recallLoc,probeLoc] = createRims(prefs,window.centerX,window.centerY,1,1);
        stim.rimLocs(t,:) = discMat; stim.faceLocs(t,:) = faceMat; stim.handLocs(t,:) = handMat;
        stim.recallLoc(t,:) = recallLoc; stim.probeLoc(t) = probeLoc;

        % DEFINE POSITION OF STIMULI
        tmpArc = prefs.arcChannels(stim.oriChan(t));
        tmpArcOff = randsample(0:prefs.arcChanVar-1,1)*randsample([-1,1],1);
        stim.arc(1,t) = tmpArc+tmpArcOff;
        stim.arc(stim.arc<0) = stim.arc(stim.arc<0)+360;
                
        % DEFINE COLOR OF STIMULI - chooses the color for the 'to be remembered' square - LK
        tmpCol = prefs.colChannels(stim.colChan(t));
        tmpColOff = randsample(0:prefs.colChanVar,1)*randsample([-1 1],1);
        stim.color(1,t) = tmpCol+tmpColOff;
        stim.color(stim.color<1) = stim.color(stim.color<1)+360; % add 360 if <= 0
        stim.col(t,:) = prefs.originalcolorwheel(stim.color(1,t),:);  % retrieve info from color wheel matrix

        % stim.color are the indices/degrees for the to-be-remembered colors
        % stim.col are the RGB vals for the to-be-remembered colors

    end % end trial loop
end

function instruct(window, mouseCondition)
    spacebar = KbName('space');

    Screen('FillRect',window.onScreen,window.bcolor); % added by LK, make background gray
    Screen('TextSize', window.onScreen, window.fontsize);
    Screen('Textcolor',window.onScreen, window.white);
    
    DrawFormattedText(window.onScreen, 'Color Recall & Recognition', 'center', window.centerY - 200);
    
    DrawFormattedText(window.onScreen, 'At the start of each trial, you will be shown three squares. Remember their colors.', 'center', window.centerY - 100);
    DrawFormattedText(window.onScreen, 'After a delay, you will be shown colored squares again, or a color wheel.', 'center', window.centerY - 50);
    
    
    % depending on how the mouse is for this participant, give different instructions
    if mouseCondition == 1 % left is match, right is mismatch judgement
        DrawFormattedText(window.onScreen, 'If you are shown colored squares again, indicate whether they are the SAME (LEFT click)', 'center', window.centerY + 50);
        DrawFormattedText(window.onScreen, 'or DIFFERENT (RIGHT click) as before.', 'center', window.centerY + 100);
    else % right is match, left is mismatch judgement
        DrawFormattedText(window.onScreen, 'If you are shown colored squares again, indicate whether they are the SAME (RIGHT click)', 'center', window.centerY + 50);
        DrawFormattedText(window.onScreen, 'or DIFFERENT (LEFT click) as before.', 'center', window.centerY + 100);
    end



    DrawFormattedText(window.onScreen, 'If you are shown a color wheel, select the color of the lightest square.', 'center', window.centerY + 150);
    
    DrawFormattedText(window.onScreen, 'Press space to begin.', 'center', window.centerY + 250);
    
    Screen('Flip', window.onScreen); % show the instructions
    
    [~,~, keyCode] = KbCheck(); % wait for spacebar to start
     % [~,~, keyCode] = KbCheck; % wait for spacebar to start

    while ~keyCode(spacebar)
        [~,~,keyCode] = KbCheck();
        % [~,~,keyCode] = KbCheck;

        keyCode(spacebar);
    end
    WaitSecs(1);

end
% 
% function breakBetween(window)
%     spacebar = KbName('space');
% 
%     Screen('FillRect',window.onScreen,window.bcolor);
%     Screen('TextSize', window.onScreen, window.fontsize);
%     Screen('Textcolor',window.onScreen, window.white);
% 
%     DrawFormattedText(window.onScreen,'Please take a break. Press space to continue the task.', 'center',window.centerY);
% 
%     Screen('Flip', window.onScreen); % show the text
% 
%     [~,~, keyCode] = KbCheck(); % wait for spacebar to continue
%     while ~keyCode(spacebar)
%         [~,~,keyCode] = KbCheck();
%         keyCode(spacebar);
%     end
% 
% end

function recogScreen(window)
    spacebar = KbName('space');

    Screen('FillRect',window.onScreen,window.bcolor);
    Screen('TextSize', window.onScreen, window.fontsize);
    Screen('Textcolor',window.onScreen, window.white);
    
    DrawFormattedText(window.onScreen,'MATCHING BLOCK', 'center',window.centerY);
    DrawFormattedText(window.onScreen,'Press space to continue when you are ready.', 'center',window.centerY+50);

    
    Screen('Flip', window.onScreen); % show the text
    % WaitSecs(1); % temp
    
    [~,~, keyCode] = KbCheck(); % wait for spacebar to continue
    while ~keyCode(spacebar)
        [~,~,keyCode] = KbCheck();
        keyCode(spacebar);
    end

end

function recallScreen(window)
    spacebar = KbName('space');

    Screen('FillRect',window.onScreen,window.bcolor);
    Screen('TextSize', window.onScreen, window.fontsize);
    Screen('Textcolor',window.onScreen, window.white);
    
    DrawFormattedText(window.onScreen,'COLOR WHEEL BLOCK', 'center',window.centerY);
    DrawFormattedText(window.onScreen,'Press space to continue when you are ready.', 'center',window.centerY+50);

    
    Screen('Flip', window.onScreen); % show the text
    % WaitSecs(1); % temp
    
    [~,~, keyCode] = KbCheck(); % wait for spacebar to continue
    while ~keyCode(spacebar)
        [~,~,keyCode] = KbCheck();
        keyCode(spacebar);
    end

end

function intermixScreen(window)
    spacebar = KbName('space');

    Screen('FillRect',window.onScreen,window.bcolor);
    Screen('TextSize', window.onScreen, window.fontsize);
    Screen('Textcolor',window.onScreen, window.white);
    
    DrawFormattedText(window.onScreen,'COLOR WHEEL AND MATCHING', 'center',window.centerY);
    DrawFormattedText(window.onScreen,'Press space to continue when you are ready.', 'center',window.centerY+50);

    
    Screen('Flip', window.onScreen); % show the text
    % WaitSecs(1); % temp
    
    [~,~, keyCode] = KbCheck(); % wait for spacebar to continue
    while ~keyCode(spacebar)
        [~,~,keyCode] = KbCheck();
        keyCode(spacebar);
    end

end

function endScreen(window)
    spacebar = KbName('space');

    Screen('FillRect',window.onScreen,window.bcolor);
    Screen('TextSize', window.onScreen, window.fontsize);
    Screen('Textcolor',window.onScreen, window.white);
    
    DrawFormattedText(window.onScreen,'Please ring the bell.', 'center',window.centerY);
    
    Screen('Flip', window.onScreen); % show the text
    
    [~,~, keyCode] = KbCheck(); % wait for spacebar to continue
    while ~keyCode(spacebar)
        [~,~,keyCode] = KbCheck();
        keyCode(spacebar);
    end

end

function conditionScreen(window)
    spacebar = KbName('space');

    Screen('FillRect',window.onScreen,window.bcolor);
    Screen('TextSize', window.onScreen, window.fontsize);
    Screen('Textcolor',window.onScreen, window.white);
    
    DrawFormattedText(window.onScreen,'Please ring the bell.', 'center',window.centerY - 50);
    DrawFormattedText(window.onScreen,'Halfway point.', 'center',window.centerY + 50);

    
    Screen('Flip', window.onScreen); % show the text
    
    [~,~, keyCode] = KbCheck(); % wait for spacebar to continue
    while ~keyCode(spacebar)
        [~,~,keyCode] = KbCheck();
        keyCode(spacebar);
    end

    WaitSecs(1); % another screen will be after this, prevent just clicking through that one

end

function [window,rect] = openWindow()
  Screen('Preference', 'SkipSyncTests', 1);

  % added by LK to fix setMouse/getMouse for retina displays (Mac OS)
  PsychImaging('PrepareConfiguration');
  PsychImaging('AddTask','General','UseRetinaResolution')
  [window.onScreen, rect] = PsychImaging('OpenWindow', 0,0, []);

  Screen('BlendFunction', window.onScreen, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
  [window.screenX, window.screenY] = Screen('WindowSize', window.onScreen); % check resolution
  window.screenRect  = [0, 0, window.screenX, window.screenY]; % screen rect
  window.centerX = window.screenX * 0.5; % center of screen in X direction
  window.centerY = window.screenY * 0.5; % center of screen in Y direction
  window.centerXL = floor(mean([0, window.centerX])); % center of left half of screen in X direction
  window.centerXR = floor(mean([window.centerX, window.screenX])); % center of right half of screen in X direction

  % Basic drawing and screen variables.
  window.black    = BlackIndex(window.onScreen);
  window.white    = WhiteIndex(window.onScreen);
  window.gray     = mean([window.black window.white]);
  window.fontsize = 45;

  % window.bcolor   = window.gray; % have this commmented out so the 'invisible squares' in recog blend in
  % issue with having the intensity value is that it's not RGB
  % if I change how the locations for probe squares are show on recog, then
  % I can add this back in b/c I won't need to draw gray squares

  window.bcolor = [127.5 127.5 127.5]; % gray

end

% This function is from Conjunctions_Exp4.m
function postpareEnvironment
  ShowCursor();
  ListenChar(0);
  Screen('CloseAll');
end

function returnToFixation(window, fixationSize)
  Screen('FillRect', window.onScreen, window.bcolor);
  Screen('DrawDots', window.onScreen, [window.centerX, window.centerY], fixationSize, 255);
  Screen('Flip', window.onScreen);
end

function prepareEnvironment
  sca; % Screen('CloseAll')
  clc; % clear command window
  clearvars;
  HideCursor;

  commandwindow; % Select the command window to avoid typing in open scripts

  % Seed the random number generator.
  RandStream.setGlobalStream(RandStream('mt19937ar', 'seed', sum(100*clock)));

  ListenChar(2); % Don't print key presses to MATLAB command window

end

% This function is from Conjunctions_Exp4.m
function [struct] = buildTrialStructure(cond,prefs)

    % Basically, this requires there to always be an even number of recall trials
    s = fieldnames(cond);
    struct = [];
    reps = 1;
    for i=1:length(s)
        % curVal = getfield(cond,s{i}); % original value was 4...
        curVal = 2;

        c = [];
        
        if iswhole(prefs.nTrialsRecall/(curVal*reps)) == 0 % total number of trials must be divisable by 2
            Screen('CloseAll');
            msgbox('Invalid Trial Number', 'modal')
            return;
        end
        
        for r = 1:reps
            for t = 1:curVal
                c = [c, repmat(t,1,prefs.nTrialsRecall/(curVal*reps))];
            end
        end
        struct = [struct; c];
        reps = reps*curVal;
    end
end

% This function is from Conjunctions_Exp4.m
function [d,ndx,varout] = iswhole(varargin)
    % ISWHOLE True for integers(whole numbers).
    %     ISWHOLE(X) is 1 for the elements of X that are integers, 0 otherwise.
    %     ISWHOLE(X1,X2,..,XN) returns a 1-by-N array with 1 for integers and 0
    %     otherwise.
    % ISWHOLE does not check for integer data type as does ISINTEGER.
    
    % Mukhtar Ullah
    % mukhtar.ullah@informatik.uni-rostock.de
    % University of Rostock
    % November 15, 2004
    %%%%%%%%%%%%%%%%%%%%%
    % Adapted SLIGHTLY from the original version of iswhole() by Mukhtar
    % Josef Lotz
    % Portland State University
    % February 276, 2007
    
    C = varargin;
    nIn = nargin;
    
    if nIn < 2
        d = C{:} == round(C{:});
    else
        d = zeros(1,nIn);
        for i = 1:nIn
            d(i) = isequal(C{i},round(C{i}));
        end
    end
    
    % Line added to retrieve the indices of the whole number
    ndx = find(d);
    temp = [C{:}];
    % As well as a vector of whole numbers found
    varout = temp(ndx);
end

% This function is from Conjunctions_Exp4.m
function [discMat,faceMat,handMat,recallLoc,probeLoc] = createRims(prefs,xPos,yPos,setSize,ss)

    discMat = []; faceMat = []; handMat = [];
    for arc = 1:setSize
        x1 = xPos(arc)-prefs.rimSize/2; y1 = yPos(arc)-prefs.rimSize/2; x2 = xPos(arc)+prefs.rimSize/2; y2 = yPos(arc)+prefs.rimSize/2;
        discPos = [x1 y1 x2 y2]';
        discMat = [discMat, discPos];
        
        x1 = xPos(arc)-prefs.faceSize/2; y1 = yPos(arc)-prefs.faceSize/2; x2 = xPos(arc)+prefs.faceSize/2; y2 = yPos(arc)+prefs.faceSize/2;
        facePos = [x1 y1 x2 y2]';
        faceMat = [faceMat, facePos];
        
        x1 = xPos(arc)-prefs.handX/2; y1 = yPos(arc)-prefs.handY/2; x2 = xPos(arc)+prefs.handX/2; y2 = yPos(arc)+prefs.handY/2;
        handPos = [x1 y1 x2 y2]';
        handMat = [handMat, handPos];
    end
    
    % defining location and shape of probe circle (recall)
    probeLoc = randsample(ss,1);
    probeX1 = xPos(probeLoc)-prefs.rimSize/2;
    probeX2 = xPos(probeLoc)+prefs.rimSize/2;
    probeY1 = yPos(probeLoc)-prefs.rimSize/2;
    probeY2 = yPos(probeLoc)+prefs.rimSize/2;

    recallLoc = [probeX1 probeY1 probeX2 probeY2];
end

% closes the screen by pressing 'Pause/Break' key - JP + JS
function [] = theCloser()
    [~, ~, keyCode] = KbCheck;
    if keyCode(KbName('Pause'))
        Screen('CloseAll');
        stop = here + please;
    end
end