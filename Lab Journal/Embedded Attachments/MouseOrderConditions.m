%% Pserudorandomize the mouse and order conditions - Lena's NURA

clear;
clc;

n = 100;
type1 = zeros(n/4, 1) + 1;
type2 = zeros(n/4, 1) + 2;
type3 = zeros(n/4, 1) + 3;
type4 = zeros(n/4, 1) + 4;

conditionList = vertcat(type1,type2, type3, type4);
conditionList = conditionList(randperm(length(conditionList)));

for i = 1:100

    switch conditionList(i)
        case 1
            bothConditions(i,:) = [0 0 0];
        case 2
            bothConditions(i,:) = [0 0 1];
        case 3
            bothConditions(i,:) = [1 1 0];
        case 4
            bothConditions(i,:) = [1 1 1];

    end   

end

disp('Done!')
save('/Users/lena/Desktop/NURA24/conditions.mat','bothConditions')