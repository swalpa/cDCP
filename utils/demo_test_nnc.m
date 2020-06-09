clc;
clear all;
load TC10_Noise05.mat
load TC10Level.mat
load TC10Catagory.mat


% Make Sure before compile
feature =feature;


numFolds=10;
VecsPerCat = getVecsPerCat(feature, Level, Catagory);

foldSizes = computeFoldSizes(VecsPerCat,numFolds);

[X_sorted, y_sorted] = randSortAndGroup(feature, Level, Catagory);
% X_sorted=feature;
% y_sorted=Level;

RoundAccMat=zeros(1,numFolds);

tic
for roundNumber= 1:numFolds
[X_train, y_train, X_val, y_val] = getFoldVectors(X_sorted, y_sorted, Catagory, VecsPerCat, foldSizes, roundNumber);
% classification test using Non Parametric classifier %
trains = X_train;
tests = X_val;
trainNum = size(trains,1);
testNum = size(tests,1);
DM = zeros(testNum,trainNum);
for i=1:testNum;
    test = tests(i,:);        
    DM(i,:) = distMATChiSquare(trains,test)';
end
CP=ClassifyOnNN(DM,y_train,y_val)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RoundAccMatNN(roundNumber)= CP;
end
RoundAccMatNN
AveAcc=sum(RoundAccMatNN)/length(RoundAccMatNN) 
Test = std(RoundAccMatNN);
toc
fprintf('Stadard Deviation of NNC Test is =  +- %f\n',Test); 


clear Train;
clear Test;




