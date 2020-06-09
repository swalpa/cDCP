clc
clear all
close all

addpath 'TC10_Classifier '/

 
 dn = '/home/swalpa/Downloads/SSELBP_Code/cDCP/images/';
 

 db=dir(strcat(dn,'*.ras'));
 k=1;
 tic
 feature=[];
 
 % Radius and Neighborhood
 R=1;
 P=8;
 

 
 for(i=1:1:length(db))
 
     fname=db(i).name;
     fname=strcat(dn,fname);
     im=imread(fname);
     
     Gray = im2double(im);
     img = (Gray-mean(Gray(:)))/std(Gray(:))*20+128; % image normalization, to remove global intensity 
    
     LBP_H = main(img);
      
                                                 
     feature = [feature; LBP_H];
     k=k+1   
 end;
 k-1
save('CDCP','feature');
toc

clc;

load CDCP.mat
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


