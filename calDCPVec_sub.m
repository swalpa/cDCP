function [DCP_S DCP_C DCP_M] = calDCPVec_sub(Input_Im, RegionRowNum, RegionColNum, PtsNum, Radius, OffSet, mapping)
%% DESCRIPTION
%   	The function is to generate an feature vector by DCP-1 or DCP-2.
%
%   REFERENCE:
%       Changxing Ding, Jonghyun Choi, Dacheng Tao, and Larry S. Davis,
%       'Multi-Directional Multi-Level Dual-Cross Patterns for Robust Face 
%       Recognition', Vol.38, No.3, pp.518-531, TPAMI 2016.
%
%   AUTHOR:
%       Changxing Ding @ University of Technology, Sydney
%
%   VERSION:
%       0.1 - 08/06/2013 first implementation
%       Matlab 2012a



%% Construct a Block to Calculate DCP code
% calculate the relative coordinates for neighbor points
Neighbors = zeros(PtsNum,4);
step = 2*pi/PtsNum;   % angle step.
for i = 1:PtsNum
    Neighbors(i,1) = -Radius(1)*cos((i-1)*step + OffSet);   % for x coordinate (vertical)
    Neighbors(i,2) = Radius(1)*sin((i-1)*step + OffSet);    % for y coordinate (horizontal)
    Neighbors(i,3) = -Radius(2)*cos((i-1)*step + OffSet);   % for x coordinate (vertical)
    Neighbors(i,4) = Radius(2)*sin((i-1)*step + OffSet);    % for y coordinate (horizontal)
end

minX = min(min(Neighbors(:,[1 3])));
maxX = max(max(Neighbors(:,[1 3])));
minY = min(min(Neighbors(:,[2 4])));
maxY = max(max(Neighbors(:,[2 4])));

% block size, each code is computed within a block of size bsizey*bsizex
bSizeX = ceil(max(maxX, 0)) - floor(min(minX, 0)) + 1;
bSizeY = ceil(max(maxY, 0)) - floor(min(minY, 0)) + 1;

% coordinates of origin (0,0) in the block
origX = 1 - floor(min(minX,0));
origY = 1 - floor(min(minY,0));


%% Preparation for the Coding
Input_Im = double(Input_Im);
[imHeight, imWidth] = size(Input_Im);

% Minimum allowed size for the input image
if(imHeight < bSizeX || imWidth < bSizeY)
  error('Too small input image. Should be at least (2*Radius+1) x (2*Radius+1)');
end

% Calculate dx and dy;
dx = imHeight - bSizeX;
dy = imWidth - bSizeY;

% Fill the center pixel matrix C.
C = Input_Im(origX:origX + dx, origY:origY + dy);
d_C = double(C);

% Initialize the Result matrix with zeros.
%Result = zeros(dx+1, dy+1);

% Initialize the result matrix with zeros.
DCP_S=zeros(dy+1,dx+1);
DCP_M=zeros(dy+1,dx+1);
DCP_C=zeros(dy+1,dx+1);



%% Compute the DCP Code Image
for i = 1:PtsNum
    %% for the first circle
    x = Neighbors(i,1) + origX;
    y = Neighbors(i,2) + origY;
    % Calculate floors, ceils and rounds for the x and y.
    fx = floor(x); cx = ceil(x); rx = round(x);
    fy = floor(y); cy = ceil(y); ry = round(y);

    % Check if interpolation is needed.
    if (abs(x - rx) < 1e-6) && (abs(y - ry) < 1e-6)
        % Interpolation is not needed, use original datatypes
        N1 = Input_Im(rx:rx+dx, ry:ry+dy);
        D1 = (N1 >= C); 
        DiffD1{i} = abs(N1 - C);
        MeanDiffD1(i) = mean(mean(DiffD1{i}));
    else
        % Interpolation needed, use double type images 
        tx = x - fx;
        ty = y - fy;
        
        % Calculate the interpolation weights.
        w1 = (1 - tx) * (1 - ty);
        w2 =      tx  * (1 - ty);
        w3 = (1 - tx) *      ty ;
        w4 =      tx  *      ty ;
        % Compute interpolated pixel values
        N1 = w1*Input_Im(fx:fx+dx, fy:fy+dy) + w2*Input_Im(cx:cx+dx, fy:fy+dy) + ...,
            w3*Input_Im(fx:fx+dx, cy:cy+dy) + w4*Input_Im(cx:cx+dx, cy:cy+dy);
        D1 = (N1 >= C); 
        DiffD1{i} = abs(N1 - C);
        MeanDiffD1(i) = mean(mean(DiffD1{i}));
    end  
    
    %% for the second circle
    x = Neighbors(i,3) + origX;
    y = Neighbors(i,4) + origY;
    % Calculate floors, ceils and rounds for the x and y.
    fx = floor(x); cx = ceil(x); rx = round(x);
    fy = floor(y); cy = ceil(y); ry = round(y);

    % Check if interpolation is needed.
    if (abs(x - rx) < 1e-6) && (abs(y - ry) < 1e-6)
        % Interpolation is not needed, use original datatypes
        N2 = Input_Im(rx:rx+dx, ry:ry+dy);
        D2 = (N2 >= N1); 
        DiffD2{i} = abs(N2 - N1);
        MeanDiffD2(i) = mean(mean(DiffD2{i}));
    else
        % Interpolation needed, use double type images 
        tx = x - fx;
        ty = y - fy;
        
        % Calculate the interpolation weights.
        w1 = (1 - tx) * (1 - ty);
        w2 =      tx  * (1 - ty);
        w3 = (1 - tx) *      ty ;
        w4 =      tx  *      ty ;
        
        % Compute interpolated pixel values
        N2 = w1*Input_Im(fx:fx+dx, fy:fy+dy) + w2*Input_Im(cx:cx+dx, fy:fy+dy) + ...,
             w3*Input_Im(fx:fx+dx, cy:cy+dy) + w4*Input_Im(cx:cx+dx, cy:cy+dy);
        D2 = (N2 >= N1); 
        DiffD2{i} = abs(N2 - N1);
        MeanDiffD2(i) = mean(mean(DiffD2{i}));
    end  
    
    % Difference threshold for DCP_M
   DiffThresholdD1 = mean(MeanDiffD1);
   DiffThresholdD2 = mean(MeanDiffD2);
   D1_M = DiffD1{i} >= DiffThresholdD1;
   D2_M = DiffD2{i} >= DiffThresholdD2;

    % Update the result matrix.
    v = 4^(i-1);
    DCP_S = DCP_S + v*(D2*2+D1);
    DCP_M = DCP_M + v*(D2_M*2+D1_M);
end

DCP_C = d_C>=mean(Input_Im(:));

%disp('size of rows and cols');
%size(Result,1)
%size(Result,2)
%figure, imshow(uint8(Result));

%Apply mapping if it is defined
if isstruct(mapping)
    bins = mapping.num;
    %sizarray = size(CLBP_S);
    for i = 1:size(DCP_S,1)
        for j = 1:size(DCP_S,2)
            DCP_S(i,j) = mapping.table(DCP_S(i,j)+1);
            DCP_M(i,j) = mapping.table(DCP_M(i,j)+1);
        end
    end
end

%disp('size of rows and cols');
%size(Result)
%size(DCP_C)


%% Divid the DCP Map and Apply Pooling
%histoDim = power(4,PtsNum);
%Result = calSpatialHistogram(Result, 0:histoDim-1, RegionRowNum, RegionColNum);

%%% This change made by Swalpa Kumar Roy
%Result = calSpatialHistogram(Result, 0:mapping.num-1, 1, 1);

end




