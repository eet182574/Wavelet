% Program that will used WaveletDecomposition() to decompose given image
% into two levels.

clc, clear all, close('all');
fprintf('*********************************************\n');
fprintf('Assignement-1 : > Multimedia Systems ELL786\n');
fprintf('Submission By:\n');
fprintf('1. Parveen Bajaj 2018EET2574\n');
fprintf('2. Abhishek Bohra 2018EET2561\n');
fprintf('*********************************************\n\n');
ImageName = input('Enter file name: ','s');
OriginalImage = imread(ImageName);
[rows,cols,dims] = size(OriginalImage);
if dims == 3
    OriginalImage = rgb2gray(OriginalImage);
end
figure,imshow(OriginalImage),title('Original Image as input');
OriginalImage = double(OriginalImage);

ReconstructionLevels = 0;

fprintf('Choose from following options\n');
fprintf('1 for Haar Wavelet\n');
fprintf('2 for DB9/7\n');
fprintf('3 for DB5/3 (GALL)\n');

choice = input('Enter your option: ','s');

switch(choice)
    
    case {'1'}
        % Haar 
        LPDFilter = [1 1]/sqrt(2);
        HPDFilter = [1 -1]/sqrt(2);
        Origin_L = 0;  % zeroth sample
        Origin_H = 0;  % zeroth sample
        choice = 'Haar';
     
    case {'2'}
        % db 9/7
        LPDFilter = [0 0.0267487574 -0.0168641184 -0.0782232665 0.2668641184 0.6029490182 0.2668641184 -0.0782232665 -0.0168641184 0.0267487574];
        HPDFilter = [0 0.0912717631 -0.0575435263  -0.5912717631 1.1150870525 -0.5912717631 -0.0575435263 0.0912717631 0 0];
        Origin_L = 6;  % zeroth sample
        Origin_H = 5;  % zeroth sample
        choice = 'DB 9/7';
    case {'3'}
        
        %db 5/3 GALL
        LPDFilter = [0 -1/8 2/8 6/8 2/8 -1/8];
        HPDFilter = [0 -1/2 1 -1/2 0 0];
        Origin_L = 4;  % zeroth sample
        Origin_H = 3;  % zeroth sample
        choice = 'DB 5/3 GALL';
end

%% Level-1 decomposition 
[Level_1DecompositionForDisplay,Level_1Decomposition,LL_1Image] = WaveletDecomposition_Fun_04(OriginalImage,LPDFilter,HPDFilter,Origin_L,Origin_H);
figure,imshow(Level_1DecompositionForDisplay,[]);
str = '1st Level Decomposition of Image '; 
title([str choice]);
ReconstructionLevels = ReconstructionLevels + 1; 
%% Level-2 decomposition 
[Level_2DecompositionForDisplay,Level_2Decomposition,LL_2Image] = WaveletDecomposition_Fun_04(LL_1Image,LPDFilter,HPDFilter,Origin_L,Origin_H);
figure,imshow(Level_2DecompositionForDisplay,[]);
str = '2nd Level Decomposition of Image '; 
title([str choice]);

Level_1Decomposition(1:size(Level_2Decomposition,1),1:size(Level_2Decomposition,2)) = Level_2Decomposition;
Level_1DecompositionForDisplay(1:size(Level_2Decomposition,1),1:size(Level_2Decomposition,2)) = Level_2DecompositionForDisplay;
figure,imshow(Level_1DecompositionForDisplay,[]);
str = '2nd Level Decomposition embedded in 1st Level Decomposition of Image';
title([str choice]);
ReconstructionLevels = ReconstructionLevels + 1;
%% Level-3 decomposition 
[Level_3DecompositionForDisplay,Level_3Decomposition,LL_3Image] = WaveletDecomposition_Fun_04(LL_2Image,LPDFilter,HPDFilter,Origin_L,Origin_H);
figure,imshow(Level_3DecompositionForDisplay,[]);
str = '3rd Level Decomposition of Image ';
title([str choice]);

Level_1Decomposition(1:size(Level_3Decomposition,1),1:size(Level_3Decomposition,2)) = Level_3Decomposition;
Level_1DecompositionForDisplay(1:size(Level_3Decomposition,1),1:size(Level_3Decomposition,2)) = Level_3DecompositionForDisplay;
figure,imshow(Level_1DecompositionForDisplay,[]);
str = '3rd Level Decomposition embedded in 1st Level Decomposition of Image';
title([str choice]);
ReconstructionLevels = ReconstructionLevels + 1;
%% Level-4 decomposition 
[Level_4DecompositionForDisplay,Level_4Decomposition,LL_4Image] = WaveletDecomposition_Fun_04(LL_3Image,LPDFilter,HPDFilter,Origin_L,Origin_H);
figure,imshow(Level_4DecompositionForDisplay,[]);
str = '4th Level Decomposition of Image ';
title([str choice]);

Level_1Decomposition(1:size(Level_4Decomposition,1),1:size(Level_4Decomposition,2)) = Level_4Decomposition;
Level_1DecompositionForDisplay(1:size(Level_4Decomposition,1),1:size(Level_4Decomposition,2)) = Level_4DecompositionForDisplay;
figure,imshow(Level_1DecompositionForDisplay,[]);
str = '4th Level Decomposition embedded in 1st Level Decomposition of Image ';
title([str choice]);
ReconstructionLevels = ReconstructionLevels + 1;
%% Quantize the coefficents using b bits
[ScaledCoefficents,ScalingFactor] = Quantize(Level_1Decomposition);
figure,imshow(ScaledCoefficents,[]);title('ScaledCoefficents');
%% Now Undo the effect of Quantization
RescaledCoefficents = ScaledCoefficents * ScalingFactor;

%% Now Reconstruct the Images 
for i = ReconstructionLevels %:-1:1
    RowRange = rows/ 2^(i-1);
    ColRange = cols/ 2^(i-1);
    Seg = Level_1Decomposition(1:RowRange,1:ColRange);
    figure, imshow(Seg,[]);
    temp = ReconstImage_ver_01(Seg,LPDFilter,HPDFilter);
    Level_1Decomposition(1:RowRange,1:ColRange) = temp;
    figure, imshow(Level_1Decomposition);
end






% %% Level-1 decomposition using Daub
% [Level_1Decomposition,LL_1Image] = WaveletDecomposition_Fun_02(OriginalImage,Daubechies4PointCoefficentsLPF);
% figure,imshow(Level_1Decomposition,[]);title('1st Level Decomposition of Image using Daub');
%  
% %% Level-2 decomposition
% [Level_2Decomposition,LL_2Image] = WaveletDecomposition_Fun_02(LL_1Image,Daubechies4PointCoefficentsLPF);
% figure,imshow(Level_2Decomposition,[]);title('2nd Level Decomposition of Image using Daub');
% 
% Level_1Decomposition(1:size(Level_2Decomposition,1),1:size(Level_2Decomposition,2)) = Level_2Decomposition;
% figure,imshow(Level_1Decomposition,[]);title('2nd Level Decomposition embedded in 1st Level Decomposition of Image using Daub');
