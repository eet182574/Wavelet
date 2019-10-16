function [ResultForDisplay,Result,ScalingScalingResultColsubsampled] = WaveletDecomposition_Fun_04(OriginalImage,ScalingCoefficents,WaveletCoefficents,Origin_L,Origin_H)
% This version 4 is different from ver03 in terms of LPDFilter and HPD
% filter and Origin of filters are passed now.
[rows,cols,dims] = size(OriginalImage);
h_s = ScalingCoefficents;
h_w = WaveletCoefficents;
FilterLength_HS = length(h_s);
FilterLength_HW = length(h_w);
%% Perform row wise scaling and wavelet transform
ScalingResultRow = zeros(rows,cols+FilterLength_HS-1);
WaveletResultRow = zeros(rows,cols+FilterLength_HW-1);
for i = 1:rows
    temp1 = padarray(OriginalImage(i,:),[0 FilterLength_HS-1],'symmetric','both');
    temp2 = conv(temp1,h_s(end:-1:1));
    ScalingResultRow(i,:) = temp2(FilterLength_HS:end-FilterLength_HS+1);
    temp1 = padarray(OriginalImage(i,:),[0 FilterLength_HW-1],'symmetric','both');
    temp2 = conv(temp1,h_w(end:-1:1));
    WaveletResultRow(i,:) = temp2(FilterLength_HW:end-FilterLength_HW+1);
end
figure,imshow(ScalingResultRow,[]),title('RowWise Scaling image');
figure,imshow(WaveletResultRow,[]),title('RowWise Wavelet image');
%% Perform subsampling along rows
% if Origin_L == 0
%     offset_hs = FilterLength_HS;  % Offset is for h(n) where origin of h(n) is at n=0. As h(n) becomes h(-n), origin of h(-n) becomes 1. And further h(-n) is flipped in conv operation, so origin goes back to 0.
%     offset_hw = FilterLength_HW; 
% else 
%     offset_hs = Origin_L - 1;
%     offset_hw = Origin_H + 1;
% end

offset_hs = FilterLength_HS;  % Offset is for h(n) where origin of h(n) is at n=0. As h(n) becomes h(-n), origin of h(-n) becomes 1. And further h(-n) is flipped in conv operation, so origin goes back to 0.
offset_hw = FilterLength_HW; 


ScalingResultRowsubsampled = zeros(rows,length(offset_hs:2:size(ScalingResultRow,2)));
WaveletResultRowsubsampled = zeros(rows,length(offset_hw:2:size(WaveletResultRow,2)));
for i = 1:rows
    ScalingResultRowsubsampled(i,:) = ScalingResultRow(i,offset_hs:2:size(ScalingResultRow,2)); 
    WaveletResultRowsubsampled(i,:) = WaveletResultRow(i,offset_hw:2:size(WaveletResultRow,2)); 
end

% figure,imshow(ScalingResultRowsubsampled,[]),title('RowWise Scaling and subsampled image');
% figure,imshow(WaveletResultRowsubsampled,[]),title('RowWise Wavelet and subsampled image');

%% Perform col wise scaling and wavelet transform
ScalingScalingResultCol = zeros(size(ScalingResultRowsubsampled,1) + FilterLength_HS-1,size(ScalingResultRowsubsampled,2));
ScalingWaveletResultCol = zeros(size(ScalingResultRowsubsampled,1) + FilterLength_HW-1,size(ScalingResultRowsubsampled,2));
WaveletScalingResultCol = zeros(size(WaveletResultRowsubsampled,1) + FilterLength_HS-1,size(WaveletResultRowsubsampled,2));
WaveletWaveletResultCol = zeros(size(WaveletResultRowsubsampled,1) + FilterLength_HW-1,size(WaveletResultRowsubsampled,2));
for i = 1:size(ScalingResultRowsubsampled,2)
    temp1 = padarray(ScalingResultRowsubsampled(:,i),[FilterLength_HS-1 0],'symmetric','both');
    temp2 = conv(temp1,h_s(end:-1:1));
    ScalingScalingResultCol(:,i) = temp2(FilterLength_HS:end-FilterLength_HS+1);
    temp2 = conv(temp1,h_w(end:-1:1));
    ScalingWaveletResultCol(:,i) = temp2(FilterLength_HS:end-FilterLength_HS+1);
    temp1 = padarray(WaveletResultRowsubsampled(:,i),[FilterLength_HS-1 0],'symmetric','both');
    temp2 = conv(temp1,h_s(end:-1:1));
    WaveletScalingResultCol(:,i) = temp2(FilterLength_HS:end-FilterLength_HS+1);
    temp2 = conv(temp1,h_w(end:-1:1));
    WaveletWaveletResultCol(:,i) = temp2(FilterLength_HS:end-FilterLength_HS+1);
end
%figure,imshow(ScalingScalingResultCol,[]),title('ColWise ScalingScaling image');
%figure,imshow(ScalingWaveletResultCol,[]),title('ColWise ScalingWavelet image');
%figure,imshow(WaveletScalingResultCol,[]),title('ColWise WaveletScaling image');
%figure,imshow(WaveletWaveletResultCol,[]),title('ColWise WaveletWavelet image');
%% Perform subsampling along cols
ScalingScalingResultColsubsampled = zeros(length(offset_hs:2:size(ScalingScalingResultCol,1)),size(ScalingScalingResultCol,2));
ScalingWaveletResultColsubsampled = zeros(length(offset_hw:2:size(ScalingWaveletResultCol,1)),size(ScalingWaveletResultCol,2));
WaveletScalingResultColsubsampled = zeros(length(offset_hs:2:size(WaveletScalingResultCol,1)),size(WaveletScalingResultCol,2));
WaveletWaveletResultColsubsampled = zeros(length(offset_hw:2:size(WaveletWaveletResultCol,1)),size(WaveletWaveletResultCol,2));

for i = 1:size(ScalingScalingResultCol,2)
    ScalingScalingResultColsubsampled(:,i) = ScalingScalingResultCol(offset_hs:2:size(ScalingScalingResultCol,1),i); 
    ScalingWaveletResultColsubsampled(:,i) = ScalingWaveletResultCol(offset_hw:2:size(ScalingWaveletResultCol,1),i); 
    WaveletScalingResultColsubsampled(:,i) = WaveletScalingResultCol(offset_hs:2:size(WaveletScalingResultCol,1),i); 
    WaveletWaveletResultColsubsampled(:,i) = WaveletWaveletResultCol(offset_hw:2:size(WaveletWaveletResultCol,1),i); 
end
% figure,imshow(ScalingScalingResultColsubsampled,[]),title('ColWise Scaling Scaling and subsampled image');
% figure,imshow(ScalingWaveletResultColsubsampled,[]),title('ColWise Scaling Wavelet and subsampled image');
% figure,imshow(WaveletScalingResultColsubsampled,[]),title('ColWise Wavelet Scaling and subsampled image');
% figure,imshow(WaveletWaveletResultColsubsampled,[]),title('ColWise Wavelet Wavelet and subsampled image');

%% pack the results in one image
Result = zeros(rows,cols);
ResultForDisplay = Result;

Result(1:(size(ScalingScalingResultColsubsampled,1)),1:(size(ScalingScalingResultColsubsampled,2))) = ScalingScalingResultColsubsampled;
Result(1:(size(ScalingScalingResultColsubsampled,1)),size(WaveletWaveletResultColsubsampled,2)+1:2*size(WaveletWaveletResultColsubsampled,2)) = ScalingWaveletResultColsubsampled;
Result(size(WaveletWaveletResultColsubsampled,1)+1:2*size(WaveletWaveletResultColsubsampled,1),1:(size(ScalingScalingResultColsubsampled,2))) = WaveletScalingResultColsubsampled;
Result(size(WaveletWaveletResultColsubsampled,1)+1:2*size(WaveletWaveletResultColsubsampled,1),size(WaveletWaveletResultColsubsampled,2)+1:2*size(WaveletWaveletResultColsubsampled,2)) = WaveletWaveletResultColsubsampled;

ResultForDisplay(1:(size(ScalingScalingResultColsubsampled,1)),1:(size(ScalingScalingResultColsubsampled,2))) = mat2gray(ScalingScalingResultColsubsampled);
ResultForDisplay(1:(size(ScalingScalingResultColsubsampled,1)),size(WaveletWaveletResultColsubsampled,2)+1:2*size(WaveletWaveletResultColsubsampled,2)) = mat2gray(ScalingWaveletResultColsubsampled);
ResultForDisplay(size(WaveletWaveletResultColsubsampled,1)+1:2*size(WaveletWaveletResultColsubsampled,1),1:(size(ScalingScalingResultColsubsampled,2))) = mat2gray(WaveletScalingResultColsubsampled);
ResultForDisplay(size(WaveletWaveletResultColsubsampled,1)+1:2*size(WaveletWaveletResultColsubsampled,1),size(WaveletWaveletResultColsubsampled,2)+1:2*size(WaveletWaveletResultColsubsampled,2)) = mat2gray(WaveletWaveletResultColsubsampled);
% figure,imshow(Result,[]),title('Scale -1 decomposition');
end