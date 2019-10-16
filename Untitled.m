clc, clear all, close('all');
ImageName = input('Enter file name: ','s');
OriginalImage = imread(ImageName);
[rows,cols,dims] = size(OriginalImage);
if dims == 3
    OriginalImage = rgb2gray(OriginalImage);
end
figure,imshow(OriginalImage),title('Original Image as input');
OriginalImage = double(OriginalImage);
h_v = [1/sqrt(2) 1/sqrt(2)];
h_w = [1/sqrt(2) -1/sqrt(2)];
%% Perform row wise scaling and wavelet transform
ScalingResultRow = zeros(rows,cols+length(h_v)-1);
WaveletResultRow = zeros(rows,cols+length(h_v)-1);
for i = 1:rows
    ScalingResultRow(i,:) = conv(OriginalImage(i,:),h_v(end:-1:1));
    WaveletResultRow(i,:) = conv(OriginalImage(i,:),h_w(end:-1:1));
end
figure,imshow(ScalingResultRow,[]),title('RowWise Scaling image');
figure,imshow(WaveletResultRow,[]),title('RowWise Wavelet image');
%% Perform subsampling along rows
offset = length(h_v) - 1;  % Offset is for h(n) where origin of h(n) is at n=0. As h(n) becomes h(-n), origin of h(-n) becomes 1. And further h(-n) is flipped in conv operation, so origin goes back to 0.
ScalingResultRowsubsampled = zeros(rows,length(1+offset:2:size(ScalingResultRow,2)));
WaveletResultRowsubsampled = ScalingResultRowsubsampled;
for i = 1:rows
    ScalingResultRowsubsampled(i,:) = ScalingResultRow(i,1+offset:2:size(ScalingResultRow,2)); 
    WaveletResultRowsubsampled(i,:) = WaveletResultRow(i,1+offset:2:size(ScalingResultRow,2)); 
end

figure,imshow(ScalingResultRowsubsampled,[]),title('RowWise Scaling and subsampled image');
figure,imshow(WaveletResultRowsubsampled,[]),title('RowWise Wavelet and subsampled image');

%% Perform col wise scaling and wavelet transform
ScalingResultCol = zeros(size(ScalingResultRowsubsampled,1) + length(h_v)-1,size(ScalingResultRowsubsampled,2));
WaveletResultCol = ScalingResultCol;

for i = 1:size(ScalingResultRowsubsampled,2)
    ScalingResultCol(:,i) = conv(ScalingResultRowsubsampled(:,i),h_v(end:-1:1));
    WaveletResultCol(:,i) = conv(WaveletResultRowsubsampled(:,i),h_w(end:-1:1));
end

figure,imshow(ScalingResultCol,[]),title('ColWise Scaling image');
figure,imshow(WaveletResultCol,[]),title('ColWise Wavelet image');
%% Perform subsampling along cols
offset = length(h_v) - 1;  % Offset is for h(n) where origin of h(n) is at n=0. As h(n) becomes h(-n), origin of h(-n) becomes 1. And further h(-n) is flipped in conv operation, so origin goes back to 0.
ScalingResultColsubsampled = zeros(length(1+offset:2:size(ScalingResultCol,1)),size(ScalingResultCol,2));
WaveletResultColsubsampled = ScalingResultColsubsampled;
for i = 1:size(ScalingResultCol,2)
    ScalingResultColsubsampled(:,i) = ScalingResultCol(1+offset:2:size(ScalingResultCol,1),i); 
    WaveletResultColsubsampled(:,i) = WaveletResultCol(1+offset:2:size(ScalingResultCol,1),i); 
end
figure,imshow(ScalingResultColsubsampled,[]),title('ColWise Scaling and subsampled image');
figure,imshow(WaveletResultColsubsampled,[]),title('ColWise Wavelet and subsampled image');


%% pack the results in one image
Result = zeros(rows,cols);
Result(1:(size(ScalingResultColsubsampled,1)),1:(size(ScalingResultColsubsampled,2))) = ScalingResultColsubsampled;

Result(size(ScalingResultColsubsampled,1)+1:2*size(ScalingResultColsubsampled,1),size(ScalingResultColsubsampled,2)+1:2*size(ScalingResultColsubsampled,2)) = WaveletResultColsubsampled;

figure,imshow(Result,[]),title('Scale -1 decomposition');