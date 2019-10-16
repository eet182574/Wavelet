function [Result] = ReconstImage_ver_01(Image,LPDFilter,HPDFilter)    % Reconstruct the image
clc, close('all');
h_s = LPDFilter; 
h_w = HPDFilter;  
FilterLength_HS = length(h_s);
FilterLength_HW = length(h_w);

%Image = Level_1Decomposition;
[rows,cols] = size(Image);
%% Extract Subband Images
SS_Image = Image(1:rows/2,1:cols/2);
SW_Image = Image(1:rows/2,cols/2 + 1:end);
WS_Image = Image(rows/2+1:end,1:cols/2);
WW_Image = Image(rows/2+1:end,cols/2+1:end);

%% Display individual subbands
figure,imshow(SS_Image,[]);title('SS Image');
figure,imshow(SW_Image,[]);title('SW Image');
figure,imshow(WS_Image,[]);title('WS Image');
figure,imshow(WW_Image,[]);title('WW Image');

%% Upsample the images along rows
SS_upsampled = zeros(rows/2,cols);
SW_upsampled = zeros(rows/2,cols);
WS_upsampled = zeros(rows/2,cols);
WW_upsampled = zeros(rows/2,cols);

for i = 1:rows/2
    SS_upsampled(i,1:2:end) = SS_Image(i,:);
    SW_upsampled(i,1:2:end) = SW_Image(i,:);
    WS_upsampled(i,1:2:end) = WS_Image(i,:);
    WW_upsampled(i,1:2:end) = WW_Image(i,:);
end
% Display zero filled individual subbands
figure,imshow(SS_upsampled,[]);title('SS UpImage');
figure,imshow(SW_upsampled,[]);title('SW UpImage');
figure,imshow(WS_upsampled,[]);title('WS UpImage');
figure,imshow(WW_upsampled,[]);title('WW UpImage');

%% Convolve the images rowwise 
SS_upsampledRowConv = zeros(size(SS_upsampled,1),size(SS_upsampled,2)+FilterLength_HS-1);
SW_upsampledRowConv = zeros(size(SW_upsampled,1),size(SW_upsampled,2)+FilterLength_HW-1);
WS_upsampledRowConv = zeros(size(WS_upsampled,1),size(WS_upsampled,2)+FilterLength_HS-1);
WW_upsampledRowConv = zeros(size(WW_upsampled,1),size(WW_upsampled,2)+FilterLength_HW-1);
for i = 1:size(SS_upsampled,1)
    temp1 = padarray(SS_upsampled(i,:),[0 FilterLength_HS-1],'symmetric','both');
    temp2 = conv(temp1,h_s);
    SS_upsampledRowConv(i,:) = temp2(FilterLength_HS:end-FilterLength_HS+1);
    temp1 = padarray(SW_upsampled(i,:),[0 FilterLength_HW-1],'symmetric','both');
    temp2 = conv(temp1,h_w);
    SW_upsampledRowConv(i,:) = temp2(FilterLength_HW:end-FilterLength_HW+1);
    temp1 = padarray(WS_upsampled(i,:),[0 FilterLength_HW-1],'symmetric','both');
    temp2 = conv(temp1,h_s);
    WS_upsampledRowConv(i,:) = temp2(FilterLength_HS:end-FilterLength_HS+1);
    temp1 = padarray(WW_upsampled(i,:),[0 FilterLength_HW-1],'symmetric','both');
    temp2 = conv(temp1,h_w);
    WW_upsampledRowConv(i,:) = temp2(FilterLength_HW:end-FilterLength_HW+1);
end

% Display convlouted results individual subbands
figure,imshow(SS_upsampledRowConv,[]);title('SS Convoluted UpImage');
figure,imshow(SW_upsampledRowConv,[]);title('SW Convoluted UpImage');
figure,imshow(WS_upsampledRowConv,[]);title('WS Convoluted UpImage');
figure,imshow(WW_upsampledRowConv,[]);title('WW Convoluted UpImage');


%% Prepare two images out of four images
S_upsampledRowConv = SS_upsampledRowConv + WS_upsampledRowConv;
W_upsampledRowConv = SW_upsampledRowConv + WW_upsampledRowConv;
figure,imshow(S_upsampledRowConv,[]);title('S = SS+SW Convoluted UpImage');
figure,imshow(W_upsampledRowConv,[]);title('W = WS+WW Convoluted UpImage');

%% Upsample along columns
S_upsampledCol = zeros(rows,cols);
W_upsampledCol = zeros(rows,cols);

for i = 1:cols
    S_upsampledCol(1:2:end,i) = S_upsampledRowConv(:,i);
    W_upsampledCol(1:2:end,i) = W_upsampledRowConv(:,i);
end
figure,imshow(S_upsampledCol,[]);title('Col UpImage');
figure,imshow(W_upsampledCol,[]);title('Col UpImage');

%% Do convolution along cols
S_upsampledColConv = zeros(size(S_upsampledCol,1)+FilterLength_HS-1,size(S_upsampledCol,2));
W_upsampledColConv = zeros(size(W_upsampledCol,1)+FilterLength_HW-1,size(W_upsampledCol,2));
for i = 1:size(S_upsampledCol,2)
    temp1 = padarray(S_upsampledCol(:,i),[FilterLength_HS-1 0],'symmetric','both');
    temp2 = conv(temp1,h_s);
    S_upsampledColConv(:,i) = temp2(FilterLength_HS:end-FilterLength_HS+1);
    temp1 = padarray(W_upsampledCol(:,i),[FilterLength_HW-1 0],'symmetric','both');
    temp2 = conv(temp1,h_w);
    W_upsampledColConv(:,i) = temp2(FilterLength_HW:end-FilterLength_HW+1);
end
H = fspecial('gaussian');
ReconstructedImage = S_upsampledColConv + W_upsampledColConv;
% ReconstructedImage = uint8(mat2gray(ReconstructedImage));
figure,imshow(imfilter(ReconstructedImage,H,'same'),[]), title('Reconstructed Image');
Result = ReconstructedImage;
end


