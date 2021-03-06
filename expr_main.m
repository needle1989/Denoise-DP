%% The initialization of the experiment
clear;
close all;
server = 0;
sub = 13;
if server == 1
    slash = '/';
elseif server == 0
    slash = '\';
end
bb1=8; % block size of KSVD
bb2=8; % block size of EBSBL
RR=4; % redundancy factor
K1=RR*bb1^2; % number of atoms in the dictionary
K2=RR*bb2^2; % number of atoms in the dictionary
sigma =10;
MM = 64; % sensor measurement number
%% Initialize original image
if sub == 1
    picture = 'Lena512.png';
elseif sub == 2
    picture = 'cameraman512.png';
elseif sub == 3
    picture = 'house512.png';
elseif sub == 4
    picture = 'barbara.png';
elseif sub == 5
    picture = 'peppers.png';
elseif sub == 6
    picture = 'lake512.png';
elseif sub == 7
    picture = 'boats.tif';
elseif sub == 8
    picture = 'pirate512.png';
elseif sub == 9
    picture = 'walkbridge512.png';
elseif sub == 10
    picture = 'woman_darkhair512.png';
elseif sub == 11
    picture = 'livingroom512.png';
elseif sub == 12
    picture = 'jetplane512.png';
elseif sub == 13
    picture = '\dataset\livingroom_b2.png';
end
[IMin0,~]=imread(['img',slash,picture]);
IMin0 = imresize(IMin0,[MM MM]);
IMin0=im2double(IMin0);
if (length(size(IMin0))>2)
    IMin0 = rgb2gray(IMin0);
end
if (max(IMin0(:))<2)
    IMin0 = IMin0*255;
end
%% Initialize image with noise
noise=zeros(size(IMin0));
way_to_noise = 3;
if way_to_noise == 1
    for m=1:size(IMin0)-bb+1
        for n=1:size(IMin0)-bb+1
            for i=m:m+7
                for j=n:n+7
                    noi = randn(bb);
                    noise(i,j)=noise(i,j)+noi(i-m+1,j-n+1);
                end
            end
        end
    end
    IMin=IMin0+noise;
elseif way_to_noise == 2
    size = noisesize(1);
    IMin=IMin0+sigma*size;
elseif way_to_noise == 3
    IMin=IMin0+sigma*randn(size(IMin0));
end
PSNRIn = 20*log10(255/sqrt(mean((IMin(:)-IMin0(:)).^2)));
SSIMIn = ssim(IMin0,IMin);
%% Implement K means cluster
% IoutAdaptive1 = expr_cluster_test(IMin,bb1);
% IoutAdaptive1 = expr_denoise_dp_test(IMin,bb1);
% IoutAdaptive1 = denoiseImageDP(IMin,bb1);
% PSNROut1 = 20*log10(255/sqrt(mean((IoutAdaptive1(:)-IMin0(:)).^2)));
% SSIMOut1 = ssim(IMin0,IoutAdaptive1);
%% KSVD to denoise the image
[IoutAdaptive1,~] = denoiseImageKSVD(IMin, sigma,K1,bb1);
% IoutAdaptive1 = denoiseImageDP(IMin,bb1);
PSNROut1 = 20*log10(255/sqrt(mean((IoutAdaptive1(:)-IMin0(:)).^2)));
SSIMOut1 = ssim(IMin0,IoutAdaptive1);
%% EBSBL_BO+KSVD to denoise the image
IoutAdaptive2 = denoiseImageDP(IMin,bb1);
% IoutAdaptive2 = denoiseImageKSVD_DP(IMin, sigma,K1,bb1);
PSNROut2 = 20*log10(255/sqrt(mean((IoutAdaptive2(:)-IMin0(:)).^2)));
SSIMOut2 = ssim(IMin0,IoutAdaptive2);
%% Paint the final result
figure;
subplot(2,2,1); 
imshow(IMin0,[]); 
title('Original clean image');

subplot(2,2,2);
imshow(IMin,[]); 
title(strcat(['Noisy image, ',num2str(PSNRIn),'dB SSIM:',num2str(SSIMIn)]));

subplot(2,2,3); 
imshow(IoutAdaptive1,[]); 
title(strcat(['Clean Image by KSVD, ',num2str(PSNROut1),'dB SSIM:',num2str(SSIMOut1)]));
% 
subplot(2,2,4); 
imshow(IoutAdaptive2,[]); 
title(strcat(['Clean Image by DP, ',num2str(PSNROut2),'dB SSIM:',num2str(SSIMOut2)]));