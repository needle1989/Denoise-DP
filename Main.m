%% The initialization of the experiment
clear;
close all;
server = 1;
if server == 1
    slash = '/';
elseif server == 0
    slash = '\';
end
bb1=8; % block size of KSVD
bb2=8; % block size of EBSBL
RR=4; % redundancy factor  ��������
K1=RR*bb1^2; % number of atoms in the dictionary
K2=RR*bb2^2; % number of atoms in the dictionary
sigma =10;
p = 13;
MM = 64; % sensor measurement number
%% Get the original image
if p == 1
    picture = 'Lena512.png';
elseif p == 2
    picture = 'cameraman512.png';
elseif p == 3
    picture = 'house512.png';
elseif p == 4
    picture = 'barbara.png';
elseif p == 5
    picture = 'peppers.png';
elseif p == 6
    picture = 'lake512.png';
elseif p == 7
    picture = 'boats.tif';
elseif p == 8
    picture = 'pirate512.png';
elseif p == 9
    picture = 'walkbridge512.png';
elseif p == 10
    picture = 'woman_darkhair512.png';
elseif p == 11
    picture = 'livingroom512.png';
elseif p == 12
    picture = 'jetplane512.png';
elseif p == 13
    picture = 'exp8.png';
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
%% Get the noisy image
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
%% KSVD to denoise the image
[IoutAdaptive1,~] = denoiseImageKSVD(IMin, sigma,K1,bb1);
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

subplot(2,2,4); 
imshow(IoutAdaptive2,[]); 
title(strcat(['Clean Image by DP, ',num2str(PSNROut2),'dB SSIM:',num2str(SSIMOut2)]));
%% Download the result of the experiment
time = datestr(datetime,'yyyymmddHHMMSS');
file_name = [time,'-',picture,'-',num2str(sigma)];
mkdir(['datas',slash],file_name);
path = ['datas',slash,file_name,slash];
saveas(gcf, [path,'all.fig']);
imwrite(uint8(IMin0), [path,'original.png']);
imwrite(uint8(IMin), [path,'noisy-',num2str(PSNRIn),'-',num2str(SSIMIn),'.png']);
imwrite(uint8(IoutAdaptive1), [path,'KSVD-',num2str(PSNROut1),'-',num2str(SSIMOut1),'.png']);
imwrite(uint8(IoutAdaptive2), [path,'KSVD+DP-',num2str(PSNROut2),'-',num2str(SSIMOut2),'.png']);