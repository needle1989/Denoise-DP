IMin0 = img{4,1};
IMin = img{3,1};
IoutAdaptive1 = img{1,1};
IoutAdaptive2 = img{2,1};
PSNRIn = 28.1763;
SSIMIn = 0.74205;
PSNROut1 = 30.903;
SSIMOut1 = 0.88317;
PSNROut2 = 31.9615;
SSIMOut2 = 0.89182;

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