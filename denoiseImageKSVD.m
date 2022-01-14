function [IOut,output] = denoiseImageKSVD(Image,sigma,K,blockSize)
%==========================================================================
%   P E R F O R M   D E N O I S I N G   U S I N G   A  D I C T I O N A R Y
%                  T R A I N E D   O N   N O I S Y   I M A G E
%==========================================================================
% function IOut = denoiseImageKSVD(Image,sigma,K,varargin)
% denoise an image by sparsely representing each block with the
% already overcomplete trained Dictionary, and averaging the represented parts.
% Detailed description can be found in "Image Denoising Via Sparse and Redundant
% representations over Learned Dictionaries", (appeared in the 
% IEEE Trans. on Image Processing, Vol. 15, no. 12, December 2006).
% This function may take some time to process. Possible factor that effect
% the processing time are:
%  1. number of KSVD iterations - the default number of iterations is 10.
%  However, fewer iterations may, in most cases, result an acceleration in
%  the process, without effecting  the result too much. Therefore, when
%  required, this parameter may be re-set.
%  2. maxBlocksToConsider - The maximal number of blocks to train on. If this 
%  number is larger the number of blocks in the image, random blocks
%  from the image will be selected for training. 
% ===================================================================
% INPUT ARGUMENTS : Image - the noisy image (gray-level scale)
%                   sigma - the s.d. of the noise (assume to be white Gaussian).
%                   K - the number of atoms in the trained dictionary.
%    Optional arguments:              
%                  'blockSize' - the size of the blocks the algorithm
%                       works. All blocks are squares, therefore the given
%                       parameter should be one number (width or height).
%                       Default value: 8.
%                  'errorFactor' - a factor that multiplies sigma in order
%                       to set the allowed representation error. In the
%                       experiments presented in the paper, it was set to 1.15
%                       (which is also the default  value here).
%                  'maxBlocksToConsider' - maximal number of blocks that
%                       can be processed. This number is dependent on the memory
%                       capabilities of the machine, and performances?
%                       considerations. If the number of available blocks in the
%                       image is larger than 'maxBlocksToConsider', the sliding
%                       distance between the blocks increases. The default value
%                       is: 250000.
%                  'slidingFactor' - the sliding distance between processed
%                       blocks. Default value is 1. However, if the image is
%                       large, this number increases automatically (because of
%                       memory requirements). Larger values result faster
%                       performances (because of fewer processed blocks).
%                  'numKSVDIters' - the number of KSVD iterations processed
%                       blocks from the noisy image. If the number of
%                       blocks in the image is larger than this number,
%                       random blocks from all available blocks will be
%                       selected. The default value for this parameter is:
%                       10 if sigma > 5, and 5 otherwise.
%                  'maxNumBlocksToTrainOn' - the maximal number of blocks
%                       to train on. The default value for this parameter is
%                       65000. However, it might not be enough for very large
%                       images
%                  'displayFlag' - if this flag is switched on,
%                       announcement after finishing each iteration will appear,
%                       as also a measure concerning the progress of the
%                       algorithm (the average number of required coefficients
%                       for representation). The default value is 1 (on).
%                  'waitBarOn' - can be set to either 1 or 0. If
%                       waitBarOn==1 a waitbar, presenting the progress of the
%                       algorithm will be displayed.
% OUTPUT ARGUMENTS : Iout - a 2-dimensional array in the same size of the
%                       input image, that contains the cleaned image.
%                    output.D - the trained dictionary.
% =========================================================================

% first, train a dictionary on the noisy image

reduceDC = 1;
[NN1,NN2] = size(Image);
waitBarOn = 1;
if (sigma > 5)%%%sigma=50   numIterOfKsvd = 10;
    numIterOfKsvd = 10;
else
    numIterOfKsvd = 5;
end
errorFactor = 1.15;
maxBlocksToConsider = 260000;
slidingDis = 1;
%blockSize = 8;
maxNumBlocksToTrainOn = 65000;
displayFlag = 1;
%hh=length(varargin);%%%%%%%%%%%����һ���ܲ��ܽ��������forѭ����ȥ��
% for argI = 1:2:length(varargin)
%     if (strcmp(varargin{argI}, 'slidingFactor'))
%         slidingDis = varargin{argI+1};
%     end
%     if (strcmp(varargin{argI}, 'errorFactor'))
%         C = varargin{argI+1};
%     end
%     if (strcmp(varargin{argI}, 'maxBlocksToConsider'))
%         maxBlocksToConsider = varargin{argI+1};
%     end
%     if (strcmp(varargin{argI}, 'numKSVDIters'))
%         numIterOfKsvd = varargin{argI+1};
%     end
%     if (strcmp(varargin{argI}, 'blockSize'))
%         bb = varargin{argI+1};
%     end
%     if (strcmp(varargin{argI}, 'maxNumBlocksToTrainOn'))
%         maxNumBlocksToTrainOn = varargin{argI+1};
%     end
%     if (strcmp(varargin{argI}, 'displayFlag'))
%         displayFlag = varargin{argI+1};
%     end
%     if (strcmp(varargin{argI}, 'waitBarOn'))
%         waitBarOn = varargin{argI+1};
%     end
% end

%if (sigma <= 5)
%    numIterOfKsvd = 5;
%end

% first, train a dictionary on blocks from the noisy image

if(prod([NN1,NN2]-blockSize+1)> maxNumBlocksToTrainOn)
    randPermutation =  randperm(prod([NN1,NN2]-blockSize+1));
    selectedBlocks = randPermutation(1:maxNumBlocksToTrainOn);

    blkMatrix = zeros(blockSize^2,maxNumBlocksToTrainOn);
    for i = 1:maxNumBlocksToTrainOn
        [row,col] = ind2sub(size(Image)-blockSize+1,selectedBlocks(i));
        currBlock = Image(row:row+blockSize-1,col:col+blockSize-1);
        blkMatrix(:,i) = currBlock(:);
    end
else
    blkMatrix = im2col(Image,[blockSize,blockSize],'sliding');%%%%%%%8*8=64   ����blkMatrix�����СΪ��64*[��NN1-bb+1��*(NN2-bb+1)]
end

param.K = K;%%%K=256  4*8*8=256
param.numIteration = numIterOfKsvd ;%sigma=50   ����numIterOfKsvd = 10;

param.errorFlag = 1; % decompose signals until a certain error is reached. do not use fix number of coefficients.
param.errorGoal = sigma*errorFactor;
param.preserveDCAtom = 0;

Pn=ceil(sqrt(K));%%Pn=16
DCT=zeros(blockSize,Pn);%%blockSize=8
for k=0:1:Pn-1
    V=cos((0:1:blockSize-1)'*k*pi/Pn);
    if k>0, V=V-mean(V); end
    DCT(:,k+1)=V/norm(V);
end
DCT=kron(DCT,DCT);%%%%%��DCT�еĴ���һ����     64*256�ľ���  ****param.numIteration

param.initialDictionary = DCT(:,1:param.K );%%%% ȡ��256�С�Ҳ����ȫ����ȡ��
param.InitializationMethod =  'GivenMatrix';

if (reduceDC)%%reduceDC=1
    vecOfMeans = mean(blkMatrix);
    blkMatrix = blkMatrix-ones(size(blkMatrix,1),1)*vecOfMeans;%%%��ȥƽ����  blkMatrix�����СΪ��64*[��NN1-bb+1��*(NN2-bb+1)]
end

if (waitBarOn)%waitBarOn=1
    counterForWaitBar = param.numIteration+1;%param.numIteration = numIterOfKsvd ;  =10
    h = waitbar(0,'Denoising In Process ...');
    param.waitBarHandle = h;
    param.counterForWaitBar = counterForWaitBar;
end


param.displayProgress = displayFlag;%displayFlag = 1;
[Dictionary,output] = KSVD(blkMatrix,param);%%%%%%%����ĵĺ���%%%%%%%%%%%%
output.D = Dictionary;

if (displayFlag)%displayFlag = 1;
    disp('finished Trainning dictionary');
end



%% denoise the image using the resulted dictionary
errT = sigma*errorFactor;
IMout=zeros(NN1,NN2);
Weight=zeros(NN1,NN2);
%blocks = im2col(Image,[NN1,NN2],[bb,bb],'sliding');
while (prod(floor((size(Image)-blockSize)/slidingDis)+1)>maxBlocksToConsider)
    slidingDis = slidingDis+1;
end
[blocks,idx] = my_im2col(Image,[blockSize,blockSize],slidingDis);

if (waitBarOn)
    newCounterForWaitBar = (param.numIteration+1)*size(blocks,2);
end


% go with jumps of 30000
for jj = 1:80:size(blocks,2)
    if (waitBarOn)
        waitbar(((param.numIteration*size(blocks,2))+jj)/newCounterForWaitBar);
    end
    jumpSize = min(jj+30000-1,size(blocks,2));
    if (reduceDC)
        vecOfMeans = mean(blocks(:,jj:jumpSize));
        blocks(:,jj:jumpSize) = blocks(:,jj:jumpSize) - repmat(vecOfMeans,size(blocks,1),1);
    end
    
    %Coefs = mexOMPerrIterative(blocks(:,jj:jumpSize),Dictionary,errT);
    Coefs = OMPerr(Dictionary,blocks(:,jj:jumpSize),errT);
    full_Coefs=full(Coefs);
    max_elem=max(max(full_Coefs));
    min_elem=min(min(full_Coefs));
    gap=(max_elem-min_elem)/100;
    full_Coefs=floor((full_Coefs-min_elem)./gap).*gap+gap/2+min_elem;
    full_Coefs=reshape(full_Coefs,[],1);
    for i = 1:size(full_Coefs,1)
        if (full_Coefs(i) >= 0.00001)
            full_Coefs(i) = 1;
        elseif(full_Coefs(i) <= -0.00001)
            full_Coefs(i) = -1;
        else
            full_Coefs(i) = 0;
        end
    end
    x=1:size(full_Coefs,1);
    %plot(x,full_Coefs);
    if (reduceDC)
        blocks(:,jj:jumpSize)= Dictionary*Coefs + ones(size(blocks,1),1) * vecOfMeans;
    else
        blocks(:,jj:jumpSize)= Dictionary*Coefs ;
    end
end

count = 1;
Weight = zeros(NN1,NN2);
IMout = zeros(NN1,NN2);
[rows,cols] = ind2sub(size(Image)-blockSize+1,idx);
for i  = 1:length(cols)
    col = cols(i); row = rows(i);        
    block =reshape(blocks(:,count),[blockSize,blockSize]);
    IMout(row:row+blockSize-1,col:col+blockSize-1)=IMout(row:row+blockSize-1,col:col+blockSize-1)+block;
    Weight(row:row+blockSize-1,col:col+blockSize-1)=Weight(row:row+blockSize-1,col:col+blockSize-1)+ones(blockSize);
    count = count+1;
end

if (waitBarOn)
    close(h);
end
IOut = (Image+0.034*sigma*IMout)./(1+0.034*sigma*Weight);

