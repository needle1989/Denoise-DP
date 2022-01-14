 function [Dictionary,output] = KSVD_DP(...
    Data,... % an nXN matrix that contins N signals (Y), each of dimension n.
    param)
% =========================================================================
%                          K-SVD algorithm
% =========================================================================
% The K-SVD algorithm finds a dictionary for linear representation of
% signals. Given a set of signals, it searches for the best dictionary that
% can sparsely represent each signal. Detailed discussion on the algorithm
% and possible applications can be found in "The K-SVD: An Algorithm for 
% Designing of Overcomplete Dictionaries for Sparse Representation", written
% by M. Aharon, M. Elad, and A.M. Bruckstein and appeared in the IEEE Trans. 
% On Signal Processing, Vol. 54, no. 11, pp. 4311-4322, November 2006. 
% =========================================================================
% INPUT ARGUMENTS:
% Data                         an nXN matrix that contins N signals (Y), each of dimension n. 
% param                        structure that includes all required
%                                 parameters for the K-SVD execution.
%                                 Required fields are:
%    K, ...                    the number of dictionary elements to train
%    numIteration,...          number of iterations to perform.
%    errorFlag...              if =0, a fix number of coefficients is
%                                 used for representation of each signal. If so, param.L must be
%                                 specified as the number of representing atom. if =1, arbitrary number
%                                 of atoms represent each signal, until a specific representation error
%                                 is reached. If so, param.errorGoal must be specified as the allowed
%                                 error.
%    preserveDCAtom...         if =1 then the first atom in the dictionary
%                                 is set to be constant, and does not ever change. This
%                                 might be useful for working with natural
%                                 images (in this case, only param.K-1
%                                 atoms are trained).
%    (optional, see errorFlag) L,...                 % maximum coefficients to use in OMP coefficient calculations.
%    (optional, see errorFlag) errorGoal, ...        % allowed representation error in representing each signal.
%    InitializationMethod,...  mehtod to initialize the dictionary, can
%                                 be one of the following arguments: 
%                                 * 'DataElements' (initialization by the signals themselves), or: 
%                                 * 'GivenMatrix' (initialization by a given matrix param.initialDictionary).
%    (optional, see InitializationMethod) initialDictionary,...      % if the initialization method 
%                                 is 'GivenMatrix', this is the matrix that will be used.
%    (optional) TrueDictionary, ...        % if specified, in each
%                                 iteration the difference between this dictionary and the trained one
%                                 is measured and displayed.
%    displayProgress, ...      if =1 progress information is displyed. If param.errorFlag==0, 
%                                 the average repersentation error (RMSE) is displayed, while if 
%                                 param.errorFlag==1, the average number of required coefficients for 
%                                 representation of each signal is displayed.
% =========================================================================
% OUTPUT ARGUMENTS:
%  Dictionary                  The extracted dictionary of size nX(param.K).
%  output                      Struct that contains information about the current run. It may include the following fields:
%    CoefMatrix                  The final coefficients matrix (it should hold that Data equals approximately Dictionary*output.CoefMatrix.
%    ratio                       If the true dictionary was defined (in
%                                synthetic experiments), this parameter holds a vector of length
%                                param.numIteration that includes the detection ratios in each
%                                iteration).
%    totalerr                    The total representation error after each
%                                iteration (defined only if
%                                param.displayProgress=1 and
%                                param.errorFlag = 0)
%    numCoef                     A vector of length param.numIteration that
%                                include the average number of coefficients required for representation
%                                of each signal (in each iteration) (defined only if
%                                param.displayProgress=1 and
%                                param.errorFlag = 1)
% =========================================================================


%isfield(param,'displayProgress'):��ʾ����param���Ƿ���displayPrograess����������򷵻�1��û���򷵻�0
if (~isfield(param,'displayProgress'))%%%ԭ���ĳ����к���param.displayProgress = displayFlag;%displayFlag = 1;  ���Դ˾�Ҳ����ִ��
    param.displayProgress = 0;
end
totalerr(1) = 99999;%������ۻ����
if (isfield(param,'errorFlag')==0)%%%param.errorFlag = 1;   �˾�Ҳ����ִ��
    param.errorFlag = 0;
end

if (isfield(param,'TrueDictionary'))%%%param��û��TrueDictionary
    displayErrorWithTrueDictionary = 1;
    ErrorBetweenDictionaries = zeros(param.numIteration+1,1);
    ratio = zeros(param.numIteration+1,1);
else
    displayErrorWithTrueDictionary = 0;%%ִ�д˾�
	ratio = 0;%����ͷ��˵��
end
if (param.preserveDCAtom>0)  %param.preserveDCAtom = 0;
    FixedDictionaryElement(1:size(Data,1),1) = 1/sqrt(size(Data,1));
else
    FixedDictionaryElement = [];%ִ�д˾�
end
% coefficient calculation method is OMP with fixed number of coefficients
if (size(Data,2) < param.K)%K=256    size(Data,2)=249*249  �˾䲻����if����
    disp('Size of data is smaller than the dictionary size. Trivial solution...');
    Dictionary = Data(:,1:size(Data,2));
    return;
elseif (strcmp(param.InitializationMethod,'DataElements'))%%�Ƚ������ַ����Ƿ����    param.InitializationMethod =  'GivenMatrix';
    Dictionary(:,1:param.K-param.preserveDCAtom) = Data(:,1:param.K-param.preserveDCAtom);
elseif (strcmp(param.InitializationMethod,'GivenMatrix'))%% param.InitializationMethod =  'GivenMatrix';   ִ�д˾�
    Dictionary(:,1:param.K-param.preserveDCAtom) = param.initialDictionary(:,1:param.K-param.preserveDCAtom);%param.initialDictionary = DCT(:,1:param.K );%%%% ȡ��256�С�Ҳ����ȫ����ȡ��
%param.preserveDCAtom=0   param.K-param.preserveDCAtom=K=256   ��ʼ���ֵ����DCT�ֵ�
end
% reduce the components in Dictionary that are spanned by the fixed
% elements
if (param.preserveDCAtom)%param.preserveDCAtom = 0;   �˾䲻ִ��
    tmpMat = FixedDictionaryElement \ Dictionary;
    Dictionary = Dictionary - FixedDictionaryElement*tmpMat;
end




%%���������ˣ���������
%normalize the dictionary.   ���ֵ���й�һ��
Dictionary = Dictionary*diag(1./sqrt(sum(Dictionary.*Dictionary)));%64*256    *256*256�����Խ��������ĵ�����diag(1./sqrt(sum(Dictionary.*Dictionary)))  ��sum(Dictionary.*Dictionary)��Ϊ�Խ�������һ���Խǵľ���
Dictionary = Dictionary.*repmat(sign(Dictionary(1,:)),size(Dictionary,1),1); % multiply in the sign of the first element.  64*256  64*256
totalErr = zeros(1,param.numIteration);%param.numIteration = numIQV  CaW    SterOfKsvd=10 ;    %sigma=50   ����numIterOfKsvd = 10;  
    
% the K-SVD algorithm starts here.
%Kerr's changes
%blockSize=8;
blockSize=ceil(sqrt(size(Data,1)));
blkMatrix = im2col(Data,[blockSize,blockSize],'distinct');%һ�б�ʾһ����,��blockSize^2��
blkNum=size(blkMatrix,2);
blkMtx=zeros(size(blkMatrix));%��������Ѿ�����õĿ�
for iterNum = 1:param.numIteration  %param.numIteration = numIterOfKsvd=10
    % find the coefficients
    if (param.errorFlag==0)  %param.errorFlag = 1;   
        %CoefMatrix = mexOMPIterative2(Data, [FixedDictionaryElement,Dictionary],param.L);
        CoefMatrix = OMP([FixedDictionaryElement,Dictionary],Data, param.L); %size(Data,2)=249*249
    else  
        %CoefMatrix = mexOMPerrIterative(Data, [FixedDictionaryElement,Dictionary],param.errorGoal);
        %CoefMatrix = OMPerr([FixedDictionaryElement,Dictionary],Data, param.errorGoal);%%%%%%%%%%param.errorGoal = sigma*C;   ϡ�����
        
        [CoefMatrix,~,~] = bmtl_DP(Dictionary,Data);
        
        param.L = 1;
    end
    
    replacedVectorCounter = 0;
	rPerm = randperm(size(Dictionary,1));%size(Dictionary,2)=256  ����һ�¾�֪���ú������÷���(����1��256�������������û���غϵ�����)
    for j = rPerm  %j��ֵΪ��1��256���������ֵ��û���ظ��ģ�
        [betterDictionaryElement,CoefMatrix,addedNewVector] = I_findBetterDictionaryElement(Data,...%%%%%%%%�ο����ڿ�ṹ���ֵ�ѧϰ    
            [FixedDictionaryElement,Dictionary],j+size(FixedDictionaryElement,2),...
            CoefMatrix,param.L);
        Dictionary(:,j) = betterDictionaryElement;%%%%%�ѿ���
        if (param.preserveDCAtom)%param.preserveDCAtom = 0;   �˾䲻ִ��
            tmpCoef = FixedDictionaryElement\betterDictionaryElement;
            Dictionary(:,j) = betterDictionaryElement - FixedDictionaryElement*tmpCoef;
            Dictionary(:,j) = Dictionary(:,j)./sqrt(Dictionary(:,j)'*Dictionary(:,j));
        end
        replacedVectorCounter = replacedVectorCounter+addedNewVector;%%%%ʵ��֤�������w.jpgͼ�񣩣�ֵ�ۼ���һ��
    end

    
    if (iterNum>1 & param.displayProgress)%param.displayProgress = 1
        if (param.errorFlag==0)%param.errorFlag = 1;
            output.totalerr(iterNum-1) = sqrt(sum(sum((Data-[FixedDictionaryElement,Dictionary]*CoefMatrix).^2))/prod(size(Data)));
            disp(['Iteration   ',num2str(iterNum),'   Total error is: ',num2str(output.totalerr(iterNum-1))]);
        else %ִ�д˾�
            output.numCoef(iterNum-1) = length(find(CoefMatrix))/size(Data,2);%%CoefMatrix�����з�0Ԫ�صĳ���/DATE������
            disp(['Iteration   ',num2str(iterNum),'   Average number of coefficients: ',num2str(output.numCoef(iterNum-1))]);
        end
    end
    if (displayErrorWithTrueDictionary ) %displayErrorWithTrueDictionary = 0;
        [ratio(iterNum+1),ErrorBetweenDictionaries(iterNum+1)] = I_findDistanseBetweenDictionaries(param.TrueDictionary,Dictionary);%%%%%%
        disp(strcat(['Iteration  ', num2str(iterNum),' ratio of restored elements: ',num2str(ratio(iterNum+1))]));
        output.ratio = ratio;
    end
    
   Dictionary = I_clearDictionary(Dictionary,CoefMatrix(size(FixedDictionaryElement,2)+1:end,:),Data);%%%%%%%%%%size(FixedDictionaryElement,2)=0  CoefMatrix��256��
       
%     h = waitbar(0,'Denoising In Process ...');
%     param.waitBarHandle = h;
    if (isfield(param,'waitBarHandle'))
        waitbar(iterNum/param.counterForWaitBar);
    end
end

output.CoefMatrix = CoefMatrix;
Dictionary = [FixedDictionaryElement,Dictionary];%% FixedDictionaryElement = [];%ִ�д˾�



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  findBetterDictionaryElement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%���ֵ�ԭ��D�Ľⶨ��ΪU�еĵ�һ�У���ϵ������CoefMatrix�Ľⶨ��ΪV�ĵ�һ����S(1��1)�ĳ˻�
function [betterDictionaryElement,CoefMatrix,NewVectorAdded] = I_findBetterDictionaryElement(Data,Dictionary,j,CoefMatrix,numCoefUsed)
if (length(who('numCoefUsed'))==0)
    numCoefUsed = 1;
end
relevantDataIndices = find(CoefMatrix(j,:)); % the data indices that uses the j'th dictionary element.relevantDataIndices = find(Coefs(3,:));
if (length(relevantDataIndices)<1) %(length(relevantDataIndices)==0)
    ErrorMat = Data-Dictionary*CoefMatrix;
    ErrorNormVec = sum(ErrorMat.^2);
    [d,i] = max(ErrorNormVec);
    betterDictionaryElement = Data(:,i);%ErrorMat(:,i); %
    betterDictionaryElement = betterDictionaryElement./sqrt(betterDictionaryElement'*betterDictionaryElement);
    betterDictionaryElement = betterDictionaryElement.*sign(betterDictionaryElement(1));
    CoefMatrix(j,:) = 0;
    NewVectorAdded = 1; %%%%%ʵ��֤�������w.jpgͼ�񣩣�ֵ�ۼ���һ��
%     liuzhe=1  û���д˾䣬˵��ϡ������ÿһ�ж��з����Ԫ��
    return;
end
%%
NewVectorAdded = 0;
tmpCoefMatrix = CoefMatrix(:,relevantDataIndices); %��ϡ������з�0 ��ȡ����  tmpCoefMatrix�ߴ�Ϊ��256*length(relevantDataIndices)
tmpCoefMatrix(j,:) = 0;% the coeffitients of the element we now improve are not relevant.
%errors =(Data(:,relevantDataIndices) - tmpCoefMatrix); % vector of errors that we want to minimize with the new element    D:64*256     tmpCoefMatrix�ߴ�Ϊ��256*length(relevantDataIndices)  Data(:,relevantDataIndices):64*relevantDataIndices
errors =(Data(:,relevantDataIndices) - Dictionary*tmpCoefMatrix);
% % the better dictionary element and the values of beta are found using svd.
% % This is because we would like to minimize || errors - beta*element ||_F^2. 
% % that is, to approximate the matrix 'errors' with a one-rank matrix. This
% % is done using the largest singular value.
[betterDictionaryElement,singularValue,betaVector] = svds(errors,1);%%%%%%%����ȡ���˵�һ������     errors�Ĵ�СΪ;64*relevantDataIndices   M=64  N=relevantDataIndices     betterDictionaryElement*singularValue*betaVector'���ƵĿ��Ա�ʾerrors
%a=[1 2 3 4;5 6 7 8;9 10 11 12;2 4 6 7.99999]; [u,s,v]=svds(a)   u*s*v'    [u,s,v]=svds(a,1):ȡ���ĵ�һ���ɷ� 
%����svds������aΪM*N�ľ�����ôu:M*M   S:M*N(��д��M*M)   V=N*M    V'=M*N
%����svd������aΪM*N�ľ��� ��ôu:M*M   S:M*N             V=N*N    V'=N*N
%���ֵ�ԭ��D�Ľⶨ��ΪU�еĵ�һ�У���ϵ������CoefMatrix�Ľⶨ��ΪV�ĵ�һ����S(1��1)�ĳ˻�    ����Ǻ���  ���� ���ģ�����������������������������
CoefMatrix(j,relevantDataIndices) = singularValue*betaVector';% *signOfFirstElem  s*v'    [u,s,v]=svds(a,1):ȡ���ĵ�һ���ɷ� �����Դ�ʱs*v'�����СΪ 1*N����CoefMatrix(j,relevantDataIndices)ҲΪ��1*N     betterDictionaryElement:M*1,��64*1������




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  findDistanseBetweenDictionaries
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ratio,totalDistances] = I_findDistanseBetweenDictionaries(original,new)
% first, all the column in oiginal starts with positive values.
catchCounter = 0;
totalDistances = 0;
for i = 1:size(new,2)
    new(:,i) = sign(new(1,i))*new(:,i);
end
for i = 1:size(original,2)
    d = sign(original(1,i))*original(:,i);
    distances =sum ( (new-repmat(d,1,size(new,2))).^2);
    [minValue,index] = min(distances);
    errorOfElement = 1-abs(new(:,index)'*d);
    totalDistances = totalDistances+errorOfElement;
    catchCounter = catchCounter+(errorOfElement<0.01);
end
ratio = 100*catchCounter/size(original,2);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  I_clearDictionary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Dictionary = I_clearDictionary(Dictionary,CoefMatrix,Data)
T2 = 0.99;
T1 = 3;
K=size(Dictionary,1); %%K=256
Er=sum((Data-Dictionary*CoefMatrix).^2,1); % remove identical atoms��ɾ����ͬ��ԭ�ӣ�  �����   CoefMatrix(j,relevantDataIndices)�Ĵ�СΪ256*relevantDataIndices
G=Dictionary'*Dictionary; %256*256
G = G-diag(diag(G));%���磺G=magic(3)     diag(diag(G))   Ҳ���ǽ��Խǵ�Ԫ�ظ�ֵΪ0
for jj=1:1:K
    if max(G(jj,:))>T2 || length(find(abs(CoefMatrix(jj,:))>1e-7))<=T1 
        [val,pos]=max(Er);
        clearDictionary=1; %%%%%%%%%%%%%%%%%%%%%%%%��������if�������ж��ٴ�
        Er(pos(1))=0;%���������ֵ��ֵΪ0
        Dictionary(:,jj)=Data(:,pos(1))/norm(Data(:,pos(1)));%%norm(Data(:,pos(1))����������ģ   �������൱�ڹ�һ��
        G=Dictionary'*Dictionary;
        G = G-diag(diag(G));
    end
end

