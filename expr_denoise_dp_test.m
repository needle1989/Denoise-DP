function X = expr_denoise_dp_test(Y,blockSize)

[m, n] = size(Y);
% Pn=ceil(sqrt(K));
% DCT=zeros(blockSize,Pn);
% for k=0:1:Pn-1
%     V=cos((0:1:blockSize-1)'*k*pi/Pn);
%     if k>0, V=V-mean(V); end
%     DCT(:,k+1)=V/norm(V);
% end
% DCT=kron(DCT,DCT);

MM = 8; NN = 16; % 其中MM是DCT字典的行数，NN是DCT字典经过分频采样后的列数
V = sqrt(2 / MM)*cos((0:MM-1)'*(0:NN-1)*pi/MM/2); 
V(1,:) = V(1,:) / sqrt(2);
DCT=kron(V,V);  
Phi = DCT;

% Phi = rand(64, 256
blkMatrix = im22col(Y,blockSize,1);
% [tmpR, tmpC] = size(blkMatrix);
% tmp = zeros(blockSize*blockSize, floor(tmpC/4)+1);
% n = 1;
% for i = 1:4:tmpC
%     tmp(:,n) = blkMatrix(:,i);
%     n = n+1;
% end
% blkMatrix = tmp;

reduceDC = 1;
if (reduceDC)%%reduceDC=1
    vecOfMeans = mean(blkMatrix);
    blkMatrix = blkMatrix-ones(size(blkMatrix,1),1)*vecOfMeans;
end
load('iter10.mat','X1','kappa','alpha');
X_reshape = Phi * X1 + ones(size(blkMatrix,1),1) * vecOfMeans; 
result = col22im(X_reshape, m, 1);
X = result;
end
% [M, N] = size(X);
% count = zeros(64, 64);
% times = zeros(64, 64);
% final = zeros(64, 64);
% for i=1:1:M
%     for j=1:1:N
%         m = 2*ceil(j/25)-1+ceil(i/16)-1;
%         n = 2*mod(j-1,25)+1+mod(i-1,16);
%         count(m, n) = count(m, n) + X(i, j);
%         times(m, n) = times(m, n) + 1;
%     end
% end
% for i=1:1:64
%     for j=1:1:64
%         final(i,j) = count(i,j)/times(i,j);
%     end
% end

