function res = expr_cluster(Y, blkSize)

[m, n] = size(Y);
% 其中MM是DCT字典的行数，NN是DCT字典经过分频采样后的列数
MM = 8; NN = 16;
V = sqrt(2 / MM)*cos((0:MM-1)'*(0:NN-1)*pi/MM/2); 
V(1,:) = V(1,:) / sqrt(2);
DCT=kron(V,V);  
Phi = DCT;
blkMatrix = im22col(Y,blkSize,1);
reduceDC = 1;
if (reduceDC)
    vecOfMeans = mean(blkMatrix);
    blkMatrix = blkMatrix-ones(size(blkMatrix,1),1)*vecOfMeans;
end
kBlkMatrix = blkMatrix';

% img: n*p n为数据个数 p为单个数据维度 k:聚类个数
% idx: 每个样本所在类别 c:k*p k个聚类质心位置

[idx,c] = kmeans(kBlkMatrix,3);
c1 = [];c2 = [];c3 = [];
cc1 = 0;cc2 = 0;cc3 = 0;
for i=1:size(idx, 1)
    if idx(i,1)==1
        cc1 = cc1+1;
        c1(:,cc1) = blkMatrix(:,i);
    elseif idx(i,1)==2
        cc2 = cc2+1;
        c2(:,cc2) = blkMatrix(:,i);
    elseif idx(i,1)==3
        cc3 = cc3+1;
        c3(:,cc3) = blkMatrix(:,i);
    end
end
[X1 , ~ , ~]=bmtl_DP(Phi,c1);
[X2 , ~ , ~]=bmtl_DP(Phi,c2);
[X3 , ~ , ~]=bmtl_DP(Phi,c3);
res = 2
res = 1;