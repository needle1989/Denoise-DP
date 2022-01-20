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
% img: n*p n为数据个数 p为单个数据维度
% idx: 每个样本所在类别 c:k*p k个聚类质心位置
[idx,c] = kmeans(blkMatrix,800);
res = 1;
[X , ~ , ~]=bmtl_DP(Phi,blkMatrix);
res = 1;