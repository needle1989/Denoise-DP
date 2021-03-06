function res = expr_cluster(Y, blkSize)

[m, n] = size(Y);
% 生成DCT字典
% 其中MM是DCT字典的行数，NN是DCT字典经过分频采样后的列数
MM = 8; NN = 16;
V = sqrt(2 / MM)*cos((0:MM-1)'*(0:NN-1)*pi/MM/2); 
V(1,:) = V(1,:) / sqrt(2);
DCT=kron(V,V);  
Phi = DCT;
% convolution
blkMatrix = im22col(Y,blkSize,1);
% reduce dc
% block matrix 64*3249
reduceDC = 1;
if (reduceDC)
    vecOfMeans = mean(blkMatrix);
    blkMatrix = blkMatrix-ones(size(blkMatrix,1),1)*vecOfMeans;
end
% kmeans聚类
% img: n*p n为数据个数 p为单个数据维度 k:聚类个数
% idx: 每个样本所在类别 c:k*p k个聚类质心位置
% 一次聚类：64*800
tBlkMatrix = blkMatrix';
[idx,c] = kmeans(tBlkMatrix,800);
tc = c';
% DP,获取二次聚类
[X1 , ~ , ~, kappa, alpha] = bmtl_DP(Phi,tc);
% result
kappa(kappa<0.1) = 0;
alpha(alpha<=1) = 0;
[row,col] = size(kappa);
[blkRow,blkCol] = size(blkMatrix);
dicKappa = [];
newMatrix = [];
% 生成kappa字典，准备还原
for i = 1:col
    [temp] = find(kappa(:,i)~=0);
    [tempRow, ~] = size(temp);
    preb = [];
    for j = 1:tempRow
        [sub] = find(alpha(:,temp(j))~=0);
        [subRow, ~] = size(sub);
        subDic = [];
        for k = 1:subRow
            subDic(:,k) = Phi(:,sub(k));
        end
        hatb = subDic*(subDic'*subDic)\subDic'*tc(:,i);
        preb = preb+hatb*temp(tempRow);
    end
    dicKappa(:,i) = preb;
end
% 还原至一次聚类前的图像
for i = 1:blkCol
    newMatrix(:,i) = dicKappa(:,idx(i));
end
% recover DC
X_reshape = newMatrix + ones(size(blkMatrix,1),1) * vecOfMeans;
result = col22im(X_reshape, m, 1);
res = result;
% c1 = [];c2 = [];c3 = [];
% cc1 = 0;cc2 = 0;cc3 = 0;
% for i=1:size(idx, 1)
%     if idx(i,1)==1
%         cc1 = cc1+1;
%         c1(:,cc1) = blkMatrix(:,i);
%     elseif idx(i,1)==2
%         cc2 = cc2+1;
%         c2(:,cc2) = blkMatrix(:,i);
%     elseif idx(i,1)==3
%         cc3 = cc3+1;
%         c3(:,cc3) = blkMatrix(:,i);
%     end
% end
% [X1 , ~ , ~]=bmtl_DP(Phi,c1);
% [X2 , ~ , ~]=bmtl_DP(Phi,c2);
% [X3 , ~ , ~]=bmtl_DP(Phi,c3);