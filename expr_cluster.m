function res = expr_cluster(img)
% img: n*p n为数据个数 p为单个数据维度
% idx: 每个样本所在类别 c:k*p k个聚类质心位置
[idx,c] = kmeans(img,3);
res = 1;