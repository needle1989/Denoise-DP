function result = col22im(x, k, stride)
[M, N] = size(x);
blksize = sqrt(M);
blkNum = sqrt(N);
total = zeros(k);
weight = zeros(k);
result = zeros(k);

for i =1:1:M
    for j = 1:1:N
%         num = (j-1)*blocksize+i-(j-1)*4;
%         m = floor(num/(k+1))+1;
%         n = num - (m-1)*k;
        yblk = floor((j-1)/blkNum)+1;
        xblk = j-blkNum*(yblk-1);
        yMin = (yblk-1)*stride+1;
        xMin = (xblk-1)*stride+1;
        xMax = floor((i-1)/(blksize))+1;
        yMax = i-blksize*(xMax-1);
        m = yMin-1+yMax;
        n = xMin-1+xMax;
%         m = floor((i-1)/blksize) + floor((j-1)/sqrt(N)) + 1;
%         n = mod(i-1,blksize) + mod(j-1,sqrt(N)) + 1;
        total(m, n) = total(m, n) + x(i, j);
        weight(m, n) = weight(m, n) + 1;
    end
end
for i =1:1:k
    for j = 1:1:k
        result(i,j) = total(i,j) / weight(i,j);
    end
end
end

