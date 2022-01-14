function result = im22col(x, blksize, stride)
[M, N] = size(x);
outM = floor((M - blksize) / stride) + 1;
outN = floor((N - blksize) / stride) + 1;
out = zeros(blksize*blksize, outM*outN);
n = 1;
for i = 1:outM
    yMin = (i-1)*stride+1;
    yMax = yMin+blksize-1;
    for j = 1:outN
        xMin = (j-1)*stride+1;
        xMax = xMin+blksize-1;
        blk_reshape = reshape(x(yMin:yMax,xMin:xMax),blksize*blksize,1);
        out(:,n) = blk_reshape;
        n = n+1;
    end
end
result = out;
end

