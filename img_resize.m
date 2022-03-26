clear;
clc;
pic = imread('img/dataset/livingroom.png');
[row,col]=size(pic);
imshow(pic);
%[x,y] = ginput(2);    %确定图像上的两点利用ginput函数，返回值是两点的坐标
x = [1;128];
y = [26;153];
pic_1 = imcrop(pic,[x(1),y(1),abs(x(1)-x(2)),abs(y(1)-y(2))]);

%利用imcrop函数对图像进行切割，输入参数是一个定点坐标，

%从该定点出发向右abs(x(1)-x(2))，向下abs(y(1)-y(2))的区域进行切割

figure,imshow(pic_1);
res = 1;
imwrite(pic_1,'img/dataset/1.png');