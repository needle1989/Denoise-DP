# 基于稀疏特性的图像去噪

1852141 李德涛

## 卷积

以64 * 64图像为例

将lena512 resize为64 * 64图像

block size设为8

经过col22im步长为1的卷积（distinct）得：

sqrt(64 - block size + 1)

即64 * 3249的block matrix，共3249个信号

## K means 聚类

无监督的聚类算法，对于给定的样本集，按照样本之间的距离大小，将样本集划分为K个簇；让簇内的点尽量紧密的连在一起，而让簇间的距离尽量的大。

### IO：

input

data: n * p n为数据个数 p为单个数据维度 

k:聚类个数

output

idx: 每个样本所在类别 

c:k * p k个聚类质心位置

### 缺陷：

* K值的选取不好把握
* 对于不是凸的数据集比较难收敛
* 如果各隐含类别的数据不平衡，比如各隐含类别的数据量严重失衡，或者各隐含类别的方差不同，则聚类效果不佳。
* 采用迭代方法，得到的结果只是局部最优。
* 对噪音和异常点比较的敏感。

## 二次聚类

kmeans聚类产生K个簇，即有k个聚类中心

通过DP产生新的P个聚类中心

观察kappa获取k在p聚类中心的概率值，已知非零位置的使用位及其概率

可获知每个信号在p组基中的概率

假设k个聚类中心对应的每一幅图都share相同概率

对于k个聚类中心每一个聚类点对应不同图像做加权/最小二乘恢复

## 恢复

## KSVD

K-SVD 算法为信号的线性表示找到字典

OMP算法：

给定一个过完备字典矩阵，其中它的每列表示一种原型信号的原子。输出一个信号y，它可以被表示成这些原子的稀疏线性组合。信号 y 可以被表达为 y = Dx ，或者。 

最初做法是使用DP代替OMP生成稀疏系数矩阵 coefficient matrix

字典矩阵中所谓过完备性，指的是原子的个数远远大于信号y的长度(其长度很显然是n)，即n<<k

function Dictionary = I_clearDictionary(Dictionary,CoefMatrix,Data)



设置阈值：kappa(kappa<0.1) = 0

寻找非零元：[row,col] = find(kappa~=0)

64图像 一次DP用时约为20 min

10次DP用时2h 此时row为长1099的矢量 表明kappa仍未完全收敛

15次DP 1026

30次DP 1044



字典 64 256 

将含噪图像矩阵以8 8为一块 将每一个块转化为64列向量 生成新的矩阵 对此矩阵和字典 （ksvd）利用omp算法求解该矩阵在字典下的稀疏系数coefs 重复更新直到所有列处理完成

Kappa 800 800 认为第n个块属于第m类的概率

认为第688块属于507 508 ... 513类 属于这些类的概率分别为：

![image-20220211151508189](/Users/ember/Library/Application Support/typora-user-images/image-20220211151508189.png)

alpha 256 800

得到kappa中的概率后 取这些类对应的alpha

![image-20220211151124772](/Users/ember/Library/Application Support/typora-user-images/image-20220211151124772.png)

并设置threshold，大于或者小于此值（取决于alpha是precision还是variance）认为这个alpha对应的字典原子被这个类用到了，把这个类用到的所有原子的set找出来

这样就找出来了一个类应该用到的原子的set（为整个字典的subset）

重复该过程就有多组subset对应不同的类

即每个收敛后的类对应的字典原子数组



最小二乘的封闭解：

A为当前类的sub dictionary

b为图像的patch x为稀疏系数

要得到去噪声的b, /hat b = Ax

计算公式：/hat b = a(at*a)^(-1)atb

32 * 32 image：

![image-20210829150033760](/Users/ember/Denoise-DP/report/8.25.assets/image-20210829150033760.png)

把某个类下的某个图片用该类对应subset重新表示

对于对应多个类的就分别表示，然后按照kappa对应的概率加权表示

将表示结果拉回到原始图像上完成一次更新，再更新字典，重复多次结束迭代

![image-20220312152424029](/Users/ember/Library/Application Support/typora-user-images/image-20220312152424029.png)

todo：

* 通过实验验证kmeans算法缺陷：

kmeans可能并不适用于卷积后图像块的聚类

1 增大kmeans聚类个数使其逼近单次DP聚类以验证可行性

2 比较二次聚类后cluster个数 与单次DP结果进行对比

3 取较小图块比较其具体特征（32 * 32）

* 寻找kmeans代替算法以实现相似的效果
* 通过实验室服务器同步对64图像进行DP实验

