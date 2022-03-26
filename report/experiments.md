### exp 1

iter 20

img 64 house

stride 2

sigma 10

block 8

![image-20210825145628106](8.25.assets/image-20210825145628106.png)

### exp 2

iter 20

img 64 lena

stride 2

sigma 10

block 8

![image-20210825175026627](8.25.assets/image-20210825175026627.png)

### exp 3

iter 20

img 64 Barbara

stride 2

sigma 10

block 8

![image-20210825203533203](8.25.assets/image-20210825203533203.png)

### exp 4

iter 10

img 64 house

stride 1

sigma 10

block 8

![image-20210827092122917](8.25.assets/image-20210827092122917.png)

### exp 5

iter 10

img 128 house

stride 4

sigma 10

block 8

<img src="8.25.assets/image-20210827112842615.png" alt="image-20210827112842615"  />

### exp 6 (DP+KSVD)

KSVD iter 10

DP iter 5

img 32 house

stride 1

sigma 10

block 8

![image-20210829100442454](8.25.assets/image-20210829100442454.png)

initial kappa in a new way

![image-20210829151921896](8.25.assets/image-20210829151921896.png)

### exp 6 (DP)

DP iter 5

img 32 house

stride 1

sigma 10

block 8

![image-20210829150033760](8.25.assets/image-20210829150033760.png)

### exp 7

DP iter 1

ksvd iter 2

img 133 taken from 512 image lena compressed to 32

stride 1

sigma 10

block 8

![image-20210911083140109](8.25.assets/image-20210911083140109.png)

### exp 8.1

DP iter 1

ksvd iter 2

img 64 taken from 512 image house compressed to 32

stride 1

sigma 10

block 8

![image-20210911102557523](8.25.assets/image-20210911102557523.png)

### exp 8.2

DP iter 1

ksvd iter 2

img 64 taken from 512 image house

stride 1

sigma 10

block 8

![image-20210918180934889](8.25.assets/image-20210918180934889.png)

### exp 9

DP iter 10

img 64 taken from 512 image house

stride 1

sigma 10

block 8

![image-20220319161945158](/Users/ember/Library/Application Support/typora-user-images/image-20220319161945158.png)

### exp 10

DP iter 10

img 64 lena

stride 1

sigma 10

block 8

![image-20220319162048614](/Users/ember/Library/Application Support/typora-user-images/image-20220319162048614.png)

### exp 11

DP iter 5

img 64 house

stride 1

sigma 10

block 8

![image-20220319162136262](/Users/ember/Library/Application Support/typora-user-images/image-20220319162136262.png)

### exp 12

DP iter 10

img 64 house

stride 1

sigma 10

block 8

### exp 13

DP iter 10

img 64 livingroom block

stride 1

sigma 10

block 8

### exp 14

DP iter 20

img 64 pepper

stride 2

sigma 10

block 8

### exp 14-2

DP iter 10

img 64 pepper

stride 1

sigma 10

block 8

### exp 15

DP iter 10

img 64 house b2

stride 2

sigma 10

block 8

### exp 16

DP iter 10

img 64 livingroom b2

stride 2

sigma 10

block 8

h=findobj(gcf,'type','image');

img=get(h,'CData');
