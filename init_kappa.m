% x = -5:0.1:5;
% y = -5:0.1:5;
% [X,Y] = meshgrid(x,y);
function init_kappa = init_kappa(si)
a = -10;b = 10;
x = a:(b-a)/(2*si-2):b;
y = (2*pi)^0.5*exp(-0.5.*x.*x);

e = y(si)*ones(1,si);
kappa = diag(e);
for i = 1:(si-1)
    vec = y(i)*ones(1,i);
    d = diag(vec,si-i);
    dm = diag(vec,i-si);
    kappa = kappa+d+dm;
end
vec = sum(kappa);
[M,N] = size(kappa);
B = repmat(vec,M,1);
kappa = kappa./B;
init_kappa = kappa;