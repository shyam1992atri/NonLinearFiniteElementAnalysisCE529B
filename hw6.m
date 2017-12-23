clear all
clc
g=[3.7825 3.31 4.665;
    3.31 4.88 5.42;
    4.66 5.42 8.53];
[n,l]=eig(g)
sq_l=sqrt(l)
u=n*sq_l*n'
f=[1.75 1 1.5;
    0.6 1.8 1.2;
    0.6 0.8 2.2]
r=f*inv(u)

syms l1 l2 
w=0.5*(l1^2+l2^2+1/(l1*l2)^2-3)-0.5*(l1^-2+l2^-2+(l1*l2)^2-3);
t1=(1/l1)*diff(w,l1)
t2=(1/l2)*diff(w,l2)

pretty(t1)
pretty(t2)

l1=sq_l(1,1)
l2=sq_l(2,2)

t=zeros(3,3);

t(1,1)=eval(t1)
t(2,2)=eval(t2)

t_old=n*t*n'