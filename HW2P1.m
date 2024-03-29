t=[5.565 6.62 9.33;
    6.62 7.76 10.84;
    9.33 10.84 15.06];
s=[1.484 1.763 2.488;
    1.763 2.0693 2.89;
    2.488 2.89 4.016];
f=[1.75 1 1.5;
    0.6 1.8 1.2;
    0.6 0.8 2.2];

tau=f*s*f';
g1=[1.75 .6 .6]';
g2=[1 1.8 .8]';
g3=[1.5 1.2 2.2]';

% pk2f=t*eye(3)*[g1 g2 g3]'
% pk1f=s*[g1 g2 g3]'*3.75
% cf=tau*eye(3)

pk2f=[g1 g2 g3]*t*eye(3)*[g1 g2 g3]'
pk1f=[g1 g2 g3]*s*[g1 g2 g3]'*3.75
cf=tau*eye(3)*3.75

% pk2f=t*eye(3)*[g1 g2 g3]'
% pk1f=s*[g1 g2 g3]'*3.75
% cf=tau*eye(3)
% 
% ti=t*[g1 g2 g3];
% tf=ti*[0 0 0;
%     0 0 0;
%     0 0 1]*(1/2.732)
% cf=tau*[0 0 0;
%     0 0 0;
%     0 0 1]*(1/2.737)

pk2f=[g1 g2 g3]*t*eye(3)*[g1 g2 g3]'/2
pk1f=[g1 g2 g3]*s*[g1 g2 g3]'*3.75/2.732
cf=tau*eye(3)*3.75/2.732
