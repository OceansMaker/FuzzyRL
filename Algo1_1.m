function []=offline_zero_sum


clear all; close all; clc;

mc = 685;
dfm = 1.35;
drm = 1.55;
mfw = 42;
mrw = 45;
pa = 1380;
sf1 = 17500;
sr1 = 21000;
sf2 = 210000;
sr2 = 210000;
sfw = 1100;
srw = 1100;
l1 = 1/mc + dfm^2/pa;
l2 = 1/mc - dfm*drm/pa;
l3 = 1/mc + drm^2/pa;
v=1e-2/pi;


A1 =  [  0         1         -1         0          0        0      0         0
    -l1*sf1+10000*l1   -l1*sfw    l1*sfw    -l2*sr1+10000*l2    -l2*srw   l2*srw   0         0
    sf1/mfw-10000/mfw    sfw/mfw   -sfw/mfw     0         0         0    -sf2/mfw    0    
       0         0         0          0         1         -1     0         0
    -l2*sf1+10000*l2    -l2*sfw   l3*sfw    -l3*sr1+10000*l3    -l3*srw    l3*srw  0         0
      0          0         0       -(10000/mrw)+sr1/mrw    srw/mrw  -srw/mrw  0       -sr2/mrw
      0          0         1         0          0         0      0        0
      0          0         0         0          0         1      0        0]
  
A2 =  [0         1         -1         0          0        0      0         0
    -l1*sf1+10000*l1*v   -l1*sfw    l1*sfw    -l2*sr1+10000*l2    -l2*srw   l2*srw   0         0
    sf1/mfw-10000*v/mfw    sfw/mfw   -sfw/mfw     0         0         0    -sf2/mfw    0    
       0         0         0          0         1         -1     0         0
    -l2*sf1+10000*l2*v    -l2*sfw   l3*sfw    -l3*sr1+10000*l3    -l3*srw    l3*srw  0         0
      0          0         0       sr1/mrw+(-10000/mrw)    srw/mrw  -srw/mrw  0       -sr2/mrw
      0          0         1         0          0         0      0        0
      0          0         0         0          0         1      0        0]


B=[  0    0
    l1    l2
   -1/mfw  0
    0      0
    l2     l3
    0     -1/mrw
    0      0
    0      0]

D=[0 0
   0 0
   0 0
   0 0
   0 0
   0 0
   -1 0
   0 -1];

Q=0.1*diag([1 1 1 1 1 1 1 1]);
R=1*eye(2);

P10_save=[10];
P20_save=[10];
P1k=zeros(8);P1=10*eye(8);
P2k=zeros(8);P2=10*eye(8);

gamma=1.2;
it=0;

while norm(P1-P1k)>1e-7  & it<200    
    
    it=it+1;
    
    P1k=P1;
    
    A11=A1-B*inv(R)*B'*P1k+gamma^(-2)*D*D'*P1k;
    Q1=Q+P1k*B*inv(R)*B'*P1k-gamma^(-2)*P1k*D*D'*P1k;
    P1=lyap(A11',Q1);
    P10_save=[P10_save;norm(P1)];
      
end

it=0;

while norm(P2-P2k)>1e-7  & it<200    
    
    it=it+1    
    P2k=P2;
    
    A22=A2-B*inv(R)*B'*P2k+gamma^(-2)*D*D'*P2k;
    Q2=Q+P2k*B*inv(R)*B'*P2k-gamma^(-2)*P2k*D*D'*P2k;
    P2=lyap(A22',Q2);
    P20_save=[P20_save;norm(P2)];
    
end

save ('Algo1_1.mat','P10_save','P20_save')
figure; box on;hold on;
plot([0:length(P10_save)-1],P10_save,'o-','LineWidth',1.5)
plot([0:length(P20_save)-1],P20_save,'o-','LineWidth',1.5)
legend('||P1||','||P2||')
xlabel('Number of iterations')
P1
P2

K1=inv(R)*B'*P1;
L1=gamma^(-2)*D'*P1;
sprintf(' %.8f',K1);
sprintf(' %.8f',L1);

K2=inv(R)*B'*P2;
L2=gamma^(-2)*D'*P2;
sprintf(' %.8f',K2);
sprintf(' %.8f',L2);















