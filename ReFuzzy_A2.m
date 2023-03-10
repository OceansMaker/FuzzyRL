function []=fuzzy_RL()
clc;
close all;
clear;
x_save=[];
t_save=[];
u1_save = [];
t1_save =[];
flag=1;  

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

A =  [0         1         -1         0          0        0      0         0
    -l1*sf1+10000*l1*v   -l1*sfw    l1*sfw    -l2*sr1+10000*l2    -l2*srw   l2*srw   0         0
    sf1/mfw-10000*v/mfw    sfw/mfw   -sfw/mfw     0         0         0    -sf2/mfw    0    
       0         0         0          0         1         -1     0         0
    -l2*sf1+10000*l2*v    -l2*sfw   l3*sfw    -l3*sr1+10000*l3    -l3*srw    l3*srw  0         0
      0          0         0       sr1/mrw+(-10000/mrw)    srw/mrw  -srw/mrw  0       -sr2/mrw
      0          0         1         0          0         0      0        0
      0          0         0         0          0         1      0        0]
  A2 = A; 

B=[  0    0
    l1    l2
   -1/mfw  0
    0      0
    l2     l3
    0     -1/mrw
    0      0
    0      0];

D=[0 0
   0 0
   0 0
   0 0
   0 0
   0 0
   -1 0
   0 -1];

[xn,un]=size(B);
[xn,wn]=size(D);

Q=0.1*diag([1 1 1 1 1 1 1 1]);
R=1*eye(2);
r=1.2;
K=zeros(un,xn); 
L=zeros(wn,xn);
N=200; 
T=.02;
x0=10*[0.5;0.8;0.5;0.8;0.8;0.8;0.8;0.8];
i1=(rand(1,100)-.5)*1000;
i2=(rand(1,100)-.5)*1000;

Dxx=[];XX=[];XU=[]; XW=[]; 
X=[x0;kron(x0',x0')';kron(x0,zeros(un,1));kron(x0,zeros(wn,1))]';

for i=1:N
    [t,X]=ode45(@mysys, [(i-1)*T,i*T],X(end,:));
    size(X);
    Dxx=[Dxx;kron(X(end,1:xn),X(end,1:xn))-kron(X(1,1:xn),X(1,1:xn))];
    XX=[XX;X(end,xn+1:xn+xn^2)-X(1,xn+1:xn+xn^2)];
    XU=[XU;X(end,xn+xn^2+1:xn+xn^2+xn*un)-X(1,xn+xn^2+1:xn+xn^2+xn*un)];
    XW=[XW;X(end,xn+xn^2+xn*un+1:end)-X(1,xn+xn^2+xn*un+1:end)];
    size(XX);
    size(XU);
    size(XW);
    x_save=[x_save;X(end,1:8)];
    t_save=[t_save;t];
end

Dxx=processing_Dxx(Dxx); 
P_old=zeros(xn);P=10*eye(xn);
it=0;
p_save=[10];
k_save=[];

while norm(P-P_old)>1e-7 && it<20  
    it=it+1;
    P_old=P;
    QK=Q+K'*R*K-r^2*L'*L;
    X2=XX*kron(eye(xn),K');
    X3=XX*kron(eye(xn),L');
    size(Dxx);
    size(-X2-XU);
    size(X3-XW);
  
    X1=[Dxx, -X2-XU, 2*r*r*(X3-XW)];
    Y=-XX*QK(:);
    pp=X1\Y;
    P=reshape_p(pp)
    p_save=[p_save,norm(P)];
    BPv1=pp(end-(xn*un+xn*wn-1):end-(xn*wn));
    BPv2=pp(end-(xn*wn-1):end);
    K=reshape(BPv1,un,xn)
    L=reshape(BPv2,wn,xn)
end

P;
sprintf(' %.8f',K)
sprintf(' %.8f',L)

K2=K;
L2=L

for i=N+1 :300
        [t,X]=ode45(@mysys, [(i-1)*T,i*T],X(end,:));
        x_save=[x_save;X(end,1:8)];

end
size(x_save);
load ('Algo1_1.mat','P20_save');
figure(1); box on; hold on;
plot([0:length(P20_save)-1],P20_save,'s--','LineWidth',1.5)
plot([0:length(p_save)-1],p_save,'o-','LineWidth',1.5)
legend('||P_2 by Algorithm 1||','||P_2 by Algorithm 2||')
xlabel('Number of iterations')
set(gca, 'FontSize', 12)

figure(2) 
plot([0:length(x_save)-1]*0.02,x_save(:,1),'r-','LineWidth',1.5);hold on;
plot([0:length(x_save)-1]*0.02,x_save(:,2),'r--','LineWidth',1.5);hold on;
plot([0:length(x_save)-1]*0.02,x_save(:,3),'g-','LineWidth',1.5);hold on;
plot([0:length(x_save)-1]*0.02,x_save(:,4),'g--','LineWidth',1.5);hold on;
plot([0:length(x_save)-1]*0.02,x_save(:,5),'b-','LineWidth',1.5);hold on;
plot([0:length(x_save)-1]*0.02,x_save(:,6),'b--','LineWidth',1.5);hold on;
plot([0:length(x_save)-1]*0.02,x_save(:,7),'m-','LineWidth',1.5);hold on;
plot([0:length(x_save)-1]*0.02,x_save(:,8),'m--','LineWidth',1.5);
grid on; box on;
legend('x_1','x_2','x_3','x_4','x_5','x_6','x_7','x_8')
xlabel('s/times')

figure();hold 
plot(t1_save, u1_save);
legend('u_2(1,1)','u_2(2,1)')

function dX=mysys(t,X)
        x=X(1:xn);
        u=zeros(1,un);
        w=zeros(1,wn);
               
        if t>=N*T   
            flag=0;
        end
        if flag==1
              u=1*[exp(-0.1*t)*sin(10*t)+0.01*sin(10*t);exp(-0.1*t)*sin(10*t)+0.1*sin(10*t)]; % constructing the
            w=1*[exp(-0.01*t)*sin(100*t)+0.01*sin(15*t);exp(-0.01*t)*sin(10*t)+0.1*sin(10*t)];  % exploration noise
        else
            u=-K*x;
            w=L*x;
        end
        
        u1_save=[u1_save u];
        t1_save=[t1_save t];
                
        dx=A*x+B*u+D*w;
        dxx=kron(x',x')';
        dux=kron(x',u')';
        dwx=kron(x',w')';
        dX=[dx;dxx;dux;dwx];
end

    function P=reshape_p(p)
        P=zeros(xn);
        ij=0;
        for i=1:xn
            
            for j=1:i
                ij=ij+1;
                P(i,j)=p(ij);
                P(j,i)=P(i,j);
            end
        end
    end

    function Dxx=processing_Dxx(Dxx)
        ij=[]; ii=[];
        for i=1:xn
            ii=[ii (i-1)*xn+i];
        end
        for i=1:xn-1
            for j=i+1:xn
                ij=[ij (i-1)*xn+j];
            end
        end
        Dxx(:,ii)=Dxx(:,ii)/2;
        Dxx(:,ij)=[];
        Dxx=Dxx*2;
    end
end





