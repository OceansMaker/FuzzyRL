function []=fuzzy_RL()

clc;clear;close all;
x_save=[];
t_save=[];

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

A =  [  0         1         -1         0          0        0      0         0
    -l1*sf1+10000*l1   -l1*sfw    l1*sfw    -l2*sr1+10000*l2    -l2*srw   l2*srw   0         0
    sf1/mfw-10000/mfw    sfw/mfw   -sfw/mfw     0         0         0    -sf2/mfw    0    
       0         0         0          0         1         -1     0         0
    -l2*sf1+10000*l2    -l2*sfw   l3*sfw    -l3*sr1+10000*l3    -l3*srw    l3*srw  0         0
      0          0         0       -(10000/mrw)+sr1/mrw    srw/mrw  -srw/mrw  0       -sr2/mrw
      0          0         1         0          0         0      0        0
      0          0         0         0          0         1      0        0];
  A1 = A

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

K=[0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1
    0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2];
L=[0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1
    0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2];


P10= [0.4588    0.0105    0.0009    0.0547    0.0016   -0.0001    0.0421    0.0439
    0.0105    0.0209   -0.0001    0.0044    0.0020   -0.0000   -0.1215   -0.0018
    0.0009   -0.0001    0.0023    0.0009   -0.0000    0.0000    0.0102   -0.0044
    0.0547    0.0044    0.0009    0.7123    0.0130    0.0016    0.0130    0.1155
    0.0016    0.0020   -0.0000    0.0130    0.0195   -0.0000   -0.0757   -0.1162
   -0.0001   -0.0000    0.0000    0.0016   -0.0000    0.0024    0.0064    0.0106
    0.0421   -0.1215    0.0102    0.0130   -0.0757    0.0064   12.1146    0.4450
    0.0439   -0.0018   -0.0044    0.1155   -0.1162    0.0106    0.4450   11.9370];

N=200; 

x0=10*[0.5;0.8;0.5;0.8;0.8;0.8;0.80;0.8];
i1=(rand(1,100)-.5)*1000;
i2=(rand(1,100)-.5)*1000;

p1_save=[];p10_save=[];Dxx=[];XX=[];XU=[]; XW=[];
X=[x0;kron(x0',x0')';kron(x0,zeros(un,1));kron(x0,zeros(wn,1))]';

u1_save = [];
t1_save =[];

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
p_save=[];
k_save=[];

while norm(P-P_old)>1e-7 && it<20 
  
    p1_save=[p1_save;norm(P)];
    p10_save=[p10_save;norm(P10)];       
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
    BPv1=pp(end-(xn*un+xn*wn-1):end-(xn*wn));
    BPv2=pp(end-(xn*wn-1):end);
    K=reshape(BPv1,un,xn)
    L=reshape(BPv2,wn,xn)
end

P;
K;

sprintf(' %.8f',K)
sprintf(' %.8f',L)

K1=K;
L1=L

for i=N+1 :300
        [t,X]=ode45(@mysys, [(i-1)*T,i*T],X(end,:));
        x_save=[x_save;X(end,1:8)];
end

figure();hold 
plot(t1_save, u1_save);
legend('u(1,1)','u(2,1)')

load ('Algo1_1.mat','P10_save');
figure; box on;hold on;
plot([0:length(P10_save)-1],P10_save,'s--','LineWidth',1.5)
plot([0:length(p1_save)-1],p1_save,'o-','LineWidth',1.5)

legend('||P_1 by Algorithm 1||','||P_1 by Algorithm 2||')
xlabel('Number of iterations')
set(gca, 'FontSize', 12)

% figure
% plot([0:length(x_save)-1]*0.02,x_save(:,1),'r-','LineWidth',1.5);hold on;
% plot([0:length(x_save)-1]*0.02,x_save(:,2),'r--','LineWidth',1.5);hold on;
% plot([0:length(x_save)-1]*0.02,x_save(:,3),'g-','LineWidth',1.5);hold on;
% plot([0:length(x_save)-1]*0.02,x_save(:,4),'g--','LineWidth',1.5);hold on;
% plot([0:length(x_save)-1]*0.02,x_save(:,5),'b-','LineWidth',1.5);hold on;
% plot([0:length(x_save)-1]*0.02,x_save(:,6),'b--','LineWidth',1.5);hold on;
% plot([0:length(x_save)-1]*0.02,x_save(:,7),'m-','LineWidth',1.5);hold on;
% plot([0:length(x_save)-1]*0.02,x_save(:,8),'m--','LineWidth',1.5);
% % axis([-0.5,it-.5,-5,15])
% grid on; box on;
% legend('x_1','x_2','x_3','x_4','x_5','x_6','x_7','x_8')
% xlabel('s/times')

function dX=mysys(t,X)
        x=X(1:xn);
        u=zeros(un,1);
        w=zeros(wn,1);       
        
        if t>=N*T  
            flag=0;
        end
        if flag==1
            u=[exp(-0.1*t)*sin(10*t)+0.01*sin(10*t);exp(-0.1*t)*sin(10*t)+0.1*sin(10*t)]; 
            w=[exp(-0.01*t)*sin(100*t)+0.01*sin(15*t);exp(-0.01*t)*sin(10*t)+0.1*sin(10*t)];
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





