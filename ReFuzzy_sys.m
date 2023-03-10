function Fuzzy_sys
clear,clc,close all;

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

x_save=[0.5;0.8;0.5;0.8;0.8;0.8;0.8;0.8];
u_save = [];
w_save = [];
t_save =[];
K_save=[];

A1 =  [  0         1         -1         0          0        0      0         0
    -l1*sf1+10000*l1   -l1*sfw    l1*sfw    -l2*sr1+10000*l2    -l2*srw   l2*srw   0         0
    sf1/mfw-10000/mfw    sfw/mfw   -sfw/mfw     0         0         0    -sf2/mfw    0    
       0         0         0          0         1         -1     0         0
    -l2*sf1+10000*l2    -l2*sfw   l3*sfw    -l3*sr1+10000*l3    -l3*srw    l3*srw  0         0
      0          0         0       -(10000/mrw)+sr1/mrw    srw/mrw  -srw/mrw  0       -sr2/mrw
      0          0         1         0          0         0      0        0
      0          0         0         0          0         1      0        0];
A2 =  [0         1         -1         0          0        0      0         0
    -l1*sf1+10000*l1*v   -l1*sfw    l1*sfw    -l2*sr1+10000*5*l2    -l2*srw   l2*srw   0         0
    sf1/mfw-10000*v/mfw    sfw/mfw   -sfw/mfw     0         0         0    -sf2/mfw    0    
       0         0         0          0         1         -1     0         0
    -l2*sf1+10000*l2*v    -l2*sfw   l3*sfw    -l3*sr1+10000*l3    -l3*srw    l3*srw  0         0
      0          0         0       sr1/mrw+(-10000/mrw)    srw/mrw  -srw/mrw  0       -sr2/mrw
      0          0         1         0          0         0      0        0
      0          0         0         0          0         1      0        0];

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

K1 =  [ -0.00547587,  0.00119340,  -0.00020583, -0.00002853, -0.00038729, -0.00000784,  -0.00921304, -0.00063489
        0.00557233, -0.00107277, 0.00009858,  0.00000019,  0.00052692,  -0.00009880, 0.00740761,             0];

K2 = [-0.05164353   0.00455362 -0.00049103 -0.00592751 -0.00651158 0.00343425 0.02348853 -0.00153019
0.05206647 -0.00445414  0.00038896  0.00606836  0.00671899  -0.00359308  -0.02581707   0.00000000];


L1 = [-0.02925566,  0.08434360, -0.00709414, -0.00908143, 0.05258369, -0.00441404,  -8.41286730, -0.30905189
 -0.03067345, 0.00129591, 0.00308370,  -0.08020773, 0.08065797, -0.00738134, -0.30931470, -8.28977223];

L2 = [  -0.35197331  0.08305894 -0.00744537 0.02263422 0.06859692 -0.00407208 -8.61759266 -0.34844429
-0.08378345  -0.01283675  0.00260714  -0.09508426  0.08148994  -0.00730563   -0.34924038 -8.31054704];

T=0.02;
t_end = 20;

xk=[0.5;0.8;0.5;0.8;0.8;0.8;0.8;0.8];

for i = 1:T:t_end
    
    x1=xk(1);
    
    x_save=[x_save, xk];
    
    
    if x1==0
        h1=1;
        h2=0;
    else
        h1=(sin(x1)-v*x1)/(x1*(1-v));
        h2=(x1-sin(x1))/(x1*(1-v));
    end
    
    K=h1*K1+h2*K2;
    L=h1*L1+h2*L2;
    
    u=-K*xk;
    w= L*xk;
    u_save=[u_save u];
    w_save=[w_save w];
     
     
    K_save=[K_save;K];
    
    A = h1 * A1 + h2 * A2;
    
    [t,X]=ode45(@mysys, [(i-1)*T,i*T],xk);
     
    xk = X(end,:)';
    
end


save ('Algo2.mat','x_save')
sprintf(' %.8f',K);
sprintf(' %.8f',L);

function dx=mysys(t,x)
    
    u= - K * x;
    w= - L * x;
     
     dx=A*x+B*u+D*w;     
    end

% figure(1);
% plot([0:T:T*(length(x_save)-1)],x_save(1,:),'Color',[1 0 0],'LineWidth',1.5);hold on;
% plot([0:T:T*(length(x_save)-1)],x_save(2,:),'Color',[0 0 1],'LineWidth',1.5);hold on;
% plot([0:T:T*(length(x_save)-1)],x_save(3,:),'Color',[0 1 0],'LineWidth',1.5);hold on;
% plot([0:T:T*(length(x_save)-1)],x_save(4,:),'Color',[0 0.39216 0],'LineWidth',1.5);hold on;
% plot([0:T:T*(length(x_save)-1)],x_save(5,:),'Color',[0.28235 0.23922 0.5451],'LineWidth',1.5);hold on;
% plot([0:T:T*(length(x_save)-1)],x_save(6,:),'Color',[1 0.5 0],'LineWidth',1.5);hold on;
% plot([0:T:T*(length(x_save)-1)],x_save(7,:),'Color',[1 0.5 0.5],'LineWidth',1.5);hold on;
% plot([0:T:T*(length(x_save)-1)],x_save(8,:),'Color',[1 0 1],'LineWidth',1.5);
% grid on; box on;
% legend('x_1','x_2','x_3','x_4','x_5','x_6','x_7','x_8')
% xlabel('Time/s')
% xlim([0, 6])

% axes('position',[0.3,0.25,0.4,0.2])
% plot([0:T:T*(length(x_save)-1)],x_save(1,:),'Color',[1 0 0],'LineWidth',1.5);hold on;
% plot([0:T:T*(length(x_save)-1)],x_save(2,:),'Color',[0 0 1],'LineWidth',1.5);hold on;
% plot([0:T:T*(length(x_save)-1)],x_save(3,:),'Color',[0 1 0],'LineWidth',1.5);hold on;
% plot([0:T:T*(length(x_save)-1)],x_save(4,:),'Color',[0 0.39216 0],'LineWidth',1.5);hold on;
% plot([0:T:T*(length(x_save)-1)],x_save(5,:),'Color',[0.28235 0.23922 0.5451],'LineWidth',1.5);hold on;
% plot([0:T:T*(length(x_save)-1)],x_save(6,:),'Color',[1 0.5 0],'LineWidth',1.5);hold on;
% plot([0:T:T*(length(x_save)-1)],x_save(7,:),'Color',[1 0.5 0.5],'LineWidth',1.5);hold on;
% plot([0:T:T*(length(x_save)-1)],x_save(8,:),'Color',[1 0 1],'LineWidth',1.5)
% grid on; box on;
% xlim([min(1.8),max(3.2)]);
% ylim([min(-0.5),max(0.5)]);

figure();hold 
plot([0:T:T*(length(u_save)-1)], u_save,'LineWidth',1.5);
legend('u(1)','u(2)')
box on;
xlabel('Time/s')
xlim([0, 5])
set(gca, 'FontSize', 12)

figure();hold 
plot([0:T:T*(length(w_save)-1)], w_save,'LineWidth',1.5);
legend('w(1)','w(2)')
box on;
xlabel('Time/s')
xlim([0, 5])
ylim([-8, 4])
set(gca, 'FontSize', 12)

end