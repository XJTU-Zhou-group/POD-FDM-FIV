% Cantilever vibration
clear
clf 
global n
n=100; % number of nodes

%space nodes
xoverl=[0:1/n:1]; 
xoverl1=xoverl(:,2:n); 


load('U_340.mat');
%fluid velocity & time delay
U=3.4;
tau=2*pi/U;
%initial condition
b=[4.730041 7.853205 10.995608 14.137166 17.278759];
y=0.005*ones(5,1);
for ii=1:1:n-1
    x=1/n*ii;
    for i=1:1:5
        gama(i)=(sinh(b(i))+sin(b(i)))/(cos(b(i))-cosh(b(i)));
        phi(i)=cos(b(i)*x)-cosh(b(i)*x)+(sin(b(i)*x)-sinh(b(i)*x))*(gama(i));
        phi2(i)=(-sin(b(i)*x)-sinh(b(i)*x)+(cos(b(i)*x)-cosh(b(i)*x))*(gama(i)))*b(i);
    end
    q(ii)=phi*y;
end
%initial displacement and velocity
w0=q; 
x0=[zeros(1,n-1) w0]; 

%FDM Stiffness Matrix
S=[7 -4 1 zeros(1,n-4)]; 
S=[S;-4 6 -4 1 zeros(1,n-5)]; 
for i=3:n-3 
    S=[S;zeros(1,i-3) 1 -4 6 -4 1 zeros(1,n-i-3)]; 
end 
S=[S;zeros(1,n-5) 1 -4 6 -4]; 
S=[S;zeros(1,n-4) 1 -4 7]; 

%calculate the POD load

snap=tuv(2:n,:);
n_snap=5000;

[PHI,D,V] = svd(snap(:,1:n_snap),'econ');

global PHI1
global r
r=2;
PHI1=PHI(:,1:r);

%POD Stiffness, Damping Matrix
SS=PHI1'*S*PHI1;
II=PHI1'*eye(n-1,n-1)*PHI1;
CC=PHI1'*eye(n-1,n-1)*PHI1;

global A
global B
global C

%parameters
a=[0.0145149308099809;0.00524150279249309;0.00151828734229666;0.0259888836042596];
beta=0.24;
K=1000;

a1=a(1)+a(2)*U;
a2=a(3)*n^4;
a3=a(4)*U^2;
a4=(1-beta)*K*n;

A=[zeros(r,r) II;-a2*SS -a1*CC]; 

B=[zeros(2*n-2,2*n-2)];
B(n-1+n/2,n/2)=-a4;

C=[zeros(r,r) zeros(r,r);-a3*II zeros(r,r)]; 

% xx0=[zeros(r,1);PHI1'*x0(n:2*n-2)'];
% % tspan=linspace(0,1.787*(n^2),300); 
% % tspan=linspace(0,5,1000);
% 
% 
% % [t,x]=ode45(@beamrhs,tspan,x0); 
% 
% tottime=500;
% tspan=linspace(0,tottime,10000);
% 
% 
% sol=dde23(@beamrhs,tau,xx0,tspan);

y=[zeros(r,1);PHI1'*x0(n:2*n-2)'];
totTime = 1000;
time1 = 5;
timeCont = 0;
tt = [];
t = 0.0;
t1 = tic;
t2 = tic;
% options=ddeset('events',@(t,y,Z,ZP) EventsL2N2L(t,y,Z,ZP,n,dG));
sol = ddensd(@beamrhs,tau,[],y,[timeCont timeCont+time1]);
tint = timeCont:0.01:timeCont+time1;
y1 = (interp1(sol.x(1,1:end),sol.y(:,1:end)',tint))';
t = [t tint(2:end)];
y = [y y1(:,2:end)];
timeCont = timeCont+time1;
while timeCont<totTime
    if rem(timeCont,10)==0
        toc(t2);
        tt = [tt toc(t2)];
        disp(['time = ',num2str(timeCont)]);
        t2 = tic;
    end
    [~,nn] = min(abs(sol.x-(timeCont-2)));
    nn1 = length(sol.x);
    sols = struct('solver','ddensd','x',sol.x(:,nn:end),'y',sol.y(:,nn:end),'yp',sol.yp(:,nn:end),'IVP',0,'history',zeros(2*n-2,1));
    sol = ddensd(@beamrhs,tau,[],sols,[timeCont timeCont+time1]);
    tint = timeCont:0.01:timeCont+time1;
    y1 = (interp1(sol.x(1,nn1-nn+2:end),sol.y(:,nn1-nn+2:end)',tint))';
    t = [t tint(2:end)];
    y = [y y1(:,2:end)];
    timeCont = timeCont+time1;
end
timeCost = toc(t1);

X_pod_all=PHI1*y(1:r,:);
V_pod_all=PHI1*y(r+1:2*r,:);

%%Post-processing
% cd('C:\Users\10645\Desktop\POD_beam vibration\ddenst\figure');

% 
figure(1)
w_real=tuv(2:n,1:25000);
imagesc(w_real);
colormap(jet);
set(gca,'XTick',[1 5000 10000 15000 20000 25000],'XTickLabel',{'0','50','100','150','200','250'});
set(gca,'YTick',[1 25 50 74 99],'YTickLabel',{'1','0.75','0.5','0.25','0'});
set(gcf,'Position',[100 100 670 200]);
h=colorbar;
% set(get(h,'title'),'string','Defletion');
xlabel('t');
ylabel('x');
set(gca,'FontName','Times New Roman','FontSize',12.5,'FontWeight','bold'); 
set(gca,'linewidth',1.5);

name1=['U_',num2str(U*10),'_',num2str(r),'_1.jpg'];
% name2=['U_',num2str(uu),'.mat'];
saveas(gcf,name1);

figure(2)
w_pod=X_pod_all(:,1:25000);
imagesc(w_pod);
colormap(jet);
set(gca,'XTick',[1 5000 10000 15000 20000 25000],'XTickLabel',{'0','50','100','150','200','250'});
set(gca,'YTick',[1 25 50 74 99],'YTickLabel',{'1','0.75','0.5','0.25','0'});
set(gcf,'Position',[100 100 670 200]);
h=colorbar;
% set(get(h,'title'),'string','Defletion');
xlabel('t');
ylabel('x');
set(gca,'FontName','Times New Roman','FontSize',12.5,'FontWeight','bold'); 
set(gca,'linewidth',1.5);

name1=['U_',num2str(U*10),'_',num2str(r),'_2.jpg'];
% name2=['U_',num2str(uu),'.mat'];
saveas(gcf,name1);


figure(3)
w_real=tuv(2:n,24800:25000);
imagesc(w_real);
colormap(jet);
set(gca,'XTick',[1 50 100 150 200],'XTickLabel',{'248','248.5','249','249.5','250'});
set(gca,'YTick',[1 25 50 74 99],'YTickLabel',{'1','0.75','0.5','0.25','0'});
set(gcf,'Position',[500 200 400 200]);
% h=colorbar;
xlabel('t');
ylabel('x');
set(gca,'FontName','Times New Roman','FontSize',12.5,'FontWeight','bold'); 
set(gca,'linewidth',1.5);

name1=['U_',num2str(U*10),'_',num2str(r),'_3.jpg'];
% name2=['U_',num2str(uu),'.mat'];
saveas(gcf,name1);

figure(4)
w_pod=X_pod_all(:,24800:25000);
imagesc(w_pod);
colormap(jet);
set(gca,'XTick',[1 50 100 150 200],'XTickLabel',{'248','248.5','249','249.5','250'});
set(gca,'YTick',[1 25 50 74 99],'YTickLabel',{'1','0.75','0.5','0.25','0'});
set(gcf,'Position',[500 200 400 200]);
% h=colorbar;
xlabel('t');
ylabel('x');
set(gca,'FontName','Times New Roman','FontSize',12.5,'FontWeight','bold'); 
set(gca,'linewidth',1.5);

name1=['U_',num2str(U*10),'_',num2str(r),'_4.jpg'];
% name2=['U_',num2str(uu),'.mat'];
saveas(gcf,name1);

L_W1=3;
L_W2=1;

figure(5)
plot(tuv(1,1:5000),tuv(n/2+1,1:5000),'g','linewidth',L_W1);
hold on
plot(tuv(1,5000:25000),tuv(n/2+1,5000:25000),'r','linewidth',L_W1);
hold on
plot(t(1:25000),X_pod_all(n/2,1:25000),'bl-.','linewidth',L_W2);
set(gcf,'Position',[0 500 1000 300]);
xlabel('t');
ylabel('w');
set(gca,'FontName','Times New Roman','FontSize',13,'FontWeight','bold'); 
set(gca,'linewidth',2);

name1=['U_',num2str(U*10),'_',num2str(r),'_5.jpg'];
% name2=['U_',num2str(uu),'.mat'];
saveas(gcf,name1);
% 
figure(6)
plot(tuv(1,24000:25000),tuv(n/2+1,24000:25000),'r','linewidth',L_W1);
hold on
plot(t(24000:25000),X_pod_all(n/2,24000:25000),'bl-.','linewidth',L_W2);
set(gcf,'Position',[1000 500 400 300]);
xlabel('t');
ylabel('w');
set(gca,'FontName','Times New Roman','FontSize',13,'FontWeight','bold'); 
set(gca,'linewidth',2);

name1=['U_',num2str(U*10),'_',num2str(r),'_6.jpg'];
% name2=['U_',num2str(uu),'.mat'];
saveas(gcf,name1);

figure(7)
plot(tuv(n/2+1,end-10000:end),tuv(n+n/2,end-10000:end),'r','linewidth',L_W1)
hold on
plot(X_pod_all(n/2,end-10000:end),V_pod_all(n/2,end-10000:end),'bl-.','linewidth',L_W2);
set(gcf,'Position',[1000 500 400 300]);
xlabel('displacement');
ylabel('velocity');
set(gca,'FontName','Times New Roman','FontSize',13,'FontWeight','bold'); 
set(gca,'linewidth',2);

name1=['U_',num2str(U*10),'_',num2str(r),'_7.jpg'];
% name2=['U_',num2str(uu),'.mat'];
saveas(gcf,name1);

% plot modes
% [PHI,B,V] = svd(snap(:,1:n_snap),'econ');
% for i=1:1:6
%     figure(i);
% plot(PHI1(:,i))
% end

function xdot=beamrhs(t,x,Z,ZP)
global A
global B
global C
global PHI1
global n
global r
KK=[Z(1:r,1); zeros(r,1)];
x_all=[PHI1*x(1:r,1);zeros(n-1,1)];
F_imp=B*x_all.^3;
FF_imp=[zeros(r,1);PHI1'*F_imp(n:2*n-2,1)];
xdot=A*x+FF_imp+C*KK; 
end