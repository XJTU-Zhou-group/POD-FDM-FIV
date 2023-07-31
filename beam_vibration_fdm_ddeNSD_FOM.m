% Cantilever vibration
clear
clf 
n=100; % number of nodes


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
%write stiffness
S=[7 -4 1 zeros(1,n-4)]; 
S=[S;-4 6 -4 1 zeros(1,n-5)]; 
for i=3:n-3 
    S=[S;zeros(1,i-3) 1 -4 6 -4 1 zeros(1,n-i-3)]; 
end 
S=[S;zeros(1,n-5) 1 -4 6 -4]; 
S=[S;zeros(1,n-4) 1 -4 7]; 
global A
global B
global C

xoverl=[0:1/n:1]; 
xoverl1=xoverl(:,2:n); 

w0=q; 
% x0=[zeros(1,n-1) w0]; 

%parameters
a=[0.0145149308099809;0.00524150279249309;0.00151828734229666;0.0259888836042596];
beta=0.24;
K=1000;

for U=3.4:0.01:4.8S
% U=5.2;
tau=2*pi/U;

a1=a(1)+a(2)*U;
a2=a(3)*n^4;
a3=a(4)*U^2;
a4=(1-beta)*K*n;

A=[zeros(n-1,n-1) eye(n-1,n-1);-a2*S -a1*eye(n-1,n-1)]; 

B=[zeros(2*n-2,2*n-2)];
B(n-1+n/2,n/2)=-a4;

C=[zeros(n-1,n-1) zeros(n-1,n-1);-a3*eye(n-1,n-1)  zeros(n-1,n-1)]; 

% tspan=linspace(0,1.787*(n^2),300); 
% tspan=linspace(0,5,1000);


% [t,x]=ode45(@beamrhs,tspan,x0); 

% tspan=linspace(0,100,10000);
% sol=dde23(@beamrhs,tau,x0,tspan);

y=[zeros(n-1,1);w0'];
totTime = 500;
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
%         toc(t2);
        tt = [tt toc(t2)];
%         disp(['time = ',num2str(timeCont)]);
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
toc(t1);
timeCost = toc(t1);

% sol=ddensd(@beamrhs,tau,[],x0,[0 100]);
% tint = 0:0.01:100;
% y = deval(sol,tint);
% tuv = [tint;y];

tuv = [t;y];
% dlmwrite('tuvDDENF.csv',tuv,'precision',10);


name1=['U_',num2str(U*100),'.mat'];
save(name1,'tuv');

end

% w=sol.y(1:n-1,:)'; 
% sz=size(sol.x); 
% w=[zeros(sz(1,2),1) w zeros(sz(1,2),1)]; 
% [T,X]=meshgrid(sol.x',xoverl); 

% figure(1) 
% clf 
% set(gca,'Box','on') 
% xp=[ 0 0 -.15 -.13 -.16 -.14 -.15 0]; 
% yp=[-.037 .037 .02 .005 -.015 -.02 -.03 -.037]; 
% patch(xp,yp,'r'); 
% hold on 
% plot([0 1.1],[0 0]) 
% L=plot(xoverl,zeros(1,n+1),'k','EraseMode','xor','LineWidth',[2.5]); 
% hold on 
% xp=[ 1 1 1.15 1.13 1.16 1.14 1.15 1]; 
% yp=[-.037 .037 .015 .005 -.015 -.015 -.025 -.037]; 
% patch(xp,yp,'r'); 
% hold on 
% text(1.02,0,'x/L') 
% hold on 
% axis([-.2 1.2 -.5 .5]) 
% thandl=text(0.2,-0.35,'Press Enter to Set Initial Condition'); 
% pause 
% for i=1:sz(1,2) 
%  set(L,'YData',w(i,:)); 
%  if i==1 
%  set(thandl,'String','Press Enter to Animate'); 
%  pause 
%  end 
%  drawnow; 
% end 
% set(thandl,'String','Press Enter to Continue'); 
% pause 
% figure(2) 
% mesh(T,X,w') 
% colormap([0 0 0]); 
% view(60,30) 
% %axis([0 45 0 1 -.3 .3])
% xlabel('Dimensionless Time, \tau','rotation',-31) 
% ylabel('Dimensionless Location, x/L','rotation',12) 
% zlabel('Displacement, w(x/L,\tau)') 


figure(3)
plot(tuv(n/2+1,end-10000:end),tuv(n+n/2,end-10000:end));
figure(4)
plot(tuv(1,:),tuv(n/2+1,:));


% % max(x(:,n/2))
% % min(x(:,n/2))
% t_out=sol.x';
% v_out=sol.y(n-1+n/2,:)';
% y_out=sol.y(n/2,:)';

function xdot=beamrhs(t,x,Z,ZP)
global A
global B
global C
nn=size(A,2)/2;
KK=[Z(1:nn,1); zeros(nn,1)];
xdot=A*x+B*x.^3+C*KK; 
end

