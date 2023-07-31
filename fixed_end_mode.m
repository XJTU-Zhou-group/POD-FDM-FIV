clear
clc


n=100;
%fixed end beam

beta=[4.730041 7.853205 10.995608 14.137166 pi*(5.5) pi*(6.5) pi*(7.5) pi*(8.5) pi*(9.5) pi*(10.5) pi*(11.5)];


x_position=[0:1/n:1];

for i=1:1:11
    Y_r(i)=(sin(beta(i))+sinh(beta(i)))/(cos(beta(i))-cosh(beta(i)));
    Y=cos(beta(i)*x_position)-cosh(beta(i)*x_position)+Y_r(i)*(sin(beta(i)*x_position)-sinh(beta(i)*x_position));
    Y_phi(:,i)=Y';
end

% for i=1:1:6
%     figure(i);
%     plot(x_position,Y_phi(:,i))
% end



% plot(snap(30,:));

%POD_mode
r_p=10;
% load('v_3_4_snap.mat');

n_snap=5000;
n_p=100;
x_pod=[0:1/n_p:1];

%U_3.4

load('U_340.mat');

snap=tuv(2:n_p,:);

[PHI,D_pod,V_pod] = svd(snap(:,1:n_snap),'econ');

PHI1=PHI(:,1:r_p);

PHI_pod=zeros(n_p+1,r_p);
PHI_pod(2:n_p,:)=PHI1;

PHI_1=PHI_pod;

%U_4.0

load('U_400.mat');

snap=tuv(2:n_p,:);

[PHI,D_pod,V_pod] = svd(snap(:,1:n_snap),'econ');

PHI1=PHI(:,1:r_p);

PHI_pod=zeros(n_p+1,r_p);
PHI_pod(2:n_p,:)=PHI1;

PHI_2=PHI_pod;

%U_4.8

load('U_480.mat');

snap=tuv(2:n_p,:);

[PHI,D_pod,V_pod] = svd(snap(:,1:n_snap),'econ');

PHI1=PHI(:,1:r_p);

PHI_pod=zeros(n_p+1,r_p);
PHI_pod(2:n_p,:)=PHI1;

PHI_3=PHI_pod;

%     plot(x_position,Y_phi(:,i)/10,'-.k','LineWidth',3)
%     hold on
%     plot(x_pod,PHI_1(:,i),'b','LineWidth',5)
%     hold on
%     plot(x_pod,PHI_2(:,i),'-ro','LineWidth',2,'MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','r')
%     hold on
%     plot(x_pod,PHI_3(:,i),'Color',[0 1 0],'LineWidth',4)


for i=[1 2]
    figure(i);
    plot(x_position,Y_phi(:,i)/10,'-.k','LineWidth',3)
    hold on
    plot(x_pod,PHI_1(:,i),'b','LineWidth',5)
    hold on
    plot(x_pod,PHI_2(:,i),'ro','LineWidth',2,'MarkerSize',6,'MarkerEdgeColor','r')
    hold on
    plot(x_pod,PHI_3(:,i),'Color',[0 1 0],'LineWidth',4)
    hold off
    set(gca,'FontName','Times New Roman','FontSize',18,'FontWeight','bold'); 
    set(gca,'linewidth',3);
    name1=['D_m_',num2str(i),'.jpg'];
    saveas(gcf,name1);
end

for i=3
    figure(i);
    plot(x_position,Y_phi(:,i)/10,'-.k','LineWidth',3)
    hold on
    plot(x_pod,PHI_1(:,i),'b','LineWidth',4)
    hold on
    plot(x_pod,-PHI_2(:,i),'ro','LineWidth',2,'MarkerSize',6,'MarkerEdgeColor','r')
    hold on
    plot(x_pod,PHI_3(:,i),'Color',[0 1 0],'LineWidth',4)
    hold off
    set(gca,'FontName','Times New Roman','FontSize',18,'FontWeight','bold'); 
    set(gca,'linewidth',3);
    name1=['D_m_',num2str(i),'.jpg'];
    saveas(gcf,name1);
end

for i=4:1:6
    figure(i);
    plot(x_position,Y_phi(:,i)/10,'-.k','LineWidth',3)
    hold on
    plot(x_pod,-PHI_1(:,i),'b','LineWidth',4)
    hold on
    plot(x_pod,-PHI_2(:,i),'ro','LineWidth',2,'MarkerSize',6,'MarkerEdgeColor','r')
    hold on
    plot(x_pod,-PHI_3(:,i),'Color',[0 1 0],'LineWidth',4)
    hold off
    set(gca,'FontName','Times New Roman','FontSize',18,'FontWeight','bold'); 
    set(gca,'linewidth',3);
    name1=['D_m_',num2str(i),'.jpg'];
    saveas(gcf,name1);
end





% 

% for i=1:1:20
%     DMD_mode(:,i)=real(PHI_DMD(:,bindex_p(i)));
% end
