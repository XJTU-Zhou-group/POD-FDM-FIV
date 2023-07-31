clear

UU=[340; 400; 480];

for i=1:1:3
    U_i=UU(i);
    U=UU(i);
    name1=['U_',num2str(U),'.mat'];
    load(name1);

    [U,S,V] = svd(tuv(2:100,1:10000),'econ');
    
    Lamd=diag(S);
    
    All=sum(Lamd);
    
    for i=1:99
        B(i,1)=sum(Lamd(1:i,1))/All;
        B1(i,1)=1-B(i,1);
        B1(i,2)=i;
    end
    figure;
    semilogy(B1(1:40,2),B1(1:40,1),'-bo','linewidth',2.5,'MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','r');
    xlim([1,40]);
    set(gcf,'Position',[300 500 400 300]);
    set(gca,'FontName','Times New Roman','FontSize',18,'FontWeight','bold'); 
    set(gca,'linewidth',3);
    name1=['U_',num2str(U_i),'_40.jpg'];
% name2=['U_',num2str(uu),'.mat'];
    saveas(gcf,name1);
end