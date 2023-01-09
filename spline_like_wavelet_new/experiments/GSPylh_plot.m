function GSPylh_plot(v,W,x)
%% Command:GSPylh_plot(v,W, x)
%% v: 顶点集, N*2的矩阵,每一行是[0,1]X[0,1]中的一个点
%% W: 无向图的邻接矩阵, 对称的方阵列, 各元素非负, 顶点元素全为0
%% x: 定义在顶点集上的实信号

n=size(W,1); 

%% 将顶点变换到单位平行四边形中
v=v(:,1)+(1+1i)*v(:,2)/2;
xv=v+1i*x;

%clf; %%清除当前图形窗口
hold on;

%%用灰色绘制顶点所在的平面
patch([0,0.5,1.5,1],[0,0.5,0.5,0],'k','FaceAlpha',0.1);
% patch([0,0.5,1.5,1],[0,1.5,1.5,0],'k','FaceAlpha',0.1);

%%用黑色绘制顶点
plot(real(v),imag(v),'k.')

%%用灰色线段绘制顶点之间的边
for k=1:n
    for j=k+1:n
        if W(k,j)~=0
            
%             plot([real(v(k)),real(v(j))],[imag(v(k)),imag(v(j))],'color',[0.8 0.8 0.8],'LineWidth',0.8);
            plot([real(v(k)),real(v(j))],[imag(v(k)),imag(v(j))],'Color','g','LineWidth',0.6);
%             plot(real(v(k)),imag(xv(k)),'.','markersize',10,'Color',[1 0 0]);
        end
    end
end

%%绘制各顶点的函数值, 用蓝色代表小于0的值, 红色代表大于0的值
for k=1:n
    tmpx=[real(v(k)),real(xv(k))];
    tmpy=[imag(v(k)),imag(xv(k))];           
    if x(k)<0
%         plot(tmpx,tmpy,'-','Color',[0.6 0.8 0.6]);
        plot(tmpx,tmpy,'-','Color',[0 0 1]);
        hold on
        plot(real(v(k)),imag(xv(k)),'.','LineWidth',0.6,'markersize',6,'Color',[1 0 0]);
        ylim([-0.5,3])
%         plot(tmpx,tmpy,'-','Color',[0.49,1,0.63],real(v(k)),imag(xv(k)),'.','LineWidth',1.2,'markersize',8,'MarkerFaceColor',[0.49,1,0.63]);
%         orif=plot(tmpx,tmpy,'r-',real(v(k)),imag(xv(k)),'r.','LineWidth',1.2,'markersize',8);   
    elseif x(k)>0
%         plot(tmpx,tmpy,'-','Color',[0.6 0.8 0.6]);
        plot(tmpx,tmpy,'-','Color',[0.4 0.8 0.4]);
        hold on
        plot(real(v(k)),imag(xv(k)),'.','LineWidth',0.6,'markersize',6,'Color',[0 0 1]);
        ylim([-0.5,3])
%         plot(tmpx,tmpy,'b-',real(v(k)),imag(xv(k)),'b.','LineWidth',1.2,'markersize',8);    
    else
        ;
    end
end
hold off;
            