function GSPylh_plot(v,W,x)
%% Command:GSPylh_plot(v,W, x)
%% v: ���㼯, N*2�ľ���,ÿһ����[0,1]X[0,1]�е�һ����
%% W: ����ͼ���ڽӾ���, �ԳƵķ�����, ��Ԫ�طǸ�, ����Ԫ��ȫΪ0
%% x: �����ڶ��㼯�ϵ�ʵ�ź�

n=size(W,1); 

%% ������任����λƽ���ı�����
v=v(:,1)+(1+1i)*v(:,2)/2;
xv=v+1i*x;

%clf; %%�����ǰͼ�δ���
hold on;

%%�û�ɫ���ƶ������ڵ�ƽ��
patch([0,0.5,1.5,1],[0,0.5,0.5,0],'k','FaceAlpha',0.1);
% patch([0,0.5,1.5,1],[0,1.5,1.5,0],'k','FaceAlpha',0.1);

%%�ú�ɫ���ƶ���
plot(real(v),imag(v),'k.')

%%�û�ɫ�߶λ��ƶ���֮��ı�
for k=1:n
    for j=k+1:n
        if W(k,j)~=0
            
%             plot([real(v(k)),real(v(j))],[imag(v(k)),imag(v(j))],'color',[0.8 0.8 0.8],'LineWidth',0.8);
            plot([real(v(k)),real(v(j))],[imag(v(k)),imag(v(j))],'Color','g','LineWidth',0.6);
%             plot(real(v(k)),imag(xv(k)),'.','markersize',10,'Color',[1 0 0]);
        end
    end
end

%%���Ƹ�����ĺ���ֵ, ����ɫ����С��0��ֵ, ��ɫ�������0��ֵ
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
            