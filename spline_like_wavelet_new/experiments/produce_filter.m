function [K,GA]=produce_filter(W,mth)
%% produce the spline-like filters
%W: the normalized adjacency
%mth:1: proposed method; 2:cited method


%% basic information
N=length(W);
W=(W+W')/2;
[U,v]=eig(W);
diagv=diag(v);
[diagv,sid]=sort(diagv,'descend');
U=U(:,sid);

%% setting
r=1; %dimension of Z0
s=4; %dimension of Z1
J=5; %length of w

Gamma=zeros(N,J);
for i=1:J
    Gamma(:,i)=diagv.^(i-1);
end
Gr=Gamma(1:r,:);
Gs=Gamma(N+1-s:N,:);
Gc=Gamma(r+1:N-s,:);
Ur=U(:,1:r);
Us=U(:,N+1-s:N);


%% produce A B
%find a full rank s ubmatrix Urr and a full rank Uss
restids=(1:N);%randperm(N,N);
tmpUr=Ur(restids,:); %shuffle the columns
eps=1e-2;
[~,rids]=rref(tmpUr',eps); %Ur'(:,rids) is a basis for the range of Ur'
rids=restids(rids);


restids(ismember(restids,rids))=[];
tmpUs=Us(restids,:);
if rank(tmpUs)<s
    %T is not full row rank
    sprinft('the rest submatrix of Us is not full rank')
else
    [~,sids]=rref(tmpUs',eps);
    sids=restids(sids);
    %restids(ismember(restids,sids))=[];
    uN=U(:,N);
    % use the polarity to partition A B
    % randomly partition the vertex set into A and B
    %A=[restids(1:floor(length(restids)/2)),rids];
    %B=[restids(floor(length(restids)/2)+1:end),sids];
    % design a sampling pattern K
    K=ones(N,1);
    K(uN<0)=-1;K(sids)=-1; %B
    K=diag(K);
    
    
    %% optimization
    tol=1e-5;
    h_ideal=zeros(N,1);
    h_ideal(diagv>median(diagv))=1;%linspace(0,1,floor(N/2)); %desired ideal filter
    h_des=h_ideal;
    tolv=sort(uniquetol(diagv,1e-4),'descend');
    Gr=zeros(r,J);
    Gs=zeros(s,J);
    Gc=zeros(length(tolv)-r-s,J);
    for j=1:J
        Gr(:,j)=tolv(1:r).^(j-1);
        Gs(:,j)=tolv(length(tolv)-s+1:end).^(j-1);
        Gc(:,j)=tolv(r+1:length(tolv)-s).^(j-1);
    end
    n=length(tolv);
    
    
    regular_term=[zeros(N,1),Gamma(:,1:J-1)]*diag(0:(J-1));
    alpha=0.01; %super peremeter of regular term
	%weights=abs(linspace(-1,1,N));
    %weights=zeros(N,1)';
    %weights(1:30)=1;weights(end-30:end)=1;

    %% opt proposed methods
    if mth==1
        cvx_begin
        variable w(J)
        %minimize(norm(weights*(h_des-1/2*(1+Gamma*w)),inf))
        minimize(norm(h_des-1/2*(1+Gamma*w),inf)+alpha*norm(regular_term*w,2))
        %minimize(norm(h_des-1/2*(1+Gamma*w),inf))
        subject to
        [Gr;Gs]*w==[ones(r,1);-ones(s,1)] %Gr*w == ones(r,1) && Gs*w == -ones(s,1);
        [Gc;-Gc]*w<=ones(2*(n-r-s),1)-tol % Gc*w<=ones(n-r-s,1)-tol && -Gc*w<=ones(n-r-s,1)-tol
        cvx_end
    
    else
        %% literature
        cvx_begin
        variable w(J)
        minimize(norm(h_des-1/2*(1+Gamma*w),inf))
        subject to
        w'*ones(size(w))==1
        w>=1e-4
        cvx_end
    end

    
    GA=U*diag(Gamma*w)*U';
%     h=figure;
%     plot(1-diagv,1+Gamma*w)

%     set(h, 'PaperPosition', [-0.75 -0.3 15.5 14]);%设置图的位置，-0.75，0.2表示图片左下角的坐标（pdf页面左下角为（0，0）），26.5，26表示图片的宽高
%     set(h, 'PaperSize', [14 14]); %Keep the same paper size，设置pdf纸张的大小，分别表示pdf页面的长宽
%     saveas(h, ['figs/hl_logo',num2str(N),'.pdf']);
end