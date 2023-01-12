function [results,approx,sigu]=spline_wav(t,x,W,G,mth)
%% decompose signal
%W: normalized W
%x:signal
%t: number of layers
N=length(x);
dxl=x;
results=cell(t,5);
coords=G.coords;
ncoord=coords;
ncoord(:,1)=(coords(:,1)-min(coords(:,1)))./(max(coords(:,1))-min(coords(:,1)));
ncoord(:,2)=(coords(:,2)-min(coords(:,2)))./(max(coords(:,2))-min(coords(:,2)));
nW=full(G.W);
rD=diag(sum(nW));
tar=x;%<----
for i=1:t
    xl=dxl;
    %% construct filterbank
    [K,GA]=produce_filter(W,mth);
    sigu=2/svds(eye(size(K))+K*GA,1,'smallest');
    HL=rD^(-1/2)*1/2*(eye(N)+GA)*rD^(1/2);
    HH=rD^(-1/2)*1/2*(eye(N)-GA)*rD^(1/2);
    HINV=2*rD^(-1/2)*pinv(eye(N)+K*GA)*rD^(1/2); %the synthesis filter
    
    %% filtering and sampling
    xh=HH*xl; %make it a zero-DC filterbank
    xl=HL*xl;
    
    
    dxl=xl(diag(K)==1);
    dxh=xh(diag(K)==-1);
    
    %compute relative error
    tar=tar(diag(K)==1);
    re=norm(dxl-tar)/norm(tar);
    disp(re)
    %% graph reduction
    nL=diag(sum(nW))-nW;
    tic
    rL=kronL(nL,K,N);
    rD=diag(diag(rL));
    nW=rD-rL;
    N=length(nW);
    W=rD^(-1/2)*nW*rD^(-1/2);
    ncoord=ncoord(diag(K)==1,:);
    toc

    %save results
    results(i,:)={dxh,K,HINV,W,dxl}; 
end
approx=dxl;
end
