function reconx=reconstruction(results,approx)
%% reconstruction
t=size(results,1); %number of layers
for i=1:t
    dxl=approx;
    tmp=results(t-i+1,:);
    dxh=tmp{1};
    K=tmp{2};
    HINV=tmp{3};
    %upsampling
    xl=zeros(length(K),1);
    xh=zeros(length(K),1);
    xl(diag(K)==1)=dxl;
    xh(diag(K)==-1)=dxh;
    
    %filtering
    approx=HINV*(xl+xh);
end
reconx=approx;
end