%% 显示重构误差具有上界
close all;clear;rng('shuffle')
ngs=10;N=100;
graphs={};Ws={};
for i=1:ngs
    G=gsp_sensor(N);
    graphs{i}=G;
    Ws{i}=full(G.W);
end


%%
num=100;
res=zeros(ngs,num);
sigs=zeros(ngs,1);
fv=zeros(N,1); %vector of frequencies

for j=1:ngs
    G=graphs{j};
    [u,v]=eig(full(G.L));
    for n=1:num
        fv(1:100)=sort(rand(100,1)*10,'descend');
        fv=fv/norm(fv,2); %normalization
        x=u*fv;
       %% decomposition
        layers=1;
        [results,approx,sigu]=spline_wav(layers,x,Ws{j},G,1);

        %% reconstruction
        tmpr=results;
        for i=1:layers
            tmpr{i,1}=zeros(size(tmpr{i,1}));
        end
        reconx=reconstruction(tmpr,approx);
        %record the reconstruction error
        res(j,n)=norm(x-reconx);
        sigs(j)=sigu;
    end
end

h=figure;plot(res','linewidth',1.5);hold on 
plot(sigs,'.-')
set(h, 'PaperPosition', [-0.25 -0.3 10 10]);%设置图的位置，-0.75，0.2表示图片左下角的坐标（pdf页面左下角为（0，0）），26.5，26表示图片的宽高
set(h, 'PaperSize', [9.5 9.5]); %Keep the same paper size，设置pdf纸张的大小，分别表示pdf页面的长宽
saveas(h, 'figs/app_error_sensor.pdf');
