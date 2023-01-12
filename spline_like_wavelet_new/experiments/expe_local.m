clear;close all
rng('shuffle')

%% graph info
N=512;
G=gsp_ring(N);
coords=G.coords;
ncoords=coords;
ncoords(:,1)=(coords(:,1)-min(coords(:,1)))./(max(coords(:,1))-min(coords(:,1)));
ncoords(:,2)=(coords(:,2)-min(coords(:,2)))./(max(coords(:,2))-min(coords(:,2)));
oricoor=ncoords;
W=full(G.W);
D=diag(sum(W));
W=D^(-1/2)*W*D^(-1/2);
[u,v]=eig(eye(N)-W);
diagv=diag(v);
[diagv,sid]=sort(diagv);
u=u(:,sid);
%% graph signal
% 
ori_signal=zeros(N,1);ori_signal(1:floor(N/2))=1; %discontinuous signal
e=randn(N,1)*0.00;
signal=ori_signal+e;
param.colorbar=1;param.climits=[min(signal),max(signal)];
figure;
gsp_plot_signal(G,signal,param);title('Original Signal')


%% decomposition
layers=2; %number of layers
[results,approx]=spline_wav(layers,signal,W,G,1);

%% reconstruction
tmpr=results;
for i=1:layers
    tmpr{i,1}=zeros(size(tmpr{i,1}));
end
reconx=reconstruction(tmpr,approx);

%%  plot 

%%
eigfig=figure;
subplot(2,2,1);
GSPylh_plot(ncoords,W,signal);
title('origianl signal')
for i=1:layers
    approx=results{i,end};
    approx(abs(approx)<1e-6)=0;
    iW=results{i,end-1};
    iK=results{i,2};
    coords=coords(diag(iK)==1,:);
    iG=gsp_graph(iW);
    iG.coords=coords;
    ncoords=iG.coords;
    ncoords(:,1)=(iG.coords(:,1)-min(iG.coords(:,1)))./(max(iG.coords(:,1))-min(iG.coords(:,1)));
    ncoords(:,2)=(iG.coords(:,2)-min(iG.coords(:,2)))./(max(iG.coords(:,2))-min(iG.coords(:,2)));
    subplot(2,2,i+1)
    GSPylh_plot(ncoords,iW,approx);
    title([num2str(i),'-layer approximation'])
end
subplot(2,2,4)
reconx(abs(reconx)<1e-6)=0;
disp(norm(reconx-signal)/norm(signal))
GSPylh_plot(oricoor,W,reconx);title('reconstruction from LP channel');
%% save

set(eigfig, 'PaperPosition', [-0.75 -0.3 16.5 14]);
set(eigfig, 'PaperSize', [15 14]);
saveas(eigfig, 'figs/ring_local.pdf');

    

