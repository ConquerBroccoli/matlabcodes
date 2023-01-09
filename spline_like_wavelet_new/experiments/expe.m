clear;close all


%% graph info
G=gsp_logo();
N=G.N;
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
Gnodes=find(coords(:,1)<200);
Snodes=find(coords(:,1)>=200 & coords(:,1)<407);
Pnodes=find(coords(:,1)>=407);
xcoords=coords(:,1);
ori_signal=(max(xcoords)-(xcoords))/(max(xcoords)-min(xcoords))*0.7+0.5;

e=randn(N,1)*0.05;
signal=ori_signal+e;
tmp=randperm(N);
signal(tmp(1:50))=0;

%%%%%
param.colorbar=1;param.climits=[min(signal),max(signal)];
figure;
subplot(121)
GSPylh_plot(ncoords,W,ori_signal);title('Original Signal')
subplot(122)
GSPylh_plot(ncoords,W,signal);title('Noisy Signal')
figure;plot(u'*signal,'linewidth',1.5);title('GFT of the Original Signal')

%% decomposition
layers=2; %number of layers
[results,approx]=spline_wav(layers,signal,W,G);

%% reconstruction
tmpr=results;
for i=1:layers
    tmpr{i,1}=zeros(size(tmpr{i,1}));
end
reconx=reconstruction(tmpr,approx);

%%  plot
%%
eigfig=figure;
coords=G.coords;
for i=1:layers
    ap=results{i,end};
    ap(abs(ap)<1e-6)=0;
    iW=results{i,end-1};
    iK=results{i,2};
    coords=coords(diag(iK)==1,:);
    iG=gsp_graph(iW);
    iG.coords=coords;
    ncoords=iG.coords;
    ncoords(:,1)=(iG.coords(:,1)-min(iG.coords(:,1)))./(max(iG.coords(:,1))-min(iG.coords(:,1)));
    ncoords(:,2)=(iG.coords(:,2)-min(iG.coords(:,2)))./(max(iG.coords(:,2))-min(iG.coords(:,2)));
    subplot(1,layers,i)
    GSPylh_plot(ncoords,iW,ap);
    title([num2str(i),'-layer approximation'])
end

%% save


set(eigfig, 'PaperPosition', [-0.75 -0.3 16.5 9]);
set(eigfig, 'PaperSize', [15 9]);
saveas(eigfig, 'figs/logo_MRAregu2liter.pdf');

    

