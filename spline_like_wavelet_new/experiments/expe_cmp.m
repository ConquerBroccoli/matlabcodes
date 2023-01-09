%% 对比多种方法
clear;close all
addpath(genpath('G:\codes\matlab\spline_like_wavelet_new\spline_like_wavelet_new\experiments'));
addpath narangcodes/

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

%%%%%
param.colorbar=1;param.climits=[min(signal),max(signal)];
figure;
subplot(121)
GSPylh_plot(ncoords,W,ori_signal);title('Original Signal')
subplot(122)
GSPylh_plot(ncoords,W,signal);title('Noisy Signal')
figure;plot(u'*signal,'linewidth',1.5);title('GFT of the Original Signal')

%% decomposition

%%%% our method
layers=1; %number of layers
[regres,regapp]=spline_wav(layers,signal,W,G,1);

%%%% spline cited method
[litres,litapp]=spline_wav(layers,signal,W,G,2);

%%%% QMF
filter_type='approximate';filtername='mayer';
[exf_hat,QMFchannel_info,QMFColorednodes] = short_QMF_filterbank_demo_2(G,signal,filter_type,filtername);
approx_narang_recon=sum(exf_hat(:,[1,2]),2); %LL,HLchannel
QMFMSE = norm(ori_signal - approx_narang_recon);
QMFRE=QMFMSE/norm(ori_signal);


%%%% Biorth
filter_type='approximate';
[polyf_hat,Biorchannel_info,BiorColorednodes] = short_Biorth_filterbank_demo2(G,signal,filter_type);
poly_narang_recon=sum(polyf_hat(:,[1,2]),2);
biorMSE = norm(ori_signal - poly_narang_recon);
biorRE=biorMSE/norm(ori_signal);

%%  plot 

%% first 2 methods
eigfig=figure;

coords=G.coords;
ap=regres{end};
ap(abs(ap)<1e-6)=0;
iW=regres{end-1};
iK=regres{2};
coords=coords(diag(iK)==1,:);
ncoords=coords;
ncoords(:,1)=(ncoords(:,1)-min(ncoords(:,1)))./(max(ncoords(:,1))-min(ncoords(:,1)));
ncoords(:,2)=(ncoords(:,2)-min(ncoords(:,2)))./(max(ncoords(:,2))-min(ncoords(:,2)));
subplot(1,2,1)
GSPylh_plot(ncoords,iW,ap);title('regOpt')

ap=litres{end};
ap(abs(ap)<1e-6)=0;
iW=litres{end-1};
iK=litres{2};
coords=G.coords(diag(iK)==1,:);
ncoords=coords;
ncoords(:,1)=(ncoords(:,1)-min(ncoords(:,1)))./(max(ncoords(:,1))-min(ncoords(:,1)));
ncoords(:,2)=(ncoords(:,2)-min(ncoords(:,2)))./(max(ncoords(:,2))-min(ncoords(:,2)));
subplot(1,2,2)
GSPylh_plot(ncoords,iW,litapp);title('literOpt')

%% qmf graph reduction
L=full(G.L);
Kqmf=-1*ones(N,1);
Kqmf(Biorchannel_info(1).nodes)=1;
Kqmf(Biorchannel_info(2).nodes)=1;
Kqmf=diag(Kqmf);
rL=kronL(L,Kqmf,N);
rD=diag(diag(rL));
rW=rD-rL;
ncoords=G.coords(diag(Kqmf)==1,:);
ncoords(:,1)=(ncoords(:,1)-min(ncoords(:,1)))./(max(ncoords(:,1))-min(ncoords(:,1)));
ncoords(:,2)=(ncoords(:,2)-min(ncoords(:,2)))./(max(ncoords(:,2))-min(ncoords(:,2)));
qmf_x=approx_narang_recon(diag(Kqmf)==1);
bior_x=poly_narang_recon(diag(Kqmf)==1);
reqmf=norm(ori_signal(diag(Kqmf)==1)-qmf_x)/norm(ori_signal(diag(Kqmf)==1));
rebior=norm(ori_signal(diag(Kqmf)==1)-bior_x)/norm(ori_signal(diag(Kqmf)==1));
narfig=figure;
subplot(1,2,1)
GSPylh_plot(ncoords,rW,approx_narang_recon(diag(Kqmf)==1));title('graphQMF')
subplot(1,2,2)
GSPylh_plot(ncoords,rW,poly_narang_recon(diag(Kqmf)==1));title('graphBior')


%% save

set(eigfig, 'PaperPosition', [-0.75 -0.3 16.5 9]);
set(eigfig, 'PaperSize', [15 9]);
saveas(eigfig, 'figs/logo_cmp.pdf');

set(narfig, 'PaperPosition', [-0.75 -0.3 16.5 9]);
set(narfig, 'PaperSize', [15 9]);
saveas(narfig, 'figs/logo_cmp2.pdf');

    

