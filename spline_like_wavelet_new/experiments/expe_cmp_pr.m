clear;close all
addpath(genpath('E:\codes\matlabcodes\spline_like_wavelet_new\spline_like_wavelet_new\experiments'));
addpath narangcodes/

%% graph info
G=gsp_logo();
N=G.N;
coords=G.coords;
ncoords=coords;
ncoords(:,1)=(coords(:,1)-min(coords(:,1)))./(max(coords(:,1))-min(coords(:,1)));
ncoords(:,2)=(coords(:,2)-min(coords(:,2)))./(max(coords(:,2))-min(coords(:,2)));
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


signal=ori_signal;

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
reg_reconx=reconstruction(regres,regapp);
%%%% spline cited method
[litres,litapp]=spline_wav(layers,signal,W,G,2);
lit_reconx=reconstruction(litres,litapp);
%%%% QMF
filter_type='approximate';filtername='mayer';
[exf_hat,QMFchannel_info,QMFColorednodes] = short_QMF_filterbank_demo_2(G,signal,filter_type,filtername);
approx_narang_recon=sum(exf_hat(:,[1,2,3,4]),2); %LL,HLchannel
QMFMSE = norm(ori_signal - approx_narang_recon);
QMFRE=QMFMSE/norm(ori_signal);


%%%% Biorth
filter_type='approximate';
[polyf_hat,Biorchannel_info,BiorColorednodes] = short_Biorth_filterbank_demo2(G,signal,filter_type);
poly_narang_recon=sum(polyf_hat(:,[1,2,3,4]),2);
biorMSE = norm(ori_signal - poly_narang_recon);
biorRE=biorMSE/norm(ori_signal);

%%  plot 

%% first 2 methods
eigfig=figure;
subplot(1,2,1)
GSPylh_plot(ncoords,W,reg_reconx);title('regOpt')

subplot(1,2,2)
GSPylh_plot(ncoords,W,lit_reconx);title('literOpt')

%% qmf graph reduction
rereg=norm(signal-reg_reconx)/norm(signal);
relit=norm(signal-lit_reconx)/norm(signal);
reqmf=norm(signal-approx_narang_recon)/norm(signal);
rebior=norm(signal-poly_narang_recon)/norm(signal);
narfig=figure;
subplot(1,2,1)
GSPylh_plot(ncoords,W,approx_narang_recon);title('graphQMF')
subplot(1,2,2)
GSPylh_plot(ncoords,W,poly_narang_recon);title('graphBior')


%% save

set(eigfig, 'PaperPosition', [-0.75 -0.3 16.5 9]);
set(eigfig, 'PaperSize', [15 9]);
saveas(eigfig, 'figs/logo_pr1.pdf');

set(narfig, 'PaperPosition', [-0.75 -0.3 16.5 9]);
set(narfig, 'PaperSize', [15 9]);
saveas(narfig, 'figs/logo_pr.pdf');

    

