close all; clear
root='E:\codes\matlabcodes\wavelets_for_J\Jcode_opt_method\FIRSTMM_DB\';
Winfo=load([root,'FIRSTMM_DB_A.txt']);
node_graphID=load([root,'FIRSTMM_DB_graph_indicator.txt']);
edge_attr=load([root,'FIRSTMM_DB_edge_attributes.txt']);
node_attr=load([root,'FIRSTMM_DB_node_attributes.txt']);
coords_all=load([root,'FIRSTMM_DB_coordinates.txt']);

%% get W and signal
Ng=max(node_graphID); %number of graphs
% sizesOFgraph=hist(node_graphID,unique(node_graphID)); %size of each graph
gid=38; %the chosen graph
nodes=find(node_graphID==gid); %ids of nodes of the chosen graph
edges=Winfo(sum(ismember(Winfo,nodes),2)==2,:); %find the edges of the graph
disp(find(sum(ismember(Winfo,nodes),2)==1));
edge_weights=edge_attr(sum(ismember(Winfo,nodes),2)==2,:);
%edgeweights contain distance of nodes and curvature of edge, add them up

N=length(nodes);
W=zeros(N,N);
for i=1:length(edges)
    W((nodes==edges(i,1)),(nodes==edges(i,2)))=1/sum(edge_weights(i,:));
end
W=(W+W')/2;
%--graph
G=gsp_graph(W);
G.coords=coords_all(nodes,:);
xcoords=G.coords(:,1);
%--signal
signal=node_attr(nodes);
signal=2*(signal-min(signal))/(max(signal)-min(signal));
ori_signal=(xcoords-min(xcoords))/(max(xcoords)-min(xcoords))*5;
e=randn(N,1)*0.5;
signal=ori_signal+e;

param.colorbar=1;param.climits=[min(signal),max(signal)]*1.2;
eigfig=figure;
subplot(121)
gsp_plot_signal(G,signal,param);title('Original Signal')
view([-37 0])
subplot(122)
[u,v]=eig(eye(N)-W);
diagv=diag(v);
[diagv,sid]=sort(diagv);
u=u(:,sid);
plot(u'*signal);title('GFT of the Original Signal')


%% decomposition
layers=2; %number of layers
[results,approx]=spline_wav(layers,signal,W);

%% reconstruction
tmpr=results;
for i=1:layers
    tmpr{i,1}=zeros(size(tmpr{i,1}));
end
reconx=reconstruction(tmpr,approx);

%%  plot 
%----plot original and reconstructed signal

figure;
subplot(121)
gsp_plot_signal(G,signal,param);title('Original Signal');
subplot(122)
gsp_plot_signal(G,reconx,param);title('Perfectly Reconstructed Signal');

%%
coords=G.coords;
eigfig=figure;
for i=1:layers
    dxh=results{i,1};
    approx=results{i,end};
    approx(abs(approx)<1e-6)=0;
    iW=results{i,end-1};
    iK=results{i,2};
    coords=coords(diag(iK)==1,:);
    iG=gsp_graph(iW);
    iG.coords=coords;
    ncoords=(iG.coords-min(iG.coords))./(max(iG.coords)-min(iG.coords));
    subplot(1,2,i)
    gsp_plot_signal(iG,approx,param)
    view([-37 0])
    title([num2str(i),'-layer approximation'])
end


    

