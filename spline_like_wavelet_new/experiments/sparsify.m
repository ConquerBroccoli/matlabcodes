function sW=sparsify(N,W,Q,krL)
%sW:sparsified W
%N:number of vertices
%W:new adjacency matrix corresponding to krL
%Q:super parameter, number of samples
%krL:Kron reduction of L,N*N

sW=zeros(N,N);
I=eye(N);
%resistance distance of vertices i and j
%pinv:pseudoinverse
pkL=pinv(krL);
drg=@(i,j)((I(:,i)-I(:,j))'*pkL*(I(:,i)-I(:,j)));
triW=triu(W);
%matrix of resistance distances of i,j
DrGW=zeros(N,N);
%calculatr the sum of drgs
for i=1:N
    cids=find(triW(i,:));
    for j=cids
        DrGW(i,j)=drg(i,j);
    end
end
DrGW=DrGW.*triW;
denom=sum(sum(DrGW));
%probability matrix
pdf=DrGW/denom;
%---probability distribution----
%choose a random edge according to the distribution
edges=randsrc(Q,1,[1:N*N;reshape(pdf',[1,N*N])]); %the id of the edge
uni_edge=unique(edges);

erid=ceil(uni_edge/N); %change the id into row and column number
ecid=mod(uni_edge,N);
ecid(ecid==0)=N;

rcids=[erid,ecid];
uni_rcids=unique(rcids,'row'); %unique indices (r,c)
char_rcids=num2str(rcids); 
count_rcids=countcats(categorical(cellstr(char_rcids))); %count the numbers of each coord

erid=uni_rcids(:,1); ecid=uni_rcids(:,2);
linearids=sub2ind(size(sW),erid,ecid);
sW(linearids)=sW(linearids)+W(linearids)./(Q*pdf(linearids)).*count_rcids;
% nsW=sW;
% sW=zeros(N,N);
% for q=1:Q
%     %choose a random edge according to the distribution
% %     edge=randsrc(1,1,[1:N*N;reshape(pdf',[1,N*N])]); %the id of the edge
% edge=edges(q);
%     erid=ceil(edge/N); %change the id into row and column number
%     ecid=mod(edge,N);
%     if ecid==0
%         ecid=N;
%     end
%     sW(erid,ecid)=sW(erid,ecid)+W(erid,ecid)/(Q*pdf(erid,ecid));
% end

sW=(sW+sW')/2;
% sum(sum(sW>0))
end