function reducedL=kronL(nL,J,N)
%L:non-normalized Laplacian
%J:sampling matrix
%N:preserved vertices
pids=find(diag(J==1));
%discarded vertices
disca=1:N;
disca(diag(J==1))=[];
% reducedL=nL(pids,pids)-nL(pids,disca)*inv(nL(disca,disca))*(nL(disca,pids));
reducedL=nL(pids,pids)-nL(pids,disca)*(nL(disca,disca)\nL(disca,pids));
end