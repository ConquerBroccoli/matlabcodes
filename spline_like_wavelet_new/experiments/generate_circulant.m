function G=generate_circulant(N,conn)
%% generate circulant graphs
%conn: connected hops,1 should be in it, lies in [1,N/2]

%%
G=gsp_ring(N);
coords=G.coords;
row1=zeros(N,1); %the first row of the circulant matrix
row1(conn+1)=1;
row1(N-conn+1)=1;
W=toeplitz(row1);
G=gsp_graph(W);
G.coords=coords;

end