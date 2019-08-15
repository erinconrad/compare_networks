function cs = configuration_similarity(adj_flat)

%{
This function takes a matrix of the form YxT, where Y = N*(N-1)/2 where N s
the number of channels and T is the number of times. Each of the T columns
is a flattened Yx1 vector representing the symmetric adjacency matrix.

It outputs a TxT matrix where each element i,j represents the Pearson
correlation coefficient between the ith and jth columns of the original
input matrix. 

The goal is to get the similarity between every time point.
%}

%% Make the matrix
% returns a matrix of the pairwise linear correlation coefficient between 
% each pair of columns in the input matrix 
cs = corr(adj_flat); 

%% Redefine the matrix such that cs(i,j) when i=j is 0
% Make diagonal matrix with ones on the diagonal and zeros elsewhere
D = diag(ones(size(cs,1),1));

% Make it a logical matrix
D = logical(D);

% Set elements of cs along D to 0
cs(D) = 0;


end