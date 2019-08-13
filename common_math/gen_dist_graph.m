function A = gen_dist_graph(pt,whichPt)

%{
This function generates an adjacency matrix weighted by 1/the distance
between electrode pairs
%}



locs = pt(whichPt).new_elecs.locs;
n = size(locs,1);

A = zeros(n,n);

% Reconstruct upper triangular matrix
for i = 1:n
    for j = 1:i-1
        A(i,j) = 1/norm(locs(i,:)-locs(j,:));
    end
end

% Reflect it over the diagonal to make symmetric matrix
A = A + A';

% Find anything that is inf and set it to zero (this happens when 2
% electrodes are mistakenly assigned the same location, as is the case in
% patient 8 'LG61' and 'MST4')
A(A==inf) = 0;

if 0
    figure
    imagesc(A)
    colorbar
    pause
    close(gcf)
end

end