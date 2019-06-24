function adj_fake_flat = find_row_flat(ch_idx)

%{ 
The purpose of this function is to take a flattened adjacency matrix of the
form nch*(nch-1)/2 rows and m columns (where each column is a time point),
as well as a logical array (nch x 1 in size) indicating which chs to
remove, and it returns a logical array nch*(nch-1)/2 x 1 in size saying
which rows in the flattened array it corresponds to
%}

nchs = length(ch_idx);

% Make a fake adjacency matrix and mark the appropriate rows and columns
adj_fake = zeros(nchs,nchs);
adj_fake(ch_idx,:) = 1;
adj_fake(:,ch_idx) = 1;

% Flatten this
adj_fake_flat = zeros(nchs*(nchs-1)/2,1);
count = 0;
for i = 1:nchs
    for j = 1:i-1
        count = count + 1;
        if adj_fake(j,i) == 1
            adj_fake_flat(count) = 1;
        end
    end
end

adj_fake_flat = logical(adj_fake_flat);


end