function M = get_community_affiliations(cs,gamma,test_gamma)

%{
This function gets community affiliations from the configuration similarity
vector cs. It calls a function from the Brain Connectivity Toolbox
community_louvain (written by Mika Rubinov, U Cambridge 2015-2016). 

If called with test_gamma == 1, it will run through a vector of gammas and
plots the modularity Q as a function of gamma
%}

gammas = [0.8:0.01:1.3];


if test_gamma == 0
    M=community_louvain(cs,gamma);
    
else
    
    % initialize ouput array of Qs
    Q_test = zeros(length(gammas),1);
    
    % Run through the gammas and get Q for each
    for l = 1:length(gammas)
        curr_gamma = gammas(l);
        [~,Q_test(l)] = community_louvain(cs,curr_gamma);
        
    end
    
    % Plot Q as a function of gamma
    figure
    plot(gammas,Q_test);
    
end


end