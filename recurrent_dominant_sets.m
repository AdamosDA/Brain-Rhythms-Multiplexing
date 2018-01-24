function [ best_groups, best_cost_function ] = recurrent_dominant_sets(WA,Iterations)
% [ best_groups, best_cost_function ] = recurrent_dominant_sets(WA,Iterations)
%      
%      Run iterative dominant sets in a recurrent fashion,
%      evaluating results based on best-cohesiveness scores.                      
%      WA -> Weighted Adjacency matrix  
%      Iterations -> Number of iterations to run
%
%      Returns the most cohessive groups, ordered by highest-achieved
%      cohesiveness scores in the suggested number of 'Iteration' loops.
%
%      Ver.20160608 
%      (C) D. Adamos, d.adamos@ieee.org, http://neuroinformatics.gr

latest_cost_function=0; 
for L=1:Iterations
rng(1); % Seeding the random generator, used for a quick demo. In order for the algorithm to realisticly converge, replace with the following line ('shuffle') and increase iterations. 
%rng('shuffle');

clear groups no_groups cost_function sorted_score sorted_cost_function;
[groups,no_groups,cost_function,f_ini]=iterative_dominant_set_extraction(WA);
[~,sorted_score]=sort(cost_function,'descend');sorted_groups=groups;for ii=1:no_groups; sorted_groups(find(groups==sorted_score(ii)))=ii; end;sorted_cost_function=cost_function(sorted_score);
%latest_cost_function % 4Debug
%sorted_cost_function % 4Debug
   if sorted_cost_function(1)>latest_cost_function(1),
        latest_groups=sorted_groups;latest_cost_function=sorted_cost_function;
        % brk=1
        % [latest_cost_function]
    end
    if (size(sorted_cost_function,1)>1 && sorted_cost_function(1)>=latest_cost_function(1) && sorted_cost_function(2)>latest_cost_function(2))
            latest_groups=sorted_groups;latest_cost_function=sorted_cost_function;
        % brk=2   
        %[latest_cost_function]
    end
    if (size(sorted_cost_function,1)>2 && sorted_cost_function(1)>=latest_cost_function(1) && sorted_cost_function(2)>=latest_cost_function(2) && sorted_cost_function(3)>latest_cost_function(3))
            latest_groups=sorted_groups;latest_cost_function=sorted_cost_function;
        % brk=3  
        % [latest_cost_function]
    end
    if (size(sorted_cost_function,1)>3 && sorted_cost_function(1)>=latest_cost_function(1) && sorted_cost_function(2)>=latest_cost_function(2) && sorted_cost_function(3)>=latest_cost_function(3) && sorted_cost_function(4)>latest_cost_function(4))
            latest_groups=sorted_groups;latest_cost_function=sorted_cost_function;
        % brk=4  
        % [latest_cost_function]
    end
end
best_groups=latest_groups;
best_cost_function=latest_cost_function;

end

