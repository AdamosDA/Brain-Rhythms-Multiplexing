function [groups,no_groups,cost_function,f_ini]=iterative_dominant_set_extraction(A)
%
%  [groups,no_groups,cost_function_values,f_ini]=iterative_dominant_set_extraction(A)
%  
%  Subtractive Clustering based on the function listed below
%
%  groups contains class labels : e.g. [1 2 1 1 3 1 0 1 1] 
%   ---> zero denotes spurious nodes in the graph
%   different graph-components have been extracted sequentially 
%
%   cost_function --> tabulates the corresponding 'cluster-quality' for each component
%
%    f_ini ---> cluster quality for the whole graph taken as single component
%
%  It is based on the function
%  [sel_list,rest_list,ordered_list,memberships,cost_function]=dominant_set_extraction(A)
%  A is a similarity (weighted adjacency) matrix 
%
%  The algorithm has been adopted from
%  PAMI, vol.29(1),2007,pp.167 ---> Dominant Sets & Pairwise Clustering
%
%  Original adaptation by N.A. Laskaris
%  Updated version here by D.A. Adamos
%
%  For citation use: 
%  Adamos et al., "In quest of the missing neuron: Spike sorting based on dominant-sets clustering",
%  Computer methods and programs in biomedicine (2012) vol.107(1), pp.28-35 
%  http://dx.doi.org/10.1016/j.cmpb.2011.10.015
% 
%  For more information see: http://neurobot.bio.auth.gr/spike-sorting/
%


A=A-diag(diag(A)); % to enfroce no self loop for each node
[N,N]=size(A);

groups=zeros(1,N); rem_list=[1:N];i=1; cost_function=[];
zn=ones(1,N);zn=zn/sum(zn); f_ini=zn*A*zn'; f_current=f_ini+0.5*f_ini;
while (f_current-f_ini)*(1-isempty(rem_list))>0        %f_current >= f_ini
[sel_list,rest_list]=dominant_set_extraction(A(rem_list,rem_list));
group_new=rem_list(sel_list); groups(group_new)=i; 
z=zeros(1,N); z(group_new)=1; z=z/sum(z); f_current=z*A*z';
cost_function(i)=f_current;
rem_list=setdiff(rem_list,group_new);
i=i+1;
end
no_groups=i-1;

