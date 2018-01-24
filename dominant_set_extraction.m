function [sel_list,rest_list,ordered_list,memberships,cost_function]=dominant_set_extraction(A)
%
%  [sel_list,rest_list,ordered_list,memberships,cost_function]=dominant_set_extraction(A)
%                  A is a similarity (weighted adjacency) matrix
%  sel_list is the set of selected nodes forming the most prominent graph component  
%  rest_list is the compement of the above list
%  cost_function tabulates the corresponding 'cluster-quality' for the detected component
%  
%
%  additional outputs ---> an ordered list of all nodes and the
%  corresponding list of memberships to the dominant setordered_list,memberships 
%
%  The algorithm has been adapted from
%  PAMI, vol.29(1),2007,pp.167 ---> Dominant Sets & Pairwise Clustering
%  and can be called recursively for hierarchical extraction of clusters    
%
%  Original adaptation by N.A. Laskaris
%  Updated version here by D.A. Adamos
%
%  For citation use: 
%  Adamos et al., "In quest of the missing neuron: Spike sorting based on dominant-sets clustering",
%  Computer methods and programs in biomedicine (2012) vol.107(1), pp.28-35 
%  http://dx.doi.org/10.1016/j.cmpb.2011.10.015
% 
%  For more information see: http://neuroinformatics.gr
%


A=A-diag(diag(A)); % no self loop 
[N,N]=size(A);

X=rand(1,N);X=X/sum(X); % initialization 
f_0=X*A*X'; % initial value of the ojective function to be maximized 
X2=X;  f1=0; f2=f_0;
while (f2-f1)/f_0 > 0.00000001 % if current iteration do not improve more than a tolerance stop  
f1=f2;X2=X2.*[(X2*A)/(X2*A*X2')];f2=X2*A*X2';
end, X=X2; 
cost_function=f2;
sel_list=setdiff([1:N],find(X<eps)); 
                   % in the current implementation we delete all nodes with meaningless membership
                   % in the future maybe we can select the most contributing instead  
                   %  i.e. list=find(X>0.01*max(X))
rest_list=setdiff([1:N],sel_list);                   

%___ ordering the nodes --> according to membership valus
[memberships,ordered_list]=sort(X2);

% ---> two alternative ways for selecting the exact number of nodes 
% first,
% for i=1:N-1;ff(i)=mean(mean(A(ordered_list([end-i:end]),ordered_list([end-i:end])))); end, plot(ff)
% second,
% plot(memberships)