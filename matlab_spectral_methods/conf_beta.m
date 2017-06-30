clear all; 
%close all; 

finvals = struct2cell(load('finalvalues.mat'));

%column vector of all ones 
%fvals = ones(5,1); 

%finvals

beta = beta_check(finvals); 

beta ;
