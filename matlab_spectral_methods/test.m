clear all; 
close all; 

%column vector of all ones 
%fvals = ones(5,1); 
fvals = rand(5,1); 

b = beta(fvals); 

b ;

options.Display = 'iter'; 
options.MaxIter = 50;

finalvalues = fminsearch(@beta,fvals,options) 

potential = chebfun(finalvalues) ; 

plot(potential) 

save('finalvalues.mat', 'finalvalues');