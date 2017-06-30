clear all;
close all;



L = 10;
dom = [-L,L];
% d2 = chebop(@(x,u) -diff(u,2) + x.^2.*u, dom, 'dirichlet');
% [V,D] = eigs(d2,1,0); plot(V(:,1)); norm(V(:,1),2);
% 
% X = chebfun('x', dom);
% 
% v1 = V(:,1);
% Xv1 = X.*v1;
% p = Xv1'*v1
% 
% 
H = @(E) chebop(@(x,u) -diff(u,2) + x.^2.*u + E*x.*u, dom, 'dirichlet');

 
 X = chebfun('x', dom);

%v1 = V(:,1);
%Xv1 = X.*v1;
%p = Xv1'*v1


numpts = 5; 
emin = -.01;
emax = .01; 
E = chebpts(numpts,[emin,emax]);

%hold on;

%create an empty vector that will be filled sequentially in the for loop
%below. 
Z=[];

for  e = E'
    [V,D] = eigs(H(e),4,0);
    %plot(V(:,1));
    
    v1 = V(:,1);
    Xv1 = X.*v1;
    p = Xv1'*v1;
    Z = [Z ; p];
    
end


%computing the polarization function. 
P = chebfun(Z,[emin,emax]); 

P1 = diff(P); 
%alpha is the the derivative at P(0) 
alpha = P1(0)

P2 = diff(P,2); 
%beta is the second derivative at zero 
beta = P2(0)

P3 = diff(P,3); 
%gamma is the third derivative at zero. 
gamma = P3(0) 


optimize = fminsearch(H,0); 

optimize 



%polarization = chebfun(p, chebpts(numpts)); 
%  
% beta  = diff(chebfun(Z),2); 
% 
% beta(0)
% 
% gamma = diff(chebfun(Z),3); 
% 
% gamma(0) 


 
%hold off;