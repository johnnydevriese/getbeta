function b = beta(fvals) 



L = 10;
dom = [-L,L];

v = chebfun(fvals, dom);

H = @(E) chebop(@(x,u) -diff(u,2) + v(x).*u + E*x.*u, dom, 'dirichlet');

numpts = 5; 
emin = -.01;
emax = .01; 
E = chebpts(numpts,[emin,emax]);

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


b = - abs(beta)