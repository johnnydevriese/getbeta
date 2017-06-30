
function beta = beta_check(finvals) 

format long e; 

%this function checks the potential function of the optimization. 

L = 1;
dom = [-L,L];


X = chebop(@(x,u) x.*u, dom, 'dirichelet');
%creating a chebyshev function using the final values found using the
%optimization
v = chebfun(finvals, dom);

%The Hamiltonian using the potential chebyshev function. 
H = @(E) chebop(@(x,u) -diff(u,2)/2 + v(x).*u + E*x.*u, dom, 'dirichlet');

%plot(v);

%defining a heuristic interval by taking band gap at zero electric field
%and multiplying by a constant (0.1 in this case)
[V,D] = eigs(H(0),2,'sr');

interval =  .1 * (D(2,2) - D(1,1));

%%%%%computing beta using the derivative method%%%%%%%
numpts = 33 ; 
emin = -interval;
emax = interval; 
E = chebpts(numpts,[emin,emax]);


%intializing an array which will be filled sequentially by the for loop 
Z=[];
%for loop that fills 'Z' with the different values corresponding to varying
%the electric field. 
for  e = E'
    [V,D] = eigs(H(e),1,'SR');
    %plot(V(:,1));
    
    v1 = V(:,1);
    Xv1 = X * v1;
    p = Xv1'* v1 / (v1'*v1);
    Z = [Z ; p]; 
    
end

%computing the polarization function. 
P = chebfun(Z,[emin,emax]); 

% P1 = diff(P); 
% %alpha is the the derivative at P(0) 
% alpha = P1(0)

P2 = diff(P,2); 
%beta is the second derivative at zero 
beta = P2(0); 

% P3 = diff(P,3); 
% %gamma is the third derivative at zero. 
% gamma = P3(0) 





%%%%%computing beta using sum over states%%%%

k = 25;

[V,D] = eigs(H(0),k,'SR');

%creating transition moment matrix

for ii=1:k; T(:,ii) = X*V(:,ii); end;
T = V'*T;

T = T -( T(1,1)* eye(k) );


%checking the sum rules. 
sum_matrix = (2 * T) *( D * T ) - ((T)^2 * D ) - (D * (T)^2) ;
 
%surf(sum_matrix);

%this is the same as V' * (X * V(:,1))
R = T(:,1) ;

L = R';

%This is the same as X_00
T(1,1); 

D = D - D(1,1) * eye(k) ;

d = diag(D(2:end,2:end));

d = d.^-1; 

D(2:end,2:end) = diag(d); 

%multiply the energy level difference matrix by the vector. The 'left side'
%is complex conjugated. 
Y = D * R ;

X = L * D  ;

%computing the final value of beta. 
Z = T * Y ;

beta_2 = X * Z * 1;

beta_3 = beta_2 * 6;

%error = abs(beta_1 - beta_2) / abs(beta_2)

%b = beta_3






%beta max's--eigenvalues must be recomputed because I've overwritten the
%ones I found earlier... 

[V,D] = eigs(H(0),2,'SR');

%computing beta max.
%this is kuzyk's beta max. 
b_kuz =  3^(1/4) * 1 / (D(2,2) - D(1,1))^(7/2) ; 

%this is jareds beta max (1 in numerator could be the length of the
%interval). 
RADICAL_max = 1 ^ 3 / ( (D(2,2) - D(1,1)) ^ (2) ) ; 

beta_intrinsic = beta_3 / b_kuz ; 

%b = - abs(beta_norm);

beta_intrinsic

%B_abs_error = abs(B_sos - 1)

%j radicals calculations in comparison with kuzyks limit.
RADICAL_int = beta_3 / RADICAL_max ;






