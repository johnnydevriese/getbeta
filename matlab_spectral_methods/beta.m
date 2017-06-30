function b = beta(fvals) 



L = 1;
dom = [-L,L];

X = chebop(@(x,u) x.*u, dom, 'dirichelet');

v = chebfun(fvals, dom);

H = @(E) chebop(@(x,u) -diff(u,2) + v(x).*u + E*x.*u, dom, 'dirichlet');

% numpts = 17; 
% emin = -.01;
% emax = .01; 
% E = chebpts(numpts,[emin,emax]);
% 
% Z=[];

% for  e = E'
%     [V,D] = eigs(H(e),1,'SR');
%     %plot(V(:,1));
%     
%     v1 = V(:,1);
%     Xv1 = X.*v1;
%     p = Xv1'*v1 / (v1'*v1);
%     Z = [Z ; p]; 
%     
% end
% 
% %computing the polarization function. 
% P = chebfun(Z,[emin,emax]); 
% 
% % P1 = diff(P); 
% % %alpha is the the derivative at P(0) 
% % alpha = P1(0)
% 
% P2 = diff(P,2); 
% %beta is the second derivative at zero 
% beta = P2(0); 
% 
% % P3 = diff(P,3); 
% % %gamma is the third derivative at zero. 
% % gamma = P3(0) 

%%%%%computing beta using sum over states%%%%

k = 3;

[V,D] = eigs(H(0),k,'SR');

%creating transition moment matrix

for ii=1:k; T(:,ii) = X*V(:,ii); end;
T = V'*T;

T = T -( T(1,1)* eye(k) );

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

Q = L * D  ;

%computing the final value of beta. 
Z = T * Y ;

beta_2 = Q * Z * 1;

beta_3 = beta_2 * 3; 

%error = abs(beta_1 - beta_2) / abs(beta_2)





% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% 
% 
% 
k = 40;

[V,D] = eigs(H(0),k,'SR');

%creating transition moment matrix

clear T; 

for ii=1:k; T(:,ii) = X*V(:,ii); end;
T = V'*T;

T = T -( T(1,1)* eye(k) );

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

beta_true = beta_2 * 3; 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


b = abs((beta_3/beta_true) - 1); 




[V,D] = eigs(H(0),2,'SR');

%computing beta max.
%this is kuzyk's beta max. 
b_kuz =  3^(1/4) * 1 / (D(2,2) - D(1,1))^(7/2) ; 

%this is jareds beta max (1 in numerator could be the length of the
%interval). 
b_max = 1 ^ 3 / ( (D(2,2) - D(1,1)) ^ (2) ) ; 

beta_norm = beta_3 / b_kuz ; 

%b = - abs(beta_norm);

RADICAL_int = beta_3 / b_max;

