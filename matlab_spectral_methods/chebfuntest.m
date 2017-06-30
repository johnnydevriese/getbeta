clear all;
close all;

format 

L = 1;
dom = [-L,L];

% base ham
H = chebop(@(x,u) -diff(u,2) +  x.*u , dom, 'dirichlet');

%the hamiltonian 
H_pert = @(E) chebop(@(x,u) H*u + E*x.*u, dom, 'dirichlet');
 
%the X operator 
X = chebop(@(x,u) x.*u, dom, 'dirichlet');


%defining the interval by taking the second eigenvalue minus the first
%eigenvalue with zero electric field.

[V,D] = eigs(H,2,'sr');

interval =  .1 * (D(2,2) - D(1,1));

%checking to make sure that it is normalized. (taking the L2 norm by
%default)
Hamnorm = norm( H * V(:,1) - D(1,1) * V(:,1) ) ;
 
%Hamnorm 

%making sure the ground state is normalized
psi_norm = norm( V(:,1) ); 

%psi_norm 

%checking the relative error
error = Hamnorm / psi_norm; 

%error

%checking the polarization. 
polarization =( (X*V(:,1))' * V(:,1) )/( V(:,1)' * V(:,1) ); 

%polarization



numpts = 9; 
emin = -interval;
emax = interval; 
E = chebpts(numpts,[emin,emax]);


% % % % % testing the 'BAND GAP' e.g. the first and second eigenvalues as a
% function of the electric field %%%%
% Z=[];
% W=[];
% 
% for  e = E'
%     [V,D] = eigs(H(e),2,'sr');
%     
%     EV1 = D(1,1);
%     EV2 = D(2,2) ;
%     
%     W = [W ; EV2] ; 
%     
%     Z = [Z ; EV1];
%     
% end
% 
% %Z 
% 
% %W
% 
% %plot(E,Z); hold on; plot(E,W); hold off;  
% 
% test_Z = chebfun(Z,[emin,emax]);
% 
% test_W = chebfun(W,[emin,emax]);
% 
% plot(test_Z); hold on; plot(test_W); hold off;


%%%%%%%% computing the polarization for a give potential by varying the
%%%%%%%% electric field. 

Q=[];

for  e = E'
    [V,D] = eigs(H_pert(e),4,'SR');
    %plot(V(:,1));
    
    v1 = V(:,1);
    Xv1 = X*v1;
    p = Xv1'*v1;
    Q = [Q ; p];
    
end

P = chebfun(Q,[emin,emax]); 

P1 = diff(P); 

%alpha is the the derivative at P(0) 
alpha = P1(0);


P2 = diff(P,2); 
%beta is the second derivative at zero 
beta_1 = P2(0)

%BETA = norm(P2, inf)

beta_realtive_error = beta_1 * (emax - emin)^2 / norm(P, inf) 

%this is kuzyks beta max. 
% b_max =  (3^(1/4) * 1 )/( (D(2,2) - D(1,1))^(7/2) )
% 
% beta_int = beta / b_max 


%%%%%computing beta using sum over states%%%%

k = 20;

[V,D] = eigs(H,k,'SR');

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

X = L * D  ;

%computing the final value of beta. 
Z = T * Y ;

beta_2 = X * Z * 6

error = abs(beta_1 - beta_2) / abs(beta_2)




%%%%computing beta in a different manner. %%%%%%


% 
% [V,D] = eigs(H,1,'sr');
% 
% P = chebop(@(x,u) u - V(:,1) * ( V(:,1)' * u ) ); 
% 
% H_shift = @(alpha) chebop(@(x,u) H*(P*u) - alpha.*(P*u), dom, 'dirichlet');
% 
% H_tilde = H_shift(D(1,1));  
% 
% % psi_hat = X * V(:,1); 
% % 
% % psi_hat = psi_hat - V(:,1) * (V(:,1)' * V(:,1) )    ; 
% 
% psi_hat = X * V(:,1); 
% 
% 
% phi = H_tilde \ psi_hat
% 
% X_00 = X * V(:,1); 
% 
% X_00 = X_00' * V(:,1); 
% 
% 
% 
% beta_ver2 = (X * phi)'*phi;% - (X_00 * phi)'*phi ;
% 
% 
% % beta_ver2 = beta_ver2 * phi ;
% 
% beta_ver2 = beta_ver2 * 3 









