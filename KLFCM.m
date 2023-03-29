function [PKL1, eta,P,C, dist, J,P2] = KLFCM(X,Xbar, k, b,m,n,alpha,beta)
iter = 0;
[N, p] = size(X);
P = rand(N, k);
P = P./(sum(P, 2)*ones(1, k));
J_prev = inf; J = []; eta_prev=inf; eta=[];
while true
    iter = iter + 1;
    t = P.^b;
    C = ((X+beta*Xbar)'*t)'./(sum(t)'*ones(1, p)*(1+beta));
    dist = sum(X.*X, 2)*ones(1, k) + (sum(C.*C, 2)*ones(1, N))'-2*X*C'+ beta*(sum(Xbar.*Xbar, 2)*ones(1, k) + (sum(C.*C, 2)*ones(1, N))'-2*Xbar*C');

    PP=reshape(P, m, n, k);
    for i=1:k
        P1=meafilt2(PP(:,:,i));
        P2(:,i)=P1(:);
    end
    Pbar=P2./(P2*ones(k,1)*ones(1,k));
    t3=Pbar.*exp(-1.*dist./alpha);
    PKL1=P./Pbar;
    KL=sum(sum(P.*PKL1));
    
    PPP=P;
    P = t3./(sum(t3, 2)*ones(1, k));
    J_cur = sum(sum((P.^b).*dist));
    J = [J J_cur];
    
    eta_cur=norm(PPP-P,2);
    eta= [eta eta_cur];
    eta_prev = eta_cur; 
    
    if eta_cur < 1e-5
        break;
    end
    
%     if norm(J_cur-J_prev, 'fro') < 1e-3,
%         break;
%     end
   display(sprintf('#iteration: %03d, objective function: %f', iter, J_cur));
   J_prev = J_cur; 
end