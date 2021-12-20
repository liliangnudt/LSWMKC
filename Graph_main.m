function [Kstar, Z, gamma, omega, obj] = Graph_main(KH,alpha)

global neibour

iter = 0;
MaxIter = 100;
flag = 1;

[num,num,numker] = size(KH); 

omega = ones(numker,1)/numker; %initial kernel weights
avgKer = mycombFun(KH,omega);
Kstarold = avgKer;             %initial K* = average kernel
alpha_initial = 0;             %initial alpha = 0
[Z_old, gamma] = update_gamma(KH,neibour,omega,Kstarold,alpha_initial); %initial Z/gamma
Z_initial = Z_old;

Kstar_fnorm = [];
term_1 = [];
term_2 = [];
term_3 = [];
obj = [];

while flag
    iter = iter+1;
    %--------update weight---------------
    omega = zeros(1,numker);
    delta = zeros(1,numker);
    for p =1: numker
        delta(1,p) = trace(KH(:,:,p)*Z_old);
    end
    omega = delta./(norm(delta,2));
    
    %--------update Z--------------
    [Z] = update_Z(Kstarold,KH,alpha,gamma,omega); 
    Z_old = Z;
    
    %--------update K*---------------
    [Kstar] = update_Kstar(Z);
    
    %--------calculate K* norm-------
    Kstar_fnorm(end+1) = norm(Kstar-Kstarold,'fro');
    Kstarold = Kstar;
    
    %--------covergence--------------
    if (iter>=20 && abs((Kstar_fnorm(iter)-Kstar_fnorm(iter-1)))/norm(Kstar,'fro')<1e-3) || iter>MaxIter
        flag =0;
    end
    
    %-----------obj------------------
    term_1_temp = 0;
    for p = 1:numker
        term_1_temp = term_1_temp + (- trace(omega(p) * KH(:,:,p) * Z));
    end
    term_1(end+1) = term_1_temp;
    
    for i = 1:num
        term_2_temp(i) = norm(Z(i,:),2)^2;
    end
    term_2(end+1) = gamma'*term_2_temp';
    
    term_3(end+1) = alpha.*norm(Kstar - Z, 'fro')^2;
    
end
obj = term_1+term_2+term_3;
