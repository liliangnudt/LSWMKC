
function [W, gamma] = update_gamma(KH, k,omega,Kstar,alpha) 

[num, num, numker] = size(KH); %

W = zeros(num);
gamma_temp = zeros(num,1);
D = zeros(1,num);

for i = 1:num
    Kv_temp = zeros(1,num);
    for p = 1:numker
        Kv_temp = Kv_temp + omega(p) * KH(i,:,p);
    end
    D = -(2*alpha*Kstar(i,:) + Kv_temp);
    [dumb, idx] = sort(D, 2); 

    id = idx(1,2:k+2); 
    di = D(1, id); 
    W(i,id) = (di(k+1)-di)/(k*di(k+1)-sum(di(1:k))+eps); 
    gamma_temp(i) = k/2*di(k+1) - 1/2*sum(di(1:k)) - alpha;
end;
gamma = gamma_temp;







