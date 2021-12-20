function [Z] = update_Z(Kstar,KH,alpha,gamma,omega)

[num,num,numker] = size(KH);

for i = 1: num  
    Kv_temp = zeros(1,num);
    for p = 1:numker
        Kv_temp = Kv_temp + omega(p) * KH(i,:,p);
    end
    
    ft = -(2*alpha*Kstar(i,:) + Kv_temp);
    
    Z_hat = -ft/2/(alpha+gamma(i));
    
    indx = 1:num;
    indx(i) = [];
    [Z(i,indx), ~] = EProjSimplex_new(Z_hat(:,indx));
end

end

