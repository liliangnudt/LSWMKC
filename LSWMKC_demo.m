clear
clc
warning off;

path = './';
addpath(genpath(path));
%%
global neibour
neibour = 5;    %initial neighbours

dataName = 'YALE';
disp(dataName);
load([path,'datasets/',dataName,'_Kmatrix'],'KH','Y');
Y(Y==-1)=2;
numclass = length(unique(Y)); %cluster
numker = size(KH,3);          %view
num = size(KH,1);             %sample number

KH = kcenter(KH);
KH = knorm(KH);
%%
alpha_range = 2.^[0:1:10];

for alpha_indx = 1:length(alpha_range)
    alpha = alpha_range(alpha_indx);
    
    [Kstar,Z,gamma,omega,obj] = Graph_main(KH,alpha);
    
    Kstar = kcenter(Kstar);
    Kstar = knorm(Kstar);
    [H, ~] = eigs(Kstar, numclass, 'la');
    [res(:,alpha_indx)] = myNMIACCV2(H,Y,numclass);
    
    fprintf('ACC:%4.4f \t NMI:%4.4f \t Pur:%4.4f \t Rand:%4.4f \n',...
        [ res(1,alpha_indx),res(2,alpha_indx),res(3,alpha_indx),res(4,alpha_indx) ]);
end
%%
[~,max_indx] = max(res(1,:,:),[],'all','linear'); %max_ACC
res_opt = res(:,max_indx);

fprintf('ACC_max:%4.4f \t NMI_max:%4.4f \t Pur_max:%4.4f \t Rand_max:%4.4f \t indx:%4.0f \n',...
    [ res(1,max_indx) res(2,max_indx) res(3,max_indx) res(4,max_indx) max_indx]);


