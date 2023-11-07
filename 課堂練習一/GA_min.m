clear;
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 參數設定
cost_function=@M1;
sol_space=nncopy([-600 600],10,1);
pop_size=30;
bits=30; % 每個變數的位元數
survival_rate=0.3; % 0.5表示10個coromosomes中，最優的5個直接進入父母池
mut_rate=0.01; % 子代chromosomes的總位元數*mut_rate即為突變位元數
max_gen=200;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%主程式
m=size(sol_space,1); % m：變數個數
r=sol_space(:,2)-sol_space(:,1); % r：解空間範圍差(m*1矩陣)
keep=ceil(survival_rate*pop_size); % 從最優的keep個直接作為父母池

% Initialization
bin=round(rand(m*bits,pop_size)); % 以亂數產生(m*bits)*pop_size的二進制位元
quant=(0.5.^[1:bits])/sum(0.5.^[1:bits]); % quantization levels normalized

% Evaluation
q=size(bin,2); % 確認pop_size
Chi_En_01=reshape(bin,bits,[]); % 將二進制位元排列成每行僅含一個變數
Chi_En_02=quant*Chi_En_01; % 將二進制位元轉換成十進制位元(介於0~1之間)
Chi_En_03=Chi_En_02.*nncopy(r',1,q)+nncopy(sol_space(:,1)',1,q); % 將十進位數轉換成unnormalized
x=reshape(Chi_En_03,m,[]); % 將候選解排列成m*pop_size
for i=1:pop_size
    f(1,i)=cost_function(x(:,i));
end
[f,I]=sort(f,'ascend');
x=x(:,I); % 將候選解依f由小至大排列
bin=bin(:,I); % 將bin依f由小至大排列
F(1)=min(f);

for i=2:max_gen
    
    % Selection
    M=ceil((pop_size-keep)/2); % 交配數(族群總數扣掉keep數後的一半)
    prob=fliplr([1:keep]/sum([1:keep])); % 依序產生 keep/(1+2+...+keep)、...、1/(1+2+...+keep)數列
    odds=[0 cumsum(prob(1:keep))]; % 產生[0 keep/(1+2+...+keep) ... 1]的機率分佈向量
    pick1=rand(1,M); % 1*M 的亂數
    pick2=rand(1,M); % 1*M 的亂數
    ic=1;
    while ic<=M
        for id=2:keep+1
            if pick1(ic)<=odds(id) & pick1(ic)>odds(id-1)
                ma(ic)=id-1; % 媽媽(1*M 指標)
            end
            if pick2(ic)<=odds(id) & pick2(ic)>odds(id-1)
                pa(ic)=id-1; % 爸爸(1*M 指標)
            end
        end
        ic=ic+1;
    end
    
    % Crossover
    ix=1:2:(M*2-1); % ix=1,3,5,..
    xp=ceil(rand(1,M)*(bits*m-1)); % 1*M 的交配點(介於1~(bits*m-1)的整數)
    bin(:,keep+ix)=[bin(1:xp,ma);bin(xp+1:bits*m,pa)]; % 接續在survival後的第一批子代
    bin(:,keep+ix+1)=[bin(1:xp,pa);bin(xp+1:bits*m,ma)]; % 接續在survival後的第二批子代
    
    % Mutation
    nmut=ceil(m*bits*(2*M)*mut_rate); % 突變位元數
    mrow=ceil(rand(1,nmut)*(bits*m)); % 選列(1~(bits*m)的整數值)
    mcol=keep+ceil(rand(1,nmut)*(2*M)); % 選行(keep+1~keep+2M，僅針對子代進行突變)
    for j=1:nmut
        bin(mrow(j),mcol(j))=abs(bin(mrow(j),mcol(j))-1);
    end
    
    % Evaluation
    q=size(bin,2); % 確認pop_size
    Chi_En_01=reshape(bin,bits,[]); % 將二進制位元排列成每行僅含一個變數
    Chi_En_02=quant*Chi_En_01; % 將二進制位元轉換成十進制位元(介於0~1之間)
    Chi_En_03=Chi_En_02.*nncopy(r',1,q)+nncopy(sol_space(:,1)',1,q); % 將十進位數轉換成unnormalized
    x=reshape(Chi_En_03,m,[]); % 將候選解排列成m*pop_size
    for j=1:pop_size
        f(1,j)=cost_function(x(:,j));
    end
    [f,I]=sort(f,'ascend');
    x=x(:,I); % 將候選解依f由小至大排列
    bin=bin(:,I); % 將bin依f由小至大排列
    F(i)=min(f);
end
X=x(:,1)
F(end)