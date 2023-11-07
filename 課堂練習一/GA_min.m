clear;
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% �ѼƳ]�w
cost_function=@M1;
sol_space=nncopy([-600 600],10,1);
pop_size=30;
bits=30; % �C���ܼƪ��줸��
survival_rate=0.3; % 0.5���10��coromosomes���A���u��5�Ӫ����i�J������
mut_rate=0.01; % �l�Nchromosomes���`�줸��*mut_rate�Y�����ܦ줸��
max_gen=200;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%�D�{��
m=size(sol_space,1); % m�G�ܼƭӼ�
r=sol_space(:,2)-sol_space(:,1); % r�G�ѪŶ��d��t(m*1�x�})
keep=ceil(survival_rate*pop_size); % �q���u��keep�Ӫ����@��������

% Initialization
bin=round(rand(m*bits,pop_size)); % �H�üƲ���(m*bits)*pop_size���G�i��줸
quant=(0.5.^[1:bits])/sum(0.5.^[1:bits]); % quantization levels normalized

% Evaluation
q=size(bin,2); % �T�{pop_size
Chi_En_01=reshape(bin,bits,[]); % �N�G�i��줸�ƦC���C��ȧt�@���ܼ�
Chi_En_02=quant*Chi_En_01; % �N�G�i��줸�ഫ���Q�i��줸(����0~1����)
Chi_En_03=Chi_En_02.*nncopy(r',1,q)+nncopy(sol_space(:,1)',1,q); % �N�Q�i����ഫ��unnormalized
x=reshape(Chi_En_03,m,[]); % �N�Կ�ѱƦC��m*pop_size
for i=1:pop_size
    f(1,i)=cost_function(x(:,i));
end
[f,I]=sort(f,'ascend');
x=x(:,I); % �N�Կ�Ѩ�f�Ѥp�ܤj�ƦC
bin=bin(:,I); % �Nbin��f�Ѥp�ܤj�ƦC
F(1)=min(f);

for i=2:max_gen
    
    % Selection
    M=ceil((pop_size-keep)/2); % ��t��(�ڸs�`�Ʀ���keep�ƫ᪺�@�b)
    prob=fliplr([1:keep]/sum([1:keep])); % �̧ǲ��� keep/(1+2+...+keep)�B...�B1/(1+2+...+keep)�ƦC
    odds=[0 cumsum(prob(1:keep))]; % ����[0 keep/(1+2+...+keep) ... 1]�����v���G�V�q
    pick1=rand(1,M); % 1*M ���ü�
    pick2=rand(1,M); % 1*M ���ü�
    ic=1;
    while ic<=M
        for id=2:keep+1
            if pick1(ic)<=odds(id) & pick1(ic)>odds(id-1)
                ma(ic)=id-1; % ����(1*M ����)
            end
            if pick2(ic)<=odds(id) & pick2(ic)>odds(id-1)
                pa(ic)=id-1; % ����(1*M ����)
            end
        end
        ic=ic+1;
    end
    
    % Crossover
    ix=1:2:(M*2-1); % ix=1,3,5,..
    xp=ceil(rand(1,M)*(bits*m-1)); % 1*M ����t�I(����1~(bits*m-1)�����)
    bin(:,keep+ix)=[bin(1:xp,ma);bin(xp+1:bits*m,pa)]; % ����bsurvival�᪺�Ĥ@��l�N
    bin(:,keep+ix+1)=[bin(1:xp,pa);bin(xp+1:bits*m,ma)]; % ����bsurvival�᪺�ĤG��l�N
    
    % Mutation
    nmut=ceil(m*bits*(2*M)*mut_rate); % ���ܦ줸��
    mrow=ceil(rand(1,nmut)*(bits*m)); % ��C(1~(bits*m)����ƭ�)
    mcol=keep+ceil(rand(1,nmut)*(2*M)); % ���(keep+1~keep+2M�A�Ȱw��l�N�i�����)
    for j=1:nmut
        bin(mrow(j),mcol(j))=abs(bin(mrow(j),mcol(j))-1);
    end
    
    % Evaluation
    q=size(bin,2); % �T�{pop_size
    Chi_En_01=reshape(bin,bits,[]); % �N�G�i��줸�ƦC���C��ȧt�@���ܼ�
    Chi_En_02=quant*Chi_En_01; % �N�G�i��줸�ഫ���Q�i��줸(����0~1����)
    Chi_En_03=Chi_En_02.*nncopy(r',1,q)+nncopy(sol_space(:,1)',1,q); % �N�Q�i����ഫ��unnormalized
    x=reshape(Chi_En_03,m,[]); % �N�Կ�ѱƦC��m*pop_size
    for j=1:pop_size
        f(1,j)=cost_function(x(:,j));
    end
    [f,I]=sort(f,'ascend');
    x=x(:,I); % �N�Կ�Ѩ�f�Ѥp�ܤj�ƦC
    bin=bin(:,I); % �Nbin��f�Ѥp�ܤj�ƦC
    F(i)=min(f);
end
X=x(:,1)
F(end)