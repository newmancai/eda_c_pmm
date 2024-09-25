%%以对S11曲线进行向量拟合为例
close all
clc;
clear;
addpath('..');
pmm_setup

opts=pmm_default;


%% LOAD DATA
% [G,W,F,H]=pmm('channel.S2P',5, opts);
% inputfile='symind.s2p';
inputfile='../data/IEEE39BUS26.s1p'
% [inputfile,outputfile,subcktname]=ParseInput(inputfile);

sobj = sparameters(inputfile);
Fr = sobj.Frequencies;
Sx = sobj.Parameters;
[Frequency,Sx]=ldstone(inputfile,opts);

Frequency=Fr';  %%1.5-6GHz(需要进行向量拟合的频段)
Omega=2*pi*Frequency;  %%rad frequency(Omega)
%%获取需要拟合的数据
% S11_Real=interp1(S11_Real_Data(:,1),S11_Real_Data(:,2),Frequency,'spline');  %%用于计算复数向量S11
% S11_Imag=interp1(S11_Imag_Data(:,1),S11_Imag_Data(:,2),Frequency,'spline');  %%用于计算复数向量S11
% %%根据实部和虚部计算S11幅值(linear)
% S11_abs_linear=abs(S11_Real+S11_Imag*i);
% %%S11幅值(dB)
% S11_abs_dB=20*log10(S11_abs_linear);
% %%S11复数(根据S11实部和虚部计算)
% S11=S11_Real+S11_Imag*i;  %%以复数形式表示S11(包括S11的幅值和相位，一般所说的为其幅值)
% f=reshape(Sx(1,1,:),[size(Sx,3), 1]);
for wwwi=1:size(Sx,3)
   f(wwwi)=Sx(:,:,wwwi); 
end
 
%%初始化极点(Initialize the poles)(传递函数的阶数为O，为传递函数极点个数), ...
%%因为设置的极点为共轭复数极点，其总是成对出现，设置为偶数
O=10;  %%传递函数的阶数
% Frequency(1)
% Frequency(size(Frequency,1))
% Frequency(size(Frequency,1))
Belta=linspace(1,Frequency(1,size(Frequency,2)),O/2+2)'  %%1.5-6GHz包括起点和终点的采样频率点(Imaginary Parts of Poles)
Belta(1)=[];  %%舍弃起始频率点
Belta(end)=[];  %%舍弃终末频率点
%%这里的传递函数阶数为O(对应O/2个谐振峰值点)
Alpha=Belta/100;  %%Real Parts of Poles(Alpha的长度为O/2)
An=zeros(2*length(Belta),1);  %%Complex Poles(起始极点)(An的长度为O)
for o=1:length(Belta)
    An(o*2-1)=-Alpha(o)+i*Belta(o);
    An(o*2)=-Alpha(o)-i*Belta(o);
end

%%计算矩阵A和向量B
A=zeros(length(Frequency),2*length(An)+1);
B=zeros(length(Frequency),1);
for k=1:length(Frequency)
    %%Calculation of Matrix A(矩阵A的列数为2*O+2)
    for o=1:length(An)
        A(k,o)=1/(i*Omega(k)-An(o));
%         Omega(k)
%         An(o)
%         break
    end
%     break
    A(k,length(An)+1)=1;
%     A(k,length(An)+2)=i*Omega(k);
    for o=1:length(An)
        A(k,length(An)+1+o)=-f(k)/(i*Omega(k)-An(o));
%         -f(k)/(i*Omega(k)-An(o))
%         break
    end
%     break
    %%Calculation of Matrix B
    B(k)=f(k);
end
X=pinv(A)*B;   %%pinv(A)为计算矩阵A的伪逆矩阵，求解得到X

%%Unknown Variables
C=X(1:length(An));
d=X(length(An)+1);
% h=X(length(An)+2);
C1=X(length(An)+2:end);
 
syms x;
sigma=@(x)sum(C1./(x-An))+1;
An1=solve(sigma(x)==0,x);
An1=double(An1);  %%sigma函数的零点(新极点)作为下一次迭代的起始极点
% An1 
%%sigma函数的零点(新极点)与起始极点的差值
Error=sum(abs(An1-An));
%%迭代求解该S11曲线的零点和极点
tau=5e-4;  %%Threshold to Stop
 
while Error>tau
    An=An1;  %%新极点作为起始极点[O,1]
    %%计算矩阵A和向量B
    A=zeros(length(Frequency),2*length(An)+1);
    for k=1:length(Frequency)
        %%Calculation of Matrix A
        for o=1:length(An)
            A(k,o)=1/(i*Omega(k)-An(o));
        end
        A(k,length(An)+1)=1;
%         A(k,length(An)+2)=i*Omega(k);
        for o=1:length(An)
            A(k,length(An)+1+o)=-f(k)/(i*Omega(k)-An(o));
        end
    end
    X=pinv(A)*B;
 
 %%Unknown Variables
    C=X(1:length(An));
    d=X(length(An)+1);
%     h=X(length(An)+2);
    C1=X(length(An)+2:end);
    
    sigma=@(x)sum(C1./(x-An))+1;
    An1=solve(sigma(x)==0,x);  
    An1=double(An1);  %%sigma函数的零点(新极点)作为下一次迭代的起始极点
    %%sigma函数的零点(新极点)与起始极点的差值
    Error=sum(abs(An1-An))
end
 
An=An1  %%新极点作为起始极点[O,1](传递函数的极点) 
%%计算矩阵A和向量B
A=zeros(length(Frequency),2*length(An)+1);
for k=1:length(Frequency)
    %%Calculation of Matrix A
    for o=1:length(An)
        A(k,o)=1/(i*Omega(k)-An(o));
    end
    A(k,length(An)+1)=1;
%     A(k,length(An)+2)=i*Omega(k);
    for o=1:length(An)
       A(k,length(An)+1+o)=-f(k)/(i*Omega(k)-An(o));
    end
end
X=pinv(A)*B;
 
%%Unknown Variables
C=X(1:length(An))
d=X(length(An)+1)
% h=X(length(An)+2);
C1=X(length(An)+2:end);
 
RS=zeros(length(Frequency),1);
for k=1:length(Frequency)
     RS(k)=sum(C./(i*Omega(k)-An))+d;
end


RS_Linear=abs(RS);


 
%%绘制(sigma*f)函数和sigma函数的幅值曲线
sigma_f=zeros(length(Frequency),1);
sigma=zeros(length(Frequency),1);
sigmaf=zeros(length(Frequency),1);
for k=1:length(Frequency)
    sigma_f(k)=sum(C./(i*Omega(k)-An))+d;
    sigma(k)=sum(C1./(i*Omega(k)-An))+1;
    sigmaf(k)=sigma(k)*f(k);
end


error1=0;
error2=0;
for wwwk=1:size(RS,1)
    error1=error1+norm(RS(wwwk)-f(wwwk),2)
    error2=error2+norm(f(wwwk),2)
end
disp(['error= ', num2str(error1/error2)]);%<0.1
disp(['error(DC) = ',num2str(norm(f(1)-RS(1),2))]);%<10^-10
singularValues = svd(d);
if all(singularValues < 1)
    disp(['passivity = ', 'passivity']);
else
    disp(['passivity = ', 'non-passivity']);
end



Error_f=sigmaf-sigma_f;
 
figure(1);
plot(Frequency,RS_Linear,'b-',Frequency,abs(f),'r--','linewidth',2);
xlabel('Frequency(GHz)','FontSize',14);
ylabel('Magnitude','FontSize',14);
legend('Vector Fitting','Simulation(CST)','FontSize',14);
 
figure(2);
plot(Frequency,abs(sigma_f),'b-',Frequency,abs(sigma),'r--', ...
    Frequency,abs(sigmaf),'k:','linewidth',2);
xlabel('Frequency(GHz)','FontSize',14);
ylabel('S11(Linear)','FontSize',14);
legend('\sigma_f_i_tf','\sigma_f_i_t','(\sigmaf)_f_i_t','FontSize',14);
 
figure(3);
plot(Frequency,abs(Error_f),'k:','linewidth',2);
xlabel('Frequency(GHz)','FontSize',14);