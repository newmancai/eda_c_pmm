clear
clc
fs=50;         %采样频率为50Hz
dt=1/fs;       %采样时间间隔为0.02s
n = 400;       %数据共400个点
t = 0:dt:(n-1)*dt;
%给定的分析信号y，信号有3个模式，信号阶数为5 (2+1+2)
y = 2.*exp(-0.1*t).*cos(2*pi*3.*t+pi/3)+1.*exp(-0.05*t)+1.*cos(2*pi*1.*t-pi);

N=fix(n/2);%扩展阶数，一般取数据点数/2

%对应于式（19）的矩阵
x0=zeros(N+1,N+1);
for ri=1:N+1
    for rj=1:N+1
        for nn=N+1:n
            x0(ri,rj)=x0(ri,rj)+conj(y(nn-ri+1))*y(nn-rj+1);
        end
    end
end    

%% 采用SVD-TLS方法计算式（19）                  
[~,S,V]=svd(x0);          %奇异值分解

%归一化奇异值法确定阶数
for i=1:1:N                
    if(S(i,i)/S(1,1)<0.001)  %0.001为判定阶数的阈值
         break;
    end
end
ord=i-1;

SP=zeros(ord+1,ord+1);
for sj=1:ord
    for si=1:(N+1-ord)
        SP=SP+(S(sj,sj)^2)*V(si:si+ord,sj)*V(si:si+ord,sj)';
    end
end
SP1=inv(SP);
xx=SP1(2:ord+1,1)/SP1(1,1);   %多项式的系数z
zz=zeros(ord+1,1);
zz(1)=1;
zz(2:ord+1)=xx;
zn=roots(zz);%得到特征方程的根

Va=zeros(n,ord);    %生成范得蒙矩阵Va，对应式（20）中的Z
for k=1:n
   Va(k,:)=zn.'.^(k-1);        
end 

bn=(Va'*Va)^-1*Va'*y';

ff=abs(atan2(imag(zn),real(zn))/(2*pi*dt));%初步计算频率
DD=log(abs(zn))/dt;%初步计算衰减因子

%将频率从小到大排序
[F2,I]=sort(ff);
%剔除重复特征值对应的结果
m=0;
for k=1:ord-1
    if F2(k) ~= F2(k+1) && F2(k)~=0
        continue;
    end 
    m=m+1;
    l=I(k);
    F(m)=ff(l);%频率
    D(m)=DD(l);%阻尼比
    if F(m)==0
        A(m)=abs(bn(l));
    else
        A(m)=2*abs(bn(l));
    end
    theta(m)=angle(bn(l));%初相
    %theta(m)=atan2(imag(bn(l)),real(bn(l)));%初相
end
Dzn=-D./sqrt(D.*D+4*pi*pi.*F.*F)*100;%阻尼比
% Dzn=-log(abs(zn))/(abs(log(zn))*dt);%阻尼比

%% 画图
nt=0:1:n-1;
yy=zeros(size(y));     %根据辨识结果拟合曲线
 for ii=1:length(A)
     yy=yy+A(ii)*exp(D(ii)*dt*(nt)).*cos(2*pi*F(ii)*dt*(nt)+theta(ii)); 
 end
t = 0:dt:(n-1)*dt;
figure(1)
plot(t,y,t,yy,'--');
xlabel('时间(s)');
ylabel('幅值');
legend('真实','拟合');
