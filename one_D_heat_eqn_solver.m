%solve 1D heat equation

%ut=mu*uxx

% physical parameter
mu=1/4;
% numerical parameter
N=50;    %一條線上分割的數量
dx=1/N;  %單位長度
x=0:dx:1; %在x軸上以dx為單位得到的點
dt=dx^2;  %單位時間
alpha=mu*dt/dx^2;  %將所有固定的係數整理起來
t_final=1;%計算總時間
nt=N^2; %總時間量,nt*dt=1

%initial value
u=zeros(N,1);
u=cos(pi*(dx/2:dx:1-dx/2)'); %初始條件t=0
uext=exp(-pi^2*mu*t_final)*cos(pi*(dx/2:dx:1-dx/2)');% 參考的解

%construct matrix
A=(1-2*alpha)*diag(ones(N,1))...
+(alpha)*diag(ones(N-1,1),1)...
+(alpha)*diag(ones(N-1,1),-1);

A(1,1)=1-alpha;
A(N,N)=1-alpha;
%矩陣迭代
for k=1:nt
    u=A*u;
end
%compute error
err=u-uext;
norm(err);
