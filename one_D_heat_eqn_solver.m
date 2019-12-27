%solve 1D heat equation

%ut=mu*uxx

% physical parameter
mu=1/4;
% numerical parameter
N=50;    %�@���u�W���Ϊ��ƶq
dx=1/N;  %������
x=0:dx:1; %�bx�b�W�Hdx�����o�쪺�I
dt=dx^2;  %���ɶ�
alpha=mu*dt/dx^2;  %�N�Ҧ��T�w���Y�ƾ�z�_��
t_final=1;%�p���`�ɶ�
nt=N^2; %�`�ɶ��q,nt*dt=1

%initial value
u=zeros(N,1);
u=cos(pi*(dx/2:dx:1-dx/2)'); %��l����t=0
uext=exp(-pi^2*mu*t_final)*cos(pi*(dx/2:dx:1-dx/2)');% �ѦҪ���

%construct matrix
A=(1-2*alpha)*diag(ones(N,1))...
+(alpha)*diag(ones(N-1,1),1)...
+(alpha)*diag(ones(N-1,1),-1);

A(1,1)=1-alpha;
A(N,N)=1-alpha;
%�x�}���N
for k=1:nt
    u=A*u;
end
%compute error
err=u-uext;
norm(err);
