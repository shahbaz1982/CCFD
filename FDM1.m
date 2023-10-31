clc
clear all
close all
% Constants and parameters
a = 0;          % Starting point of the interval
b = 1;          % Ending point of the interval
h = 1/8;        % Step size
N = (b - a) / h; % Number of steps

% Parameters
phi=0.01;
ros=0.001;
rof=0.02;
co2=1.5;
neta=0.6;
co3=0.5;
kn=0.1;
kf=0.2;
lambda=0.5;

l = 1/(((1-phi)^2.5)*(1-phi+phi*(ros/rof)));          % Parameter l
m = (2*co2)/((1-phi+phi*(ros/rof))*(neta + co3)^4);          % Parameter m
n = (kn/kf)/(1-phi+phi*(ros/rof));          % Parameter n
p = 0.1;          % Parameter perantle
q = (2*lambda*co2)/((neta + co3)^3);          % Parameter q
r = 1;          % Parameter epsilon
s = 4*lambda;          % Parameter s

% Initialize arrays
t = a:h:b;
f = zeros(N+1,1);
g = zeros(N+1,1);
g(1,1)=1;





% Create coefficient matrices
A1 = zeros(N-1, N-1);
A2 = zeros(N-1, N-1);
A3 = zeros(N-1, N-1);

I=eye(N-1, N-1);

X1=zeros(N-1, 1);
X2=zeros(N-1, 1);


B1 = zeros(N-1, N-1);
B2 = zeros(N-1, N-1);
B3 = zeros(N-1, N-1);
B4 = zeros(N-1, N-1);

for ii = 1:5
for i = 1:N-1
    % Calculate coefficients for the equations
    co1 = l / (2*h^3);
    co2 = 1 / h^2;
    co3 = 1 / (4*h^2);
    co4 = 1 / (2 * h);
    co5 =  1 / h;
    
    
    % Fill the coefficient matrices
    A1(i,i)=-2*co1;A1(i, i+1) = -2*co1; A1(i+1,i)=2*co1;
    A1(i,i+2)=co1; A1(i+2,i)=-co1;
    
        
    A2(i,i)=-2*co2*f(i+1);A2(i,i+1) = co2*f(i+1);A2(i+1,i)=co2*f(i+2);
    
    A3(i, i+1) = co3*(f(i+2)-f(i)); A3(i+1,i)=-co3*(f(i+2)-f(i));
   % A3(1,2)=co3*f(i+1); A3(N-1,N-2)=co3*f(i-1);
   
   X1(1,1)= 2*(l/h)+((h^4)/4)+f(2,1);
   X1(2,1)=l/h;
       
    B1(i,i)=-2*co2;B1(i, i+1) = co2; B1(i+1,i)=co2;
    
    B2(i, i+1) = co4*f(i+1); B2(i+1,i)=-co4*f(i+2);
    
    B3(i,i)=f(i+1);

    B4(i, i+1) = co5*g(i+1); B4(i+1,i)=-co5*g(i+2);
    
    X2(1,1)= -(n/(h^2))*g(2,1) + (p/(2*h))*f(2,1)*g(1,1) + p*h*g(2,1) + s*(h^2)/4;     
    X2(N-1,1)= -(n/(h^2))*g(N+1,1) + (p/(2*h))*f(N-1,1)*g(N,1);
   
    
       
end

A1 = A1(1:N-1,1:N-1)
A2 = A2(1:N-1,1:N-1)
A3 = A3(1:N-1,1:N-1)


B1 = B1(1:N-1,1:N-1)
B2 = B2(1:N-1,1:N-1)
B3 = B3(1:N-1,1:N-1)
B4 = B4(1:N-1,1:N-1)

    J  = l*A1 + A2 - A3;
    K = m*I;
    L = n*B1 + p*B2 +q*B3;
    M = p*B4 + q*r*I + s*A3;
    
    F = inv(J+K*inv(L)*M)*(X1 - K*inv(L)*X2)
    G = inv(L)*(M*F + X2)
end    

t1 = t(2: N)
DF = diff(F);
DFF = [1; DF]

% Plotting the solution
figure;
plot(t1(:), DFF, '-o', t1(:), G, '-x');
xlabel('t');
ylabel('df(t), g(t)');
legend('df(t)', 'g(t)');
grid on;
