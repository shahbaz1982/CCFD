clear all 
close all
 clc
 format short
 a=0
 b=1
 N=100
 dx=(b-a)/N
 X=0:dx:1;
 XR=[X(2:N)];
 
 %RHS
  LW=1; %Left
  RW=0; %Right
  R=zeros(N-1,1);
  R(1,1)=R(1,1)-LW*((1/(dx^2))-(1/(2*XR(1)*dx)));
  R(N-1,1)=R(N-1,1)-RW*((1/(dx^2))+(1/(2*XR(N-1)*dx)));

 %LHS
 O=ones(N-1,1);
 DW=spdiags([-1*O 0*O 1*O], -1:1, N-1,N-1);
  
 DDW=spdiags([1*O -2*O 1*O], -1:1, N-1,N-1);
  
  
  for i=1:N-1
     for j=1:N-1
  TW(i,j)=XR(i)*DW(i,j);
     end
  end
  
  
  
  W=inv((1/(dx^2))*DDW+(1/(2*dx))*TW)*R;
  
  %RHS Theta
  Br=1;
  LO=1; %Left
  ROO=0; %Right
  
  RO = -Br*((1/(2*dx))*DW*W).^2;
  RO(1,1)=RO(1,1)-LO*((1/(dx^2))-(1/(2*XR(1)*dx)));
  RO(N-1,1)=RO(N-1,1)-ROO*((1/(dx^2))+(1/(2*XR(N-1)*dx)));

  
  THETA=inv((1/(dx^2))*DDW+(1/(2*dx))*TW)*RO;
  
  
  
  %Graph
  plot(X,[LW; W; RW]','-*',X,[LO; THETA; ROO]','-*')
  legend('W','Theta')
  
  
  
  
  
  
 
 