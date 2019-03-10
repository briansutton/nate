function kappa = ivpl2conduni(al,be,ga,ab,c,m,n)
%IVPL2CONDUNI   Estimate the absolute condition number of a second-order
%collocation problem on a piecewise-uniform grid.
%
%   KAPPA = IVPL2CONDUNI(AL,BE,GA,[A B],[ C10 C11 C12; C20 C21 C22 ],M,N)
%   estimates the absolute condition number of the collocation problem
%   associated with the IVP
%
%     AL(T)*U(T) + BE(T)*U'(T) + GA(T)*U''(T) = G(T),
%     C10*U(A) + C11*U'(A) + C12*U''(A) = D1,
%     C20*U(A) + C21*U'(A) + C22*U''(A) = D2,
%
%   on a piecewise-uniform grid.
%
%   Inputs:
%
%     AL: coefficient function of U(T) in the differential equation
%
%     BE: coefficient function of U'(T) in the differential equation
%
%     GA: coefficient function of U''(T) in the differential equation
%
%     AB = [A B]: endpoints of the problem domain
%
%     C = [ C10 C11 C12; C20 C21 C22 ]: matrix of coefficients in the
%     initial conditions
%
%     M: degree of the grid of collocation nodes on each subinterval
%
%     N: number of subintervals
%
%   Outputs:
%
%     KAPPA: an estimate for the absolute condition number of the finite
%     collocation system
%
%   References:
%
%     William W. Hager. Condition estimates. SIAM J. Sci. Statist. Comput.
%     5(2):311-316, 1984.
%
%   Example: A damped and potentially forced oscillator
%   2*U + 0.1*U' + U'' = G(T).
%
%     al = @(t) 2;
%     be = @(t) 0.1;
%     ga = @(t) 1;
%     a = 0; b = 90;
%     c = [ 1 0 0; 0 1 0 ];
%     m = 3; n = 1000;
%     ivpl2conduni(al,be,ga,[a b],c,m,n)
%
%   Copyright 2019 Brian Sutton

narginchk(7,7);
natecheck('ivpl2conduni',al,be,ga,ab,c,m,n);
s = warning('query','MATLAB:nearlySingularMatrix');
warning('off','MATLAB:nearlySingularMatrix');
a = ab(1);
b = ab(2);
als = sampleuni(al,[a b],m,n);
bes = sampleuni(be,[a b],m,n);
gas = sampleuni(ga,[a b],m,n);
l = (b-a)/n;
J1 = [ -ones(m+2,1) eye(m+2) ];
J2 = [ -ones(m+1,1) eye(m+1) ];
K1 = intmatrixuni(m+1);
K2 = intmatrixuni(m);
E1 = resamplematrixuni(m,m+2);
E2 = resamplematrixuni(m,m+1);
v = zeros((3*m+6)*n,1);
kappaabs = nan(2,1);
x = 1/((3*m+6)*n)*ones((3*m+6)*n,1);
for trial = 1:3
  for k = 1:5
    y = Ltinvtimes(als,bes,gas,c,m,n,l,J1,J2,K1,K2,E1,E2,x);
    %ycheck = L'\x; %[[temp]]
    %max(abs(y-ycheck)) %[[temp]]
    xi = nzsign(y);
    z = Linvtimes(als,bes,gas,c,m,n,l,J1,J2,K1,K2,E1,E2,xi);
    %zcheck = L\xi; %[[temp]]
    %max(abs(z-zcheck)) %[[temp]]
    [M,j] = max(abs(z));
    if M<=z'*x
      break;
    end
    e = [ zeros(j-1,1); 1; zeros((3*m+6)*n-j,1) ];
    v(j) = 1;
    x = e;
  end
  kappaabs(trial) = norm(Ltinvtimes(als,bes,gas,c,m,n,l,J1,J2,K1,K2,E1,E2,x),1);
  x = 1/((3*m+6)*n-sum(v))*(ones((3*m+6)*n,1)-v);
end
kappa = max(kappaabs);
warning(s);

end

function y = nzsign(x)

y = sign(x);
y(y==0) = 1;

end

function y = Linvtimes(als,bes,gas,c,m,n,l,J1,J2,K1,K2,E1,E2,x)

u = nan(m+3,n);
v = nan(m+2,n);
w = nan(m+1,n);
e = reshape([x(end-1); x(1:(m+3)*n-1)],m+3,n);
f = reshape([x(end); x((m+3)*n:(m+3)*n+(m+2)*n-2)],m+2,n);
g = reshape(x((m+3)*n+(m+2)*n-1:end-2),m+1,n);
lambdas1 = K1(end,:)';
lambdas2 = K2(end,:)';
L1 = [ c(1,1)*eye(1,m+3) c(1,2)*eye(1,m+2) c(1,3)*eye(1,m+1) ;
       J1                -l*K1(2:end,:)    zeros(m+2,m+1)    ;
       c(2,1)*eye(1,m+3) c(2,2)*eye(1,m+2) c(2,3)*eye(1,m+1) ;
       zeros(m+1,m+3)    J2                -l*K2(2:end,:)    ;
       diag(als(:,1))*E1 diag(bes(:,1))*E2 diag(gas(:,1))    ];
rhs1 = [e(:,1);f(:,1);g(:,1)];
sol1 = L1\rhs1;
u(:,1) = sol1(1:m+3);
v(:,1) = sol1(m+4:2*m+5);
w(:,1) = sol1(2*m+6:end);
cumsumv = zeros(m+2,1);
cumsumw = zeros(m+1,1);
for j = 2:n
  Lj = [ eye(m+3)          -l*K1             zeros(m+3,m+1) ;
         zeros(m+2,m+3)    eye(m+2)          -l*K2 ;
         diag(als(:,j))*E1 diag(bes(:,j))*E2 diag(gas(:,j)) ];
  cumsumv = cumsumv+v(:,j-1);
  cumsumw = cumsumw+w(:,j-1);
  rhsj = [ e(:,j)+u(1)*ones(m+3,1)+ones(m+3,1)*l*lambdas1'*cumsumv ;
           f(:,j)+v(1)*ones(m+2,1)+ones(m+2,1)*l*lambdas2'*cumsumw ;
           g(:,j) ];
  solj = Lj\rhsj;
  u(:,j) = solj(1:m+3);
  v(:,j) = solj(m+4:2*m+5);
  w(:,j) = solj(2*m+6:end);
end
y = [ reshape(u,(m+3)*n,1); reshape(v,(m+2)*n,1); reshape(w,(m+1)*n,1) ];

end

function y = Ltinvtimes(als,bes,gas,c,m,n,l,J1,J2,K1,K2,E1,E2,x)

e = nan(m+3,n);
f = nan(m+2,n);
g = nan(m+1,n);
u = reshape(x(1:(m+3)*n),m+3,n);
v = reshape(x((m+3)*n+1:(2*m+5)*n),m+2,n);
w = reshape(x((2*m+5)*n+1:end),m+1,n);
lambdas1 = K1(end,:)';
lambdas2 = K2(end,:)';
cumsume = zeros(m+3,1);
cumsumf = zeros(m+2,1);
for j = n:-1:2
  Lj = [ eye(m+3)       zeros(m+3,m+2) (diag(als(:,j))*E1)' ;
         -l*K1'         eye(m+2)       (diag(bes(:,j))*E2)' ;
         zeros(m+1,m+3) -l*K2'         diag(gas(:,j))       ];
  rhsj = [ u(:,j)                                ;
           v(:,j)+l*lambdas1*ones(1,m+3)*cumsume ;
           w(:,j)+l*lambdas2*ones(1,m+2)*cumsumf ];
  solj = Lj\rhsj;
  e(:,j) = solj(1:m+3);
  f(:,j) = solj(m+4:2*m+5);
  g(:,j) = solj(2*m+6:end);
  cumsume = cumsume+e(:,j);
  cumsumf = cumsumf+f(:,j);
end
L1 = [ J1'             zeros(m+3,m+1)  (diag(als(:,1))*E1)' c(1,1)*eye(m+3,1) c(2,1)*eye(m+3,1) ;
       -l*K1(2:end,:)' J2'             (diag(bes(:,1))*E2)' c(1,2)*eye(m+2,1) c(2,2)*eye(m+2,1) ;
       zeros(m+1,m+2)  -l*K2(2:end,:)' diag(gas(:,1))       c(1,3)*eye(m+1,1) c(2,3)*eye(m+1,1) ];
rhs1 = [ u(:,1)+sum(sum(e(:,2:n)))*eye(m+3,1) ;
         v(:,1)+l*lambdas1*ones(1,m+3)*cumsume+sum(sum(f(:,2:n)))*eye(m+2,1) ;
         w(:,1)+l*lambdas2*ones(1,m+2)*cumsumf ];
sol1 = L1\rhs1;
e(:,1) = sol1([end-1 1:m+2]);
f(:,1) = sol1([end m+3:2*m+3]);
g(:,1) = sol1(2*m+4:end-2);
e = reshape(e,(m+3)*n,1);
f = reshape(f,(m+2)*n,1);
g = reshape(g,(m+1)*n,1);
y = [ e(2:end) ; f(2:end) ; g ; e(1) ; f(1) ];

end

