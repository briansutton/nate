function kappa = ivpl1condchebpw(al,be,asbs,c,m)
%IVPL1CONDCHEBPW   Absolute condition number of a first-order linear
%collocation problem on a piecewise-Chebyshev grid.
%
%   KAPPA = IVPL1CONDCHEBPW(AL,BE,ASBS,[C0 C1],M) estimates the absolute
%   condition number of the collocation problem associated with
%
%     AL(T)*U(T) + BE(T)*U'(T) = G(T), C0*U(A) + C1*U'(A) = D
%
%   on a piecewise-Chebyshev grid.
%
%   Inputs:
%
%     AL: coefficient function of U(T) in the differential equation
%
%     BE: coefficient function of U'(T) in the differential equation
%
%     ASBS = [ A1 A2 ... AN; B1 B2 ... BN ]: endpoints of the subintervals
%     [AJ,BJ] in the piecewise-defined grid
%
%     C = [C0 C1]: vector of coefficients in the initial condition
%
%     M: degree of the grid of collocation nodes
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
%   Copyright 2019 Brian Sutton

narginchk(5,5);
natecheck('ivpl1condchebpw',al,be,asbs,c,m)
s = warning('query','MATLAB:nearlySingularMatrix');
warning('off','MATLAB:nearlySingularMatrix');
n = size(asbs,2);
ls = asbs(2,:)-asbs(1,:);
als = samplechebpw(al,asbs,m);
bes = samplechebpw(be,asbs,m);
J = [ -ones(m+1,1) eye(m+1) ];
K = intmatrixcheb(m);
E = resamplematrixcheb(m,m+1);
v = zeros((2*m+3)*n,1);
kappaabs = nan(2,1);
x = 1/((2*m+3)*n)*ones((2*m+3)*n,1);
for trial = 1:3
  for k = 1:5
    y = Ltinvtimes(als,bes,c,m,n,ls,J,K,E,x);
    xi = nzsign(y);
    z = Linvtimes(als,bes,c,m,n,ls,J,K,E,xi);
    [M,j] = max(abs(z));
    if M<=z'*x
      break;
    end
    e = [ zeros(j-1,1); 1; zeros((2*m+3)*n-j,1) ];
    v(j) = 1;
    x = e;
  end
  kappaabs(trial) = norm(Ltinvtimes(als,bes,c,m,n,ls,J,K,E,x),1);
  x = 1/((2*m+3)*n-sum(v))*(ones((2*m+3)*n,1)-v);
end
kappa = max(kappaabs);
warning(s);

end

function y = nzsign(x)

y = sign(x);
y(y==0) = 1;

end

function y = Linvtimes(als,bes,c,m,n,ls,J,K,E,x)

u = nan(m+2,n);
v = nan(m+1,n);
e = reshape([x(end); x(1:(m+2)*n-1)],m+2,n);
f = reshape(x((m+2)*n:end-1),m+1,n);
lambdas = K(end,:)';
L1 = [ c(1)*eye(1,m+2)  c(2)*eye(1,m+1)     ;
       J                -ls(1)/2*K(2:end,:) ;
       diag(als(:,1))*E diag(bes(:,1))      ];
rhs1 = [ e(:,1); f(:,1) ];
sol1 = L1\rhs1;
u(:,1) = sol1(1:m+2);
v(:,1) = sol1(m+3:end);
cumsumv = zeros(m+1,1);
for j = 2:n
  Lj = [ eye(m+2)         -ls(j)/2*K     ;
         diag(als(:,j))*E diag(bes(:,j)) ];
  cumsumv = cumsumv+ls(j-1)/2*v(:,j-1);
  rhsj = [ e(:,j)+u(1)*ones(m+2,1)+ones(m+2,1)*lambdas'*cumsumv ;
           f(:,j)                                               ];
  solj = Lj\rhsj;
  u(:,j) = solj(1:m+2);
  v(:,j) = solj(m+3:end);
end
y = [ reshape(u,(m+2)*n,1); reshape(v,(m+1)*n,1) ];

end

function y = Ltinvtimes(als,bes,c,m,n,ls,J,K,E,x)

e = nan(m+2,n);
f = nan(m+1,n);
u = reshape(x(1:(m+2)*n),m+2,n);
v = reshape(x((m+2)*n+1:end),m+1,n);
lambdas = K(end,:)';
cumsume = zeros(m+2,1);
for j = n:-1:2
  Lj = [ eye(m+2)    (diag(als(:,j))*E)' ;
         -ls(j)/2*K' diag(bes(:,j))      ];
  rhsj = [ u(:,j)                                     ;
           v(:,j)+ls(j)/2*lambdas*ones(1,m+2)*cumsume ];
  solj = Lj\rhsj;
  e(:,j) = solj(1:m+2);
  f(:,j) = solj(m+3:end);
  cumsume = cumsume+e(:,j);
end
L1 = [ J'                   (diag(als(:,1))*E)' c(1)*eye(m+2,1) ;
       -ls(1)/2*K(2:end,:)' diag(bes(:,1))      c(2)*eye(m+1,1) ];
rhs1 = [ u(:,1)+sum(sum(e(:,2:n)))*eye(m+2,1)       ;
         v(:,1)+ls(1)/2*lambdas*ones(1,m+2)*cumsume ];
sol1 = L1\rhs1;
e(:,1) = sol1([end 1:m+1]);
f(:,1) = sol1(m+2:end-1);
e = reshape(e,(m+2)*n,1);
f = reshape(f,(m+1)*n,1);
y = [ e(2:end) ; f ; e(1) ];

end

