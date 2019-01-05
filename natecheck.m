function natecheck(fn,varargin)
%NATECHECK   Check inputs.
%
%   NATECHECK(FN,...) checks the inputs supplied to the function FN in the
%   NATE library.
%
%   Inputs:
%
%     FN: name of the function
%
%     ...: inputs to the function
%
%   An error is raised if any input is found to be invalid.
%
%   Copyright 2019 Brian Sutton

global natemsgid

natemsgid = 'NATE:natecheck';

try
  
  fn = lower(fn);
  switch fn
    case 'antiderivcheb'
      [pa,qs,ab] = deal(varargin{:});
      validateattributes(pa,{'numeric'},{'scalar'},fn,'PA',1);
      validateattributes(qs,{'numeric'},{'column'},fn,'QS',2);
      validateattributes(ab,{'numeric'},{'numel',2},fn,'AB',3);
    case 'antiderivgen'
      [pa,qs,xs,xs_,a] = deal(varargin{:});
      validateattributes(pa,{'numeric'},{'scalar'},fn,'PA',1);
      validateattributes(qs,{'numeric'},{'column'},fn,'QS',2);
      validateattributes(xs,{'numeric'},{'column'},fn,'XS',3);
      if length(xs)~=length(qs)+1
        error(natemsgid,'Input 3 (XS) must have one more entry than input 2 (QS).');
      end
      validateattributes(xs_,{'numeric'},{'column'},fn,'XS_',4);
      if length(xs_)~=length(qs)
        error(natemsgid,'Input 4, XS_, must have the same number of entries as input 2, QS.');
      end
      validateattributes(a,{'numeric'},{'scalar'},fn,'A',5);
    case 'antiderivuni'
      [pa,qs,ab] = deal(varargin{:});
      validateattributes(pa,{'numeric'},{'scalar'},fn,'PA',1);
      validateattributes(qs,{'numeric'},{'2d','nonempty'},fn,'QS',2);
      validateattributes(ab,{'numeric'},{'vector','numel',2},fn,'AB',3);
    case 'backsubstitute'
      [U,b] = deal(varargin{:});
      validateattributes(U,{'numeric'},{'square'},fn,'U',1);
      if norm(tril(U,-1),1)~=0
        error(natemsgid,'Input 1, U, must be an upper-triangular matrix.')
      end
      validateattributes(b,{'numeric'},{'column'},fn,'B',2);
    case 'baryweights'
      xs = varargin{1};
      validateattributes(xs,{'numeric'},{'column'},fn,'XS',1);
    case 'bernstein'
      [ab,rho] = deal(varargin{:});
      validateattributes(ab,{'numeric'},{'numel',2},fn,'AB',1);
      validateattributes(rho,{'numeric'},{'scalar','>=',1},fn,'RHO',2);
    case 'definitecheb'
      [qs,ab] = deal(varargin{:});
      validateattributes(qs,{'numeric'},{'column'},fn,'QS',1);
      validateattributes(ab,{'numeric'},{'vector','numel',2},fn,'AB',2);
    case 'definitegen'
      [qs,ab,xs,xs_] = deal(varargin{:});
      validateattributes(qs,{'numeric'},{'column'},fn,'QS',1);
      validateattributes(ab,{'numeric'},{'vector','numel',2},fn,'AB',2);
      validateattributes(xs,{'numeric'},{'column'},fn,'XS',3);
      if length(xs)~=length(qs)+1
        error('Input 3, XS, must have one more entry than input 1, QS.');
      end
      validateattributes(xs_,{'numeric'},{'column'},fn,'XS_',4);
      if length(xs_)~=length(qs)
        error(natemsgid,'Input 4, XS_, must have the same number of entries as input 1, QS.');
      end
    case 'definiteuni'
      [qs,ab] = deal(varargin{:});
      validateattributes(qs,{'numeric'},{'2d','nonempty'},fn,'QS',1);
      validateattributes(ab,{'numeric'},{'vector','numel',2},fn,'AB',2);
    case 'diffmatrix_'
      [xs_,xs,ws] = deal(varargin{:});
      validateattributes(xs_,{'numeric'},{'column'},fn,'XS_',1);
      validateattributes(xs,{'numeric'},{'column'},fn,'XS',2);
      validateattributes(ws,{'numeric'},{'column'},fn,'WS',3);
    case 'diffmatrixgen'
      [xs_,xs] = deal(varargin{:});
      validateattributes(xs_,{'numeric'},{'column'},fn,'XS_',1);
      validateattributes(xs,{'numeric'},{'column'},fn,'XS',2);
    case 'diffmatrixsquare_'
      [xs,ws] = deal(varargin{:});
      validateattributes(xs,{'numeric'},{'column'},fn,'XS',1);
      validateattributes(ws,{'numeric'},{'column'},fn,'WS',2);
      if length(xs)~=length(ws)
        error(natemsgid,'Input 1, XS, and input 2, WS, must have the same number of entries.');
      end
    case 'disptable'
      args = varargin;
      if isempty(args)
        error(natemsgid,'At least one column must be provided.');
      end
      if mod(length(args),3)~=0
        error(natemsgid,'Three inputs must be supplied for every column of the table: a heading, a vector of values, and a format string.');
      end
      n = length(args)/3;
      for j = 1:n
        validateattributes(args{3*(j-1)+1},{'char'},{},fn,sprintf('heading for column %d',j),3*(j-1)+1);
        validateattributes(args{3*(j-1)+2},{'numeric'},{'vector'},fn,sprintf('values for column %d',j),3*(j-1)+2);
        validateattributes(args{3*(j-1)+3},{'char'},{},fn,sprintf('format string for column %d',j),3*(j-1)+3);
      end
    case 'evalmatrix_'
      [xs_,xs,ws] = deal(varargin{:});
      validateattributes(xs_,{'numeric'},{'column'},fn,'XS_',1);
      validateattributes(xs,{'numeric'},{'column'},fn,'XS',2);
      validateattributes(ws,{'numeric'},{'column'},fn,'WS',3);
      if length(xs)~=length(ws)
        error(natemsgid,'Input 2, XS, and input 3, WS, must have the same number of entries.');
      end
    case 'evalmatrixgen'
      [xs_,xs] = deal(varargin{:});
      validateattributes(xs_,{'numeric'},{'column'},fn,'XS_',1);
      validateattributes(xs,{'numeric'},{'column'},fn,'XS',2);
    case 'genopivot'
      [A,b] = deal(varargin{:});
      validateattributes(A,{'numeric'},{'square'},fn,'A',1);
      validateattributes(b,{'numeric'},{'column'},fn,'B',2);
      if size(b,1)~=size(A,1)
        error(natemsgid,'Input 1, A, and input 2, B, must have the same number of rows.');
      end
    case 'gepp'
      [A,b] = deal(varargin{:});
      validateattributes(A,{'numeric'},{'square'},fn,'A',1);
      validateattributes(b,{'numeric'},{'column'},fn,'B',2);
      if size(b,1)~=size(A,1)
        error(natemsgid,'Input 1, A, and input 2, B, must have the same number of rows.');
      end
    case 'gridcheb'
      [ab,m] = deal(varargin{:});
      validateattributes(ab,{'numeric'},{'vector','numel',2},fn,'AB',1);
      validateattributes(m,{'numeric'},{'scalar','positive','integer'},fn,'M',2);
    case 'gridchebpw'
      [asbs,m] = deal(varargin{:});
      validateattributes(asbs,{'numeric'},{'nrows',2},fn,'ASBS',1);
      validateattributes(m,{'numeric'},{'scalar','positive','integer'},fn,'M',2);
    case 'griduni'
      [ab,m,n] = deal(varargin{:});
      validateattributes(ab,{'numeric'},{'numel',2},fn,'AB',1);
      validateattributes(m,{'numeric'},{'scalar','nonnegative','integer'},fn,'M',2);
      validateattributes(n,{'numeric'},{'scalar','nonnegative','integer'},fn,'N',3);
    case 'ieee754'
      x = varargin{1};
      validateattributes(x,{'numeric'},{'scalar'},fn,'X',1);
    case 'indefinitecheb'
      m = varargin{1};
      validateattributes(m,{'numeric'},{'scalar','positive','integer'},fn,'M',1);
    case 'indefinitegen'
      [xs,xs_] = deal(varargin{:});
      validateattributes(xs,{'numeric'},{'column'},fn,'XS',1);
      validateattributes(xs_,{'numeric'},{'column'},fn,'XS_',2);
      if length(xs)~=length(xs_)+1
        error(natemsgid,'Input 1, XS, must have one more entry than input 2, XS_.');
      end
    case 'indefiniteuni'
      m = varargin{1};
      validateattributes(m,{'numeric'},{'scalar','nonnegative','integer'},fn,'M',1);
    case 'infnorm'
      [f,ab,n] = deal(varargin{:});
      validateattributes(f,{'function_handle'},{},fn,'F',1);
      validateattributes(ab,{'numeric'},{'numel',2},fn,'AB',2);
      validateattributes(n,{'numeric'},{'scalar','positive','integer'},fn,'N',3);
    case 'interpcheb'
      [ps,ab] = deal(varargin{:});
      validateattributes(ps,{'numeric'},{'column'},fn,'PS',1);
      validateattributes(ab,{'numeric'},{'numel',2},fn,'AB',2);
    case 'interpchebpw'
      [ps,asbs] = deal(varargin{:});
      validateattributes(ps,{'numeric'},{'2d'},fn,'PS',1);
      validateattributes(asbs,{'numeric'},{'nrows',2},fn,'ASBS',2);
      if size(ps,2)~=size(asbs,2)
        error(natemsgid,'Input 1, PS, and input 2, ASBS, must have the same number of columns.');
      end
    case 'interpgen'
      [xs,ps] = deal(varargin{:});
      validateattributes(xs,{'numeric'},{'column'},fn,'XS',1);
      validateattributes(ps,{'numeric'},{'column'},fn,'PS',2);
      if length(xs)~=length(ps)
        error(natemsgid,'Input 1, XS, and input 2, PS, must have the same number of entries.');
      end
    case 'interpuni'
      [ps,ab] = deal(varargin{:});
      validateattributes(ps,{'numeric'},{'2d','nonempty'},fn,'PS',1);
      validateattributes(ab,{'numeric'},{'numel',2},fn,'AB',2);
    case 'intmatrix_'
      [xs,ws,xs_,a] = deal(varargin{:});
      validateattributes(xs,{'numeric'},{'column'},fn,'XS',1);
      validateattributes(xs,{'numeric'},{'column'},fn,'WS',2);
      if length(xs)~=length(ws)
        error(natemsgid,'Input 1, XS, and input 2, WS, must have the same number of entries.');
      end
      validateattributes(xs,{'numeric'},{'column'},fn,'XS_',3);
      if length(xs)~=length(xs_)+1
        error(natemsgid,'Input 1, XS, must have one more entry than input 3, XS_.');
      end
      validateattributes(a,{'numeric'},{'scalar'},fn,'A',4);
    case 'intmatrixcheb'
      m = varargin{1};
      validateattributes(m,{'numeric'},{'scalar','positive'},fn,'M',1);
    case 'intmatrixgen'
      [xs,xs_,a] = deal(varargin{:});
      validateattributes(xs,{'numeric'},{'column'},fn,'XS',1);
      validateattributes(xs_,{'numeric'},{'column'},fn,'XS_',2);
      if length(xs)~=length(xs_)+1
        error(natemsgid,'Input 1, XS, must have one more entry than input 2, XS_.');
      end
      validateattributes(a,{'numeric'},{'scalar'},fn,'A',3);
    case 'intmatrixuni'
      m = varargin{1};
      validateattributes(m,{'numeric'},{'scalar','nonnegative','integer'},fn,'M',1);
    case 'ivp1matcheb'
      [m,l] = deal(varargin{:});
      validateattributes(m,{'numeric'},{'scalar','positive','integer'},fn,'M',1);
      validateattributes(l,{'numeric'},{'scalar'},fn,'L',2);
    case 'ivp1matuni'
      [m,l] = deal(varargin{:});
      validateattributes(m,{'numeric'},{'scalar','nonnegative','integer'},fn,'M',1);
      validateattributes(l,{'numeric'},{'scalar'},fn,'L',2);
    case 'ivp2matcheb'
      [m,l] = deal(varargin{:});
      validateattributes(m,{'numeric'},{'scalar','positive','integer'},fn,'M',1);
      validateattributes(l,{'numeric'},{'scalar'},fn,'L',2);
    case 'ivp2matuni'
      [m,l] = deal(varargin{:});
      validateattributes(m,{'numeric'},{'scalar','nonnegative','integer'},fn,'M',1);
      validateattributes(l,{'numeric'},{'scalar'},fn,'L',2);
    case 'ivpl1cheb'
      [al,be,g,ab,c,d,m] = deal(varargin{:});
      validateattributes(al,{'function_handle'},{},fn,'AL',1);
      validateattributes(be,{'function_handle'},{},fn,'BE',2);
      validateattributes(g,{'function_handle'},{},fn,'G',3);
      validateattributes(ab,{'numeric'},{'numel',2},fn,'AB',4);
      validateattributes(c,{'numeric'},{'numel',2},fn,'C',5);
      validateattributes(d,{'numeric'},{'scalar'},fn,'D',6);
      validateattributes(m,{'numeric'},{'scalar','positive','integer'},fn,'M',7);
    case 'ivpl1condcheb'
      [al,be,ab,c,m] = deal(varargin{:});
      validateattributes(al,{'function_handle'},{},fn,'AL',1);
      validateattributes(be,{'function_handle'},{},fn,'BE',2);
      validateattributes(ab,{'numeric'},{'numel',2},fn,'AB',3);
      validateattributes(c,{'numeric'},{'numel',2},fn,'C',4);
      validateattributes(m,{'numeric'},{'scalar','positive','integer'},fn,'M',5);
    case 'ivpl1condchebpw'
      [al,be,asbs,c,m] = deal(varargin{:});
      validateattributes(al,{'function_handle'},{},fn,'AL',1);
      validateattributes(be,{'function_handle'},{},fn,'BE',2);
      validateattributes(asbs,{'numeric'},{'nrows',2},fn,'ASBS',3);
      validateattributes(c,{'numeric'},{'numel',2},fn,'C',4);
      validateattributes(m,{'numeric'},{'scalar','positive','integer'},fn,'M',5);
    case 'ivpl1conduni'
      [al,be,ab,c,m,n] = deal(varargin{:});
      validateattributes(al,{'function_handle'},{},fn,'AL',1);
      validateattributes(be,{'function_handle'},{},fn,'BE',2);
      validateattributes(ab,{'numeric'},{'numel',2},fn,'AB',3);
      validateattributes(c,{'numeric'},{'numel',2},fn,'C',4);
      validateattributes(m,{'numeric'},{'scalar','nonnegative','integer'},fn,'M',5);
      validateattributes(n,{'numeric'},{'scalar','positive','integer'},fn,'N',6);
    case 'ivpl1gen'
      [al,be,g,a,c,d,ss,ts,ts_] = deal(varargin{:});
      validateattributes(al,{'function_handle'},{},fn,'AL',1);
      validateattributes(be,{'function_handle'},{},fn,'BE',2);
      validateattributes(g,{'function_handle'},{},fn,'G',3);
      validateattributes(a,{'numeric'},{'scalar'},fn,'A',4);
      validateattributes(c,{'numeric'},{'numel',2},fn,'C',5);
      validateattributes(d,{'numeric'},{'scalar'},fn,'D',6);
      validateattributes(ss,{'numeric'},{'column'},fn,'SS',7);
      validateattributes(ts,{'numeric'},{'column'},fn,'TS',8);
      if length(ts)~=length(ss)+1
        error(natemsgid,'Input 8, TS, must have one more entry than input 7, SS.');
      end
      validateattributes(ts_,{'numeric'},{'column'},fn,'TS_',9);
      if length(ts_)~=length(ss)
        error(natemsgid,'Input 7, SS, and input 9, TS_, must have the same number of entries.');
      end
    case 'ivpl1uni'
      [al,be,g,ab,c,d,m,n] = deal(varargin{:});
      validateattributes(al,{'function_handle'},{},fn,'AL',1);
      validateattributes(be,{'function_handle'},{},fn,'BE',2);
      validateattributes(g,{'function_handle'},{},fn,'G',3);
      validateattributes(ab,{'numeric'},{'numel',2},fn,'AB',4);
      validateattributes(c,{'numeric'},{'numel',2},fn,'C',5);
      validateattributes(d,{'numeric'},{'scalar'},fn,'D',6);
      validateattributes(m,{'numeric'},{'scalar','nonnegative','integer'},fn,'M',7);
      validateattributes(n,{'numeric'},{'scalar','positive','integer'},fn,'N',8);
    case 'ivpl2cheb'
      [al,be,ga,g,ab,c,d,m] = deal(varargin{:});
      validateattributes(al,{'function_handle'},{},fn,'AL',1);
      validateattributes(be,{'function_handle'},{},fn,'BE',2);
      validateattributes(ga,{'function_handle'},{},fn,'GA',3);
      validateattributes(g,{'function_handle'},{},fn,'G',4);
      validateattributes(ab,{'numeric'},{'numel',2},fn,'AB',5);
      validateattributes(c,{'numeric'},{'nrows',2,'ncols',3},fn,'C',6);
      validateattributes(d,{'numeric'},{'column','numel',2},fn,'D',7);
      validateattributes(m,{'numeric'},{'scalar','positive','integer'},fn,'M',8);
    case 'ivpl2condcheb'
      [al,be,ga,ab,c,m] = deal(varargin{:});
      validateattributes(al,{'function_handle'},{},fn,'AL',1);
      validateattributes(be,{'function_handle'},{},fn,'BE',2);
      validateattributes(ga,{'function_handle'},{},fn,'GA',3);
      validateattributes(ab,{'numeric'},{'numel',2},fn,'AB',4);
      validateattributes(c,{'numeric'},{'nrows',2,'ncols',3},fn,'C',5);
      validateattributes(m,{'numeric'},{'scalar','positive','integer'},fn,'M',6);
    case 'ivpl2condchebpw'
      [al,be,ga,asbs,c,m] = deal(varargin{:});
      validateattributes(al,{'function_handle'},{},fn,'AL',1);
      validateattributes(be,{'function_handle'},{},fn,'BE',2);
      validateattributes(ga,{'function_handle'},{},fn,'GA',3);
      validateattributes(asbs,{'numeric'},{'nrows',2},fn,'ASBS',4);
      validateattributes(c,{'numeric'},{'nrows',2,'ncols',3},fn,'C',5);
      validateattributes(m,{'numeric'},{'scalar','positive','integer'},fn,'M',6);
    case 'ivpl2conduni'
      [al,be,ga,ab,c,m,n] = deal(varargin{:});
      validateattributes(al,{'function_handle'},{},fn,'AL',1);
      validateattributes(be,{'function_handle'},{},fn,'BE',2);
      validateattributes(ga,{'function_handle'},{},fn,'GA',3);
      validateattributes(ab,{'numeric'},{'numel',2},fn,'AB',4);
      validateattributes(c,{'numeric'},{'nrows',2,'ncols',3},fn,'C',5);
      validateattributes(m,{'numeric'},{'scalar','nonnegative','integer'},fn,'M',6);
      validateattributes(n,{'numeric'},{'scalar','positive','integer'},fn,'N',7);
    case 'ivpl2gen'
      [al,be,ga,g,a,c,d,ss,ts,ts_,ts__] = deal(varargin{:});
      validateattributes(al,{'function_handle'},{},fn,'AL',1);
      validateattributes(be,{'function_handle'},{},fn,'BE',2);
      validateattributes(ga,{'function_handle'},{},fn,'GA',3);
      validateattributes(g,{'function_handle'},{},fn,'G',4);
      validateattributes(a,{'numeric'},{'scalar'},fn,'A',5);
      validateattributes(c,{'numeric'},{'nrows',2,'ncols',3},fn,'C',6);
      validateattributes(d,{'numeric'},{'column','numel',2},fn,'D',7);
      validateattributes(ss,{'numeric'},{'column'},fn,'SS',8);
      validateattributes(ts,{'numeric'},{'column'},fn,'TS',9);
      validateattributes(ts_,{'numeric'},{'column'},fn,'TS_',10);
      validateattributes(ts__,{'numeric'},{'column'},fn,'TS__',11);
      m = length(ts__)-1;
      if length(ss)~=m+1||length(ts)~=m+3||length(ts_)~=m+2
        error(natemsgid,'SS, TS, TS_, and TS__ must have degrees M, M+2, M+1, and M, respectively');
      end
    case 'ivpl2uni'
      [al,be,ga,g,ab,c,d,m,n] = deal(varargin{:});
      validateattributes(al,{'function_handle'},{},fn,'AL',1);
      validateattributes(be,{'function_handle'},{},fn,'BE',2);
      validateattributes(ga,{'function_handle'},{},fn,'GA',3);
      validateattributes(g,{'function_handle'},{},fn,'G',4);
      validateattributes(ab,{'numeric'},{'numel',2},fn,'AB',5);
      validateattributes(c,{'numeric'},{'nrows',2,'ncols',3},fn,'C',6);
      validateattributes(d,{'numeric'},{'column','numel',2},fn,'D',7);
      validateattributes(m,{'numeric'},{'scalar','nonnegative','integer'},fn,'M',8);
      validateattributes(n,{'numeric'},{'scalar','positive','integer'},fn,'N',9);
    case 'ivpnl1chebpw'
      [f,dfdu,ab,ua,m,n,kmax,tol] = deal(varargin{:});
      validateattributes(f,{'function_handle'},{},fn,'F',1);
      validateattributes(dfdu,{'function_handle'},{},fn,'DFDU',2);
      validateattributes(ab,{'numeric'},{'numel',2},fn,'AB',3);
      validateattributes(ua,{'numeric'},{'scalar'},fn,'UA',4);
      validateattributes(m,{'numeric'},{'scalar','positive','integer'},fn,'M',5);
      validateattributes(n,{'numeric'},{'scalar','positive','integer'},fn,'N',6);
      validateattributes(kmax,{'numeric'},{'scalar','nonnegative','integer'},fn,'KMAX',7);
      validateattributes(tol,{'numeric'},{'scalar','nonnegative'},fn,'TOL',8);
    case 'ivpnl1condcheb'
      [dfdu,ab,m,p] = deal(varargin{:});
      validateattributes(dfdu,{'function_handle'},{},fn,'DFDU',1);
      validateattributes(ab,{'numeric'},{'numel',2},fn,'AB',2);
      validateattributes(m,{'numeric'},{'scalar','positive','integer'},fn,'M',3);
      validateattributes(p,{'function_handle'},{},fn,'P',4);
    case 'ivpnl1condchebpw'
      [dfdu,asbs,m,p] = deal(varargin{:});
      validateattributes(dfdu,{'function_handle'},{},fn,'DFDU',1);
      validateattributes(asbs,{'numeric'},{'nrows',2},fn,'ASBS',2);
      validateattributes(m,{'numeric'},{'scalar','positive','integer'},fn,'M',3);
      validateattributes(p,{'function_handle'},{},fn,'P',4);
    case 'ivpnl1conduni'
      [dfdu,ab,m,n,p] = deal(varargin{:});
      validateattributes(dfdu,{'function_handle'},{},fn,'DFDU',1);
      validateattributes(ab,{'numeric'},{'numel',2},fn,'AB',2);
      validateattributes(m,{'numeric'},{'scalar','nonnegative','integer'},fn,'M',3);
      validateattributes(n,{'numeric'},{'scalar','positive','integer'},fn,'N',4);
      validateattributes(p,{'function_handle'},{},fn,'P',5);
    case 'ivpnl1gen'
      [f,dfdu,a,ua,qstart,ss,ts,ts_,kmax,tol] = deal(varargin{:});
      validateattributes(f,{'function_handle'},{},fn,'F',1);
      validateattributes(dfdu,{'function_handle'},{},fn,'DFDU',2);
      validateattributes(a,{'numeric'},{'scalar'},fn,'A',3);
      validateattributes(ua,{'numeric'},{'scalar'},fn,'UA',4);
      validateattributes(qstart,{'function_handle'},{},fn,'QSTART',5);
      validateattributes(ss,{'numeric'},{'column'},fn,'SS',6);
      validateattributes(ts,{'numeric'},{'column'},fn,'TS',7);
      validateattributes(ts_,{'numeric'},{'column'},fn,'TS_',8);
      validateattributes(kmax,{'numeric'},{'scalar','nonnegative','integer'},fn,'KMAX',9);
      validateattributes(tol,{'numeric'},{'scalar','nonnegative'},fn,'TOL',10);
      m = length(ss)-1;
      if length(ts)~=m+2||length(ts_)~=m+1
        error(natemsgid,'SS, TS, and TS_ must have degrees M, M+1, and M, respectively.');
      end
    case 'ivpnl1uni'
      [f,dfdu,ab,ua,m,n,kmax,tol] = deal(varargin{:});
      validateattributes(f,{'function_handle'},{},fn,'F',1);
      validateattributes(dfdu,{'function_handle'},{},fn,'DFDU',2);
      validateattributes(ab,{'numeric'},{'numel',2},fn,'AB',3);
      validateattributes(ua,{'numeric'},{'scalar'},fn,'UA',4);
      validateattributes(m,{'numeric'},{'scalar','nonnegative','integer'},fn,'M',5);
      validateattributes(n,{'numeric'},{'scalar','positive','integer'},fn,'N',6);
      validateattributes(kmax,{'numeric'},{'scalar','nonnegative','integer'},fn,'KMAX',7);
      validateattributes(tol,{'numeric'},{'scalar','nonnegative'},fn,'TOL',8);
    case 'ivpnl2chebpw'
      [f,dfdu,dfdv,ab,ua,va,m,n,kmax,tol] = deal(varargin{:});
      validateattributes(f,{'function_handle'},{},fn,'F',1);
      validateattributes(dfdu,{'function_handle'},{},fn,'DFDU',2);
      validateattributes(dfdv,{'function_handle'},{},fn,'DFDV',3);
      validateattributes(ab,{'numeric'},{'numel',2},fn,'AB',4);
      validateattributes(ua,{'numeric'},{'scalar'},fn,'UA',5);
      validateattributes(va,{'numeric'},{'scalar'},fn,'VA',6);
      validateattributes(m,{'numeric'},{'scalar','positive','integer'},fn,'M',7);
      validateattributes(n,{'numeric'},{'scalar','positive','integer'},fn,'N',8);
      validateattributes(kmax,{'numeric'},{'scalar','nonnegative','integer'},fn,'KMAX',9);
      validateattributes(tol,{'numeric'},{'scalar','nonnegative'},fn,'TOL',10);
    case 'ivpnl2condcheb'
      [dfdu,dfdv,ab,m,p,q] = deal(varargin{:});
      validateattributes(dfdu,{'function_handle'},{},fn,'DFDU',1);
      validateattributes(dfdv,{'function_handle'},{},fn,'DFDV',2);
      validateattributes(ab,{'numeric'},{'numel',2},fn,'AB',3);
      validateattributes(m,{'numeric'},{'scalar','positive','integer'},fn,'M',4);
      validateattributes(p,{'function_handle'},{},fn,'P',5);
      validateattributes(q,{'function_handle'},{},fn,'Q',6);
    case 'ivpnl2condchebpw'
      [dfdu,dfdv,asbs,m,p,q] = deal(varargin{:});
      validateattributes(dfdu,{'function_handle'},{},fn,'DFDU',1);
      validateattributes(dfdv,{'function_handle'},{},fn,'DFDV',2);
      validateattributes(asbs,{'numeric'},{'nrows',2},fn,'ASBS',3);
      validateattributes(m,{'numeric'},{'scalar','positive','integer'},fn,'M',4);
      validateattributes(p,{'function_handle'},{},fn,'P',5);
      validateattributes(q,{'function_handle'},{},fn,'Q',6);
    case 'ivpnl2conduni'
      [dfdu,dfdv,ab,m,n,p,q] = deal(varargin{:});
      validateattributes(dfdu,{'function_handle'},{},fn,'DFDU',1);
      validateattributes(dfdv,{'function_handle'},{},fn,'DFDV',2);
      validateattributes(ab,{'numeric'},{'numel',2},fn,'AB',3);
      validateattributes(m,{'numeric'},{'scalar','nonnegative','integer'},fn,'M',4);
      validateattributes(n,{'numeric'},{'scalar','positive','integer'},fn,'N',5);
      validateattributes(p,{'function_handle'},{},fn,'P',6);
      validateattributes(q,{'function_handle'},{},fn,'Q',7);
    case 'ivpnl2gen'
      [f,dfdu,dfdv,a,ua,va,rstart,ss,ts,ts_,ts__,kmax,tol] = deal(varargin{:});
      validateattributes(f,{'function_handle'},{},fn,'F',1);
      validateattributes(dfdu,{'function_handle'},{},fn,'DFDU',2);
      validateattributes(dfdv,{'function_handle'},{},fn,'DFDV',3);
      validateattributes(a,{'numeric'},{'scalar'},fn,'A',4);
      validateattributes(ua,{'numeric'},{'scalar'},fn,'UA',5);
      validateattributes(va,{'numeric'},{'scalar'},fn,'VA',6);
      validateattributes(rstart,{'function_handle'},{},fn,'RSTART',7);
      validateattributes(ss,{'numeric'},{'column'},fn,'SS',8);
      validateattributes(ts,{'numeric'},{'column'},fn,'TS',9);
      validateattributes(ts_,{'numeric'},{'column'},fn,'TS_',10);
      validateattributes(ts__,{'numeric'},{'column'},fn,'TS__',11);
      validateattributes(kmax,{'numeric'},{'scalar','nonnegative','integer'},fn,'KMAX',12);
      validateattributes(tol,{'numeric'},{'scalar','nonnegative'},fn,'TOL',13);
      m = length(ss)-1;
      if length(ts)~=m+3||length(ts_)~=m+2||length(ts__)~=m+1
        error(natemsgid,'SS, TS, TS_, and TS__ must have degrees M, M+2, M+1, and M, respectively.');
      end
    case 'ivpnl2uni'
      [f,dfdu,dfdv,ab,ua,va,m,n,kmax,tol] = deal(varargin{:});
      validateattributes(f,{'function_handle'},{},fn,'F',1);
      validateattributes(dfdu,{'function_handle'},{},fn,'DFDU',2);
      validateattributes(dfdv,{'function_handle'},{},fn,'DFDV',3);
      validateattributes(ab,{'numeric'},{'numel',2},fn,'AB',4);
      validateattributes(ua,{'numeric'},{'scalar'},fn,'UA',5);
      validateattributes(va,{'numeric'},{'scalar'},fn,'VA',6);
      validateattributes(m,{'numeric'},{'scalar','nonnegative','integer'},fn,'M',7);
      validateattributes(n,{'numeric'},{'scalar','positive','integer'},fn,'N',8);
      validateattributes(kmax,{'numeric'},{'scalar','nonnegative','integer'},fn,'KMAX',9);
      validateattributes(tol,{'numeric'},{'scalar','nonnegative'},fn,'TOL',10);
    case 'lunopivot'
      A = deal(varargin{:});
      validateattributes(A,{'numeric'},{'square'},fn,'A',1);
    case 'lupp'
      A = deal(varargin{:});
      validateattributes(A,{'numeric'},{'square'},fn,'A',1);
    case 'newfig'
      [m,n] = deal(varargin{:});
      validateattributes(m,{'numeric'},{'scalar','positive','integer'},fn,'M',1);
      validateattributes(n,{'numeric'},{'scalar','positive','integer'},fn,'N',2);
    case 'newton'
      [f,fprime,x0,kmax,tol] = deal(varargin{:});
      validateattributes(f,{'function_handle'},{},fn,'F',1);
      validateattributes(fprime,{'function_handle'},{},fn,'FPRIME',2);
      validateattributes(x0,{'numeric'},{'scalar'},fn,'X0',3);
      validateattributes(kmax,{'numeric'},{'scalar','nonnegative','integer'},fn,'KMAX',4);
      validateattributes(tol,{'numeric'},{'scalar','nonnegative'},fn,'TOL',5);
    case 'partition_'
      [ab,n] = deal(varargin{:});
      validateattributes(ab,{'numeric'},{'numel',2},fn,'AB',1);
      validateattributes(n,{'numeric'},{'scalar','positive','integer'},fn,'N',2);
    case 'plotfun'
      [f,ab] = deal(varargin{:});
      validateattributes(f,{'function_handle'},{},fn,'F',1);
      validateattributes(ab,{'numeric'},{'numel',2},fn,'AB',2);
    case 'plotfun3'
      [fx,fy,fz,ab] = deal(varargin{:});
      validateattributes(fx,{'function_handle'},{},fn,'FX',1);
      validateattributes(fy,{'function_handle'},{},fn,'FY',2);
      validateattributes(fz,{'function_handle'},{},fn,'FZ',3);
      validateattributes(ab,{'numeric'},{'numel',2},fn,'AB',4);
    case 'plotimagpart'
      [f,lims] = deal(varargin{:});
      validateattributes(f,{'function_handle'},{},fn,'F',1);
      validateattributes(lims,{'numeric'},{'vector','numel',6},fn,'LIMS',2);
    case 'plotpartition'
      if nargin==2
        asbs = deal(varargin{:});
        validateattributes(asbs,{'numeric'},{'nrows',2},fn,'ASBS',1);
      else
        [ab,n] = deal(varargin{:});
        validateattributes(ab,{'numeric'},{'numel',2},fn,'AB',1);
        validateattributes(n,{'numeric'},{'scalar','positive','integer'},fn,'N',2);
      end
    case 'plotrealpart'
      [f,lims] = deal(varargin{:});
      validateattributes(f,{'function_handle'},{},fn,'F',1);
      validateattributes(lims,{'numeric'},{'vector','numel',6},fn,'LIMS',2);
    case 'plotsample'
      [xs,ys] = deal(varargin{:});
      validateattributes(xs,{'numeric'},{'2d'},fn,'XS',1);
      validateattributes(ys,{'numeric'},{'2d'},fn,'YS',2);
      if numel(xs)~=numel(ys)
        error(natemsgid,'Input 1, XS, and input 2, YS, must have the same number of entries.');
      end
    case 'polynomialeval'
      [c,x] = deal(varargin{:});
      validateattributes(c,{'numeric'},{'vector'},fn,'C',1);
      validateattributes(x,{'numeric'},{'scalar'},fn,'X',2);
    case 'randpoints_'
      n = deal(varargin{:});
      validateattributes(n,{'numeric'},{'scalar','nonnegative','integer'},fn,'N',1);
    case 'resamplematrixcheb'
      [l,m] = deal(varargin{:});
      validateattributes(l,{'numeric'},{'scalar','positive','integer'},fn,'L',1);
      validateattributes(m,{'numeric'},{'scalar','positive','integer'},fn,'M',2);      
    case 'resamplematrixuni'
      [l,m] = deal(varargin{:});
      validateattributes(l,{'numeric'},{'scalar','nonnegative','integer'},fn,'L',1);
      validateattributes(m,{'numeric'},{'scalar','nonnegative','integer'},fn,'M',2);
    case 'rungecurve'
      ab = deal(varargin{:});
      validateattributes(ab,{'numeric'},{'numel',2},fn,'AB',1);
    case 'samplecheb'
      [f,ab,m] = deal(varargin{:});
      validateattributes(f,{'function_handle'},{},fn,'F',1);
      validateattributes(ab,{'numeric'},{'numel',2},fn,'AB',2);
      validateattributes(m,{'numeric'},{'scalar','positive','integer'},fn,'M',3);      
    case 'samplechebpw'
      [f,asbs,m] = deal(varargin{:});
      validateattributes(f,{'function_handle'},{},fn,'F',1);
      validateattributes(asbs,{'numeric'},{'nrows',2},fn,'ASBS',2);
      validateattributes(m,{'numeric'},{'scalar','positive','integer'},fn,'M',3);      
    case 'sampleuni'
      [f,ab,m,n] = deal(varargin{:});
      validateattributes(f,{'function_handle'},{},fn,'F',1);
      validateattributes(ab,{'numeric'},{'numel',2},fn,'AB',2);
      validateattributes(m,{'numeric'},{'scalar','nonnegative','integer'},fn,'M',3);      
      validateattributes(n,{'numeric'},{'scalar','positive','integer'},fn,'N',4);      
    case 'zeros_'
      [xs,ws,ps] = deal(varargin{:});
      validateattributes(xs,{'numeric'},{'column'},fn,'XS',1);
      validateattributes(ws,{'numeric'},{'column'},fn,'WS',2);
      validateattributes(ps,{'numeric'},{'column'},fn,'PS',3);
      m = length(xs)-1;
      if length(ws)~=m+1||length(ps)~=m+1
        error(natemsgid,'Inputs 1, 2, and 3 (XS, WS, and PS) must have the same length.');
      end
    case 'zeroscheb'
      [ps,ab,tol] = deal(varargin{:});
      validateattributes(ps,{'numeric'},{'column'},fn,'PS',1);
      validateattributes(ab,{'numeric'},{'numel',2},fn,'AB',2);
      validateattributes(tol,{'numeric'},{'scalar','nonnegative'},fn,'TOL',3);
    case 'zerosgen'
      [xs,ps] = deal(varargin{:});
      validateattributes(xs,{'numeric'},{'column'},fn,'XS',1);
      validateattributes(ps,{'numeric'},{'column'},fn,'PS',2);
      if length(xs)~=length(ps)
        error(natemsgid,'Input 1, XS, and input 2, PS, must have the same length.');
      end
    case 'zerosuni'
      [ps,ab,tol] = deal(varargin{:});
      validateattributes(ps,{'numeric'},{'2d'},fn,'PS',1);
      validateattributes(ab,{'numeric'},{'numel',2},fn,'AB',2);
      validateattributes(tol,{'numeric'},{'scalar','nonnegative'},fn,'TOL',3);
    otherwise
      error(natemsgid,'natecheck not implemented for %s',fn);
  end

catch e
  
  throwAsCaller(e);
  
end

