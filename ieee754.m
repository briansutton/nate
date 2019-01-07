function [sbin,ebin,fbin] = ieee754(x)
%IEEE754   IEEE 754 double-precision representation.
%
%   [SBIN,EBIN,FBIN] = IEEE754(X) reveals the IEEE 754 representation of a
%   double-precision number.
%
%   Inputs:
%
%     X: a double-precision number
%
%   Outputs:
%
%     SBIN: sign bit for X in binary notation
%
%     EBIN: shifted exponent for X in binary notation
%
%     FBIN: significand for X in binary notation, except for the implicit
%     leading 1
%
%   Example:
%
%     [sbin,ebin,fbin] = ieee754(-9.25)
%     s = bin2dec(sbin)
%     be = bin2dec(ebin)-1023
%     f = 2^(-52)*bin2dec(['1' fbin])
%     (-1)^s*f*2^be
%
%   Copyright 2019 Brian Sutton

narginchk(1,1);
natecheck('ieee754',x);
if imag(x)~=0, warning('Imaginary part is ignored.'); end
x = real(x);
if x==0
  if 1/x<0, sbin = '1'; else sbin = '0'; end
else
  if x<0, sbin = '1'; else sbin = '0'; end
end
[fbin,ebin] = log2(x);
fbin = abs(fbin);
if ebin<-1022
  error([ 'Denormalized numbers are not supported. ' ...
          '(x is too close to zero.)' ]);
end
if isinf(x)
  fbin = 0; ebin = 1024;
elseif isnan(x)
  fbin = 1-2^(-53); ebin = 1024;
elseif x==0
  fbin = 0; ebin = -1023;
else
  fbin = 2*fbin;
  ebin = ebin-1;
end
fbin = dec2bin(pow2(fbin-1,52),52);
ebin = dec2bin(ebin+1023,11);

