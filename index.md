# *Numerical Analysis: Theory and Experiments*<br/>Brian Sutton

[*Numerical Analysis: Theory and Experiments*](http://bookstore.siam.org/ot161/) is a textbook on numerical analysis. Numerical methods are designed and implemented, and then they are analyzed through a combination of mathematical theory and numerical experimentation. Problem areas include interpolation, integration, linear systems, zero finding, and differential equations.

This package is a library of MATLAB routines that accompany the book.

## Installation

The library is installed as follows:
1. Download the package [nate-master.zip](https://github.com/briansutton/nate/archive/master.zip).
1. Expand the zip file if it was not expanded automatically.
1. Start MATLAB.
1. Use the `pathtool` dialog to add the new package folder to the MATLAB path:
   ```
   >> pathtool
   ```
1. Verify that the library routines are available as follows:
   ```
   >> nate
   Numerical Analysis: Theory and Experiments
   ```
1. Set the Command Window preferences in MATLAB to `long` numeric format and `compact` numeric display.

## Demonstration

The following code computes and plots the arctan function.
```
>> g = @(x) 1/(1+x^2); a = 0; b = 10;
>> fa = 0;
>> qs = samplecheb(g,[a b],60);
>> ps = antiderivcheb(fa,qs,[a b]);
>> p = interpcheb(ps,[a b]);
>> newfig;
>> plotfun(p,[a b],'displayname','p');
>> ylim([0 pi/2]);
>> legend('location','southeast');
```
![arctan](https://github.com/briansutton/nate/raw/master/arctan.png "arctan")
