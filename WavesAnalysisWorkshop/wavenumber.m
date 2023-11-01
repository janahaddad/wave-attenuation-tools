function [L, k] = wavenumber(omega,h)
%k = wavenumber(omega,h)
%
%k is a vector of same size as omega and depth containing the calculated wave numbers
%omega is the wave frequencies in rad/s
%h is the water depth
%sigma and h must be scalars or vectors of the same dimensions
%
%modified from R.Dalrymple's java code
%then modified from J.Rosman's m-file (J.Haddad) 

g=9.81;

a0=(omega.*omega.*h)./g; 
b1=1.0./(tanh(a0.^0.75));
a1=a0.*(b1.^0.666);
da1=1000.0;
   
d1=ones(size(h));
while(max(d1)==1)  
d1 = (abs(da1./a1) > .00000001);
	th=tanh(a1);
	ch=cosh(a1);
	f1=a0-(a1.*th);
	f2= - a1.*((1.0./ch).^2) -th;
	da1= -f1./f2;
	a1=a1+da1;
end

k=a1./h;
L = 2*pi./k;
   
  
