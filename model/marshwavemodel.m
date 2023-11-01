function [H,h,Dv,Db,Df,gamma,KC,Re,Ca,Cd,sigma,L,k,Ur,S,Er,Uc,a] = ...
    marshwavemodel(dx,n,Hrms0,h0,z,T,veg,ls,g,rho,Nv,bv,E,...
    gammac,Cf,BRK,FRIC,deplim,DRAG,wavemode)
%% Initialize & Preallocate 
h  = zeros(1,n);    % local water depth
H2 = zeros(1,n);    % H2 = wave height squared 
H2(1) = Hrms0^2;    % H2 at x = 0
H = zeros(1,n);
H(1) = Hrms0;       % rms wave height at x = 0
cg = zeros(1,n);    % cg = group velocity = c*(1/2)*(1 + 2*k*h/sinh(2*k*h))
V  = zeros(1,n);    % vegetation coefficient (k, h, ah dependent)
R  = zeros(1,n);    % breaking coefficient (Tp, h dependent) 
F  = zeros(1,n);    % friction coefficient
Dv = zeros(1,n); 
Db = zeros(1,n); 
Df = zeros(1,n); 
ah = zeros(1,n);
a = zeros(1,n);     % a = alpha = relative vegetation height
I = zeros(1,n);     % second moment of area for circ stem x-section 
EI = zeros(1,n);    % flexural rigidity (Pa * m^4) 
Er = zeros(1,n);
Ur = zeros(1,n);   % ursell number = HL^2/h^2
S = zeros(1,n);    % wave steepness = H/L
Re  = zeros(1,n);
Ca = zeros(1,n);
KC = zeros(1,n); 
Cd = zeros(1,n);
beta = zeros(1,n);
Uc = zeros(1,n); 

%% bottom slope -- this is used for BRK = 1 or 2 
for i = 1:n-1  
beta(i+1) = (z(i+1) - z(i))/(dx);
end 
%% Turn off model output for water depths < deplim by making it dry :
for i = 1:n 
    h(i) = h0 - z(i);
    if h(i) < deplim
        h(i) = NaN;
    end
end
%% Generate vectors that can be computed outside solver loop
sigma  = (2*pi)./T;
[L,k] = wavenumber(sigma,h);
gamma = 0.23 +1.42.*tan(beta);

for i = 1:n  
    
 % wave group velocity (m/s):
 cg(i) = sqrt((g/k(i))*tanh(k(i)*h(i)))...   
     *0.5*(1 + 2*k(i)*h(i)/sinh(2*k(i)*h(i)));  
 
 % compute V coeff for veg dissipation:
 if veg(i) == 0  % if there's no veg there, make all the veg variables 0. 
     ah(i) = 0;
     a(i) = 0;
     bv(i) = 0;
     Nv(i) = 0; 
     V(i) = 0;
 else            % if there is veg there, compute veg variables. 
     bv(i)= bv(i);
     Nv(i) = Nv(i)+1;
     ah(i) = ls(i)+0.001; 
     Er(i) = ah(i)/h(i);
     a(i) = min(Er(i),1);
     r = bv(i)/2; I(i) = (pi*r^4)/4; %m^4, for circular section 
     EI(i) = E*I(i); 
       if strcmp(wavemode,'mono')==1  
       V(i) = -(2/(3*pi)) * rho * bv(i) * Nv(i) * (k(i)*g/(2*sigma))^3 * ...
      ((sinh(k(i)*(a(i)*h(i)))^3 + 3*sinh(k(i)*(a(i)*h(i)))) / (3*k(i)*(cosh(k(i)*h(i))^3)));
       elseif strcmp(wavemode,'rand')==1
       V(i) = -(1/(2*sqrt(pi))) * rho * bv(i) * Nv(i) * (k(i)*g/(2*sigma))^3 * ...
      ((sinh(k(i)*(a(i)*h(i)))^3 + 3*sinh(k(i)*(a(i)*h(i)))) / (3*k(i)*(cosh(k(i)*h(i))^3)));
       end 
 end
 
 % compute F coeff for bottom friction dissipation:
 if FRIC == 1  
     F(i) = - (rho*Cf/(16*sqrt(pi))) * (2*pi*(1/T)/sinh(k(i)*h(i))).^3;
 end 
 
 % compute R coeff for wave breaking dissipation:
B = 1;
 if BRK == 0 % no breaking
     R(i) = 0;
 elseif BRK == 1 %TG83 eq 25, gamma =f(slope)
     R(i) = - (3*sqrt(pi)/16)*B*rho*g*(1/T)  / ((gamma(i)^4)*(h(i)^5));
 elseif BRK == 2 %TG83 eq 26, gamma =f(slope)
%      R(i) = - (3*sqrt(pi)/16)*B*rho*g*(1/T)  / ((gamma(i)^2)*(h(i)^3));
     R(i) = - (3*sqrt(pi)/16)*B*rho*g*(1/T) ;
 elseif BRK == 3 % TG83 eq 25 constant gamma
     R(i) = - (3*sqrt(pi)/16)*B*rho*g*(1/T)  / ((gammac^4)*(h(i)^5));
 elseif BRK == 4 %TG83 eq 26 constant gamma
     R(i) = - (3*sqrt(pi)/16)*B*rho*g*(1/T) ; 
 end  
 
end 

%% Solve by forward difference
for j = 1:n-1  
    Ur(j) = H(j)*L(j)^2/h(j)^3; % compute Ursell number
    S(j)  = H(j)/L(j);          % compute wave steepness 
    UCoeff = cosh(k(j)* a(j)*h(j))/sinh(k(j)*h(j));
     Uc(j) =(pi*H(j)/T) * UCoeff;
     if veg(j) == 1 
     KC(j) = (Uc(j)*T)/bv(j);
     Re(j) = (Uc(j)*bv(j))/1e-6;
     Ca(j) = (rho.*bv(j).*(Uc(j).^2).*(a(j).*h(j)).^3)./EI(j);
     end 
 % compute the drag coefficients Cd that will be used to compute
 % veg-induced dissipation, based on the DRAG input: 
 if DRAG > 0        % Cd = CONSTANT VALUE
     if veg(j) == 1
%      UCoeff = cosh(k(j)* a(j)*h(j))/sinh(k(j)*h(j));
%      Uc =(pi*H(j)/T) * UCoeff;

     Cd(j) = DRAG; 
        if KC(j) > 20 
         Cd(j) = DRAG;
        else
            Cd(j) = NaN;
        end
     else 
     Cd(j) = 0;  % needs to be 0 (not NaN) for wave heights to compute

     end 
 end 
 
 if DRAG == -1   % Cd(KC) = (DF18 fit) 
     if veg(j) == 1 %&& KC(j) > 20
        
     Cd(j) = 43.5*KC(j)^-0.55; 
     else 
     Cd(j) = 0;

     end 
 end 

if DRAG == -2   % Cd(Re) = DF18 fit 
  if veg(j) == 1  %&& KC(j) > 20 

     Cd(j) = 149.1*Re(j)^-0.55;

 else 
     Cd(j) = 0;

  end 
end 

if DRAG == -3
    if veg(j) == 1 
%         Cd(j)  =6*Ca(j)^-0.3;
        Cd(j)  =8.3*Ca(j)^-0.33;
    else 
        Cd(j) = 0;
    end 
end 

if DRAG == -3.1   % Cd(Ca) for KC > 100, Cd(KC) for KC <100 
  if veg(j) == 1 && KC(j) > 80
    Cd(j)  = 8.3*Ca(j)^-0.332;
  elseif veg(j) == 1 && KC(j) < 80
    Cd(j)  = 43.5*KC(j)^-0.55; 
  else 
     Cd(j) = 0;
  end 
end 

if DRAG == -3.2   % Cd(Ca) for Ca > 10, Cd(KC) for Ca<10
  if veg(j) == 1 && Ca(j) > 20
    Cd(j)  = 6.88*Ca(j)^-0.28;
  elseif veg(j) == 1 && Ca(j)<=20
    Cd(j)  = 43.5*KC(j)^-0.55; 
  else 
     Cd(j) = 0;
  end 
end 

if DRAG == -3.3 % Cd(Ca), with E input varying with emergence (ie. if Er > 1), 
% ... didn't really work... 
end 

if DRAG == -3.4 % Cd(Ca) with fitted a, b, & E (=2.7067e7 Pa)
    if veg(j) == 1 
%         Cd(j)  =6*Ca(j)^-0.3;
        Cd(j)  =8.916*Ca(j)^-0.2587;
    else 
        Cd(j) = 0;
    end 
end 

if DRAG == -3.5 % Cd(Ca*Nv) with fitted a, b, & E (=3.351e8 Pa)
    if veg(j) == 1 
%         Cd(j)  =6*Ca(j)^-0.3;
%         Cd(j)  =8.916*(Ca(j).*Nv(j))^-0.2587;
%          Cd(j)  =20.2*(Ca(j).*Nv(j))^-0.24;
%     Cd(j)  =41.79*(Ca(j).*Nv(j))^-0.2587; %   E (=2.7067e7 Pa)
    Cd(j)  =40*(Ca(j).*Nv(j))^-0.4; 
    else 
        Cd(j) = 0;
    end 
end 

if DRAG == -4   % Cd(KC,Er) = (DF18 fit) 
     if veg(j) == 1 %&& KC(j) > 20
     Cd(j) = -1.94+ 1.34*Er(j) + 43.5*KC(j)^-0.55;  
     else 
     Cd(j) = 0;
     end 
 end 

if DRAG == -5   % Cd(Re,Er) = DF18 fit 
  if veg(j) == 1 && KC(j) > 20
     Cd(j) = -2.19+ 1.42*Er(j) + 149.1*Re(j)^-0.55;
  else 
     Cd(j) = 0;
  end 
end 

if DRAG == -6  % Cd(Ca,Er) = DF18 fit 
  if veg(j) == 1 && KC(j) > 20
     Cd(j) = -0.01 - 0.004*Er(j) + 8.3*Ca(j)^-0.33;
  else 
     Cd(j) = 0;
  end 
end 

if DRAG == -7   % Cd(psiKC) = (DF18 fit)              % Appraoch 3a
     if veg(j) == 1 && KC(j) > 20
     Cd(j) = (22 + 10*Er(j))*KC(j)^-0.53; 
     else 
     Cd(j) = 0;
     end 
 end 
if DRAG == -8   % Cd(psiRe) = DF18 fit                % Approach 3b
  if veg(j) == 1 && KC(j) > 20
     Cd(j) = (113*Er(j)+ 343)*Re(j)^-0.77; 
  else 
     Cd(j) = 0;
  end 
end 
if DRAG == -9  % Cd(psiCa) = DF18 fit                  % Approach 3c
  if veg(j) == 1 && KC(j) > 20
      Cd(j) = (11.7*Er(j)-1.95)*Ca(j)^-0.35; 
  else 
     Cd(j) = 0;
  end 
end 


if  DRAG == -10  % Cd(Re) from AS2014
    if veg(j) == 1
     UCoeff = cosh(k(j)* a(j)*h(j))/sinh(k(j)*h(j)); % z here = -h + alpha*h
     Uc =(pi*H(j)/T) * UCoeff;
     KC(j) = (Uc*T)/bv(j);
     Re(j) = (Uc*bv(j))/1e-6;
     Cd(j) = 0.76 + (744.2/Re(j))^1.27;  % **need to add limits on KC
    else 
     Cd(j) = 0;
    end 
end 
     
if DRAG == -11   % Cd(KC) from AS2014
     if veg(j) == 1
     UCoeff = cosh(k(j)* a(j)*h(j))/sinh(k(j)*h(j));
     Uc =(pi*H(j)/T) * UCoeff;
     KC(j) = (Uc*T)/bv(j);
     Re(j) = (Uc*bv(j))/1e-6;
     Cd(j) = 1.1 + (27.4/KC(j))^3.08;  % **need to add limits on KC
    else 
     Cd(j) = 0;
     end 
end 
     
if  DRAG == -12   % Cd(KC) from Jadhav et al 2013 (need to specify STEM ht as ls to use this) 
     if veg(j) == 1
     UCoeff = 1/(sinh(k(j)*h(j)));
     Ub =(pi*H(j)/T) * UCoeff;  % KC in Jetal2013 uses Ub
     KC(j) = (Ub*T)/bv(j);
       if KC(j)<135 && KC(j)>25
         Cd(j) =70*(KC(j)^-0.86); 
        else
         Cd(j) = NaN;
       end 
     else 
     Cd(j) = 0;
     KC(j) = NaN;
     Re(j) = NaN;
     end 
end 



% compute period-averaged rates of dissipation for veg & bottom friction:
Dv(j) = (V(j)*Cd(j)*(H(j))^3);
Df(j) = (F(j)*H(j)^3);
 
% compute period-averaged rate of dissipation for breaking:
if BRK == 1 
    Db(j) = (R(j)*(H(j)^7)); 
elseif BRK == 2   
    M = ( H(j)/(gamma(j)*h(j)) )^2 ; 
    bracket = 1 - (1 / ((1 + M)^(5/2))); 
    Db(j) = (R(j)*(H(j)^5)*bracket); %% check this 
elseif BRK == 3 
    Db(j) = (R(j)*(H(j)^7)); 
elseif BRK == 4 
    M = ( H(j)/(gammac*h(j)) )^2 ; 
    bracket = 1 - (1 / ((1 + M)^(5/2))); 
    Db(j) = (R(j)*(H(j)^3/h(j))*M*bracket); 
elseif BRK == 0
    Db(j) = 0;
end

%compute wave height at next grid node: 
H2(j+1) = (H(j)^2*(cg(j)/cg(j+1))) + ...
        dx*( (Dv(j) ) + ( Db(j) ) + ( Df(j) ) )...
        /(0.125*rho*g*cg(j+1));
H(j+1) = ( sqrt(H2(j+1)));
% H2(j+1) = (H(j)^2) + ...
%     dx*( (Dv(j) ) + ( Db(j) ) + ( Df(j) ) )...
%         /(0.125*rho*g*cg(j));
%     
% H(j+1) = ( sqrt(H2(j+1)));
%    
end  % end solver loop 
end  % end function 

%% Solve wave dispersion function
function [L, k] = wavenumber(omega,h)
%k = wavenumber(omega,h)
%k is a vector of same size as omega and depth containing the calculated wave numbers
%omega is the wave frequencies in rad/s
%h is the water depth
%sigma and h must be scalars or vectors of the same dimensions
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
end 