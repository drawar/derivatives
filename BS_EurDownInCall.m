%Group 19 , A0098071 Khuong Bich Ngoc and A0098100 Le Hoang Van
% Black Scholes program for European down and in call  options
% call syntax: OptVal=BS_EurDownInCall(S0,X,r,T,H,sigma,q)

function value= BS_EurDownInCall( S0,X,r,T,H,sigma,q)

lambda = (r - q + sigma^2/2)/(sigma^2);
x1 = log(S0/H)/(sigma*sqrt(T)) + lambda*sigma*sqrt(T);
y1 = log(H./S0)/(sigma*sqrt(T)) + lambda*sigma*sqrt(T);

C1 = S0.*normcdf(x1)*exp(-q*T);
C2 = X*exp(-r*T)*normcdf(x1 - sigma*sqrt(T));
C3 = S0.*exp(-q*T).*(H./S0).^(2*lambda).*normcdf(y1);
C4 = X*exp(-r*T)*(H./S0).^(2*lambda - 2).*normcdf(y1 - sigma*sqrt(T));

value = BS_EurVanillaCall(S0,X,r,T,sigma,q) - C1 + C2 + C3 - C4;

end