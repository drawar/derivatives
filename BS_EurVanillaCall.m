%Group 19, A0098071 Khuong Bich Ngoc and A0098100 Le Hoang Van
% Black-Scholes formulae for European vanilla call
% call syntax: c = BS_EurVanillaCall( S0, X, r, T, sigma, q)

function [ value ] = BS_EurVanillaCall( S0, X, r, T, sigma, q)

d1 = (log(S0/X) + (r - q + sigma^2/2)*T)/(sigma*sqrt(T));
d2 = d1 - sigma*sqrt(T);

value = S0.*exp(-q*T).*normcdf(d1) - X*exp(-r*T).*normcdf(d2);

end