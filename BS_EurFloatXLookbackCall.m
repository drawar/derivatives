%Group 19, A0098071 Khuong Bich Ngoc and A0098100 Le Hoang Van
% Black-Scholes formulae for floating strike European lookback call
% call syntax: c = BS_FloatLookbackCall(S, m, r, T, sigma, q)

function val = BS_FloatLookbackCall( S, m, r, q, T, sigma )

d1 = (log(S/m)+(r-q+sigma^2/2)*T)/(sigma*sqrt(T));
d2 = d1 - sigma*sqrt(T);
d3 = -d1 + 2*(r-q)*sqrt(T)/sigma;

val = S*exp(-q*T)*normcdf(d1) - m*exp(-r*T)*normcdf(d2)...
    + S*exp(-r*T)*sigma^2/(2*(r-q))*(S/m)^(-2*(r-q)/sigma^2)*normcdf(d3)...
    - S*exp(-q*T)*sigma^2/(2*(r-q))*normcdf(-d1);
end

