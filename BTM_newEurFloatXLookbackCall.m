%Group 19, A0098071 Khuong Bich Ngoc and A0098100 Le Hoang Van
%BTM for newly issue floating strike European lookback call
%Call syntax OptVal=BTM_newEurFloatXLookbackCall( S0, r, q, T, sigma, N )

function val = BTM_newEurFloatXLookbackCall( S, r, q, T, sigma, N )

dt=T/N;
dx=sigma*sqrt(dt);
u=exp(dx);
d=1/u;
df=exp(-r*dt);
p=(exp((r-q)*dt)-d)/(u-d);

%initialization 

j=0:N;
W=1-d.^j;
W=[W(1),W];

%backward iteration

for n=N:-1:1
    j=2:n+1;
    W(j) = df*(p*u*W(j+1)+(1-p)*d*W(j-1));
    W(1)=W(2);
end

val=S*W(1);

end

