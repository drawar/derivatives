%Group 19, A0098071 Khuong Bich Ngoc and A0098100 Le Hoang Van
%BTM for not newly issue floating strike European lookback call
%Call syntax OptVal=BTM_runningEurFloatXLookbackCall( S, r, q, T, sigma, N, m)

function val = BTM_runningEurFloatXLookbackCall( S, r, q, T, sigma, N, m)

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

%In case m/S is not on the grid point, find the nearest gridpoint

k=log(m/S)/log(d)+2;
k_f=floor(k);

%backward iteration

for n=N:-1:1
    j=2:n+1;
    W(j)=df*(p*u*W(j+1)+(1-p)*d*W(j-1));
    W(1)=W(2);
end

val=S*( (d^k-d^k_f)*W(k_f+1) + (d^(k_f+1)-d^k)*W(k_f))/(d^(k_f+1)-d^k_f);

end

