%Group 19 , A0098071 Khuong Bich Ngoc and A0098100 Le Hoang Van
% BTM program for European down and out call options
% call syntax: OptVal=BTM_EurDownOutCall(S0,X,r,T,H,sigma,q,N)

function value = BTM_EurDownOutCall( S, X, r, T, H, sigma, q, N )

dt=T/N;
dx=sigma*sqrt(dt);
u=exp(dx);
d=1/u;
df=exp(-r*dt);
p=(exp((r-q)*dt)-d)/(u-d);

% initialization

j=1:N+1;
S_N = S*u.^(2*j-N-2);
V = (S_N-X).*(S_N>H);

% backward recursive through time
for n=N-1:-1:0
    j=1:(n+1);
    S_n=S*u.^(2*j-n-2);
    V=df*(p*V(j+1)+ (1-p)*V(j)).*(S_n>H);
end

value=V(1);
    
end

