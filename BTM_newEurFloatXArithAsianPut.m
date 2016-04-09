%Group 19, A0098071 Khuong Bich Ngoc and A0098100 Le Hoang Van
%BTM for newly issue floating strike Asian arithmetic-average put 
%Call syntax OptVal=BTM_newEurFloatXArithAsianPut( S, r, q, T, sigma, N )

function val = BTM_newEurFloatXArithAsianPut( S, r, q, T, sigma, N )

dt=T/N;
dx=sigma*sqrt(dt);
u=exp(dx);
d=1/u;
df=exp(-r*dt);
p=(exp((r-q)*dt)-d)/(u-d);


% cell array containing all possible average values

A=cell(1,N+1);
A{1}=S;

%j is the index mapping funtion tau(k) in lectured
j=(sum((dec2bin((1:2^N)-1)-48),2)+1)';

for n=1:N;
    k=1:2^n;
    l=ceil(k/2);   
    ST=S*u.^(2*j(k)-n-2);
    A{n+1}(k)=(n*A{n}(l)+ST)/(n+1);
end

%values at maturity

A{N+1}=max(A{N+1}-ST,0);

%backward iterations

for n=N:-1:1
    k=1:2^(n-1);
    A{n}(k)=df*((1-p)*A{n+1}(2*k-1)+p*A{n+1}(2*k));
end

val=A{1};

end

