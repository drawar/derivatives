%Group 19, A0098071 Khuong Bich Ngoc and A0098100 Le Hoang Van
%BTM for running fixed strike Asian lookback  put 
%Call syntax OptVal=FSGM_runningAmeFixedXLookbackPut(S, r, q, T, sigma, X, m, N)

function val = FSGM_runningAmeFixedXLookbackPut( S, r, q, T, sigma, X, m ,N)

% initialize & setting up parameters

dt=T/N;
dx=sigma*sqrt(dt);
u=exp(dx);
d=1/u;
df=exp(-r*dt);
p=(exp((r-q)*dt)-d)/(u-d);

%grid values
A=S*u.^(-N:N);
A(A>m)=m;

% teminal values
V=repmat(A',1,2*N+1);
V(:,1:2:2*N+1)=max(-V(:,1:2:2*N+1)+X,0);

%backward recursive fsg
   
for n=N-1:-1:0
    for j=-n+N+1:2:n+N+1
        k=-n+N+1:n+N+1;

        k_u=k;
        k_d=min(k,(j-N-2)+N+1);
		    
        V(k,j)=max(df*(p*V(k_u,j+1)+(1-p)*V(k_d,j-1)), X-A(k)');
        
    end    
    V(:,-n+N+2:2:n+N)=repmat(A',1,n);  
end

val = V(N+1,N+1);

end

