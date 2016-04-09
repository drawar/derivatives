%Group 19, A0098071 Khuong Bich Ngoc and A0098100 Le Hoang Van
%BTM for not newly issue floating strike Asian arithmetic-average put 
%Call syntax OptVal=FSGM_runningAsianGeomFixedXPut( S, r, q, T, sigma, N, rho, X, Nr, a, i)
%Nr time steps of the running average
%a is current geometric average
%i is the type of interpolation  
%   i=1 for nearest point interpolation
%   i=2 for linear interpolation
%   i=3 for quadratic interpolation

function val = FSGM_runningAsianGeomFixedXPut( S, r, q, T, sigma, N , rho, X, Nr, a, i)

dt=T/N;
dx=sigma*sqrt(dt);
u=exp(dx);
d=1/u;
df=exp(-r*dt);
p=(exp((r-q)*dt)-d)/(u-d);
dY=rho*dx;
m=1/rho;

% Initialize, finding a k such that a falls to the normal grid built
% We would have  S*exp(k_af*dY)< a < S*exp((k_af+1)*dY)

k_a=log(a/S)/dY;
k_af=floor(k_a);

%Build the first tree and grid with S0=S, and A0=S*exp((k_af)*dY);

A=S*exp(k_af*dY)*exp(dY*(-N*m:N*m));
V=repmat(A',1,2*N+1);
V(:,1:2:2*N+1)=max(-V(:,1:2:2*N+1)+X,0);

%backward recursive fsg
for n=N-1:-1:0
    k=(N-n)*m+1:(N+n)*m+1;
    for j=-n+N+1:2:n+N+1
        
        %incremental formula
        A_u=exp((log(S*u^(j-N))+((n+Nr+1)*log(V(k,j)))) /(n+Nr+2));
        A_d=exp((log(S*u^(j-N-2))+((n+Nr+1)*log(V(k,j)))) /(n+Nr+2));
        
        %calculate the nearest grid points for interpolation
        k_u=floor(log(A_u/S)/dY)+N*m+1-k_af;
        k_d=floor(log(A_d/S)/dY)+N*m+1-k_af;

        %Different interpolation based on i
        if (i==1) 
        
           V_u=V(k_u,j+1).*( (A_u-A(k_u)') <= (0.5*(A(k_u+1)-A(k_u))'))...
              +V(k_u+1,j+1).*( (A_u-A(k_u)') > (0.5*(A(k_u+1)-A(k_u))'));
          
           V_d=V(k_d,j-1).*( (A_d-A(k_d)') <= (0.5*(A(k_d+1)-A(k_d))'))...
              +V(k_d+1,j-1).*( (A_d-A(k_d)') > (0.5*(A(k_d+1)-A(k_d))'));
          
        elseif (i==2)
            
            V_u=(V(k_u+1,j+1).*(A_u-(A(k_u))')...
                 +V(k_u,j+1).*((A(k_u+1))'-A_u))./(A(k_u+1)-A(k_u))';
            V_d=(V(k_d+1,j-1).*(A_d-(A(k_d))')...
                 +V(k_d,j-1).*((A(k_d+1))'-A_d))./(A(k_d+1)-A(k_d))';
             
        elseif (i==3)
            
            V_u= (A(k_u)'-A_u).*(A(k_u+1)'-A_u)./((A(k_u)-A(k_u-1)).*(A(k_u+1)-A(k_u-1)))'.* V(k_u-1,j+1)...
                +(A(k_u-1)'-A_u).*(A(k_u+1)'-A_u)./((A(k_u-1)-A(k_u)).*(A(k_u+1)-A(k_u)))'.* V(k_u,j+1)...
                +(A(k_u-1)'-A_u).*(A(k_u)'-A_u)./((A(k_u-1)-A(k_u+1)).*(A(k_u)-A(k_u+1)))'.* V(k_u+1,j+1);
            
            V_d= (A(k_d)'-A_d).*(A(k_d+1)'-A_d)./((A(k_d)-A(k_d-1)).*(A(k_d+1)-A(k_d-1)))'.* V(k_d-1,j-1)...
                +(A(k_d-1)'-A_d).*(A(k_d+1)'-A_d)./((A(k_d-1)-A(k_d)).*(A(k_d+1)-A(k_d)))'.* V(k_d,j-1)...
                +(A(k_d-1)'-A_d).*(A(k_d)'-A_d)./((A(k_d-1)-A(k_d+1)).*(A(k_d)-A(k_d+1)))'.* V(k_d+1,j-1);
        end    
        V(k,j)=max(df*(p*V_u+(1-p)*V_d),X-A(k)')  ;
    end   
    V(k,-n+N+2:2:n+N)=repmat(A(k)',1,n);
end


%Build a second tree and grid with S0=S, and A0=S*exp((k_af+1)*dY);

A2=S*exp((k_af+1)*dY)*exp(dY*(-N*m:N*m));
V2=repmat(A2',1,2*N+1);
V2(:,1:2:2*N+1)=max(-V2(:,1:2:2*N+1)+X,0);

%backward recursive fsg
for n=N-1:-1:0
    k=(N-n)*m+1:(N+n)*m+1;
    for j=-n+N+1:2:n+N+1
        
        A_u=exp((log(S*u^(j-N))+((n+Nr+1)*log(V2(k,j)))) /(n+Nr+2));
        A_d=exp((log(S*u^(j-N-2))+((n+Nr+1)*log(V2(k,j)))) /(n+Nr+2));

        k_u=floor(log(A_u/S)/dY)+N*m+1-(k_af+1);
        k_d=floor(log(A_d/S)/dY)+N*m+1-(k_af+1);

        %Different interpolation based on i
        if (i==1) 
        
           V_u=V2(k_u,j+1).*( (A_u-A2(k_u)') <= (0.5*(A2(k_u+1)-A2(k_u))'))...
              +V2(k_u+1,j+1).*( (A_u-A2(k_u)') > (0.5*(A2(k_u+1)-A2(k_u))'));
          
           V_d=V2(k_d,j-1).*( (A_d-A2(k_d)') <= (0.5*(A2(k_d+1)-A2(k_d))'))...
              +V2(k_d+1,j-1).*( (A_d-A2(k_d)') > (0.5*(A2(k_d+1)-A2(k_d))'));
          
        elseif (i==2)
            
            V_u=(V2(k_u+1,j+1).*(A_u-(A2(k_u))')...
                 +V2(k_u,j+1).*((A2(k_u+1))'-A_u))./(A2(k_u+1)-A2(k_u))';
            V_d=(V2(k_d+1,j-1).*(A_d-(A2(k_d))')...
                 +V2(k_d,j-1).*((A2(k_d+1))'-A_d))./(A2(k_d+1)-A2(k_d))';
             
        elseif (i==3)
            
            V_u= (A2(k_u)'-A_u).*(A2(k_u+1)'-A_u)./((A2(k_u)-A2(k_u-1)).*(A2(k_u+1)-A2(k_u-1)))'.* V2(k_u-1,j+1)...
                +(A2(k_u-1)'-A_u).*(A2(k_u+1)'-A_u)./((A2(k_u-1)-A2(k_u)).*(A2(k_u+1)-A2(k_u)))'.* V2(k_u,j+1)...
                +(A2(k_u-1)'-A_u).*(A2(k_u)'-A_u)./((A2(k_u-1)-A2(k_u+1)).*(A2(k_u)-A2(k_u+1)))'.* V2(k_u+1,j+1);
            
            V_d= (A2(k_d)'-A_d).*(A2(k_d+1)'-A_d)./((A2(k_d)-A2(k_d-1)).*(A2(k_d+1)-A2(k_d-1)))'.* V2(k_d-1,j-1)...
                +(A2(k_d-1)'-A_d).*(A2(k_d+1)'-A_d)./((A2(k_d-1)-A2(k_d)).*(A2(k_d+1)-A2(k_d)))'.* V2(k_d,j-1)...
                +(A2(k_d-1)'-A_d).*(A2(k_d)'-A_d)./((A2(k_d-1)-A2(k_d+1)).*(A2(k_d)-A2(k_d+1)))'.* V2(k_d+1,j-1);
        end    
        V2(k,j)=max(df*(p*V_u+(1-p)*V_d),X-A2(k)')  ;
    end   
    V2(k,-n+N+2:2:n+N)=repmat(A2(k)',1,n);
end


%The final option value is linearly interpolated between the two value
%calculated in V and V2

val=( V2(N*m+1,N+1)*(a-A(N*m+1)) ...
      +V(N*m+1,N+1)*(A2(N*m+1)-a) ) /( A2(N*m+1)-A(N*m+1) );

end