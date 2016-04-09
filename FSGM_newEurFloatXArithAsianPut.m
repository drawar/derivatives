%Group 19, A0098071 Khuong Bich Ngoc and A0098100 Le Hoang Van
%BTM for newly issue floating strike Asian arithmetic-average put 
%Call syntax OptVal=FSGM_newEurFloatXArithAsianPut(( S, r, q, T, sigma, rho, N, i)
%i is the type of interpolation  
%   i=1 for nearest point interpolation
%   i=2 for linear interpolation
%   i=3 for quadratic interpolation


function val = FSGM_newEurFloatXArithAsianPut( S, r, q, T, sigma,rho,N, i) 

%Initialize & set up parameters
dt=T/N;
dx=sigma*sqrt(dt);
u=exp(dx);
d=1/u;
df=exp(-r*dt);
p=(exp((r-q)*dt)-d)/(u-d);
dY=rho*dx;
m=1/rho;

%Set up grid 
A=S*exp(dY*(-N*m:N*m));
SN=S*u.^(-N:N);

% terminal values
V=repmat(A',1,2*N+1);
V(:,1:2:2*N+1)=max(V(:,1:2:2*N+1)-repmat(SN(1:2:2*N+1),2*N*m+1,1),0);

%backward recursive fsg
for n=N-1:-1:0
    for j=-n+N+1:2:n+N+1
        k=-n*m+N*m+1:n*m+N*m+1;
        
        A_u=(V(k,j)*(n+1)+S*u^(j-N))/(n+2);
        A_d=(V(k,j)*(n+1)+S*u^(j-N-2))/(n+2);     

        k_u=floor(log(A_u/S)/dY)+N*m+1;
        k_d=floor(log(A_d/S)/dY)+N*m+1;
        
        
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
        
        V(k,j)=df*(p*V_u+(1-p)*V_d);
        
    end
    
    V(:,-n+N+2:2:n+N)=repmat(A',1,n);
    
end

val = V(N*m+1,N+1);

end

