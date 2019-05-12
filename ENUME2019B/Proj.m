clear all
close all

% initial setup
step=100;
approx = zeros(1,step);
x = linspace(-1,1,step);
y = zeros(1,10);
x_n = zeros(1,10);
y_t = zeros(1,10);

% approximation loop
iter=1;

for N=10:10:30
    y = zeros(1,N);
    x_n = zeros(1,N);
    y_t = zeros(1,N);
    index=N/10;
    for i=8:9
        K=N*0.1*i;
        K=round(K);
        
        % performing approximation
        for s=1:step
            approx(s)=approximate(x(s),N,K);
        end
        figure(iter)
        iter=iter+1;
        hold on;
        
        % original function
        plot(x, f(x), 'k');
        
        % nodes
        x_n=generate_xn(N);
        y=generate_y(x_n);
        plot(x_n, y,'ko')
        
        % approximated function
        plot(x, approx,'m')
        title(strcat('N=', num2str(N),', K=', num2str(K)));
        legend('function','nodes', 'approximation');
        hold off
    end
end

%Task3
[N,K]=meshgrid(5:50,5:50);
for n=5:50
    for k=5:50
        if(k<n)
            rms(n-4,k-4)=RMS(1,n,k);
            me(n-4,k-4)=ME(1,n,k);
        else
            rms(n-4,k-4)=NaN;
            me(n-4,k-4)=NaN;
        end
    end
end
figure(6)
surf(K, N, rms);
figure(7)
surf(K, N, me);

% Root-mean-square error
function y = RMS(x, N, K)
    y = norm(approximate(x, N, K) - f(x));
    y = y / norm(f(x));
end

% Maximum error
function y = ME(x, N, K)
    y = norm(approximate(x, N, K) - f(x), Inf);
    y = y / norm(f(x), Inf);
end

function y = generate_fi(N,K)
    Fi = zeros(N,K);
    for n=1:N
        x_n = -1+2*(n-1)/(N-1);
        for k=1:K
            Fi(n,k )= Bsk(x_n,k,K);
        end
    end
    y = Fi;
end

function y = generate_y(x)
    N = length(x);
    yn = zeros(1,N);
     for n=1:N
       x_n = -1+2*(n-1)/(N-1);
       yn(n) = f(x_n);
     end
    y = yn;
end

function y = generate_xn(N)
    x_n=zeros(1,N);
    for n=1:N
       x_n(n)=-1+2*(n-1)/(N-1);
    end
    y=x_n;
end

function y = Bsk(x,k,K)
    xk = -1+2*((k-1)/(K-1));
    y = Bs(2*(x-xk)+2);
end

% approximating function
function y = approximate(x,N,K)
    Fi = generate_fi(N,K);
    x_n = generate_xn(N);
    y = generate_y(x_n);
    p = Fi.'*Fi\Fi.'*y.';
    y=0;
    K=length(p);
    for i=1:K
        y=y+p(i)*Bsk(x,i,K);
    end
end

function x = x_k(k,K)
    x=2*((k-1)/(K-1))-1;
end

% initial function
function y = f(x)
    y = (x+(1/3)).^2 + exp(-x-2);
end

% B_spline functions definition
function y=Bs(x)
    % x?[0, 1)
    if (x>=0 && x<1)
        y = x^3;
    % x?[1, 2)
    elseif (x>=1 && x<2)
        y = -3*((x-1)^3) + 3*((x-1)^2) + 3*(x-1) + 1;
    % x?[2, 3)
    elseif (x>=2 && x<3)
        y = 3*((x-2)^3) - 6*((x-2)^2) + 4;
    % x?[3, 4]
    elseif (x>=3 && x<=4 )
        y = -((x-3)^3) + 3*((x-3)^2) - 3*(x-3) + 1;
    % x?[0, 4]
    else
        y=0;
    end
end


