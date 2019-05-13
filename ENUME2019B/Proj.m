clear
close all

% Tasks 1 and 2 - approximation of function

% initial setup
step=100;
approx = zeros(1,step);
x = linspace(-1,1,step);

i = 1;
for N=10:10:30
    for K=10:10:30
        y_t = zeros(1,N);
        index=N/10;

        % performing approximation
        for s=1:step
            approx(s)=approximate(x(s),N,K);
        end
        figure(i)
        i = i+1;
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

% Task 3 - Systematic investigation of the dependence of the accuracy
% of approximation on the values of N and K

[N,K]=meshgrid(5:50,5:50);
rms = nan(length(N));
me = nan(length(N));

for k=4:50
    % Satisfying the given K < N condition
    for n=k+1:50
        rms(n-4,k-3)= RMS(n,k);
        me(n-4,k-3)= ME(n,k);
    end
end
figure(10)
surf(K, N, rms);
xlabel("K")
ylabel("N")
zlabel("RMS")
figure(11)
surf(K, N, me);
xlabel("K")
ylabel("N")
zlabel("ME")


% Task 4 - systematic investigation of the dependence of norms on the
% standard deviation of random errors
spc = 10;
o = logspace(-5, -1, spc);
polyo = logspace(-5, -1);
RMSEc = zeros(50,49);
MEc = zeros(50,49);
RMSEm = zeros(1,spc);
MEm = zeros(1, spc);
pts = 30;
    
for s =1:spc  
for N = 10:50
    for K = 4:N-1
        xns = generate_xn(N);
        yns = f(xns);
        xks = generate_xn(K);
        errand = randn(1,N) * o(s);
        yc = yns + errand;
        FI = zeros(N, K);
for n = 1:N
    for k = 1:K
        FI(n,k)=Bsk(xns(n), k, K);
    end
end
p = FI.' * FI \ FI.' * yc';
p = p.';
ys = zeros(1,pts);

for i = 1:pts
    fun = 0;
    for j = 1:K
        fun = fun + p(j)*Bsk(x(i), j, K);
    end
    ys(i)=fun;
end
RMSEc(n,k) = norm(ys-y) / norm(y);

    end
end
RMSEm(1,s) = min(RMSEc(RMSEc>0));

end


p = polyfit(o, RMSEm, 3);
RMSEmapp = polyval(p, polyo);
figure
loglog(o, RMSEm, "ko");
hold on
loglog(polyo, RMSEmapp, 'm');
hold off
des = strcat("RMS on sigma random");
title(des)
legend("RMS nodes", "approximation");

figure
p = polyfit(o, MEm, 3);
MEmapp = polyval(p, polyo);
loglog(o, MEm, "ko");
hold on
loglog(polyo, MEmapp, 'm');
hold off
des = strcat("ME on sigma random");
title(des)
legend("ME nodes", "approximation");



% ----------- FUNCTION DEFINITIONS -------------

% Root-mean-square error
function y = RMS(N, K)
    nom = zeros(1, N);
    denom = zeros(1, N);
    x_n = generate_xn(N);
    
    for i=1:N
        nom(1, i) = approximate(x_n(i), N, K)-f(x_n(i));
        denom(1, i) = f(x_n(i));
    end
    
    y=norm(nom) / norm(denom);
end

% Maximum error
function y = ME(N, K)
    nom = zeros(1, N);
    denom = zeros(1, N);
    x_n = generate_xn(N);
    
    for i=1:N
        nom(1, i) = approximate(x_n(i), N, K) - f(x_n(i));
        denom(1, i) = f(x_n(i));
    end
    
    y= norm(nom, Inf) / norm(denom, Inf);
end

function y = generate_phi(N,K)
    phi = zeros(N,K);
    for n=1:N
        x_n = -1+2*(n-1)/(N-1);
        for k=1:K
            phi(n, k) = Bsk(x_n,k,K);
        end
    end
    y = phi;
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
    Fi = generate_phi(N,K);
    x_n = generate_xn(N);
    y = generate_y(x_n);
    p = Fi.' * Fi \ Fi.' * y.';
    y=0;
    K=length(p);
    for i=1:K
        y=y+p(i)*Bsk(x,i,K);
    end
end

% initial function
function y = f(x)
    y = (x+(1/3)).^2 + exp(-x-2);
end

% B_spline functions definition
function y=Bs(x)
    % x [0, 1)
    if (x>=0 && x<1)
        y = x^3;
    % x [1, 2)
    elseif (x>=1 && x<2)
        y = -3*((x-1)^3) + 3*((x-1)^2) + 3*(x-1) + 1;
    % x [2, 3)
    elseif (x>=2 && x<3)
        y = 3*((x-2)^3) - 6*((x-2)^2) + 4;
    % x [3, 4]
    elseif (x>=3 && x<=4 )
        y = -((x-3)^3) + 3*((x-3)^2) - 3*(x-3) + 1;
    % x doesn't belong to [0, 4]
    else
        y=0;
    end
end



