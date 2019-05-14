clear
close all

% Tasks 1 and 2 - approximation of function

% initial setup
nodes_num=100;
approx = zeros(1,nodes_num);
x = linspace(-1,1,nodes_num);

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

% initial setup
N_const = 20;
nodes_num = 10;

n_min_rms = zeros(1,nodes_num);
k_min_rms = zeros(1,nodes_num);

n_min_me = zeros(1, nodes_num);
k_min_me = zeros(1, nodes_num);

[N,K] = meshgrid(1:N_const,1:N_const);
space = logspace(-5,-1,nodes_num);

roms = nan(N_const);
roms_min = ones(1, nodes_num);
mes = nan(N_const);
mes_min = ones(1, nodes_num);

% performing the investigation
for i=1:nodes_num
    for n=5:N_const
        for k=5:N_const
            if(k<n)
                % computing errors on corrupted data
                roms(n-4, k-4) = RMS_corrupt(n, k, space(i));
                mes(n-4, k-4) = ME_corrupt(n, k, space(i));
                % selecting pairs
                if(roms(n-4,k-4)<roms_min(i))   % RMS
                    roms_min(i)=roms(n-4,k-4);
                    n_min_rms(i)=n;
                    k_min_rms(i)=k;
                end
                if(mes(n-4,k-4)<mes_min(i))     % ME
                    mes_min(i)=mes(n-4,k-4);
                    n_min_me(i)=n;
                    k_min_me(i)=k;
                end
                
            else
                continue;
            end
        end
    end
end

% fitting the data
p_rms = polyfit(space, roms_min, 3);    % RMS
RMSE_map = polyval(p_rms, space);

p_me = polyfit(space, mes_min, 3);      % ME
ME_map = polyval(p_me, space);

% plotting
figure % RMS
loglog(space, roms_min, "ko");
hold on
loglog(space, RMSE_map,'m');
legend("RMS nodes", "approximation")
hold off

figure % ME
loglog(space, mes_min, "ko");
hold on
loglog(space, ME_map,'m');
legend("ME nodes", "approximation")
hold off



% ----------- FUNCTION DEFINITIONS -------------


% approximation of corrupted data
function y = approx_corrupt(x, N, K, sigma)
    phi = generate_phi(N,K);
    x_n = generate_xn(N);
    y = generate_y_corrupted(x_n,sigma);
    p = phi.' * phi\phi.' * y.';
    y = 0;
    K = length(p);
    for i=1:K
        y = y + p(i) * Bsk(x,i,K);
    end
end

% root mean square error of corrupted data
function y = RMS_corrupt(N,K,sigma)
    nom = zeros(1,N);
    denom = zeros(1,N);
    x_n=generate_xn(N);
    for i=1:N
        nom(1, i) = approx_corrupt(x_n(i), N, K, sigma) - f(x_n(i));
        denom(1, i) = f(x_n(i));
    end
    y= norm(nom) / norm(denom);
end

% maximum error of corrupted data
function y = ME_corrupt(N,K,sigma)
    nom = zeros(1,N);
    denom = zeros(1,N);
    x_n = generate_xn(N);
    for i=1:N
        nom(1, i) = approx_corrupt(x_n(i), N, K, sigma) - f(x_n(i));
        denom(1,i) = f(x_n(i));
    end
    y= norm(nom,Inf) / norm(denom, Inf);
end

% generation of corrupted y 
function y=generate_y_corrupted(x, sigma)
    N = length(x);
    yn = zeros(1,N);
     for n=1:N
       x_n = -1+ 2*(n-1) / (N-1);
       yn(n) = f(x_n) + randn()*sigma^2;
     end
    y=yn;
end

% Root-mean-square error (2 norm)
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

% Maximum error (infinity norm)
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

% generation of phi for approximation
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

% generate y coordinates of the function
function y = generate_y(x)
    N = length(x);
    yn = zeros(1, N);
     for n=1:N
       x_n = -1+2*(n-1)/(N-1);
       yn(n) = f(x_n);
     end
    y = yn;
end

% generate x coordinates of the function
function y = generate_xn(N)
    x_n = zeros(1,N);
    for n=1:N
       x_n(n) = -1 + 2*(n-1) / (N-1);
    end
    y=x_n;
end

% calculating xk points in given B-splines
function y = Bsk(x, k, K)
    xk = -1 + 2*((k-1) / (K-1));
    y = Bs(2 * (x - xk) + 2);
end

% approximating function
function y = approximate(x, N, K)
    Fi = generate_phi(N, K);
    x_n = generate_xn(N);
    y = generate_y(x_n);
    p = Fi.' * Fi \ Fi.' * y.';
    y=0;
    K=length(p);
    for i=1:K
        y=y + p(i) * Bsk(x, i, K);
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



