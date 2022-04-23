clear;
close all;
clc;

H = [1/10,1/20,1/40];
lambda = 1;

Table = [H',zeros(3,1),zeros(3,1),zeros(3,1)];

for p = 1:length(H)
    h = H(p);
    k = lambda*h;
    t = 0:k:1;
    x = -1:h:1;
    M = length(x);
    N = length(t);
    a = lambda/4;

    %Initial condition
   U = sin(pi.*x);
%   U = U.';
   u_old = U;
   u_new = zeros(M,N);
   u_new(:,1) = U;
   
   % Tridiagonal system Setup
   
   % Build LHS A Matrix
        e = ones(M,1);
        A = spdiags([(-a)*e e (a)*e], -1:1, M,M);
        A(1,2) = 0;
        A(M,M-1) = -lambda;
        A(M,M) = 1 + lambda + k;
        Aog = A;
        
   % Time Loop : Need to setup b for each time step
   for i = 2:N
       
       % Time iteration for TriDiagonal system
       
       % Build RHS b matrix
       
       b = zeros(M,1);
       for bb = 1:M-1
           if bb == 1
               b(bb) = exp(-t(i)).*sin(pi.*(-1 - t(i))); %Left Boundary Condish
           else
               b(bb) = -a*u_old(bb+1) + a*u_old(bb-1)...
                   + u_old(bb) - k*u_old(bb); %Crank Nicolson
            end
               b(end) = u_old(end); %Right boundary condish
       end
       bog = b;
       
       % Thomas Algorithm Forward Sweep
       for X = 1:M
           if X == 1      
               g = A(X,1);
               A(X,1) = 1;
               A(X,2) = A(X,2)/g;
               b(X) = b(X)/g;
           elseif X < M & X~=1
               A(X,X+1) = A(X,X+1)/(A(X,X)-A(X,X-1)*A(X-1,X)); %gamma
               b(X) = (b(X)-A(X,X-1)*b(X-1))/(A(X,X)-A(X,X-1)*A(X-1,X)); %rho
               A(X,X-1) = 0;
               A(X,X) = 1; %middle row
           else
               b(X) = (b(X) - A(X,X-1)*b(X-1))/(A(X,X)-A(X,X-1)*A(X-1,X));
               A(X,X) = 1;
               A(X,X-1) = 0;
           end
       end

       % Thomas Algo Backward Sweep
       for P = flip(1:M)
           if P == M
               u_new(P,i) = b(P);
           else
               u_new(P,i) = b(P) - A(P,P+1)*u_new(P+1,i);
           end
       end
      
        %Plot u_new movie
%         figure(1)
%         plot(x,u_new(:,i),'-^')
%         pause(0.2)

        uexact = exp(-t(i)).*sin(pi.*(x - t(i)));
        %Errors
        if i == N
            nrm = norm(uexact.' - u_new(:,i),2);
            error = (h*nrm)^(1/2);
            Table(p,2) = error;
            infnrm = norm(uexact.' - u_new(:,i),Inf);
            inferror = infnrm;
            Table(p,3) = inferror;
        end
        %Orders of Accuracy
        if p > 1
            r = log(Table(p-1,2)/Table(p,2))/log(2);
            Table(p,4) = r;
        end
      
        A = Aog;
        b = bog;
        u_old = u_new(:,i);

   end
%        uplot = u_new(:,N);
%        %Plot Exact with Approximate
%        ue = exp(-1).*sin(pi.*(x-1));
%        figure(p)
%        plot(x,ue,x,uplot,'o');
%        title("3.5.2(a) : Exact and approximate at $t=1$ for $h=$ "+h,'interpreter','latex')
%        xlabel('x')
%        ylabel('u(x)')
%        legend('Exact','Approximate')
% %      title("Approximate for $h=$"+h,'interpreter','latex')

end
Table