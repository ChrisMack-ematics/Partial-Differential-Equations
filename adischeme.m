close all;
clear all;
clc;

H = [1/10, 1/20, 1/40];

Table = [H.',zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1)];

for p = 1:length(H)
    h = H(p);

    %Space Discreditization, counters
    xstep = h;
    ystep = h;
    x = 0:h:1;
    y = 0:ystep:1;
    xcounter = 1;
    ycounter = 1;
    M = length(x);
    N = length(y);
    
    %Call exact
    ue =@(x,y,t) exp(1.68.*t).*sin(1.2.*(x-y)).*cosh(x+2.*y);
    
    %Exact Soln
    [XX,YY] = meshgrid(x,y); 
    UE =@(TT) exp(1.68.*TT).*sin(1.2.*(XX-YY)).*cosh(XX+2.*YY);
    
    % Tilde Matrices
    A = zeros(M,N);
    u1 = A;
    ff = A;
    f2 = A;
    v2 = A;
    
    %Time Discretization, counters
    t = 0:h:1;
    tcount = 1;
    tc = 1;
    L = length(t);

    b1 = 2;
    b2 = 1;
    mux = 1/h;
    muy = 1/h;    

    Tria = zeros(N);
    Trib = Tria;
    
    for Y = 1:N
          if Y == 1 
              Trib(1,1) = 1;
          elseif Y < N
              Trib(Y,Y+1) = -b2*muy/2;
              Trib(Y,Y-1) = -b2*muy/2; 
              Trib(Y,Y) = 1+b2*muy;
          else
              Trib(end,end) = 1;                     
          end
    end
    
    for Y = 1:N
          if Y == 1 
              Tria(1,1) = 1;
          elseif Y < N
              Tria(Y,Y+1) = -b1*mux/2;
              Tria(Y,Y-1) = -b1*mux/2;
              Tria(Y,Y) = 1+b1*mux; 
          else 
              Tria(end,end) = 1;     
          end
    end  
    
    Triainv = inv(Tria);
    Tribinv = inv(Trib);    
    
uexact_break = zeros(M,N,L);

    for T = t
       for Y = y
          for X = x             
              uexact_break(xcounter,ycounter,tcount) = exp(1.68*(T+h/2))*sin(1.2*(X-Y))*cosh(X+2*Y);
              uexact_p(xcounter,ycounter,tcount) = exp(1.68*T)*sin(1.2*(X-Y))*cosh(X+2*Y);
              xcounter = xcounter + 1;
          end
          xcounter = 1;
          ycounter = ycounter + 1;
       end
       ycounter = 1;
       tcount = tcount + 1;
    end

    % Time loop
    for n = 1:L-1
    tc = n;
    if n == 1
        
        u1(:,:) = uexact_p(:,:,n); %First Time Step

    elseif n > 1              
        for m = 2:N-1                   
            %Boundary Conditions
            ff(1,m) = ue(x(1), y(m), (t(n)+h/2));
            ff(end,m) = ue(x(end), y(m), (t(n)+h/2)); 
            
            for l = 2:(M-1)              
                ff(l,m) = (muy/2)*u1(m-1,l) + (1-muy)*u1(m,l) + (muy/2)*u1(m+1,l);                   
            end
            
            F = ff(:,m);
            
            v2(:,m) = Triainv*F;          
        end     
        for m = 1:N %new time step
            u1(m,end) = ue(x(end), y(m), t(n)+h);
            u1(m,1) = ue(x(1), y(m), t(n)+h);            
        end
        for l = 1:M
            u1(1,l) = ue(x(l), y(1),t(n)+h);
            u1(end,l) = ue(x(l), y(end), t(n)+h);
        end               
        for l = 2:M-1     
            
            f2(1,l) = ue(x(l), y(1), (t(n)+h));
            f2(end,l) = ue(x(l), y(end), (t(n)+h));
            
            for m = 2:N-1
                f2(m,l) = (mux).*v2(l-1,m) + (1-2*mux).*v2(l,m) + (mux).*v2(l+1,m);
            end
            FF = f2(:,l);
            u1(:,l) = Tribinv*FF;
        end
    end
    end
    
uexact_max = max(max(max(uexact_p)));
uexact_min = min(min(min(uexact_p)));

figure(p)
subplot(1,2,1)
surf(XX,YY,UE(t(n)+h))                
title("Exact at $t=1$ for $h=$ "+h,'interpreter','latex')
xlabel('x')
ylabel('y')
zlabel('u')
zlim([uexact_min,uexact_max])

subplot(1,2,2)
surf(XX,YY,u1)
title("Approximate at $t=1$ for $h =$ "+h,'interpreter','latex')
xlabel('x')
ylabel('y')
zlabel('u')
zlim([uexact_min,uexact_max])
pause(.05)

u1 = u1';
    
%Errors, Accuracy
sum = 0;
for Y = 1:length(y)
    for X = 1:length(x)
        sum = sum + abs(u1(X,Y)-uexact_p(X,Y,end))^2;
    end
end

Table(p,2) = H(p)*sqrt(sum);
Table(p,3) = max(max(u1(:,:)-uexact_p(:,:,end)));

if p > 1
    Table(p,4) = abs(log(Table(p,2)/Table(p-1,2)))/log(2);
    Table(p,5) = abs(log(Table(p,2)/Table(p-1,2)))/log(2);
end
end
Table