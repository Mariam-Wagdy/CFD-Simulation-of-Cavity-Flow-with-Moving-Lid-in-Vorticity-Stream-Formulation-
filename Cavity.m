clc;
clear;

K=10^4;



x0= 0;
xf= 1;

% dx= (xf-x0)/(n-1);

dx=0.025;
n=(xf-x0)/dx; %number of nodes

x= x0:dx:xf;

L= length(x);

omg=zeros(n+1,n+1,1);

epsi=ones(n+1,n+1,K);
dumb=ones(n+1,n+1,500);

u=zeros(n+1,n+1,1);
v=u;

U0=1; %%m/s
Re=100; % 1000 10000

row=length(x);

err=100;
dt=0.001;

k=1;
ex=0;
ey=ex;

tic;
while err(k)> 2 %%percent relative error
    %% left & right boundries, 2nd order forward & backward
    for j=2:1:row-1
        omg(1,j,k)= (-epsi(3,j,k)+8*epsi(2,j,k)-7*epsi(1,j,k))/(2*dx^2);
        omg(row,j,k)=(-7*epsi(row-2,j,k)+8*epsi(row-1,j,k)-epsi(row,j,k))/(2*dx^2);
    end
    
    %% upper and lower boundries, 2nd order forward & backward
    for i=1:1:row
        omg(i,1,k)= (-epsi(i,3,k)+8*epsi(i,2,k)-7*epsi(i,1,k))/(2*dx^2);
        omg(i,row,k)= (-7*epsi(i,row-2,k)+8*epsi(i,row-1,k)-epsi(i,row,k))/(2*dx^2)-3*U0/dx;
    end
    
    %% Solve for omega
    
    % near left boundry, 1st order upwind
    for i=2
        % corner////////////////////////
        for j =2:1:row-1
            [ex,ey]=upwind(u(i,j,k),v(i,j,k));
            
            A= (1-ex)/(2*dx)*(u(i+1,j,k)*omg(i+1,j,k)-u(i,j,k)*omg(i,j,k));
            B= (1+ex)/(2*dx)*(u(i,j,k)*omg(i,j,k)-u(i-1,j,k)*omg(i-1,j,k));
            
            C= (1-ey)/(2*dx)*(v(i,j+1,k)*omg(i,j+1,k)-v(i,j,k)*omg(i,j,k));
            D= (1+ey)/(2*dx)*(v(i,j,k)*omg(i,j,k)-v(i,j-1,k)*omg(i,j-1,k));
            
            omg(i,j,k+1)= dt*((-1/2)*(A+B+C+D)+1/(Re*dx^2)*(omg(i+1,j,k)-4*omg(i,j,k)+omg(i-1,j,k)+omg(i,j+1,k)+omg(i,j-1,k)))+omg(i,j,k);
        end
    end
    
    % near right boundry, 1st order upwind
    for i=row-1
        % corner////////////////////////
        for j =2:1:row-1
            [ex,ey]=upwind(u(i,j,k),v(i,j,k));
            
            A= (1-ex)/(2*dx)*(u(i+1,j,k)*omg(i+1,j,k)-u(i,j,k)*omg(i,j,k));
            B= (1+ex)/(2*dx)*(u(i,j,k)*omg(i,j,k)-u(i-1,j,k)*omg(i-1,j,k));
            
            C= (1-ey)/(2*dx)*(v(i,j+1,k)*omg(i,j+1,k)-v(i,j,k)*omg(i,j,k));
            D= (1+ey)/(2*dx)*(v(i,j,k)*omg(i,j,k)-v(i,j-1,k)*omg(i,j-1,k));
            
            omg(i,j,k+1)= dt*((-1/2)*(A+B+C+D)+1/(Re*dx^2)*(omg(i+1,j,k)-4*omg(i,j,k)+omg(i-1,j,k)+omg(i,j+1,k)+omg(i,j-1,k)))+omg(i,j,k);
        end
    end
    
    
        % near top boundry, 1st order upwind
    for j=2
        % corner////////////////////////
        for i =3:1:row-2
            [ex,ey]=upwind(u(i,j,k),v(i,j,k));
            
            A= (1-ex)/(2*dx)*(u(i+1,j,k)*omg(i+1,j,k)-u(i,j,k)*omg(i,j,k));
            B= (1+ex)/(2*dx)*(u(i,j,k)*omg(i,j,k)-u(i-1,j,k)*omg(i-1,j,k));
            
            C= (1-ey)/(2*dx)*(v(i,j+1,k)*omg(i,j+1,k)-v(i,j,k)*omg(i,j,k));
            D= (1+ey)/(2*dx)*(v(i,j,k)*omg(i,j,k)-v(i,j-1,k)*omg(i,j-1,k));
            
            omg(i,j,k+1)= dt*((-1/2)*(A+B+C+D)+1/(Re*dx^2)*(omg(i+1,j,k)-4*omg(i,j,k)+omg(i-1,j,k)+omg(i,j+1,k)+omg(i,j-1,k)))+omg(i,j,k);
        end
    end
    
    % near bottom boundry, 1st order upwind
    for j=row-1
        % corner////////////////////////
        for i =3:1:row-2
            [ex,ey]=upwind(u(i,j,k),v(i,j,k));
            
            A= (1-ex)/(2*dx)*(u(i+1,j,k)*omg(i+1,j,k)-u(i,j,k)*omg(i,j,k));
            B= (1+ex)/(2*dx)*(u(i,j,k)*omg(i,j,k)-u(i-1,j,k)*omg(i-1,j,k));
            
            C= (1-ey)/(2*dx)*(v(i,j+1,k)*omg(i,j+1,k)-v(i,j,k)*omg(i,j,k));
            D= (1+ey)/(2*dx)*(v(i,j,k)*omg(i,j,k)-v(i,j-1,k)*omg(i,j-1,k));
            
            omg(i,j,k+1)= dt*((-1/2)*(A+B+C+D)+1/(Re*dx^2)*(omg(i+1,j,k)-4*omg(i,j,k)+omg(i-1,j,k)+omg(i,j+1,k)+omg(i,j-1,k)))+omg(i,j,k);
        end
    end
    
    
    % interior points, 2nd order upwind
    for i=3:1:row-2
        for j =3:1:row-2
            [ex,ey]=upwind(u(i,j,k),v(i,j,k));
            
            A= (1-ex)/(2*dx)*(-3*u(i,j,k)*omg(i,j,k)+4*u(i+1,j,k)*omg(i+1,j,k)-u(i+2,j,k)*omg(i+2,j,k));
            B= (1+ex)/(2*dx)*(3*u(i,j,k)*omg(i,j,k)-4*u(i-1,j,k)*omg(i-1,j,k)+u(i-2,j,k)*omg(i-2,j,k));
            
            C= (1-ey)/(2*dx)*(-3*v(i,j,k)*omg(i,j,k)+4*v(i,j+1,k)*omg(i,j+1,k)-v(i,j+2,k)*omg(i,j+2,k));
            D= (1+ey)/(2*dx)*(3*v(i,j,k)*omg(i,j,k)-4*v(i,j-1,k)*omg(i,j-1,k)+v(i,j-2,k)*omg(i,j-2,k));
            
            omg(i,j,k+1)= dt*((-1/2)*(A+B+C+D)+1/(Re*dx^2)*(omg(i+1,j,k)-4*omg(i,j,k)+omg(i-1,j,k)+omg(i,j+1,k)+omg(i,j-1,k)))+omg(i,j,k);
        end
    end
    
    
    %% solve for epsi
    err2=100;
    k2=1;
    dumb(:,:,1)=epsi(:,:,k);
    
    while err2>0.1 %%relative error
        for i=2:1:row-1
            for j =2:1:row-1
                dumb(i,j,k2+1)= 1/4*(dumb(i+1,j,k2)+dumb(i-1,j,k2+1)+dumb(i,j+1,k2)+dumb(i,j-1,k2+1)+dx^2*omg(i,j,k+1));
            end
        end
        
        err2 = 100*sum(abs(epsi(:,:,k2+1)-epsi(:,:,k2)),'all')/sum(abs(epsi(:,:,k2+1)),'all');
        k2=k2+1;
        
    end
    
    epsi(:,:,k+1)=dumb(:,:,k2);
   
%     surf(x,x,omg(:,:,k)');
%     tit="vortcity Function epsi. Iteration #"+string(k);
%     title(tit);
%     xlabel('x');
%     ylabel('y');
%     drawnow;
    
     relerr(k+1)= 100*sum(abs(epsi(:,:,k+1)-epsi(:,:,k)),'all')/sum(abs(epsi(:,:,k+1)),'all');
%     plot(relerr);
%     drawnow;
    
    %% BC    
    for j=2:row-1
        u(1,j,k+1)= 0;
        u(row,j,k+1)= 0;
        v(1,j,k+1)= 0;
        v(row,j,k+1)= 0;
   end
   
    for i=1:1: row
        u(i,1,k+1)=0;
        u(i,row,k+1)=U0;
        v(i,1,k+1)=0;
        v(i,row,k+1)=0;
    
    end
    
    %% update interior points
    for i=2:1:row-1
        for j =2:1:row-1
            u(i,j,k+1)=1/(2*dx)*(epsi(i,j+1,k+1)-epsi(i,j-1,k+1));
            v(i,j,k+1)=-1/(2*dx)*(epsi(i+1,j,k+1)-epsi(i-1,j,k+1));
        end
    end
    
    %% convergence criteria
    err(k+1) = sum(abs(omg(:,:,k+1)-omg(:,:,k)),'all')/sum(abs(omg(:,:,k+1)),'all'); %% mean-squared error (MSE)
    
%     plot(err);
%     surf(x,x,epsi(:,:,k+1))

% streamslice(x,x,u(:,:,k)',v(:,:,k)',2)
% xlabel("x");
% ylabel("y");
% 
% drawnow;
    
    k=k+1;
    
end
toc;

figure;
tit= "Stream Vorticity Function Omega for Re="+Re+" dt="+dt; 
surf(x,x,omg(:,:,k)');
title(tit);
xlabel("x");
ylabel("y");

figure;
tit= "Stream Function Epsi for Re="+Re+" dt="+dt; 
surf(x,x,epsi(:,:,k)');
title(tit);
xlabel("x");
ylabel("y");

figure;
axis tight 
tit= "Streamlines for Re="+Re+" dt="+dt; 
streamslice(x,x,u(:,:,k)',v(:,:,k)',2)
title(tit);
xlabel("x");
ylabel("y");

figure;
axis tight 
tit= "Vector Field for Re="+Re+" dt="+dt; 
quiver(x,x,u(:,:,k)',v(:,:,k)')
title(tit);
xlabel("x");
ylabel("y");

figure;
tit= "Omega Error for Re="+Re+" dt="+dt; 
plot(err, '-x')
title(tit);
xlabel("Iteration");
ylabel("Error");

figure;
tit= "Epsi Error for Re="+Re+" dt="+dt; 
plot(relerr, '-x')
title(tit);
xlabel("Iteration");
ylabel("Error");


