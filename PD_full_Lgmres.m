%<<<<<<< Updated upstream
%modificado por jccf julio 2025

%clear all
%clc
%This is a basic gmres program with restart
%Inputs include A,x0,b,m,tol

%function [logres,tiempoC,ciclos]= PD_full_Lgmres(A,b,mPD, alpha, delta,itermax)
%=======
function [logres,tiempoC,ciclos]= PD_full_Lgmres(A,b,mPD, lL, alpha, delta,itermax)
%>>>>>>> Stashed changes

tic;         %Time Control
tol=1e-9;
%load bcsstk13.mat;
%load sherman5.mat;
%load add20.mat;
%load raefsky1.mat;
maxit=itermax;
%maxit=1000;
%A=Problem.A;
%b=Problem.b;
%color='r'; Name_Matrix='Sherman 4';
%b=ones(size(A,1),1);
x0=zeros(size(A,1),1);
m=mPD;
%m=28;
minitial=m;
% d=lL;
 lL=1;
d=lL;
%s = m + d;
n=size (A,2);
flag=0;
%restart=1;
%r=b-A*x0;
%res(1,:)=norm(r);
%logres(1,:)= norm(r0)/res(1,1);
%iter(1,:)=1;
%miteracion(1,1)=m;
%Preallocating for speed
w=zeros(n,m+d);
z=zeros(n,1);
ij=1; %Para matriz z
%%%%%%%
restart=1;
r=b-A*x0;
res(1,:)=norm(r);
%beta0=res(1,:);
logres(1,:)=(norm(r)/res(1,1));
iter(1,:)=restart;
miteracion(1,1)=minitial;
mmin=1;
mmax=n-1; %se puede considerar que no tiene cota superior, antes de las 1000 iteraciones  no logra alcanzar mmax con alpha=-3 y delta=5
mstep=1;
alpha0=alpha;
delta0=delta;

%%%%%

while flag==0
  %%%%%%%%%%%%%%

     if iter(size(iter,1),:) ~=1

        [miter]=pdrule(m,minitial,mmin,res,iter(size(iter,1),:),mstep, mmax,alpha0, delta0); %cab
        m=miter(1,1);
        minitial=miter(1,2);
    else
        m=minitial;
    end

    miteracion(iter(size(iter,1),:)+1,1)=m;
    beta=norm(r);
    v(:,1)=r/beta;
    h=zeros(m+1,m);
    if size(logres,1)==1
        for j=1:m                       %modified gram schmidt--Arnoldi
            w(:,j)=A*v(:,j);
            for i=1:j
                h(i,j)=w(:,j)'*v(:,i);
                w(:,j)=w(:,j)-h(i,j)*v(:,i);
            end
            h(j+1,j)=norm(w(:,j));
            if h(j+1,j)==0
                m=j;
                h2=zeros(m+1,m);
                for k=1:m
                    h2(:,k)=h(:,k);
                end
                h=h2;
            else
                v(:,j+1)=w(:,j)/h(j+1,j);
            end
        end
        g=zeros(m+1,1);
        g(1,1)=beta;
        for j=1:m                       %plane rotations (QR decompostion)
            P=eye(m+1);
            sin=h(j+1,j)/(sqrt(h(j+1,j)^2 + h(j,j)^2));
            cos=h(j,j)/(sqrt(h(j+1,j)^2 + h(j,j)^2));
            P(j,j)=cos;
            P(j+1,j+1)=cos;
            P(j,j+1)=sin;
            P(j+1,j)=-sin;
            h=P*h;
            g=P*g;
        end
        R=zeros(m,m);
        G=zeros(m,1);
        V=zeros(n,m);
        for k=1:m
            G(k)=g(k);
            V(:,k)=v(:,k);
            for i=1:m
                R(k,i)=h(k,i);
            end
        end
        minimizer=R\G;
        Z=V*minimizer;
        xm=x0 + Z;
        r=b-A*xm;
        miteracion(size(miteracion,1)+1,1)=m;
        res(restart+1,:)=norm(r);
        iter(restart+1,:)=restart+1;
        %logres(size(logres,1)+1,:)=abs(g(m+1,1)/res(1,1));
        logres(size(logres,1)+1,:)=norm(r)/res(1,1);

        if logres(size(logres,1)) <tol
            flag=1;
        else
            x0=xm;                        %update and restart
            restart=restart+1;
        end

%% Calculo de z(k)
        z(:,ij)= Z;
        %ij=ij+1;
    else
        if ij<=lL
            d=ij;
            ij=ij+1;
        end
        s=m+d;
        for j=1:s                       %modified gram schmidt--Arnoldi
                if j<=m
                    w(:,j)=A*v(:,j);
                else
                    w(:,j)=A*z(:,d-(j-m-1));
                end
            for i=1:j
                h(i,j)=w(:,j)'*v(:,i);
                w(:,j)=w(:,j)-h(i,j)*v(:,i);
            end
            h(j+1,j)=norm(w(:,j));
            if h(j+1,j)==0
                s=j;
                h2=zeros(s+1,s);
                for k=1:s
                    h2(:,k)=h(:,k);
                end
                h=h2;
            else
            v(:,j+1)=w(:,j)/h(j+1,j);
            end
        end
        g=zeros(s+1,1);
        g(1,1)=beta;
        for j=1:s                       %plane rotations (QR decompostion)
            P=eye(s+1);
            sin=h(j+1,j)/(sqrt(h(j+1,j)^2 + h(j,j)^2));
            cos=h(j,j)/(sqrt(h(j+1,j)^2 + h(j,j)^2));
            P(j,j)=cos;
            P(j+1,j+1)=cos;
            P(j,j+1)=sin;
            P(j+1,j)=-sin;
            h=P*h;
            g=P*g;
        end
        R=zeros(s,s);
        G=zeros(s,1);
        V=zeros(n,s);
        for k=1:s
            G(k)=g(k);
            V(:,k)=v(:,k);
            for i=1:s
                R(k,i)=h(k,i);
            end
        end
        for k=m+1:s
            V(:,k)=z(:,d-(k-m-1));
        end
        minimizer=R\G;
        xm=x0+V*minimizer;
        iter(restart+1,:)=restart+1;
        r=b-A*xm;
        miteracion(size(miteracion,1)+1,1)=m;
        res(restart+1,:)=norm(r);
        iter(restart+1,:)=restart+1;
        logres(size(logres,1)+1,:)=norm(r)/res(1,1);
        %logres(size(logres,1)+1,:)=abs(g(s+1,1)/res(1,1));
        aux=V*minimizer;
        Z=z;
            z(:,size(z,2)+1)=aux;
            l=size(z,2);


        %if abs(g(s+1,1))/res(1,1) <tol || size(logres,1)==maxit
         if logres(size(logres,1),1) <tol || size(logres,1)==maxit
            flag=1;
         else
            x0=xm;                        %update and restart
            restart=restart+1;
         end


    end %if restarted

end  %while flag

tiempo=toc     %Imprime tiempo de ejecucion
l
%subplot(1,1,1);
%%semilogy(logres,'r')
%%hold on
%%%final_value = logres(end);
%%%semilogy(length(logres), final_value, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
%%
%%xlabel ('Número de Ciclos');
%%ylabel ('Norma relativa del residuo');
%%title('Convergencia del metodo GMRES(m),LGMRES(m,k) y LGMRES ADAPTATIVO')
%semilogy(logres,'b')
%title(Name_Matrix);
%title(color);
%xlabel('Number of Restart Cycle');ylabel('|rj|/|r0|');
% legend(['PD-GMRES(27,alpha_{P}=', num2str(alpha0),',alpha_{D}=', num2str(delta0),'), t= ', num2str(tiempo)],'Location','Best');
%title(['Example 2.2 - Complementary cycles of GMRES. Nl=', num2str(Nl),'; delta=', num2str(dl)])
 % hold on
%  subplot(2,1,2);
%  plot(miteracion,color)
%  xlabel('Number of restart cycles');ylabel('m, restart parameters');
%   hold on
lastcycle=size(logres,1);
%tiempoC= [lastcycle tiempo];
tiempoC= tiempo
ciclos= lastcycle
NormaResidual=logres;


function miter = pdrule(m,minitial,mmin,res,iter,mstep,mmax, alpha, delta) %implementacion Cabral

 ap=0;
 ad=0;

if iter >3
    
mj= m + ceil(alpha*(res(iter,:)/res(iter-1,:)) + delta*((res(iter,:) - res(iter-2,:))/(2*res(iter-1,:))));
%mj= m + round(alpha*(res(iter,:)/res(iter-1,:)) + delta*((res(iter,:) - res(iter-2,:))/(2*res(iter-1,:))));
    
 ap=res(iter,:)/res(iter-1,:);
 ad=(res(iter,:) - res(iter-2,:))/(2*res(iter-1,:));

elseif iter >2
    mj= m + ceil(alpha*(res(iter,:)/res(iter-1,:)));
    %    mj= m + round(alpha*(res(iter,:)/res(iter-1,:)));

else
    mj=minitial;
end

% ap=res(iter,:)/res(iter-1,:)
% if iter >3
% ad=(res(iter,:) - res(iter-2,:))/(2*res(iter-1,:))
%end


% if mj <= 0
%      mj=abs(mj) + mstep;
% end
if mj < mmin
    minitial = minitial + mstep;
    mj=minitial;
end

if mj > mmax
%    miter=minitial;
%     miter=miter-mstep;
    %mj=mmax-mstep;
    mj=mmax;
end

miter= [mj minitial] ;
miter= [mj minitial ap ad] ;

    

    




