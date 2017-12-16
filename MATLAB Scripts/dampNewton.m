%Chang Vang
%MAT226A Homework 2
%Algorithm for Newton's Method with Damping
%Solves F(x) = 0 for x.  F(x) can be a system of equations (SOEs)

%Input: 
%   xo - Initial guess
%   Df - Derivative or Jacobian (in the case of a SOEs)
%        NOTE: Df must be passed in as an anonymous function 
%        using the '@' symbol
%   fx - The function itself passed as an anonymous function
%   FTOL - max norm for functional evaluations
%   XTOL - max norm for successive iterates
%   kmax - max number of iterations
%   lmax - minimum damping factor

function [xo, k, EuNormFx, DX, dampfactor] = ...
         dampNewton(xo,Df,fx,FTOL,XTOL,kmax,lmax)


n = length(fx(xo));            %check if F is in R or in Rn
DX = 1;                        %initialize ||xk - xk-1||
k = 0;                         %initilize iteration counter

if fx(xo) == zeros(1,n)';      %check initial guess
    disp('Stopped. Initial Guess Produces F(x)=0 exactly.')
    disp('Initial guess = ')
    disp(xo)
    return
else
    EuNormFx = sqrt(sum(fx(xo).^2));
    while ( EuNormFx >= FTOL)&(DX >= XTOL)&(k < kmax)
        
        %solve for DeltaXk
        if n == 1       %if F is in R
            DxK = -fx(xo) / Df(xo);
        else            %if F is in Rn
            DxK = Df(xo) \ ( -fx(xo) );
        end
        
        
        %Monotinicity Test
        for l=0:lmax          %minimum lambda is 1/1024
            if k == 0
                lamda = 1/(2^l);
            else
                lamda = min( 1,2*lamda/(2^l) );
            end
            
            %solve for DeltaXkBar
            if n==1         %if F is in R
                DxKbar  = -fx(xo+lamda*DxK) / Df(xo);
            else            %if F is in Rn
                DxKbar  = Df(xo) \ ( -fx(xo+lamda*DxK) );
            end
            
            theta = 1 - (lamda/2);
            EunormDxKbar = sqrt( sum(DxKbar.^2) );
            EunormDxK  = sqrt(sum(DxK.^2));
            
            %perform monoticity test
            if (EunormDxKbar <= theta*EunormDxK)
                break
            end
            
        end %l
        
        xold = xo;
        xo = xo + lamda*DxK;
        DX = sqrt(sum( (xo - xold).^2 ));
        EuNormFx = sqrt(sum(fx(xo).^2));
        dampfactor(k+1) = lamda;    %store damping factor
        k = k + 1;                  %update iteration counter
    end %while
    
end %if fx(xo)