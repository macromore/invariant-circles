% For use and purpose: see documentation file
% 
% Dependencies:
% N/A
classdef Fourier
% For use and purpose: see documentation file
% 
% Dependencies:
% N/A
    properties % define the properties of the class
        coeffs % a list of coefficents
        len % number of modes in the fourier series
    end % end properties
    methods % define the public methods
        function obj = Fourier(c,N)
        % initialize the Fourier series
            if nargin == 1
                if size(c,1) == 2
                   N = size(c,2)-1;
                   obj.coeffs = zeros(1,2*N+1);
                   obj.coeffs(N+1) = c(1,1);
                   for i = 1:N
                      obj.coeffs(N-i+1) = c(1,i+1)/2;
                      obj.coeffs(N-i+1) = obj.coeffs(N-i+1)+1i*c(2,i+1)/2;
                      obj.coeffs(N+i+1) = c(1,i+1)/2;
                      obj.coeffs(N+i+1) = obj.coeffs(N+i+1)-1i*c(2,i+1)/2;
                   end % end for loop
                elseif mod(size(c,2),2)==1
                    obj.coeffs = c;
                else % Catch error
                    error('You need to fill in zeros so the index is -N to N');
                end % end if
            elseif nargin == 2
                fcn = inline(c); %#ok<DINLN>
                c = zeros(1,2*N+1);
                for n = -N:N
                    c(N+n+1) = fcn(n);
                end % end for loop
                obj = Fourier(c);
            end % end if
            obj.len = (length(obj.coeffs)-1)/2;
        end % end Fourier
        function m = mat(a)
        % Return the column vector of coefficents
            if size(a.coeffs,1) == 1
                m = transpose(a.coeffs); 
            else % accept vector as is
                m = a.coeffs;
            end % end if
        end % end mat
        function m = matrix(a)
        % Return the column vector of coefficents
            if size(a.coeffs,1) == 1
                m = transpose(a.coeffs); 
            else % accept vector as is
                m = a.coeffs;
            end % end if
        end % end matrix        
        function m = spMat(a)
        % Return the column vector of coefficents
            if size(a.coeffs,1) == 1
                m = sparse(transpose(a.coeffs)); 
            else % accept vector as is
                m = sparse(a.coeffs);
            end % end if
        end % end matrix
        function leng = length(obj)
        % Return the number of modes of the Fourier series
            leng = obj.len;
        end % end function length
        function disp(obj)
        % display the coefficents
            disp(obj.coeffs);
        end % end disp
        function sref = subsref(obj,n)
        % get the n-th Fourier coefficent
            n = n.subs{1};
            if abs(n) <= obj.len
                sref = obj.coeffs(obj.len + 1 + n);
            else % out of range
                sref = 0;
            end % end if
        end % end subsref
        function coeff = getCoeff(obj,n)
        % get the n-th Fourier coefficent
            if abs(n) <= obj.len
                coeff = obj.coeffs(obj.len+ 1 + n);
            else % out of range
                coeff = 0;
            end % end if
        end % end getCoeff
        function sum = plus(a,b)
        % addition of Fourier series
            [a, b] = cast(a,b);
            N = max(length(a), length(b));
            c = zeros(1,2*N+1);
            for i = -N:N
                c(i+N+1) = getCoeff(a,i) + getCoeff(b,i);
            end % end for loop
            sum = Fourier(c);
        end % end plus
        function nega = uminus(a)
        % Negation
           N = length(a);
           c = zeros(1,2*N+1);
           for i = -N:N
               c(i+N+1) = - getCoeff(a,i);
           end % end for loop
           nega = Fourier(c);
        end % end uminus
        function sum = minus(a,b)
        % subtraction
            [a, b] = cast(a,b);
            sum = a+(-b);
        end % end minus
        function prod = mtimes(a,b)
        % Multiplication with truncation
            [a, b] = cast(a,b);
            prod = truncate(Fourier(conv(mat(a),mat(b))),max(length(a),length(b)));
        end % end mtimes
        function quo = fourierinverse(a)
        % Compute the formal inverse
           if (getCoeff(a,0) ~= 0) 
               N = length(a);
               b = FourierOperator(a);
               c = truncate(Fourier(1),N);
               quo = transpose(linsolve(mat(b),mat(c)));
               quo = Fourier(quo);
           else % not invertible
               error('Non-invertible'); 
           end % end if
        end % end fourierinverse
        function quo = mrdivide(a,b)
        % division a/b
            [a, b] = cast(a,b);
            quo = a*fourierinverse(b);
        end % end mrdivide
        function power = mpower(a,n)
        % basic exponentiation
            power = 1;
            for i = 1:n
                power = power*a;
            end % end for loop
        end % end mpower
        function r = rotation(a,rho)
        % apply a rotation rho
            N = length(a);
            r = zeros(1,2*N+1);
            for k = -N:N
               r(N+k+1) = getCoeff(a,k)*exp(2*pi*1i*k*rho); 
            end % end for loop
            r = Fourier(r);
        end % end rotation 
        function new = truncate(obj,N)
        % change number of modes
            new = zeros(1,2*N+1);
            for i = -N:N
                new(N+i+1) = getCoeff(obj,i);
            end % end for loop
            new = Fourier(new);
        end % end truncate
        function sum = l1Norm(a,nu)
        % take the l^1_nu norm
            sum = 0;
            N = length(a);
            for i = -N:N
               sum = sum + abs(getCoeff(a,i)*nu^abs(i)); 
            end % end for loop
        end % end norm
        function coeffs = realcoeffs(a)
        % compute the sine and cosine series
           N = length(a);
           coeffs = zeros(2,N+1);
           coeffs(1,1) = getCoeff(a,0);
           for i = 2:N+1
              coeffs(1,i) = (getCoeff(a,i-1))+(getCoeff(a,-i+1));
              coeffs(2,i) = -1i*((getCoeff(a,i-1))-(getCoeff(a,-i+1)));
           end % end for loop
        end % end realcoeffs
        function val = endVals(a)
        % compute the average of the two tail values
           val = (abs(a.coeffs(end))+abs(a.coeffs(1)))/2;
        end % end endVals
        function derv = diff(a)
        % take the formal derivative
           N = length(a);
           derv = zeros(1,2*N+1);
           for i = -N:N
              derv(N+i+1) = 1i*2*pi*i*getCoeff(a,i);
           end % end for loop
           derv = Fourier(derv);
        end % end diff
        function inte = int(a,c)
        % Integrate with c as constant of integration
           N = length(a);
           inte = zeros(1,2*N+1);
           inte(N+1) = c;
           for i = 1:N
              inte(N+i+1) = getCoeff(a,i)/(1i*2*pi*i);
              inte(N-i+1) = getCoeff(a,-i)/(1i*2*pi*(-i));
           end % end for loop
           inte = Fourier(inte);
        end % end int
        function val = evaluate(obj,x)
        % evaluate the fourier series at x
          N = length(obj);
          val = 0;
          for n = -N:N
              val = val + getCoeff(obj,n)*exp(n*2*pi*1i*x);
          end % end for loop
        end % end evaluate
        function graph = plot(a)
        % graph the series as a function of one variable
            x = 0:1/100:1;
            y = evaluate(a,x);
            graph = plot(x,y);
        end % end plot
        function xnplus1 = newtonstep(a, xn)
        % make a newton step given a suspected zero
            xnplus1 = xn - evaluate(a,xn)/evaluate(diff(a),xn);
        end % end newtonstep
        function xn = newton(a,x0)
        % preform a newton given a suspected zero
           count = 0;
           tol = 1;
           xn = x0;
           while (count < 100000 && tol > 10^-14)
               tol = xn;
               xn = newtonstep(a,xn);
               tol = abs(tol-xn);
               count = count + 1;
           end % end while loop
        end % end newton
    end % end methods 
    methods (Access = protected)
        function [a, b] = cast(a,b)
        % convert / verify two objects are Fourier class
           if ~isa(a, 'Fourier')
                a = Fourier(a);
           end % end if
           if ~isa(b, 'Fourier')
                b = Fourier(b);
           end % end if
        end % end cast
    end % end private methods
end % end class