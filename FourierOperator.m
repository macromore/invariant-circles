classdef FourierOperator
% For use and purpose: see documentation file
%   
% Dependencies:
% N/A
    properties
        matrix
        len
    end % end properites
    methods
        function obj = FourierOperator(m, N, type)
        % initialize a FourierOperator
            if nargin == 1
                if isa(m,'Fourier')
                    N = length(m);
                    obj.len = N;
                    obj.matrix = zeros(2*N+1,2*N+1);
                    M = matrix(m);
                    obj.matrix(:,N + 1) = M;
                    for i = 1:N
                        obj.matrix(:,N + 1 + i) = [zeros(i,1); M(1:2*N+1-i)];
                        obj.matrix(:,N + 1 - i) = [M(i+1:end); zeros(i,1)];
                    end % end for loop
                else % given a matrix not Fourier series
                    if size(m,1) == size(m,2)
                        obj.matrix = m;
                        obj.len = (size(m,1)-1)/2;
                    else % convert matrix to Fourier series
                        obj = FourierOperator(Fourier(m));
                    end % end if
                end % end if
            elseif strcmp(type, 'diagonal') % given a function
                fcn = m; % Use an anonymous function% inline(m); %#ok<DINLN>
                mat = zeros(2*N+1, 2*N+1);
                for i = -N:N
                   mat(N+i+1,N+i+1) = fcn(i); 
                end % end for loop
                obj = FourierOperator(mat);
            end % end if
        end % end FourierOperator
        function disp(obj)
        % Show the coefficent matrix
           disp(obj.matrix); 
        end % end disp
        function m = mat(a)
        % return the matrix of a
           m = a.matrix; 
        end % end mat        
        function m = spMat(a)
        % return the matrix of a
           m = sparse(a.matrix); 
        end % end mat
        function leng = length(obj)
        % return the number of modes appropriate
           leng = obj.len;
        end % end length
        function entry = getEntry(obj, i , j)
        % get the i,j entry as if it was indexed as a fourier series
            if ((abs(i) <= obj.len)&&(abs(j) <= obj.len))
                entry = obj.matrix(obj.len + i + 1, obj.len + j + 1);
            else % out of bounds
                entry = 0;
            end % end if
        end % end getEntry
        function sref = subsref(obj, n)
        % get the entry like in normal matlab just centered at 0
            i = cell2mat(n.subs);
            j = i(2);
            i = i(1);
            sref = getEntry(obj, i,j);
        end % end subsref
        function sum = plus(a,b)
        % addition
            [a, b] = cast(a,b);
            sum = FourierOperator(a.matrix + b.matrix);
        end % end plus
        function nega = uminus(a)
        % unitary negation
            nega = FourierOperator(-a.matrix); 
        end % end uminus
        function sum = minus(a,b)
        % subtraction
           [a, b] = cast(a,b);
           sum = a+(-b); 
        end % end minus
        function new = truncate(obj,N)
        % change numer of modes
            if N ~= obj.len
                new = zeros(2*N+1, 2*N+1);
                for i = -N:N
                    for j = -N:N
                        new(N+i+1,N+j+1) = getEntry(obj,i,j);
                    end % end for loop
                end % end for loop
                new = FourierOperator(new);
            else % return the same size 
                new = obj;
            end % end if
        end % end turncate
        function prod = mtimes(a,b)
        % multiplication
            [a, b] = cast(a,b);
%             if isa(b,'FourierOperator')
            prod = FourierOperator(a.matrix*b.matrix);
% This part should never run since we used cast
%             elseif isa(b,'Fourier')
%                 bpc = zeros(2*N+1,1);
%                 for i = -N:N
%                    bpc(N+i+1,1) = b(i); 
%                 end % end for loop
%                 prod = Fourier(transpose(a.matrix*bpc));
%             end % end if
        end % end mtimes
        function mat = horzcat(a,b)
        % concatenate
           [a,b] = cast(a,b,1);
           mat = [a.matrix b.matrix];
        end % end horzcat
        function mat = vertcat(a,b)
        % vertical concatenate
           [a,b] = cast(a,b,1);
           mat = [a.matrix; b.matrix];
        end % end vertcat
        function quo = operatorinverse(a)
        % numerical inverse
           quo = FourierOperator(inv(a.matrix)); 
        end % end operatorinverse
        function quo = mrdivide(a,b)
        % division
           [a, b] = cast(a,b);
           quo = a*operatorinverse(b);  
        end % end mrdivide
        %possibly faster than taking the inverse of the operator?
        function preim = inverseimage(a,b)
        % compute the preimage under the matrix... numerically
            N = max(length(a), length(b));
            at = truncate(a,N);
            bt = truncate(b,N);
            preim = linsolve(at.matrix,mat(bt));
            preim = Fourier(transpose(preim));
        end % end inverseimage
    end % end methods
    methods (Access = protected)
        function [a, b] = cast(a,b)
        % convert to fourier operators
           if ~isa(a,'FourierOperator')
              a = FourierOperator(a); 
           end % end if
           if ~isa(b,'FourierOperator')
              b = FourierOperator(b);
           end % end if
            if nargin == 2
               N = max(a.len, b.len);
               a = truncate(a,N);
               b = truncate(b,N);
            end % end if
        end % end cast
    end % end private methods
end % end class