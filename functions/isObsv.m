function bool = isObsv(A,C)
% isObsv Function
%
% The 'isObsv' function determines whether the pair of matrices A and C is 
% observable. It checks if the rank of the observability matrix matches the
% size of the system's state matrix.
%
% Documentation written by ChatGPT.
%
% Syntax:
% -------
% bool = isObsv(A, C)
%
% Inputs:
% -------
% - 'A': The state transition matrix of the system, which must be a square 
%   matrix of size (n x n).
% - 'C': The output matrix of the system, which must have dimensions (m x n),
%   where m is the number of outputs and n is the number of states.
%
% Outputs:
% --------
% - 'bool': A logical value indicating whether the system is observable:
%     - 'true': If the pair (A, C) is observable.
%     - 'false': If the pair (A, C) is not observable.
%
% Description:
% ------------
% This function performs the Observability test (PBH test) on the system 
% defined by matrices 'A' (state transition matrix) and 'C' (output matrix).
% It checks the rank of the matrix formed by combining 'C' and each 
% eigenvalue of 'A'. If the rank of the matrix does not match the size of 
% the system's state matrix 'A', then the system is considered unobservable.
%
% Example:
% --------
% Observable system:
% A = [0 1; 0 0];
% C = [1 0];
% bool = isObsv(A, C);
% % bool will be true (system is observable).
% 
% Unobservable system:
% A = [0 1; 0 0];
% C = [0 1];
% bool = isObsv(A, C);
% % bool will be false (system is unobservable).
% 
% Notes:
% ------
% - This function uses the PBH (Popov-Belevitch-Hautus) test for observability, 
%   which checks if the rank of a specific matrix formed by the system matrices is full.
% - 'A' must be a square matrix, and 'C' should have the same number of columns as 'A'.
% - Eigenvalues of 'A' are used to form the matrices in the test.
%
% See also:
% ---------
% eig, rank

    n = size(A,1);
    ev = eig(A);
    I = eye(n);
    % PBH test
    bool = true;
    for i = 1:1:n
        if rank([ev(i)*I - A; C]) ~= n
            bool = false;
            break
        end
    end

end