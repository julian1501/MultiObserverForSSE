function stableBool = isMatrixStable(A)
% isMatrixStable Function
%
% The 'isMatrixStable' function checks whether the given matrix 'A' is 
% stable, which means all eigenvalues of 'A' should have strictly negative 
% real parts.
%
% Documentation written by ChatGPT.
%
% Syntax:
% -------
% stableBool = isMatrixStable(A)
%
% Inputs:
% -------
% - 'A': A square matrix for which stability is to be checked.
%     - 'A' should be a square matrix (i.e., number of rows equals the 
%       number of columns).
%
% Outputs:
% --------
% - 'stableBool': A logical value indicating if the matrix 'A' is stable:
%     - 'true': All eigenvalues of 'A' have a strictly negative real part 
%       (the matrix is stable).
%     - 'false': At least one eigenvalue of 'A' has a non-negative real 
%       part (the matrix is unstable).
%
% Description:
% ------------
% This function calculates the eigenvalues of the matrix 'A' and checks whether 
% all of them have strictly negative real parts. If any eigenvalue has a real 
% part greater than or equal to 0, the matrix is considered unstable, and the 
% function returns 'false'. Otherwise, it returns 'true', indicating that the 
% matrix is stable.
%
% Example:
% --------
% Check matrix stability:
% 
% A = [0 1; -2 -3];  % Example system matrix
% stableBool = isMatrixStable(A);
% 
% If 'A' is stable, 'stableBool' will be 'true', otherwise 'false'.
%
% Notes:
% ------
% - The function assumes the input matrix 'A' is square.
%
% See also:
% ---------
% eig, real

    % Initialize stableBool as true
    stableBool = true;
    eigenvalues = eig(A);
    n = size(eigenvalues,1);
    % Loop through eigenvalues and compare each eigenvalue to 0
    for i = 1:1:n
        if eigenvalues(i) > 0
            stableBool = false;
        end
    end
end