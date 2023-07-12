function [x, res, it, resvec] = gmres_arnoldi(A, b, eps)
	n = size(A, 1);
	maxit = 100;
	beta = norm(b);
	resvec = [beta];
	v1 = b/beta;
	V = zeros(n, maxit);
	H = zeros(maxit, maxit);
	V(:, 1) = v1;

	% Arnoldi iteration
	conv = 0;
	j = 0;
	while (conv == 0 && j < maxit)
		j = j + 1;
		w = A * V(:, j);

		% Ortogonalization of A*v_j
		H(1:j, j) = V(:, 1:j)' * w;
		w = w - V(:,1:j) * H(1:j, j);

		% Normalization of the orthogonalized vector
		H(j + 1, j) = norm(w);
		V(:, j + 1) = w / H(j + 1, j);

		% Solving the least-squares problem
		e1 = zeros(j + 1, 1);
		e1(1) = beta;
		y = H(1:j+1,1:j) \ e1;

		% Computing the residual norm
		res = norm(H(1:j+1,1:j) * y - e1);
		resvec = [resvec; res];

		% Convergence check
		if res < eps
			x = V(:,1:j) * y;
			conv = 1;
			it = j;
		end
	end
end
