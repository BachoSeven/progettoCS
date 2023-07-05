function [x, relres, it, resvec] = gmres_arnoldi(A, b, eps)
	n = size(A, 1);
	beta = norm(b);
	v1 = b/beta;
	resvec = [beta];
	j = 0;
	conv = 0;
	itmax = 100;
	V = zeros(n, itmax);
	H = zeros(itmax, itmax);
	V(:, 1) = v1;

	% Iterazione di Arnoldi
	while (conv == 0 && j < itmax)
		j = j + 1;
		w = A * V(:, j);

		% Ortogonalizzazione di A*v_j
		for i = 1:j
			H(i, j) = V(:, i)' * w;
		end
		w = w - V(:,1:j) * H(1:j, j);

		H(j + 1, j) = norm(w);
		% Normalizzazione del vettore ortogonalizzato
		V(:, j + 1) = w / H(j + 1, j);

		% Risoluzione del problema ai minimi quadrati
		e1 = zeros(j + 1, 1);
		e1(1) = beta;
		y = H(1:j+1,1:j) \ e1;

		% Calcolo del residuo
		relres = norm(H(1:j+1,1:j) * y - e1);
		resvec = [resvec; relres];

		% Controllo della convergenza
		if relres < eps
			x = V(:,1:j) * y;
			conv = 1;
			it = j;
		end
	end
end
