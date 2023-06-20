function [x, relres, it, resvec] = gmres_arnoldi(A, b, tol)
    n = size(A, 1);
    x = zeros(n,1);
    r = b - A * x;
    beta = norm(r);
    resvec = [beta];
    it=0;
    conv=0;

    while conv == 0
        it = it + 1;
        Q = zeros(n, it + 1);
        H = zeros(it + 1, it);
        Q(:, 1) = r / beta;

        for j = 1:it
            w = A * Q(:, j);

            for i = 1:j
                H(i, j) = Q(:, i)' * w;
                w = w - H(i, j) * Q(:, i);
            end

            H(j + 1, j) = norm(w);
            Q(:, j + 1) = w / H(j + 1, j);
        end

        e1 = zeros(it + 1, 1);
        e1(1) = beta;
        y = H \ e1;
        x = x + Q(:, 1:it) * y;
        relres = abs(H(it + 1, it) * y(it)) / beta;
        resvec = [resvec; relres];

        if relres < tol
            conv = 1;
        end
    end
end