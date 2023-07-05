% Carica la matrice
load coAuthorsDBLP
A = Problem.A;

% Parametri
gamma = [0.5, 0.7, 0.85, 0.99];
tol = 1e-8;
n = length(A);
e = ones(n, 1);
D = spdiags(A * e, 0, n, n);
resvecs = cell(1, length(gamma));

% Risoluzione del sistema lineare tramite GMRES
for i = 1:length(gamma)
	b = ((1-gamma(i))/n) * e;
	M = speye(n) - gamma(i) * A * D^(-1);
	tic;
	[x, relres, it, resvec] = gmres_arnoldi(M, b, tol);
	elapsed = toc;
	resvecs{i} = resvec;
	disp(['gamma = ', num2str(gamma(i)), ': ', num2str(it), ' iterazioni effettuate in ', num2str(elapsed), ' secondi']);
end

% Plot dei residui per analizzare la convergenza al variare di gamma
figure;
for i = 1:length(gamma)
	semilogy(resvecs{i}, 'o-', 'DisplayName', ['\gamma = ', num2str(gamma(i))]);
	hold on;
end
grid on;
xlabel('Numero di iterazioni');
ylabel('Residuo');
title('Convergenza di GMRES al variare di gamma');
legend('show');
print("convergenza.png");
