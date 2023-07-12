% Load the matrix
load coAuthorsDBLP
H = Problem.A;

% Parameters
gamma = [0.5, 0.7, 0.85, 0.99];
tol = 1e-8;
n = length(H);
e = ones(n, 1);
D = spdiags(H * e, 0, n, n);
resvecs = cell(1, length(gamma));

% Solving the linear system with GMRES
for i = 1:length(gamma)
	b = ((1-gamma(i))/n) * e;
	M = speye(n) - gamma(i) * H * D^(-1);
	tic;
	[x, res, it, resvec] = gmres_arnoldi(M, b, tol);
	elapsed = toc;
	resvecs{i} = resvec;
	disp(['gamma = ', num2str(gamma(i)), ': ', num2str(it), ' iterazioni effettuate in ', num2str(elapsed), ' secondi']);
end

% Plotting the residual norms in order to analyse convergence as gamma varies
figure;
for i = 1:length(gamma)
	semilogy(resvecs{i}, 'o-', 'DisplayName', ['\gamma = ', num2str(gamma(i))]);
	hold on;
end
grid on;
xlabel('Number of Iterations');
ylabel('Residual norms');
legend('show');
print("convergence.png");
