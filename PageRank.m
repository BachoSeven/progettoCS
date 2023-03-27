function [y, it] = PageRank(H, v, gamma, itmax)
% function [y, it] = PageRank(H, v, gamma, itmax)
% Calcola il vettore y del PageRank applicato alla matrice di adiacenza H
% con vettore di personalizzazione v (vettore riga) e parametro gamma
% applicando al piu' itmax iterazioni del metodo delle potenze

n = size(H,1);
%  usn = 1/n;
e = ones(n,1);
d = H*e;
d = d';
dang = d==0;
dh = d + dang*n;
dh = 1./dh;
x = rand(1,n);
x = x/sum(x);
v = v/sum(v);
if size(v,2)==1
   v = v';
end

for it=1:itmax
     y = x.*dh;
%     y = y*H + usn*sum(dang.*x);
     y = y*H + sum(dang.*y);
     y = y*gamma+(1-gamma)*v;
     err = max(abs(x-y));
     x = y;
     fprintf('iterazione = %d, errore = %d \n', it, err);
     if err<1.e-13*max(y)
          break
     end
end
fprintf('Numero totale di iterazioni: %d\n',it);
if it == itmax
     disp('raggiunto il numero massimo di iterazioni')
end

