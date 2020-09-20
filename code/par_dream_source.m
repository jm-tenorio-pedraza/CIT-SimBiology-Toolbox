function [x, p_x] = par_dream(prior, pdf, N, T,d)
[delta, c, c_star, n_CR, p_g] = deal(3, 0.1, 1e-12, 3, 0.2);
x = nan(T,d,N); p_x = nan(T,N);
CR = [1:n_CR]/n_CR; p_CR = ones(1,n_CR)/n_CR;
[J, n_id] = deal(zeros(1,n_CR));
for i=1:N, R(i, 1:N-1) = setdiff(1:N,i); end

X = prior(N,d);
for i=1:N, p_X(i, 1) = pdf(X(i, 1:d)); end
x(1, 1:d, 1:N) = reshape(X', 1, d,N); p_x(1,1:N) = p_X';

for t = 2:T
    [~, draw] = sort(rand(N-1, N));
    dX = zeros(N, d);
    lambda=(unifrnd(-c, c, N,1));
    std_X = std(X);
    for i=1:N
        D = randsample([1:delta], 1, true);
        a = R(i, draw(1:D, i)); b = R(i, draw(D+1:2*D, i));
        id(i) = randsample(1:n_CR, 1, true, p_CR);
        z = rand(1, d);
        A = find(z < CR(id(i)));
        d_star = numel(A);
        if (d_star == 0), [~, A] = min(z); d_star = 1; end
        gamma_d = 2.38/sqrt(2*D*d_star);
        g = randsample([gamma_d 1], 1, true, [1-p_g, p_g]);
        dX(i,A) = c_star*rand(1, d_star) + (1+lambda(i))*gsum(X(a,A)-X(b,A),1);
        X_p(i,1:d) = X(i, 1:d) = dX(i,1:d);   
    end
    
    parfor i=1:N
        p_Xp(i,1) = pdf(Xp(i, 1:d));
        p_acc = min(1, p_Xp(i,1)./p_X(i,1));
      
        if p_acc > rand
            X(i, :) = Xp(i, 1:d); p_X(i,1) = p_Xp(i,1);
        else
            dX(i,:) = 0;
        end
      
    end
    for i=1:N
        J(id(i)) = J(id(i)) + sum((dX(i,1:d)./std_X).^2);
        n_id(id(i)) = n_id(id(i)) + 1;
    end
    x(t,1:d,1:N) = reshape(X', 1, d, N); p_x(t, 1:N) = p_X';
    if t<T/10, p_CR = J./n_id; p_CR = p_CR/sum(p_CR); end
    [X, p_X] = check(X, mean(log(p_x(ceil(t/2):t,1:N)))); 
end