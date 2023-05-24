function dist = Dist2(x, R_orig)

[d, ~, n] = size(x);
x_err = zeros(d);
for i = 1 : n
    xi = x(:,:,i);
    x_err = x_err + xi' * R_orig(:,:, i); % compute the error in each step
end
[Ur,Lr,Vr] = svd(x_err);
S0 = diag([ones(1,d-1),det(Ur*Vr')]);
dist = real(sqrt(2*d - 2*sum(sum(Lr.*S0)) / n));

end