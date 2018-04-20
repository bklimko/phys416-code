clear, clc
B = zeros(51); C = zeros(51);
eig1 = zeros(1,20);
eig2 = zeros(1,20);
eig3 = zeros(1,20);
eig4 = zeros(1,20);
x  = randi(5,51,1);
for idx = 2:50
    B(idx, idx-1) = -1;
    B(idx, idx+1) = 1;
    C(idx, idx-1) = 1;
    C(idx, idx+1) = 1;
end

B(1,2) = 1;
B(1, end) = -1;
B(end,1) = 1;
B(end,end-1) = -1;
C(1,2) = 1;
C(1, end) = 1;
C(end,1) = 1;
C(end,end-1) = 1;

idx = 1;
nums = 0.1:0.1:2;
for th = nums
    A = (1/2)*C - (th/2*B);
    e = eig(A);
    eig1(idx) = max(abs(e));
    [l,~] = power1(A,x, 1e-3, 20);
    eig2(idx) = l;
    eig3(idx) = norm(A, 1);
    eig4(idx) = norm(A, inf);
    idx = idx + 1;
end

plot(nums, eig1, nums, eig2, nums, eig3, nums, eig4)
legend('Matlab eig function', 'Power method', '1-norm', 'inf-norm',...
    'Location', 'northwest')
xlabel('\tau/h')
ylabel('Spectral radius (|\lambda_1|)')

