function [flux, dt] = weno3(U)
    % WENO3 scheme
    vari
    flux = zeros(3, N); % flux vector​

    W = U2W(U);
    f_1 = zeros(3, N); % 正通量
    f_2 = zeros(3, N); % 负通量
    h_1 = zeros(3, N); % 插值后的f_1
    h_2 = zeros(3, N);
    f_1x = zeros(3, N); %f1x为正通量导数 f2x为负通量导数
    f_2x = zeros(3, N);

    R = ones(3, 3); % 右特征向量
    L = ones(3, 3); % 左特征向量
    D_1 = zeros(3, 3); % 正特征值
    D_2 = zeros(3, 3); % 负特征值

    A = zeros(3, 3); % 通量的jacobi矩阵
    big = zeros(1, 2); % 记录通量导数最大值
    A(1, 2) = 1;
    A(2, 3) = 0.4;

    for i = 1:N
        p = 0.4 * (U(i, 3) - 0.5 * U(i, 2)^2 / U(i, 1)); % 压强
        c = sqrt(1.4 * p / U(i, 1)); % 声速
        H = (U(i, 3) + W(i, 3)) / U(i, 1); % 焓
        u = W(i, 2);

        % 右特征向量矩阵
        R(2, 1) = u - c;
        R(2, 2) = u;
        R(2, 3) = u + c;
        R(3, 1) = H - u * c;
        R(3, 2) = 0.5 * u^2;
        R(3, 3) = H + u * c;

        % 左特征向量矩阵
        temp1 = 0.2 / c^2;
        L(1, 1) = temp1 * (0.5 * u^2 + u * c / 0.4);
        L(2, 1) = 2 * temp1 * (-0.5 * u^2 + c^2/0.4);
        L(3, 1) = temp1 * (0.5 * u^2 - u * c / 0.4);
        L(1, 2) = temp1 * (-u - c / 0.4);
        L(2, 2) = 2 * temp1 * u;
        L(3, 2) = temp1 * (-u + c / 0.4);
        L(2, 3) = -2 * temp1;
        L(1, 3) = temp1;
        L(3, 3) = temp1;

        L * R;
        % 对角矩阵
        Gamma = diag([u - c, u, u + c]);
        f_1(:, i) = R * ((Gamma + abs(Gamma)) / 2) * L * U(i, :)';
        f_2(:, i) = R * ((Gamma - abs(Gamma)) / 2) * L * U(i, :)';

        % Jacobi 矩阵
        % A = [0, 1, 0;
        %     -0.8 * u^2, 1.6 * u, 0.4;
        %     -0.3 * u^3 - c^2 * u / 0.4, 0.1 * u^2 + c^2/0.4, 1.4 * u];

    end

    dt = CFL * dx / max(abs(W(:, 2)) + sqrt(1.4 * W(:, 3) ./ W(:, 1)));

    ep = 1e-6;
    gam = zeros(1, 2);
    gam(1) = 1/3;
    gam(2) = 2/3;
    beta = zeros(3, 2);
    omg = zeros(3, 2);
    omgm = zeros(3, 1);
    q = zeros(3, 2);

    for i = 2:N - 1

        for j = 1:3
            beta(j, 1) = (f_1(j, i) - f_1(j, i - 1))^2;
            beta(j, 2) = (f_1(j, i + 1) - f_1(j, i))^2;
        end

        for j = 1:3

            for k = 1:2
                omg(j, k) = gam(k) / (ep + beta(j, k))^2;
            end

            omgm(j) = omg(j, 1) + omg(j, 2);
        end

        for j = 1:3
            q(j, 1) = -0.5 * f_1(j, i - 1) + 1.5 * f_1(j, i);
            q(j, 2) = 0.5 * f_1(j, i) + 0.5 * f_1(j, i + 1);
        end

        for j = 1:3
            h_1(j, i) = (q(j, 1) * omg(j, 1) + q(j, 2) * omg(j, 2)) / omgm(j);
        end

    end

    for i = 2:N - 1

        for j = 1:3
            beta(j, 1) = (f_2(j, i) - f_2(j, i + 1)) * (f_2(j, i) - f_2(j, i + 1));
            beta(j, 2) = (f_2(j, i - 1) - f_2(j, i)) * (f_2(j, i - 1) - f_2(j, i));
        end

        for j = 1:3

            for k = 1:2
                omg(j, k) = gam(k) / ((ep + beta(j, k)) * (ep + beta(j, k)));
            end

            omgm(j) = omg(j, 1) + omg(j, 2);
        end

        for j = 1:3
            q(j, 1) = -0.5 * f_2(j, i + 1) + 1.5 * f_2(j, i);
            q(j, 2) = 0.5 * f_2(j, i) + 0.5 * f_2(j, i - 1);
        end

        for j = 1:3
            h_2(j, i - 1) = (q(j, 1) * omg(j, 1) + q(j, 2) * omg(j, 2)) / omgm(j);
        end

    end

    for i = 3:N - 1
        f_1x(:, i) = (h_1(:, i) - h_1(:, i - 1)) / dx;
    end

    for i = 2:N - 2
        f_2x(:, i) = (h_2(:, i) - h_2(:, i - 1)) / dx;
    end

    for i = 3:N - 2
        flux(:, i) = f_1x(:, i) + f_2x(:, i);
    end

end
