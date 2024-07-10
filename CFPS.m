function phi1_est = CFPS(an_m, alpha,std_deviation_range, M,sigma)

    ti=0:0.005:0.005*(M-1);
   
    p_estimate = an_m;
    h = zeros(M-1, 1);
    G = zeros(M-1, 8);

    for ii = 1:M-1
        h(ii) = norm(p_estimate(:, ii+1))^2 - norm(p_estimate(:, 1))^2 - (alpha(ii+1)^2 - alpha(1)^2);
        G(ii, :) = [2*(p_estimate(:, ii+1) - p_estimate(:, 1))', 2*(ti(ii+1)*p_estimate(:, ii+1) - ti(1)*p_estimate(:, 1))', ...
            2*(alpha(1) - alpha(ii+1)), 2*(ti(1)*alpha(1) - ti(ii+1)*alpha(ii+1)), -(ti(1)^2 - ti(ii+1)^2), -2*(ti(1) - ti(ii+1))];
    end

    B = zeros(M-1, M);
    C = zeros(M-1, 2*M);

    for ii = 1:M-1
        B(ii, 1) = 1;
        B(ii, ii+1) = -1;
        C(ii, 1:2) = -[1,1];
        C(ii, 2*ii+1:2*ii+2) = [1,1];
    end

    G1 = G(:, 1:6);
    G2 = G(:, 7:8);
    V = null(G2');
    A1 = V' * G1;
    b1 = V' * h;
    D = V' * B;
    F = V' * C;
    Qn=std_deviation_range*eye(M);
    Qp=sigma^2*eye(2*M);
    W1 = inv(D * Qn * D' + F * Qp * F');
    phi1 = inv(A1' * W1 * A1) * A1' * W1 * b1;
    iter=0;
    while iter<2
        for ii = 1:M-1
            B(ii, 1) = 2 * (alpha(1) - phi1(5) - phi1(6) * ti(1));
            B(ii, ii+1) = -2 * (alpha(ii+1) - phi1(5) - phi1(6) * ti(ii+1));
            C(ii, 1:2) = -2 * (phi1(1:2) + ti(1) * phi1(3:4) - p_estimate(:, 1))';
            C(ii, 2*ii+1:2*ii+2) = 2 * (phi1(1:2) + ti(ii+1) * phi1(3:4) - p_estimate(:, ii+1))';
        end
        D = V' * B;
        F = V' * C;
        W1 = inv(D * Qn * D' + F * Qp * F');
        phi1 = inv(A1' * W1 * A1) * A1' * W1 * b1;
        iter=iter+1;
    end

    % stage_2
    phi = [phi1;
        phi1(6)^2 - norm(phi1(3:4))^2;
        phi1(5) * phi1(6) - phi1(1:2)' * phi1(3:4)];
    H = [eye(6);
        0, 0, -2 * phi1(3:4)', 0, 2 * phi1(6);
        -phi1(3:4)', -phi1(1:2)', phi1(6), phi1(5)];
    W2 = inv(B * Qn * B' + C * Qp * C');
    delta_phi1 = -inv(H' * G' * W2 * G * H) * H' * G' * W2 * (h - G * phi);
    phi1_est = phi1-delta_phi1;
end
