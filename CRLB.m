function crlb_matrix = CRLB(an, un, v, std_deviation_range,M,sigma)
    %ti
    ti=0:0.005:0.005*(M-1);
    for jj = 1:length(std_deviation_range)
        sigma_big=eye(2*M)*sigma^2;
        for ii = 1:M
            C(ii,ii) = 1 / std_deviation_range(jj);
            l_t(:,ii) = (an(:,ii) - un - v * ti(ii)) / norm(an(:,ii) - un - v * ti(ii));
            J_t(ii,:) = [-l_t(:,ii)', -l_t(:,ii)'*ti(ii), 1, ti(ii)];
            S_t(ii,2*ii-1:2*ii) = l_t(:,ii)';
        end
        crlb_matrix = inv(J_t' * C * J_t - J_t' * C * S_t /(S_t' * C * S_t + inv(sigma_big)) * S_t' * C * J_t);
    end
end
