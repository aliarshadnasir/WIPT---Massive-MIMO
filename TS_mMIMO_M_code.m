clc; clear;
M_vec = [100:10:150];
N = 20;
N1 = 5;

simulations = 100;
M = 500;
N = 20;
h_NLOS_store = zeros(M,N,simulations);
for sim = 1:simulations
    h_NLOS_store(:,:,sim) = 1/sqrt(2)*(randn(M,N)+1i*randn(M,N));
end

dispi = 0;
fc = 2e9;
lambda = 3e8/fc;
dist = lambda/2;
L = 6; sig_phi = 5*pi/180;
phi_store_deg = [-60+120/(2*N):120/N:60-120/(2*N)];
phi_store = phi_store_deg*pi/180; % angle of arrival
phi_ln = [phi_store_deg-25 ; phi_store_deg-15 ; phi_store_deg-5; phi_store_deg+5; phi_store_deg+15; phi_store_deg+25]*pi/180; % assuming L = 6 clusters
% phi_ln_deg = zeros(L,N);
% for l = 1:L
%     for n = 1:N
%         phi_ln_deg(l,n) = phi_store_deg(n)+80*rand;
%     end
% end
% phi_ln = phi_ln_deg*pi/180;
theta = 2*pi*dist*sin(phi_store(1:N1))/lambda;
BW = 20e6;

KR = 10^(10/10); % (Rician factor) 0 dB
P = 10.^(46/10)*1e-3; % transmit power vector (46 dBm)

noise_power_desnity = 10^(-174/10)*1e-3; % -174 dBm/Hz
s2 = noise_power_desnity*BW;
% scal_fac = 1/sqrt(s2);
scal_fac = 1e5;
s2 = s2*scal_fac^2;



g_tol = 1e-3; % convergence tolerence



e0 = 10.^((-20)/10)*1e-3; % EH threshold
e0 = e0*scal_fac^2;

eta = 0.5;
d_near_st = 20; d_near_end = 50; diff_near = d_near_end - d_near_st;
d_far_st = 65; d_far_end = 300; diff_far = d_far_end - d_far_st;

d = [d_near_st+diff_near/(2*N1):diff_near/N1:d_near_end-diff_near/(2*N1) d_far_st+diff_far/(2*(N-N1)):diff_far/(N-N1):d_far_end-diff_far/(2*(N-N1))]; % N = 10;

rate_vec = zeros(length(M_vec),1);
for M = M_vec
    display(['# of Antennas = ' num2str(M)]);
    rate_s = 0;
    iter_s = 0;
    tI_s = 0;
    mu = N*s2/P;
    for sim = 1:simulations

        h1_NLOS = h_NLOS_store(1:M,1:N1,sim); % M x N1
        h2_NLOS = h_NLOS_store(1:M,N1+1:N,sim); % M x (N-N1)
        h_LOS = exp(1i*[0:M-1]'*theta); % M x N1
        
        h_tild = zeros(M,N);
        h_tild(:,1:N1) = sqrt(KR/(1+KR))*h_LOS+sqrt(1/(1+KR))*h1_NLOS; % composite channel
        h_tild(:,N1+1:1:N) = h2_NLOS;
        
        %     beta_dB = zeros(M,N);
        %     beta_dB(:,1:N1) = 30+30*log10(repmat(d(1:N1),M,1)); % in dB
        %     beta_dB(:,N1+1:N) = 30+37.6*log10(repmat(d(N1+1:N),M,1)); % in dB
        beta_dB = zeros(N,1);
        beta_dB(1:N1) = 30+30*log10(d(1:N1)); % in dB
        beta_dB(N1+1:N) = 30+37.6*log10(d(N1+1:N)); % in dB
        beta = 10.^(-beta_dB/10);
        H = zeros(M,N);
        Theta = zeros(M,M,N);
        for n = 1:N
            for p = 1:M
                for q = 1:M
                    Theta(p,q,n) = 0;
                    for l = 1:L
                        Theta(p,q,n) = Theta(p,q,n) + exp(1i*pi*(p-q)*sin(phi_ln(l,n)))*exp(-sig_phi^2/2*(pi*(p-q)*cos(phi_ln(l,n)))^2);
                    end
                end
            end
            Theta(:,:,n) = Theta(:,:,n)/L;
            H(:,n) = sqrtm(Theta(:,:,n))*h_tild(:,n);
        end
        
        %     h = sqrt(beta).*h;
        
        H = H*scal_fac;
        
        %% initialization
        
        tE = 0.5; tI = 1-tE;
        
        H_old = H*diag(sqrt(beta));
        %     WI = H_old*inv(H_old'*H_old + mu*eye(N));
        %     WI = norm(H_old,'fro')*WI; % check
        
        %         WI = H*inv(H'*H + mu*eye(N)); % RZF
        WI = H*inv(H'*H); % ZF
        WI = norm(H,'fro')*WI; % better to do this to get nominal range values,
        % in fact if you don't do this normalization, the code may fail to work
        
        
%             WI = H_old; % conjugate beamforming (doesn't perform better)
        
        cvx_begin quiet
        variables pEE(N1,1) pI(N,1)
        expressions fac1(N1,1) fac2 fac3
        
        fac2 = 0; fac3 = 0;
        for n = 1:N
            if (n <= N1)
                fac1(n) = 0;
                for n1 = 1:N1
                    fac1(n) = fac1(n) + pEE(n1)*beta(n)*beta(n1)*abs(H(:,n)'*H(:,n1))^2;
                end
                fac2 = fac2 + beta(n)*norm(H(:,n))^2*pEE(n); % /scal_fac^2
            end
            fac3 = fac3 + norm(WI(:,n))^2*pI(n);
        end
        
        subject to
        
        for n = 1:N1
            pEE(n) >= 0;
            eta*fac1(n) >= e0/tE;
        end
        pI >= 0;
        tI*fac3+tE*fac2 <= P;
        cvx_end
        
        r = zeros(N,1);
        for n = 1:N
            int = 0;
            for n1 = 1:N
                if (n1 == n)
                else
                    int = int + abs(H(:,n)'*WI(:,n1))^2*pI(n1);
                end
            end
            r(n) = tI*log(1 + (beta(n)*abs(H(:,n)'*WI(:,n))^2*pI(n))/( beta(n)*int + s2));
        end
        
        if (dispi == 1)
            display(num2str(min(r)/log(2)));
        end
        
        %% optimization
        
        if (dispi == 1)
            display([num2str(sim) ' ' num2str('optimization running')]);
        end
        loop_exit = 0;
        iter = 0;
        obj_old = 1000;
        while (loop_exit == 0)
            iter = iter+1;
            pI_k = pI; pE_k = pEE; tI_k = tI; tE_k = tE;
            
            cvx_begin quiet
            variables tE tI inv_tI inv_tE
            variables pEE(N1,1) pI(N,1) inv_pI(N,1)
            variable min_IR
            expressions fac_EH(N1,1) r_k(N,1)
            expressions wIpI hpE LHS_P
            %         inv_tI = 1/tI; inv_pI = 1./pI;
            %         inv_tI = 1/tI; inv_tE = 1/tE;
            %         inv_pI = 1./pI;
            
            r_at_k = zeros(N,1);
            g = zeros(N,N);
            wIpI_k = 0; wIpI = 0; hpE = 0; hpE_k = 0;
            for n = 1:N
                wIpI_k = wIpI_k + norm(WI(:,n))^2*pI_k(n);
                wIpI = wIpI + norm(WI(:,n))^2*pI(n);
                if (n<=N1)
                    hpE = hpE + beta(n)*norm(H(:,n))^2*pEE(n);
                    hpE_k = hpE_k + beta(n)*norm(H(:,n))^2*pE_k(n);
                end
                int = 0;
                for n1 = 1:N
                    if (n1 == n)
                    else
                        int = int + abs(H(:,n)'*WI(:,n1))^2*pI_k(n1);
                    end
                    g(n,n1) = beta(n)*abs(H(:,n)'*WI(:,n1))^2;
                end
                r_at_k(n) = tI_k*log(1 + (beta(n)*abs(H(:,n)'*WI(:,n))^2*pI_k(n))/( beta(n)*int + s2));
                r_k(n) = 2*r_at_k(n) - r_at_k(n)*tI_k*inv_tI + (g(n,n)*pI_k(n)*tI_k)/(g(n,:)*pI_k + s2) ...
                    * (2 - pI_k(n)*inv_pI(n) -  (g(n,:)*pI-g(n,n)*pI(n)+s2)/(g(n,:)*pI_k-g(n,n)*pI_k(n)+s2) );
            end
            
            
            for n = 1:N1
                fac_EH(n) = 0;
                for n1 = 1:N1
                    fac_EH(n) = fac_EH(n) + pEE(n1)*beta(n)*beta(n1)*abs(H(:,n)'*H(:,n1))^2;
                end
            end
            LHS_P = tI_k*wIpI_k/4*(tI/tI_k + wIpI/wIpI_k)^2 + tE_k*hpE_k/(4)*(tE/tE_k + hpE/hpE_k)^2;  % tE_K*hpE_k/(4*scal_fac^2)
            
            maximize (min_IR)
            subject to
            tE >= 0; tE <= 1; tI >= 0; tI <= 1;
            tE + tI <= 1;
            inv_tE >= 1; inv_tI >= 1;
            min_IR >= 0;
            [tI 1;1 inv_tI] == semidefinite(2);
            [tE 1;1 inv_tE] == semidefinite(2);
            for n = 1:N
                r_k(n) >= min_IR;
                [pI(n) 1;1 inv_pI(n)] == semidefinite(2);
                pI(n) >= 0;
                inv_pI(n) >= 0;
            end
            LHS_P <= P;
            for n = 1:N1
                pEE(n) >= 0;
                eta*fac_EH(n) >= e0*inv_tE;
            end
            cvx_end
            
            r = zeros(N,1);
            for n = 1:N
                int = 0;
                for n1 = 1:N
                    if (n1 == n)
                    else
                        int = int + abs(H(:,n)'*WI(:,n1))^2*pI(n1);
                    end
                end
                r(n) = tI*log(1 + (beta(n)*abs(H(:,n)'*WI(:,n))^2*pI(n))/( beta(n)*int + s2));
            end
            r = r/log(2);
            
            
            
            obj_fun = min_IR;
            if ( abs((obj_fun - obj_old)/obj_fun) < g_tol )
                loop_exit = 1;
            end
            obj_old = obj_fun;
            if (dispi == 1)
                display([ num2str(iter) ' obj_exp: ' num2str(min(r)) ' obj_algo: ' num2str(min_IR/log(2))  ' '  num2str(cvx_status) ]);
            end
            
            
        end
        
        display(['SIM # ' num2str(sim) '. rate_exp = ' num2str(min(r)) ' rate_alg = ' num2str(min_IR/log(2)) ' iter ' num2str(iter) ', status ' num2str(cvx_status) ' tI: ' num2str(tI) ' tE: ' num2str(tE)]);
        rate_s = rate_s + min(r);
        iter_s = iter_s + iter;
        tI_s = tI_s + tI;
        
        if (dispi==1)
            display(num2str('************************************************************************'))
        end
        
    end
    
    display(num2str('************************************************************************'))
    
    display([' Antennas ' num2str(M) ' avg_rate ' num2str(rate_s/simulations)  ' avg_iter ' num2str(iter_s/simulations) ' avg_tI ' num2str(tI_s/simulations) ]);
    display(num2str('************************************************************************'))
    display(num2str('************************************************************************'))
    rate_vec(M==M_vec) = rate_s/simulations;
end
plot(M_vec,rate_vec);
rate_vec



