function rate = training_ddbs_e2e(Nt,B,fc,f,M,d,h,w_train,focus_loc,SNR_t,SNR_dB,Q,K,lsb_sec,sigma_tau,use_hybrid)
% TRAINING_DDBS_E2E  Baseline DDBS training, but with the SERVING beam impaired too.
%
% use_hybrid selects the SERVING-beam architecture; it must match the training
% architecture, otherwise a pure-TD serving beam becomes the bottleneck and hides
% the benefit of a hybrid training beam.
%
% Mirrors baseline training_near_rainbow_2d (argmax subcarrier -> its intended
% focus -> serve), with one correction: the data beam is built with the same
% impaired hardware (ttd_beam_impaired) instead of an ideal TTD_beam. Without
% this, impairments that degrade the beam but preserve the argmax look harmless.

if nargin < 16 || isempty(use_hybrid), use_hybrid = false; end
rate = zeros(K,1);
for idx = 1:K
    y = zeros(M,Q);
    for m = 1:M
        y(m,:) = awgn( repmat( h(m,:)*w_train(:,idx,m), [1,Q] ), SNR_dB*2/sqrt(3) );
    end
    y = abs(sum(y,2)).^2;
    [~,i] = max(y);
    w_serve = ttd_beam_impaired(Nt,B,fc,M,d, focus_loc(i,1,idx), focus_loc(i,2,idx), lsb_sec, sigma_tau, use_hybrid);
    t = 0;
    for m = 1:M
        t = t + log2(1 + SNR_t*abs(h(m,:)*w_serve(:,m))^2)/M;
    end
    if idx==1, rate(idx)=t; else, rate(idx)=max(rate(idx-1),t); end
end
end
