function y = add_awgn(x, snr_db)
% ADD_AWGN  Vectorised replacement for the baseline's per-subcarrier awgn call.
%
% The baseline (and every probe here) writes
%       for m = 1:M,  y(m,:) = awgn(repmat(h(m,:)*w(:,t,m),[1,Q]), snr_db);  end
% which is ~M*K*N_iter scalar awgn calls -- 2.5M of them at K=12, N_iter=200 --
% and the function-call overhead, not the arithmetic, is what makes N_iter=200
% impractical.
%
% The semantics are simpler than they look. MATLAB's two-argument awgn(x, snr)
% does NOT measure the signal: it assumes the power of x is 0 dBW. So the noise
% variance is a CONSTANT 10^(-snr/10) independent of the sample, and the whole
% loop collapses to one vectorised draw with an identical distribution. (Had the
% baseline written awgn(..., 'measured') this would be wrong, because the variance
% would then track each subcarrier's own power -- that is why this is spelled out
% rather than assumed.)
%
% Only the RNG draw ORDER changes, so numbers move by Monte-Carlo scatter, not by
% construction. The self-check below confirms the variance once per session.
persistent checked
nv = 10^(-snr_db/10);                       % awgn's 0-dBW-assumed noise variance
if isempty(checked)
    t = zeros(1,200000);
    a = awgn(t, snr_db);                    % awgn on a zero signal -> pure noise
    b = sqrt(nv/2)*(randn(size(t))+1j*randn(size(t)));
    fprintf('[add_awgn] self-check: var(awgn) = %.4e, var(ours) = %.4e, target %.4e\n', ...
            var(a(:)), var(b(:)), nv);
    if abs(var(a(:))/nv - 1) > 0.05
        warning('add_awgn:selftest','awgn variance does not match the 0 dBW assumption');
    end
    checked = true;
end
y = x + sqrt(nv/2)*(randn(size(x)) + 1j*randn(size(x)));
end
