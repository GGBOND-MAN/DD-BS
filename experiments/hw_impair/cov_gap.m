function [gap, cov] = cov_gap(fl, Nt, lim)
% COV_GAP  Largest coverage hole of a recalibrated focus map, in array beamwidths.
%
%   [gap, cov] = cov_gap(fl, Nt)          % lim defaults to sin(60 deg)
%
% BUG FIXED HERE (found by pilot_spacing_map): the earlier inline version took
% max(diff(sorted_foci)) over the foci alone, so a set that is DENSE but NARROW --
% a tight clump covering only part of the region -- reported a tiny gap and passed
% a criterion it should have failed. Two configurations exposed it:
%     K=4, Delta=0.125 : gap 0.6 "pass"  but rate 41%
%     K=6, Delta=0.125 : gap 0.3 "pass"  but rate 51%
% Both have comb span K*Delta = 0.5-0.75, far short of the 2 needed to tile, so
% their foci pile up in a narrow band with nothing outside it -- and a gap taken
% only BETWEEN foci cannot see the empty region beyond the outermost one.
%
% The fix is to bound the interval with the region edges before differencing, so
% the distance from the extreme focus to the edge counts as the hole it is.
if nargin < 3, lim = 0.866; end
sn = sort(reshape(fl(:,1,:), [], 1));
sn = sn(sn >= -lim & sn <= lim);
gap = max(diff([-lim; sn; lim])) / (2/Nt);      % sentinels: see above
if isempty(sn), cov = 0; else, cov = (max(sn)-min(sn))/(2*lim); end
end
