function ensure_baseline(user_dir)
% ENSURE_BASELINE  Put the DD-BS baseline code on the MATLAB path, wherever it lives.
%
%   ensure_baseline()            % auto-locate
%   ensure_baseline('C:\path\to\code_nf_distance_dependent_rainbow')
%
% The probes need the baseline helpers (near_field_channel, cal_loc, TTD_beam,
% delay_polar_2d, ...). Hard-coding a relative path breaks as soon as the
% baseline lives somewhere else, so resolve it properly instead.

if exist('near_field_channel','file')==2, return; end          % already usable

if nargin>=1 && ~isempty(user_dir)
    if exist(fullfile(user_dir,'near_field_channel.m'),'file')
        addpath(genpath(user_dir)); return;
    end
    error('ensure_baseline:notThere','near_field_channel.m is not in:\n  %s', user_dir);
end

here = fileparts(mfilename('fullpath'));
cand = { fullfile(here,'..','..','baseline_code','code_nf_distance_dependent_rainbow'), ...
         fullfile(here,'..','..','code_nf_distance_dependent_rainbow'), ...
         fullfile(here,'..','baseline_code','code_nf_distance_dependent_rainbow'), ...
         fullfile(here,'code_nf_distance_dependent_rainbow') };
for i=1:numel(cand)
    if exist(fullfile(cand{i},'near_field_channel.m'),'file')
        addpath(genpath(cand{i}));
        fprintf('[ensure_baseline] using %s\n', cand{i}); return;
    end
end

% bounded search upward (3 parents) for a folder containing the helper
p = here;
for up = 1:3
    p = fileparts(p);
    if isempty(p), break; end
    hits = dir(fullfile(p,'**','near_field_channel.m'));
    if ~isempty(hits)
        addpath(genpath(hits(1).folder));
        fprintf('[ensure_baseline] using %s\n', hits(1).folder); return;
    end
end

error('ensure_baseline:notFound', ...
 ['Could not locate the DD-BS baseline code (near_field_channel.m).\n' ...
  'Either run this once in MATLAB:\n' ...
  '    addpath(genpath(''<path to code_nf_distance_dependent_rainbow>''))\n' ...
  'or call the probe helper directly:\n' ...
  '    ensure_baseline(''<path to code_nf_distance_dependent_rainbow>'')']);
end
