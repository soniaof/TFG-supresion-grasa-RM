function [k, kpcm] = g2k(g, fov, dt, gambar)
% kTraj from Gradient, modified from Hao's g2k_hao().
% INPUT
%  g,   [Nt, Nd], G/cm
%  fov, [Nd]    , cm
%  dt,  [1]     , optional, mSec
%  gambar, [1]  , optional, 2pi*kHz/T
% OUTPUT
%  k,   [Nt, Nd], cycle/fov
%  kpcm,[Nt, Nd], cycle/cm

if ~exist('dt',     'var'), dt = 10e-3; end % unit in mSec
if ~exist('gambar', 'var'), gambar = 2*pi*42.576e3; end % 2pi * kHz/T

if sum(isnan(g(:)))>0
    g(isnan(g))=0;
end
dtSec = dt*1e-3;        % convert to Sec
gamma = gambar/2/pi/10; % convert to Hz/Gauss

kpcm = -flipud(cumsum(flipud(g),1) * dtSec * gamma);
k = bsxfun(@times, kpcm, fov);

end