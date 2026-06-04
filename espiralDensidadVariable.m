function [x,y,theta,r] = espiralDensidadVariable(rmax,nVueltas,N,alpha)

theta = linspace(0, 2*pi*nVueltas, N);
thetaMax = max(theta);

r = rmax * (theta/thetaMax).^alpha;

x = r .* cos(theta);
y = r .* sin(theta);

end