function [x,y] = espiralArquimedes(a,b,nVueltas,N)
% Espiral de Arquímedes: r = a + b*theta
%
% a        = radio inicial
% b        = separación de la espiral
% nVueltas = número de vueltas
% N        = número de puntos

theta = linspace(0, 2*pi*nVueltas, N);

r = a + b*theta;

x = r .* cos(theta);
y = r .* sin(theta);

end

