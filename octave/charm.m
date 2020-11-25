% very simple script to invert matrix octave

A = load('matrix.mat');
b = load('vector.mat');
X = linsolve(A,b);

fprintf('Solution:');
fprintf('\n');
fprintf('%d ',X);
fprintf('\n');

close();
