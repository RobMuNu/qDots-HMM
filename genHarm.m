function h=genHarm(n,a)
% Generates nth Harmonic number of degree a
% Surprisingly matlabs harmonic function doesn't allow for degrees ~= 1

if a<1
	error('Generalised harmonic number undefined for a < 1')
	return;
else
	h=0;
for w =1:n
	h = h + w^(-a);
end
end




% If you're running MATLAB 2020a you can try this
% H(n,a) = hurwitzZeta(a,1) - hurwitzZeta(a,n+1)
% Although this script is very fast so probably wont need it