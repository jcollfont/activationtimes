%% HELP:
%
%		This function estiamtes the point of minimum dvdt of the input
%		signal (assumed to be a matrix with columns representing time).
%		
%		It does so by locally approximating the signal with a deg^th order
%		polynomial and computing the dvdt on it.
%
%		INPUT:
%			- sig			- <M,T>double - signal of interest. Time coordinate are the
%										columns.
%			- winLength		- int - number of samples to be used to compute the
%								local dth order polynomial.
%			- deg			- int - order of the polynomial approximation.
%

function [tau, dy] = findMinDVDT(sig, winLength, deg)

	[M,T] = size(sig);

	tau = zeros(1,M);
	dy = zeros(M,T);
	
	for mm = 1:M
		cen = ceil(winLength/2);
		X = zeros(winLength,(deg+1));
		L = [-(cen-1):(cen-1)]'; for p=1:(deg+1), X(:,p) = L.^((deg+1)-p); end
		E = inv(X'*X)*X';
		temp = [sig(mm,:) sig(mm,end)*ones(1,cen-1)];
		a = filter(E(deg,[winLength:-1:1]),1,temp );
		dy(mm,:) = a(cen:end);

		tau(mm) = min(dy(mm,:));
		
	end
	
end