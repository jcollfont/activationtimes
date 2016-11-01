%% HELP:
%		
%		This function is a wrapper for the detection of activation times
%		from heart surface potentials.
%		Underneath it is using the spatio-temporal activation times method
%		and a further smoothing.
%
%
%		INPUT:
%			- EGM - <M,T>double - heart surface potentials from which to
%								obtain activaion times.
%			- Dtan - <C1,M>double - tangential derivative operator matrix.
%			- Ltan - <C2,M>double - tangential laplacian operator matrix.
%
%		OUTPUT:
%			- acttimes - <M,1>double - estimated timing of deploarization at
%									each node on the heart.
%			
%		AUTHOR:
%			Jaume Coll-Font <jcollfont@gmail.com>
%			Using code from Burak Erem
%
%

function [acttimes] = activationTimes_wrapper(EGM,Dtan,Ltan, alpha)
	
	% parameters
	twindow = 9;
	smoothingLambda = 10.^linspace(-10,6,1500);

	% activation times estimation
	acttimes = spatiotemporalActtimes(EGM,Dtan,twindow, alpha);

	% smooth activaation times
% 	acttimes = smoothactivationtimes( Ltan, acttimes, smoothingLambda);
	
end