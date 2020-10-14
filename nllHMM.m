function nLL=nllHMM(dat,tH0,tH1,piH0,piH1)

	nLL = double(-2*log(exp(vpa(fa_log(dat,tH1,piH1)))/...
			exp(vpa(fa_log(dat,tH0,piH0)))));