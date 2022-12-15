library(coga,quietly=TRUE)

deltaParams = function()
{
	msus = 1000

	list (	Psymp = 0.5,
		Psevere = 0.1,
		Msus = msus,
		Vsus = msus^2,
		Mlat = 5.1,
		Vlat = 7.5,
		Masy = 7,
		Vasy = 7.5,
		Mpre = 3,
		Vpre = 1.0,
		Mmld = 5,
		Vmld = 45.5,
		Msev = 11,
		Vsev = 45.5
	)
}

omicronParams = function()
{
	ommod = 0.6
	ommod2 = 0.9
	msus = 1000

	list( Psymp = 0.5,
		Psevere = 0.1,
		Msus = msus,
		Vsus = msus^2,
		Mlat = 5.1 * ommod,
		Vlat = 7.5 * ommod^2,
		Masy = 7 * ommod2,
		Vasy = 7.5 * ommod2^2,
		Mpre = 3 * ommod2,
		Vpre = 1.0 * ommod2^2,
		Mmld = 5 * ommod2,
		Vmld = 45.5 * ommod2^2,
		Msev = 11 * ommod2,
		Vsev = 45.5 * ommod2^2,
		ommod = ommod,
		ommod2 = ommod2 
	)
}

makeModel = function(Q)
{
	means = c(Q$Msus,Q$Mlat,Q$Mpre,Q$Masy,Q$Msev,Q$Mmld)
	vars = c(Q$Vsus,Q$Vlat,Q$Vpre,Q$Vasy,Q$Vsev,Q$Vmld)
	list(q=Q$Psymp, r=Q$Psevere, shape=means^2/vars, rate=means/vars) 
}

gammas = function(Q)
{
	Name = c("Sus","Lat","Asy","Pre","Mld","Sev")
	Parameter = c("U","V","W","X","Y","Z")
	Shape = Q$shape
	Rate = Q$rate
	Mean = Shape/Rate
	Variance = Shape/Rate^2
	tail = 0.025
	Lower = qgamma(tail,Shape,Rate)
	Upper = qgamma(1-tail,Shape,Rate)
	data.frame(Parameter,Name,Shape,Rate,Mean,Variance,Lower,Upper)
}

probs = function(Q)
{
	Name = c("Symptomatic","Severe")
	Parameter = c("q","r")
	Value = c(Q$q,Q$r)
	data.frame(Parameter,Name,Value)
}

FF = function(x,P,use)
{
	pcoga(x,P$shape[use],P$rate[use])
}

StateProbs01 = function(x,M,transmission)
{
	Sus = 1; Lat = 2; Asy = 3; Pre = 4; Mld = 5; Sev = 6; RecS = 7; RecA = 8
	NStates = RecA
	U = 1; V = 2; W = 3; X = 4; Y = 5; Z = 6

	s = matrix(0,nrow=NStates,ncol=length(x))
	if (transmission)
	{
		s[Sus,] = 0
		s[Lat,] = 1 - FF(x,M,c(V))
		s[Asy,] = (1-M$q) * (FF(x,M,c(V)) - FF(x,M,c(V,X)))
		s[Pre,] = M$q * (FF(x,M,c(V)) - FF(x,M,c(V,W)))
		s[Mld,] = M$q*(1-M$r) * (FF(x,M,c(V,W)) - FF(x,M,c(V,W,Z)))
		s[Sev,] = M$q*M$r * (FF(x,M,c(V,W)) - FF(x,M,c(V,W,Y)))
		s[RecS,] = M$q*M$r * FF(x,M,c(V,W,Y)) + M$q*(1-M$r)*FF(x,M,c(V,W,Z))
		s[RecA,] = (1-M$q) * FF(x,M,c(V,X))
	}
	else
	{
		s[Sus,] = 1-FF(x,M,c(U))
		s[Lat,] = FF(x,M,c(U)) - FF(x,M,c(U,V))
		s[Asy,] = (1-M$q) * (FF(x,M,c(U,V)) - FF(x,M,c(U,V,X)))
		s[Pre,] = M$q * (FF(x,M,c(U,V)) - FF(x,M,c(U,V,W)))
		s[Mld,] = M$q*(1-M$r) * (FF(x,M,c(U,V,W)) - FF(x,M,c(U,V,W,Z)))
		s[Sev,] = M$q*M$r * (FF(x,M,c(U,V,W)) - FF(x,M,c(U,V,W,Y)))
		s[RecS,] = M$q*M$r * FF(x,M,c(U,V,W,Y)) + M$q*(1-M$r)*FF(x,M,c(U,V,W,Z))
		s[RecA,] = (1-M$q) * FF(x,M,c(U,V,X))
	}
	s
}

StateProbs = function(x,p,M)
{
	p * StateProbs01(x,M,TRUE) + (1-p) * StateProbs01(x,M,FALSE)
}

NoSymptoms = c(1,1,1,1,0,0,0,1)

Infected = c(0,1,1,1,1,1,0,0)

PcrTest = function(sens,spec=1)
{
	c(spec,spec,1-sens,1-sens,1-sens,1-sens,1-sens,1-sens)
}

AntigenTest = function(sens,spec=1)
{
	c(spec,spec,1-sens,1-sens,1-sens,1-sens,spec,spec)
}

ProbInfected = function(x,p,t,M)
{
	probs = StateProbs(x,p,M)
	bot = (NoSymptoms * t) %*% probs
	top = (Infected * NoSymptoms * t) %*% probs
	top/bot
}

PostProbTransmission = function(x,p,t,M)
{
	q1 = (NoSymptoms * t) %*% StateProbs01(x,M,T)
	q0 = (NoSymptoms * t) %*% StateProbs01(x,M,F)
	q1 * p / (q1 * p + q0 * (1-p))
}

ConditionalStateProbs = function(x,p,M,t)
{
	y = p * StateProbs01(x,M,TRUE) + (1-p) * StateProbs01(x,M,FALSE)
	for (j in 1:length(y[1,]))
	{
		y[,j] = y[,j]*t
		y[,j] = y[,j]/sum(y[,j])
	}
	y
}