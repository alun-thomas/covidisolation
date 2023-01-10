library(coga,quietly=TRUE)

.Sus = 1
.Lat = 2
.Asy = 3
.Pre = 4
.Mld = 5
.Svr = 6
.RecS = 7
.RecA = 8
.NStates = .RecA

.U = 1
.V = 2
.W = 3
.X = 4
.Y = 5
.Z = 6

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

FF = function(x,P,use)
{
	pcoga(x,P$shape[use],P$rate[use])
}

StateProbs01 = function(x,M,transmission)
{
	s = matrix(0,nrow=.NStates,ncol=length(x))
	if (transmission)
	{
		s[.Sus,] = 0
		fv = FF(x,M,c(.V))
		s[.Lat,] = 1-fv
		fvx = FF(x,M,c(.V,.X))
		s[.Asy,] = (1-M$q)*(fv-fvx)
		fvw = FF(x,M,c(.V,.W))
		s[.Pre,] = M$q*(fv-fvw)
		fvwz = FF(x,M,c(.V,.W,.Z))
		s[.Mld,] = M$q*(1-M$r)*(fvw-fvwz)
		fvwy = FF(x,M,c(.V,.W,.Y))
		s[.Svr,] = M$q*M$r*(fvw-fvwy)
		s[.RecS,] = M$q*M$r*fvwy+M$q*(1-M$r)*fvwz
		s[.RecA,] = (1-M$q)*fvx
	}
	else
	{
		fu = FF(x,M,c(.U))
		s[.Sus,] = 1-fu
		fuv = FF(x,M,c(.U,.V))
		s[.Lat,] = fu-fuv
		fuvx = FF(x,M,c(.U,.V,.X))
		s[.Asy,] = (1-M$q)*(fuv-fuvx)
		fuvw = FF(x,M,c(.U,.V,.W))
		s[.Pre,] = M$q * (fuv - fuvw)
		fuvwz = FF(x,M,c(.U,.V,.W,.Z))
		s[.Mld,] = M$q*(1-M$r) * (fuvw - fuvwz)
		fuvwy = FF(x,M,c(.U,.V,.W,.Y))
		s[.Svr,] = M$q*M$r * (fuvw - fuvwy)
		s[.RecS,] = M$q*M$r * fuvwy + M$q*(1-M$r) * fuvwz
		s[.RecA,] = (1-M$q) * fuvx
	}
	s
}

StateProbs = function(x,p,M)
{
	p * StateProbs01(x,M,TRUE) + (1-p) * StateProbs01(x,M,FALSE)
}

Unconditional = c(1,1,1,1,1,1,1,1)
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

ProbInfected = function(x,p,observation,M)
{
	probs = StateProbs(x,p,M)
	bot = observation %*% probs
	top = (Infected * observation) %*% probs
	top/bot
}

PostProbTransmission = function(x,p,observation,M)
{
	q1 = observation %*% StateProbs01(x,M,T)
	q0 = observation %*% StateProbs01(x,M,F)
	q1 * p / (q1 * p + q0 * (1-p))
}

ConditionalStateProbs = function(x,p,M,observation)
{
	y = p * StateProbs01(x,M,TRUE) + (1-p) * StateProbs01(x,M,FALSE)
	for (j in 1:length(y[1,]))
	{
		y[,j] = y[,j]*observation
		y[,j] = y[,j]/sum(y[,j])
	}
	y
}

frame = function(x,xl="Days from exposure",yl="Probability")
{
	plot(x,x,type="n",ylim=c(0,1), ylab=yl, xlab=xl)
}

maximizer = function(x)
{
	(1:length(x))[x==max(x)][1]
}

skeletonPileup = function(x,ss,stateCols,stateNames)
{
	frame(x,yl="Cumulative probability")
	s = apply(ss,2,cumsum)
	lns = rev(1:length(s[,1]))
	for (i in lns)
	{
		lines(x,s[i,])
		polygon(c(x,max(x),0),c(s[i,],0,0),col=stateCols[i],border=NA)
	}

	for (i in lns)
	{
		zz = ss[i,]
		if (i > 1)
			z = s[i-1,]
		else
			z = rep(0,length(zz))

		if (max(zz) > 0)
		{
			whj = maximizer(zz)
			whx = x[whj]
			whx = min(c(whx,0.90*max(x)))
			whx = max(c(whx,0.05*max(x)))
			why = z[whj] + zz[whj]/2
			text(whx,why,stateNames[i])
		}
	}
}

pileUp = function(x,ss)
{
	stateCols = c("green","yellow","pink","orange","magenta","red","cyan","cyan")
	stateNames = c("Sus","Lat","Asy","Pre","Mld","Svr","RecS","RecA")
	skeletonPileup(x,ss,stateCols,stateNames)
}

pileup = function(x,sss,opt=2)
{
	ss = matrix(0,ncol=ncol(sss),nrow=nrow(sss)-2)

	if (opt == 2)
	{
		stateCols = c("green","yellow","pink","orange","red","cyan")
		stateNames = c("Sus","Lat","Asy","Pre","Sym","Rec")
		ss[1:4,] = sss[1:4,]
		ss[5,] = sss[5,] + sss[6,]
		ss[6,] = sss[7,] + sss[8,]
	}

	if (opt == 1)
	{
		stateCols = c("yellow","pink","orange","red","cyan","green")
		stateNames = c("Lat","Asy","Pre","Sym","Rec","Sus")
		ss[1:3,] = sss[2:4,]
		ss[4,] = sss[5,] + sss[6,]
		ss[5,] = sss[7,] + sss[8,]
		ss[6,] = sss[1,]
	}

	skeletonPileup(x,ss,stateCols,stateNames)
}

