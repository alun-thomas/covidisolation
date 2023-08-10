source("covid.R")
library(xtable,quietly=TRUE)
DMod = makeModel(deltaParams())
OPar = omicronParams()
OMod = makeModel(OPar)
PCRsens = 0.95
ANTIsens = 0.90
x = (0:100)/4

pc = function(x,d=0)
{
	round(x*100,d)
}

# Function to make table of Gamma parameters.
gammas = function(Q)
{
	Name = c("Sus","Lat","Asy","Pre","Mld","Svr")
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

# Function to make table of path splitting parameters.
probs = function(Q)
{
	Name = c("Symptomatic","Severe")
	Parameter = c("q","r")
	Value = c(Q$q,Q$r)
	data.frame(Parameter,Name,Value)
}

print(xtable(gammas(DMod),
	caption="Gamma time-in-state distribution parameters and summary
		 values for the original \\covid variants.\\bigskip",
	digits=c(0,0,0,2,3,2,2,2,2),
	label="deltagammas"),
	include.rownames=FALSE,
	caption.placement="top"
	)

print(xtable(probs(DMod),
	caption="Probability parameter values specifying the path through the 
		course of infection for the original \\covid variants.\\bigskip",
	label="deltaprobs"),
	include.rownames=FALSE,
	caption.placement="top"
	)

print(xtable(gammas(OMod),
	caption="Gamma time-in-state distribution parameters and summary
		values for more recent \\covid variants.\\bigskip",
	digits=c(0,0,0,2,3,2,2,2,2),
	label="omicrongammas"),
	include.rownames=FALSE,
	caption.placement="top"
	)

#pinf = 0.7
pinf = 0.3
pacc = 0.1
pacc2 = 0.01

clearpoint = function(x,probs,q=0.1)
{
	p = apply(probs[c(.Lat,.Asy,.Pre,.Mld,.Svr),],2,sum)
	t = min(x[p<=q])
	abline(v=t)
	text(t,0,paste(t))
	t
}


x = (0:100)/4
par(mfrow=c(2,2))
probs = StateProbs(x,pinf,DMod)
pileup(x,probs)
zz = clearpoint(x,probs,pacc)
title("(a) Unconditional")
pileup(x,ConditionalStateProbs(probs,NoSymptoms))
zz = clearpoint(x,ConditionalStateProbs(probs,NoSymptoms),pacc)
title("(b) No symptoms")
pileup(x,ConditionalStateProbs(probs,NoSymptoms*AntigenTest(ANTIsens)))
whendanti = clearpoint(x,ConditionalStateProbs(probs,NoSymptoms*AntigenTest(ANTIsens)),pacc)
title("(c) No symptoms and neg antigen test")
pileup(x,ConditionalStateProbs(probs,NoSymptoms*PcrTest(PCRsens)))
whendpcr = clearpoint(x,ConditionalStateProbs(probs,NoSymptoms*PcrTest(PCRsens)),pacc)
title("(d) No symptoms and neg PCR test")

x = (0:100)/4
par(mfrow=c(2,2))
probs = StateProbs(x,pinf,OMod)
pileup(x,probs)
zz = clearpoint(x,probs,pacc)
title("(a) Unconditional")
pileup(x,ConditionalStateProbs(probs,NoSymptoms))
zz = clearpoint(x,ConditionalStateProbs(probs,NoSymptoms),pacc)
title("(b) No symptoms")
pileup(x,ConditionalStateProbs(probs,NoSymptoms*AntigenTest(ANTIsens)))
whenoanti = clearpoint(x,ConditionalStateProbs(probs,NoSymptoms*AntigenTest(ANTIsens)),pacc)
title("(c) No symptoms and neg antigen test")
pileup(x,ConditionalStateProbs(probs,NoSymptoms*PcrTest(PCRsens)))
whenopcr = clearpoint(x,ConditionalStateProbs(probs,NoSymptoms*PcrTest(PCRsens)),pacc)
title("(d) No symptoms and neg PCR test")

x = (0:100)/4
#prior = c(1,9)/10
prior = (1:9)/10
cols = colorRampPalette(c("blue","red"))(length(prior))
par(mfcol=c(3,2))

mod = DMod
P0 = StateProbs(x,0,mod)
P1 = StateProbs(x,1,mod)

frame(x)
title("(a) Delta. No symptoms. No test.")
for (i in 1:length(prior))
{
	q = PostProbTransmission(P0,P1,NoSymptoms,prior[i])
	lines(x,q,col=cols[i])
	z = x[q==min(q)]
	abline(v=z)
	text(z,0,paste(z))
}

frame(x)
title("(c) Delta. No symptoms. Antigen neg.")
for (i in 1:length(prior))
{
	q = PostProbTransmission(P0,P1,NoSymptoms*AntigenTest(ANTIsens),prior[i])
	lines(x,q,col=cols[i])
	z = x[q==min(q)]
	abline(v=z)
	text(z,0,paste(z))
}

frame(x)
title("(e) Delta. No symptoms. PCR neg.")
for (i in 1:length(prior))
{
	q = PostProbTransmission(P0,P1,NoSymptoms*PcrTest(PCRsens),prior[i])
	lines(x,q,col=cols[i])
	z = x[q==min(q)]
	abline(v=z)
	text(z,0,paste(z))
}

mod = OMod
P0 = StateProbs(x,0,mod)
P1 = StateProbs(x,1,mod)

frame(x)
title("(b) Omicron. No symptoms. No test.")
for (i in 1:length(prior))
{
	q = PostProbTransmission(P0,P1,NoSymptoms,prior[i])
	lines(x,q,col=cols[i])
	z = x[q==min(q)]
	abline(v=z)
	text(z,0,paste(z))
}

frame(x)
title("(d) Omicron. No symptoms. Antigen neg.")
for (i in 1:length(prior))
{
	q = PostProbTransmission(P0,P1,NoSymptoms*AntigenTest(ANTIsens),prior[i])
	lines(x,q,col=cols[i])
	z = x[q==min(q)]
	abline(v=z)
	text(z,0,paste(z))
}

frame(x)
title("(f) Omicron. No symptoms. PCR neg.")
for (i in 1:length(prior))
{
	q = PostProbTransmission(P0,P1,NoSymptoms*PcrTest(PCRsens),prior[i])
	lines(x,q,col=cols[i])
	z = x[q==min(q)]
	abline(v=z)
	text(z,0,paste(z))
}

x = (0:140)/4
par(mfrow=c(2,1))

mod = DMod
frame(x)
title("(a) Delta variant")

P0 = StateProbs(x,0,mod)
P1 = StateProbs(x,1,mod)

base=ProbInfected(P0,Unconditional)
lines(x,base,col="green")

y=ProbInfected(P1,Unconditional)
lines(x,y,col="black")
z=min(x[y<base])
abline(v=z,col="black")
text(z,0,paste(z))

y=ProbInfected(P1,NoSymptoms)
lines(x,y,col="blue")
z=min(x[y<base])
abline(v=z,col="blue")
text(z,0,paste(z))

y=ProbInfected(P1,NoSymptoms*AntigenTest(ANTIsens))
lines(x,y,col="red")
z=min(x[y<base])
abline(v=z,col="red")
text(z,0,paste(z))

y=ProbInfected(P1,NoSymptoms*PcrTest(PCRsens))
lines(x,y,col="orange")
z=min(x[y<base])
abline(v=z,col="orange")
#text(z,0,paste(z))

mod = OMod
frame(x)
title("(b) Omicron variant")
P0 = StateProbs(x,0,mod)
P1 = StateProbs(x,1,mod)

base=ProbInfected(P0,Unconditional)
lines(x,base,col="green")

y=ProbInfected(P1,Unconditional)
lines(x,y,col="black")
z=min(x[y<base])
abline(v=z,col="black")
text(z,0,paste(z))

y=ProbInfected(P1,NoSymptoms)
lines(x,y,col="blue")
z=min(x[y<base])
abline(v=z,col="blue")
text(z,0,paste(z))

y=ProbInfected(P1,NoSymptoms*AntigenTest(ANTIsens))
lines(x,y,col="red")
z=min(x[y<base])
abline(v=z,col="red")
text(z,0,paste(z))

y=ProbInfected(P1,NoSymptoms*PcrTest(PCRsens))
lines(x,y,col="orange")
z=min(x[y<base])
abline(v=z,col="orange")
#text(z,0,paste(z))

CondSimStateProbs = function(x,p,M,t,sel,sims=100000)
{
        s = matrix(0,nrow=.NStates,ncol=length(x))
	wh = min((1:length(x))[x>=t])
	denom = 0
        for (i in 1:sims)
	{
		rc = asStateMatrix(rCovid(x,p,M))
		pacc = sum(sel*rc[,wh])
#		if (runif(1) < pacc)
#		{
#			denom = denom + 1
#                	s = s+rc
#		}
		denom = denom + pacc
		s = s + pacc*rc
	}
        s/denom
}

clearpoint2 = function(x,p,cl=0,q=pacc2)
{
	t = min(x[p<=q])
	abline(v=t,col=cl)
#	text(t,0,paste(t))
}

plotem = function(x,t,PP)
{
	wh = x>=t
	QQ=ProbInfected(PP,Unconditional)
	lines(x[wh],QQ[wh])
	clearpoint2(x,QQ,cl="black")
	QQ=ProbInfected(PP,NoSymptoms)
	lines(x[wh],QQ[wh],col="blue")
	clearpoint2(x,QQ,cl="blue")
	QQ=ProbInfected(PP,NoSymptoms*PcrTest(PCRsens))
	lines(x[wh],QQ[wh],col="orange")
	clearpoint2(x,QQ,cl="orange")
	QQ=ProbInfected(PP,NoSymptoms*AntigenTest(ANTIsens))
	lines(x[wh],QQ[wh],col="red")
	clearpoint2(x,QQ,cl="red")
}

set.seed(1001)
par(mfrow=c(2,2))
x = (0:100)/4
ymx = 0.10

frame(x,ymax=ymx)
t = whendanti
PP=CondSimStateProbs(x,pinf,DMod,t,NoSymptoms*AntigenTest(ANTIsens))
plotem(x,t,PP)
title("(a) Delta. Antigen test first.")

frame(x,ymax=ymx)
t = whendpcr
PP = CondSimStateProbs(x,pinf,DMod,t,NoSymptoms*PcrTest(PCRsens))
plotem(x,t,PP)
title("(b) Delta. PCR test first.")

frame(x,ymax=ymx)
t = whenoanti
PP = CondSimStateProbs(x,pinf,OMod,t,NoSymptoms*AntigenTest(ANTIsens))
plotem(x,t,PP)
title("(c) Omicron. Antigen test first.")

frame(x,ymax=ymx)
t = whenopcr
PP = CondSimStateProbs(x,pinf,OMod,t,NoSymptoms*PcrTest(PCRsens))
plotem(x,t,PP)
title("(d) Omicron. PCR test first.")

