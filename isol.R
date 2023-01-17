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

pinf = 0.7

x = (0:100)/4
par(mfrow=c(2,2))
pileup(x,StateProbs(x,pinf,DMod))
title("Unconditional")
pileup(x,ConditionalStateProbs(x,pinf,DMod,NoSymptoms))
title("No symptoms")
pileup(x,ConditionalStateProbs(x,pinf,DMod,NoSymptoms*AntigenTest(ANTIsens)))
title("No symptoms and negative antigen test")
pileup(x,ConditionalStateProbs(x,pinf,DMod,NoSymptoms*PcrTest(PCRsens)))
title("No symptoms and negative PCR test")

x = (0:100)/4
par(mfrow=c(2,2))
pileup(x,StateProbs(x,pinf,OMod))
title("Unconditional")
pileup(x,ConditionalStateProbs(x,pinf,OMod,NoSymptoms))
title("No symptoms")
pileup(x,ConditionalStateProbs(x,pinf,OMod,NoSymptoms*AntigenTest(ANTIsens)))
title("No symptoms and negative antigen test")
pileup(x,ConditionalStateProbs(x,pinf,OMod,NoSymptoms*PcrTest(PCRsens)))
title("No symptoms and negative PCR test")

x = (0:100)/4
#prior = c(1,9)/10
prior = (1:9)/10
cols = colorRampPalette(c("blue","red"))(length(prior))
par(mfcol=c(2,2))

mod = DMod
frame(x)
title("Delta variant antigen test")
for (i in 1:length(prior))
{
	p = prior[i]
	q = PostProbTransmission(x,p,AntigenTest(ANTIsens),mod)
	lines(x,q,col=cols[i])
	z = x[q==min(q)]
	abline(v=z)
	text(z,0,paste(z))
}

frame(x)
title("Delta variant PCR test")
for (i in 1:length(prior))
{
	p = prior[i]
	q = PostProbTransmission(x,p,PcrTest(PCRsens),mod)
	lines(x,q,col=cols[i])
	z = x[q==min(q)]
	abline(v=z)
	text(z,0,paste(z))
}

mod = OMod

frame(x)
title("Omicron variant antigen test")
for (i in 1:length(prior))
{
	p = prior[i]
	q = PostProbTransmission(x,p,AntigenTest(ANTIsens),mod)
	lines(x,q,col=cols[i])
	z = x[q==min(q)]
	abline(v=z)
	text(z,0,paste(z))
}

frame(x)
title("Omicron variant PCR test")
for (i in 1:length(prior))
{
	p = prior[i]
	q = PostProbTransmission(x,p,PcrTest(PCRsens),mod)
	lines(x,q,col=cols[i])
	z = x[q==min(q)]
	abline(v=z)
	text(z,0,paste(z))
}

x = (0:100)/4
par(mfrow=c(2,1))

mod = DMod
frame(x)
title("Delta variant")
base=ProbInfected(x,0,Unconditional,mod)
lines(x,base,col="green")

y=ProbInfected(x,1,Unconditional,mod)
lines(x,y,col="black")
#z=min(x[y<base])
#abline(v=z,col="black")
#text(z,0,paste(z))

y=ProbInfected(x,1,NoSymptoms,mod)
lines(x,y,col="blue")
z=min(x[y<base])
abline(v=z,col="blue")
text(z,0,paste(z))

y=ProbInfected(x,1,NoSymptoms*AntigenTest(ANTIsens),mod)
lines(x,y,col="red")
z=min(x[y<base])
abline(v=z,col="red")
text(z,0,paste(z))

y=ProbInfected(x,1,NoSymptoms*PcrTest(PCRsens),mod)
lines(x,y,col="orange")
z=min(x[y<base])
abline(v=z,col="orange")
#text(z,0,paste(z))

mod = OMod
frame(x)
title("Omicron variant")
base=ProbInfected(x,0,Unconditional,mod)
lines(x,base,col="green")

y=ProbInfected(x,1,Unconditional,mod)
lines(x,y,col="black")
#z=min(x[y<base])
#abline(v=z,col="black")
#text(z,0,paste(z))

y=ProbInfected(x,1,NoSymptoms,mod)
lines(x,y,col="blue")
z=min(x[y<base])
abline(v=z,col="blue")
text(z,0,paste(z))

y=ProbInfected(x,1,NoSymptoms*AntigenTest(ANTIsens),mod)
lines(x,y,col="red")
z=min(x[y<base])
abline(v=z,col="red")
text(z,0,paste(z))

y=ProbInfected(x,1,NoSymptoms*PcrTest(PCRsens),mod)
lines(x,y,col="orange")
z=min(x[y<base])
abline(v=z,col="orange")
#text(z,0,paste(z))

