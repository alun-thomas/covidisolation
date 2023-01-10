source("covid.R")
library(xtable,quietly=TRUE)
DMod = makeModel(deltaParams())
OPar = omicronParams()
OMod = makeModel(OPar)
PCRsens = 0.95
ANTIsens = 0.90
x = (0:100)/4

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
	Name = c("Symptomatic","Svrere")
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

x = (0:100)/4
pinf = 0.7
par(mfrow=c(2,2))
pileup(x,StateProbs(x,pinf,DMod))
title("Unconditional")
pileup(x,ConditionalStateProbs(x,pinf,DMod,NoSymptoms))
title("No symptoms")
pileup(x,ConditionalStateProbs(x,pinf,DMod,NoSymptoms*PcrTest(PCRsens)))
title("No symptoms and negative PCR test")
pileup(x,ConditionalStateProbs(x,pinf,DMod,NoSymptoms*AntigenTest(ANTIsens)))
title("No symptoms and negative antigen test")

par(mfrow=c(2,1))
x = (0:100)/4
#prior = c(1,9)/10
#prior = 5/10
prior = (1:9)/10
cols = colorRampPalette(c("blue","red"))(length(prior))

frame(x)
title("Probability that an infection occurred given negative antigen test")
for (i in 1:length(prior))
{
	p = prior[i]
	q = PostProbTransmission(x,p,AntigenTest(ANTIsens),DMod)
	lines(x,q,col=cols[i])
	z = x[q==min(q)]
	abline(v=z)
	text(z,0,paste(z))
}

frame(x)
title("Probability that an infection occurred given negative PCR test")
for (i in 1:length(prior))
{
	p = prior[i]
	q = PostProbTransmission(x,p,PcrTest(PCRsens),DMod)
	lines(x,q,col=cols[i])
	z = x[q==min(q)]
	abline(v=z)
	text(z,0,paste(z))
}

