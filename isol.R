source("covid.R")

library(xtable,quietly=TRUE)

frame=function(x)
{
	plot(x,x,type="n",ylim=c(0,1),
		ylab="Probability",
		xlab="Days from exposure")
}

DMod = makeModel(deltaParams())
OPar = omicronParams()
OMod = makeModel(OPar)

stateCols = c("cyan","yellow","pink","orange","red","red","cyan","cyan")

pileup = function(x,s)
{
	frame(x)
	for (i in 1:length(s[1,]))
		s[,i] = cumsum(s[,i])
	lns = rev((1:length(s[,1]))[-c(5,7)])
	for (i in lns)
	{
		lines(x,s[i,])
		polygon(c(x,max(x),0),c(s[i,],0,0),col=stateCols[i],border=NA)
	}
}

print(xtable(gammas(DMod),
	caption="Gamma time-in-state distribution parameters and summary
		 values for original \\covid variants.\\bigskip",
	digits=c(0,0,0,2,3,2,2,2,2),
	label="deltagammas"),
	include.rownames=FALSE,
	caption.placement="top"
	)

print(xtable(probs(DMod),
	caption="Path probability parameters for original \\covid variants.",
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

pileup(x,StateProbs(x,1,DMod))

pileup(x,ConditionalStateProbs(x,1,DMod,NoSymptoms))

pileup(x,ConditionalStateProbs(x,1,DMod,NoSymptoms*PcrTest(0.99)))

pileup(x,ConditionalStateProbs(x,1,DMod,NoSymptoms*AntigenTest(0.9)))

x = (0:100)/4
frame(x)
#for (p in (1:9)/10)
for (p in c(1,9)/10)
{
	q = PostProbTransmission(x,p,AntigenTest(0.9),DMod)
	lines(x,q)
	z = x[q==min(q)]
	abline(v=z)
	text(z,0,paste(z))
}

x = (0:100)/4
frame(x)
#for (p in (1:9)/10)
for (p in c(1,9)/10)
{
	q = PostProbTransmission(x,p,PcrTest(0.99),DMod)
	lines(x,q)
	z = x[q==min(q)]
	abline(v=z)
	text(z,0,paste(z))
}

