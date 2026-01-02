##############################
# calculating PM
######################
#important variables
#exposure (X)
#outcomes (Y)
#mediator (M)

TE.X2Y=AAM2CIMT.selected$b
var.X2Y=AAM2CIMT.selected$se^2

E.X2M=AAM2BMI.selected$b
var.X2M=AAM2BMI.selected$se^2

E.M2Y=BMI2CIMT.selected$b
var.M2Y=BMI2CIMT.selected$se^2

PM=(E.X2M*E.M2Y)/TE.X2Y

library(car)
var.matrix=matrix(data=c(var.X2Y,0,0,
                         0,var.X2M,0,
                         0,0,var.M2Y),
                  ncol=3,nrow=3)

input=c(TE.X2Y,E.X2M,E.M2Y)
names(input)=c("TE.X2Y","E.X2M","E.M2Y")

summary.PM=deltaMethod(object=input,"(E.X2M*E.M2Y)/TE.X2Y",vcov.=var.matrix)
summary.PM$variable=rownames(summary.PM)
summary.PM=data.table(summary.PM)

write.csv(summary.PM,file="proportion.csv",row.names = F)
##############################
############################################################

