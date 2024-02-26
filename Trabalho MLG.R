# Trabalho de MLG

library(survival)
library(tidyverse)
library(AdequacyModel)
library(readr)

# LEITURA DO BANCO DE DADOS ----

adesao <- read_table("adesao.txt")
adesao


# Análise descritiva ----

tempo = adesao$`"tempo"` %>% as.numeric()
censura = adesao$`"status"` %>% as.numeric()
propatraso = adesao$`"propatraso"`
comprimdia = adesao$`"comprimdia"`

km <- survfit(Surv(tempo,censura) ~ 1) # Estimação de KM 
summary(km)

par(mfrow=c(1,2))

plot(km, conf.int=F, mark.time=F, lwd=2, col="#0350a8",
     xlab= "Tempo", ylab="Função de Sobrevivência", 
     ylim = c(0.7,1))

plot(km, conf.int=F, mark.time=T, lwd=2, col="#0350a8",
     xlab= "Tempo", ylab="Função de Sobrevivência", 
     ylim = c(0.7,1))

par(mfrow=c(1,1))



# Função de Risco Acumulada - Nelson Aalen
par(mfrow=c(1,2))

ena <- survfit(coxph(Surv(tempo,censura) ~ 1,method="breslow"))
summary(ena)
plot(ena,conf.int=F,fun="cumhaz", lwd=2, col="black",
     xlab= "Tempo", ylab="Função de Risco Acumulado",)

# Gráfico TTT
TTT(tempo, col="red",lwd=2.5,grid=T,lty=2)

par(mfrow=c(1,1))



# - Análise descritiva por grupos ----

linha = adesao$`"linha"` %>% as.factor()

# Estime separadamente por grupo
km2 <- survfit(Surv(tempo,censura) ~ linha)
summary(km2)

plot(km2, conf.int=F, mark.time=F, lwd=2, col=c("blue","pink","green"),
     xlab= "Tempo", ylab="Função de Sobrevivência")
legend("bottomright", legend = c("Linha 1","Linha 2","Linha 3"), fill = c("blue", "pink", "green"))


# Gráfico TTT
dad2<-data.frame(tempo,linha)
tempo1<-dad2[dad2$linha == "1",]
TTT(tempo1$tempo, col="red", lwd=2.5, grid=TRUE, lty=2)
tempo2<-dad2[dad2$linha == "2",]
TTT(tempo2$tempo, col="red", lwd=2.5, grid=TRUE, lty=2)
tempo3<-dad2[dad2$linha == "3",]
TTT(tempo3$tempo, col="red", lwd=2.5, grid=TRUE, lty=2)

# Teste de Comparação

survdiff(Surv(tempo,censura) ~ linha, rho = 0) # Log rank - Riscos proporcionais
survdiff(Surv(tempo,censura) ~ linha, rho = 1) # Wilcoxon


ena <- survfit(coxph(Surv(tempo,censura) ~ linha,method="breslow"))
summary(ena)
plot(ena,conf.int=F,fun="cumhaz" )



# Modelo de regressão paramétrico ----


# Distribuição weibull, log normal e log logistica

km <- survfit(Surv(tempo,censura) ~ 1)
summary(km)

skm = km$surv

# Distribuição weibull
mwweb = survreg(Surv(tempo,censura) ~ 1, dist = "weibull")
summary(mwweb)

alfa_web = exp(mwweb$coefficients)
gama_web = 1 / mwweb$scale

st_weibull = exp(-(km$time/alfa_web)^gama_web)
st_weibull

# Distribuição log normal
mwlognorm = survreg(Surv(tempo,censura) ~ 1, dist = "lognormal")
summary(mwlognorm)

mu = mwlognorm$coefficients
sdev = mwlognorm$scale

st_lognorm = 1 - pnorm((log(km$time) - mu) / sdev)
st_lognorm

# Distribuição log logistica
mwloglogist = survreg(Surv(tempo,censura) ~ 1, dist = "loglogistic")
summary(mwloglogist)

alfa_loglogist = exp(mwloglogist$coefficients)
gama_loglogist = 1 / mwloglogist$scale

st_mwloglogist = 1 / (1 + (km$time/alfa_loglogist)^gama_loglogist)
st_mwloglogist

# Método 1A
plot(km, conf.int=F, mark.time=F, lwd=2, col="black",
     xlab= "Tempo", ylab="Função de Sobrevivência", ylim = c(0.7,1))
lines(c(0,km$time),c(1,st_weibull), lwd=2, col="red")
lines(c(0,km$time),c(1,st_lognorm), lwd=2, col="#0350a8")
lines(c(0,km$time),c(1,st_mwloglogist), lwd=2, col="#02752b")
legend("topright", bty = "n", legend = c("Kaplan Meier","Weibull","Log-Normal","Log-Logística"), fill = c("black", "red", "#0350a8","#02752b"))



# Método 1B
plot(skm, st_weibull, lwd=2, col="blue",
     xlab= "KM", ylab="S(t) Weibull")
lines(skm,skm, lwd=2, col="black")

plot(skm, st_lognorm, lwd=2, col="blue",
     xlab= "KM", ylab="S(t) Log Normal")
lines(skm,skm, lwd=2, col="black")

plot(skm, st_mwloglogist, lwd=2, col="blue",
     xlab= "KM", ylab="S(t) Log Logística")
lines(skm,skm, lwd=2, col="black")


extractAIC(mwweb)
extractAIC(mwlognorm)
extractAIC(mwloglogist)

BIC(mwweb)
BIC(mwlognorm)
BIC(mwloglogist)

MuMIn::AICc(mwweb)
MuMIn::AICc(mwlognorm)
MuMIn::AICc(mwloglogist)

# Modelo Weibull pode ser selecionado


# ------------------------------------------------ # 
# ------------------------------------------------ #
# ------------------------------------------------ #
# ------------------------------------------------ #
# ------------------------------------------------ #
# ------------------------------------------------ #


# Ajuste de um modelo de regressão 

modelo_01_Survival = survreg(Surv(tempo,censura) ~ comprimdia + propatraso + linha, dist = "weibull")

modelo_01_Poisson = glm(censura ~ comprimdia + propatraso + linha + offset(log(tempo)),
                        family = poisson(link="log"))

summary(modelo_01_Survival)
summary(modelo_01_Poisson)

# Surveg ajusta o modelo do valor extremo, vamos transformar
alfa_web = exp(mwweb$coefficients)
alfa_web = exp(modelo_01_Survival$coefficients[1])
gama_web = 1 / modelo_01_Survival$scale

uj = modelo_01_Poisson$fitted.values
uj

gama_est = sum(censura) / sum((uj - censura)*log(tempo))
gama_est # Estimado por GLM Poisson
gama_web # Estimado por Survival

alfa_est = exp(-modelo_01_Poisson$coefficients[1])
alfa_est 

st_weibull1 = exp(-(km$time/4000)^gama_est)
st_weibull1  

st_weibull2 = exp(-(km$time/4000)^gama_web)
st_weibull2 


# Método 1A
par(mfrow=c(1,2))

plot(km, conf.int=F, mark.time=F, lwd=2, col="black",
     xlab= "Tempo", ylab="Função de Sobrevivência", ylim = c(0.7,1))
lines(c(0,km$time),c(1,st_weibull1), lwd=2, col="red")
legend("topright", bty = "n", legend = c("Kaplan Meier","Weibull: shape = 0,949"), fill = c("black", "red"))

plot(km, conf.int=F, mark.time=F, lwd=2, col="black",
     xlab= "Tempo", ylab="Função de Sobrevivência", ylim = c(0.7,1))
lines(c(0,km$time),c(1,st_weibull2), lwd=2, col="red")
legend("topright", bty = "n", legend = c("Kaplan Meier","Weibull: shape = 0,957"), fill = c("black", "red"))

par(mfrow=c(1,1))


# --------------------------------------- #
# --------------------------------------- #
# --------------------------------------- #
# --------------------------------------- #
# --------------------------------------- #
# --------------------------------------- #
# --------------------------------------- #
# --------------------------------------- #
# --------------------------------------- #




install.packages('My.stepwise')
library(survival)
library(tidyverse)
library(AdequacyModel)
library(readr)
library(caret)
library(MASS)
library(My.stepwise)

adesao <- read_table("adesao.txt")
str(adesao)
tempo = adesao$`"tempo"` %>% as.numeric()
censura = adesao$`"status"` %>% as.numeric()
propatraso = adesao$`"propatraso"`
comprimdia = adesao$`"comprimdia"`
linha = adesao$`"linha"` %>% as.factor()

data <- data.frame(censura,log(tempo),propatraso,comprimdia,linha)

# Ajuste de um modelo de regressão 
####
modelo_nulo = glm(censura ~ 1 + offset(log(tempo)),
                  family = poisson(link="log"))
summary(modelo_nulo)

modelo_saturado = glm(censura ~ comprimdia + propatraso + linha + offset(log(tempo)),
                      family = poisson(link="log"))
summary(modelo_saturado)

modelo_step <- step(modelo_saturado,modelo_nulo, direction = "both", 
                    trace = FALSE)

summary(modelo_step)

modelo_01_Poisson = glm(censura ~ propatraso + offset(log(tempo)),
                        family = poisson(link="log"))
summary(modelo_01_Poisson)

uj = modelo_01_Poisson$fitted.values
uj

gama_est = sum(censura) / sum((uj - censura)*log(tempo))
gama_est # Estimado por GLM Poisson

odds<-exp(coef(modelo_01_Poisson))
odds
exp(confint(modelo_01_Poisson))
##------------------------------------------##
##    ANALYSING THE MODEL GOODNESS OF FIT   ##
##------------------------------------------##
## H0: The specified link function adjust   ##
## correctly the data                       ##
## H1: The specified link function does not ##
## adjust correctly the data                ##
##------------------------------------------##
pchisq(deviance(modelo_01_Poisson), df.residual(modelo_01_Poisson), lower.tail=FALSE)

##----------------------##
## Analysis of Deviance ##
##----------------------##
xtable::xtable(anova(modelo_01_Poisson, test="Chisq"))


##------------------------------##
## 2.Naglekerke (1991) R-square ##
##------------------------------##
Num<- 1-exp((modelo_01_Poisson$dev-modelo_01_Poisson$null)/nrow(adesao))
Den<- 1-exp(-modelo_01_Poisson$null/nrow(adesao))
R2_Nag<- Num/Den; R2_Nag

##--------------------##
## An alternative R^2 ##
##--------------------##
Null.dev<- modelo_01_Poisson$null.deviance
Resi.dev<- deviance(modelo_01_Poisson)
##
R2.dev<- (Null.dev-Resi.dev)/Null.dev; R2.dev



# --------------------------------------- #
# --------------------------------------- #
# --------------------------------------- #
# --------------------------------------- #
# --------------------------------------- #
# --------------------------------------- #
# --------------------------------------- #
# --------------------------------------- #
# --------------------------------------- #
# --------------------------------------- #



# ANÁLISE DE RESÍDUOS E INFLUÊNCIA #
# --------------------------------------- #
  
  modelo_final = glm(censura ~ propatraso + offset(log(tempo)),
                     family = poisson(link="log"))
summary(modelo_final)

eta.w<- as.matrix(predict(modelo_final, type="link"))

##------------------------------------##
## RESIDUAL ANALYSIS: Usual Residuals ##
##------------------------------------##
res.w_ordinar<- resid(modelo_final, type="response")
res.w_working<- resid(modelo_final, type="working")
res.w_pearson<- resid(modelo_final, type="pearson")
res.w_devianc<- resid(modelo_final, type="deviance")

x11()
par(mfrow=c(2,2))
plot(eta.w,res.w_ordinar, xlab = expression(hat(eta)), ylab ="Resíduos", main="Resíduos Ordinários")
plot(eta.w,res.w_working, xlab = expression(hat(eta)), ylab ="Resíduos", main="Resíduos Workings")
plot(eta.w,res.w_pearson, xlab = expression(hat(eta)), ylab ="Resíduos", main="Resíduos de Pearson")
plot(eta.w,res.w_devianc, xlab = expression(hat(eta)), ylab ="Resíduos", main="Resíduos Desvio")

# ------------------------------------------------ #
x11()
par(mfrow=c(2,1))
plot(res.w_devianc, xlab = "Índice", ylab ="Resíduos", main="Resíduos Deviance x Índice")
plot(propatraso, res.w_devianc, xlab = "Covariável" , ylab ="Resíduos", main="Resíduos Deviance x Covariável")

##------------------------------------------##
## RESIDUAL ANALYSIS: Studentized Residuals ##
##------------------------------------------##

## For Pearson Residuals
x11()
par(mfrow=c(2,2))
plot(modelo_final)

## For deviance residuals

st.res.dev<- rstudent(modelo_final, type="deviance")
st.res.per<- rstudent(modelo_final, type="pearson")

rstudent<- rstudent(modelo_final)
x11()
plot(rstudent, ylim=c(-2,2), xlab = "Índice", ylab ="Resíduos Studentizados")
abline(h=c(-2,2),col="red", lty=1)

##-----------------------------##
## INFLUENCE MEASURES IN GLM'S ##
##-----------------------------##

betas<-as.vector(coef(modelo_final))
X.mat<- model.matrix(modelo_final)
eta.hat<- as.matrix(predict(modelo_final, type="link"))

## Hat values
hatvalues(modelo_final)

## Cook Distance
cooks.distance(modelo_final)

## Aggregated Measures
influence.measures(modelo_final)

## Some residual graphs
x11()
par(mfrow=c(3,2))
plot(cooks.distance(modelo_final), lwd = 2, col = "blue", type="h", ylab = "Distância de Cook", xlab = "Índice")
abline(h = 0.8, col = "red", lty = 1)
plot(hatvalues(modelo_final), type = "h", lwd = 2, col = "blue", ylab= "hii", xlab = "Índice")
plot(hatvalues(modelo_final), eta.hat, lwd = 2, col = "blue", xlab = expression(hat(eta)), ylab= "hii")
plot(dffits(modelo_final), lwd = 2, col = "blue", type="h", ylab = "dffits")
plot(dfbetas(modelo_final), type = "h", lwd = 2, col = "blue", ylab= "dfbetas", xlab = "Intercepto")
plot(covratio(modelo_final), lwd = 2, col = "blue", type="h", ylab= "covratio", xlab = "Índice")


valores_ajustados <- modelo_final$fitted.values

plot(eta.w, valores_ajustados,
     xlab=expression(hat(eta)), ylab="Variável Dependente Ajustada")