library(TMB) 
# Call TMB function value
compile("prodlikd_spatial_simulations.cpp")#,"&> C:/Users/mcdonaldra/Documents/errors.txt")
# Dynamically link the C++ code
dyn.load(dynlib("prodlikd_spatial_simulations"))

# Logit function
logitp=function(p){log(p/(1-p))}
# Inverse logist function
logitpi=function(t){exp(t)/(1+exp(t))}

library(sp)
#Create square:
x_coord<-c(-2.27,2.15,2.15,-2.27)
y_coord<-c(-1.19,-1.19,0.92,0.92)
xy<-cbind(x_coord,y_coord)
p = Polygon(xy)
ps = Polygons(list(p),1)
sps = SpatialPolygons(list(ps))

yearly_tows<-150
set.seed(30)
random_loc<-spsample(sps,yearly_tows,"random")

dismat=cbind(random_loc@coords[,1],random_loc@coords[,2])

# Distance matrix
Dist= as.matrix(dist(dismat))

data<-list()
data$H<-matrix(nrow=150,ncol=4)
data$A<-matrix(nrow=150,ncol=3)
data$X<-matrix(rep(1,150),nrow=150,ncol=1)
data$s<-rep(NA,150)
# data$D<-matrix(nrow=150,ncol=150)
data$D<-Dist

par<-list()
par$betat<-log(1.878483e-05)
par$betant<-log(0.001941084)
# par$betat<-log(0.002)
# par$betant<-log(0.001)
par$theta<-logitp(0.8795768)

par$omegat<-rep(0,150)
par$omegant<-rep(0,150)

par$lognut<-0
par$lognunt<-0

par$logPhit<-log(0.07)
par$logPhint<-log(0.1)

par$logSigmat<-log(sqrt(3))
par$logSigmant<-log(sqrt(1.5))

random<-c("omegat","omegant")

map<-list(lognut=factor(NA),lognunt=factor(NA))

obj<-MakeADFun(data=data,parameters=par,map=map,random=random)
# check<-checkConsistency(obj,n=150)
# barp<-obj$simulate()
# 
# simdata<-barp
# 
# obj2<-MakeADFun(data=simdata,parameters=par,map=map,random=random)
# Opt<-nlminb(obj2$par,obj2$fn,obj2$gr)
# rep<-sdreport(obj2)

#Keep in mind that the model fitting can sometimes get lost in NAs, so will slow down a lot

simpar<-list()
simpar$betat<--5
simpar$betant<--5
simpar$theta<-logitp(0.8)

simpar$omegat<-rep(0,150)
simpar$omegant<-rep(0,150)

simpar$lognut<-0
simpar$lognunt<-0

simpar$logPhit<--2
simpar$logPhint<--2

simpar$logSigmat<-0
simpar$logSigmant<-0

time1<-Sys.time()
nrep<-500
Report_list<-list()
sim_save<-list()
diff_list<-list()
diff_list2<-list()
rep_list<-list()
msg_list<-list()
for (i in 1:nrep) {
  tryCatch({
    simdata <- obj $ simulate(complete=T)
    sim_save[[i]]<-simdata
    obj2 <- MakeADFun ( simdata , simpar , random=random, map=map )
    Opt2 <- try(nlminb(start=obj2$par,obj=obj2$fn,gr=obj2$gr,
                       control=list(eval.max=1000,iter.max=1000),silent=T),T)
    rep_list[[i]] <- sdreport(obj2,bias.correct=F)
    if (!is.null(rep_list[[i]])){
      Report_list[[i]] = obj2$report()
      msg_list[[i]] <- Opt2$message
    }
  }, error=function(e) {})
}
time2<-Sys.time()
time2-time1

fun_rep<-500

converge<-0
false<-0
x_conv<-0
singular<-0
null_vec<-nrep
chosen_ones<-c()
for (i in 1:fun_rep){
  if (!is.null(msg_list[[i]])){null_vec<-null_vec-1
  if (msg_list[[i]]=="relative convergence (4)") {converge<-converge+1
  chosen_ones<-c(chosen_ones,i)}
  if (msg_list[[i]]=="false convergence (8)") {false<-false+1}
  if(msg_list[[i]]=="singular convergence (7)"){singular<-singular+1}
  if(msg_list[[i]]=="both X-convergence and relative convergence (5)"){x_conv<-x_conv+1
  chosen_ones<-c(chosen_ones,i)}
  }
}
conv_frame_low<-data.frame(converge=converge,false=false,singular=singular,x=x_conv)

betat<-rep(NA,fun_rep)
betant<-rep(NA,fun_rep)
theta<-rep(NA,fun_rep)
phit<-rep(NA,fun_rep)
phint<-rep(NA,fun_rep)
sigmat<-rep(NA,fun_rep)
sigmant<-rep(NA,fun_rep)
for (i in chosen_ones){
  betat[i]<-summary(rep_list[[i]])[row.names(summary(rep_list[[i]]))=="betat",][1]
  betant[i]<-summary(rep_list[[i]])[row.names(summary(rep_list[[i]]))=="betant",][1]
  theta[i]<-summary(rep_list[[i]])[row.names(summary(rep_list[[i]]))=="theta",][1]
  phit[i]<-exp(summary(rep_list[[i]])[row.names(summary(rep_list[[i]]))=="logPhit",][1])
  phint[i]<-exp(summary(rep_list[[i]])[row.names(summary(rep_list[[i]]))=="logPhint",][1])
  sigmat[i]<-exp(summary(rep_list[[i]])[row.names(summary(rep_list[[i]]))=="logSigmat",][1])
  sigmant[i]<-exp(summary(rep_list[[i]])[row.names(summary(rep_list[[i]]))=="logSigmant",][1])
}

hist_plot_low<-data.frame(value=c(betat,betant,logitpi(theta),phit,phint,
                                  sigmat,sigmant),
                          parameter=rep(c("betat","betant","theta","phit",
                                          "phint","sigmat","sigmant"),each=fun_rep),
                          true=rep(c(par$betat,par$betant,logitpi(par$theta),
                                     exp(par$logPhit),exp(par$logPhint),
                                     exp(par$logSigmat),exp(par$logSigmant))
                                   ,each=fun_rep))
hist_plot_low$run<-rep("low",length(hist_plot_low$value))
hist_plot_low$parameter<-factor(hist_plot_low$parameter,
                                label=c(expression(beta [nt]),
                                        expression(beta [t]),
                                        expression(phi [nt]),
                                        expression(phi [t]),
                                        expression(sigma [nt]),
                                        expression(sigma [t]),
                                        expression(p [nt])))

library(ggplot2)

par_hist_low<-ggplot(hist_plot_low)+
  geom_histogram(aes(x=value),col="black",fill="grey")+
  geom_vline(aes(xintercept = true),col="red")+
  theme(strip.background = element_blank(),legend.position = "none")+
  facet_wrap(~parameter,scales="free",labeller=label_parsed)+
  xlab("Estimated Value")+ylab("Frequency")+
  theme_bw()

#Look if the predicted ldat and ldant are good

diff_ldat<-matrix(nrow=150,ncol=nrep)
diff_ldant<-matrix(nrow=150,ncol=nrep)
for (i in chosen_ones) {
  diff_ldat[,i]<-(sim_save[[i]]$ldat-rep_list[[i]]$value[which(names(rep_list[[i]]$value)=="ldat")])/sim_save[[i]]$ldat
  diff_ldant[,i]<-(sim_save[[i]]$ldant-rep_list[[i]]$value[which(names(rep_list[[i]]$value)=="ldant")])/sim_save[[i]]$ldant
}

hist_low<-ggplot()+
  geom_histogram(aes(x=as.vector(diff_ldat*100)),
                 col="red",fill="red",alpha=0.3,binwidth = 10)+
  geom_histogram(aes(x=as.vector(diff_ldant*100)),
                 col="blue",fill="blue",alpha=0.3,binwidth = 10)+
  coord_cartesian(xlim=c(-300,110))+
  xlab("Percent Difference Between Simulated and Estimated Indices")+
  ylab("Count")+
  theme_bw()
ggsave(filename="hist_low.png",plot=hist_low,width=6,height=4)

perc_ldat_rem_low<-1-((length(which(!(is.na(diff_ldat)))) - length(which(diff_ldat < (-4))))/length(which(!(is.na(diff_ldat)))))
perc_ldant_rem_low<-1-((length(which(!(is.na(diff_ldant)))) - length(which(diff_ldant < (-4))))/length(which(!(is.na(diff_ldant)))))

par<-list()
# par$betat<-log(1.878483e-05)
# par$betant<-log(0.001941084)
par$betat<-log(0.002)
par$betant<-log(0.001)
par$theta<-logitp(0.8795768)

par$omegat<-rep(0,150)
par$omegant<-rep(0,150)

par$lognut<-0
par$lognunt<-0

par$logPhit<-log(0.07)
par$logPhint<-log(0.1)

par$logSigmat<-log(sqrt(3))
par$logSigmant<-log(sqrt(1.5))

random<-c("omegat","omegant")

map<-list(lognut=factor(NA),lognunt=factor(NA))

obj<-MakeADFun(data=data,parameters=par,map=map,random=random)
# check<-checkConsistency(obj,n=150)
# barp<-obj$simulate()
# 
# simdata<-barp
# 
# obj2<-MakeADFun(data=simdata,parameters=par,map=map,random=random)
# Opt<-nlminb(obj2$par,obj2$fn,obj2$gr)
# rep<-sdreport(obj2)

#Keep in mind that the model fitting can sometimes get lost in NAs, so will slow down a lot

simpar<-list()
simpar$betat<--5
simpar$betant<--5
simpar$theta<-logitp(0.8)

simpar$omegat<-rep(0,150)
simpar$omegant<-rep(0,150)

simpar$lognut<-0
simpar$lognunt<-0

simpar$logPhit<--2
simpar$logPhint<--2

simpar$logSigmat<-0
simpar$logSigmant<-0

time1<-Sys.time()
nrep<-500
Report_list<-list()
sim_save<-list()
diff_list<-list()
diff_list2<-list()
rep_list<-list()
msg_list<-list()
for (i in 1:nrep) {
  tryCatch({
    simdata <- obj $ simulate(complete=T)
    sim_save[[i]]<-simdata
    obj2 <- MakeADFun ( simdata , simpar , random=random, map=map )
    Opt2 <- try(nlminb(start=obj2$par,obj=obj2$fn,gr=obj2$gr,
                       control=list(eval.max=1000,iter.max=1000),silent=T),T)
    rep_list[[i]] <- sdreport(obj2,bias.correct=F)
    if (!is.null(rep_list[[i]])){
      Report_list[[i]] = obj2$report()
      msg_list[[i]] <- Opt2$message
    }
  }, error=function(e) {})
}
time2<-Sys.time()
time2-time1

fun_rep<-500

converge<-0
false<-0
x_conv<-0
singular<-0
null_vec<-nrep
chosen_ones<-c()
for (i in 1:fun_rep){
  if (!is.null(msg_list[[i]])){null_vec<-null_vec-1
  if (msg_list[[i]]=="relative convergence (4)") {converge<-converge+1
  chosen_ones<-c(chosen_ones,i)}
  if (msg_list[[i]]=="false convergence (8)") {false<-false+1}
  if(msg_list[[i]]=="singular convergence (7)"){singular<-singular+1}
  if(msg_list[[i]]=="both X-convergence and relative convergence (5)"){x_conv<-x_conv+1
  chosen_ones<-c(chosen_ones,i)}
  }
}
conv_frame_high<-data.frame(converge=converge,false=false,singular=singular,x=x_conv)

betat<-rep(NA,fun_rep)
betant<-rep(NA,fun_rep)
theta<-rep(NA,fun_rep)
phit<-rep(NA,fun_rep)
phint<-rep(NA,fun_rep)
sigmat<-rep(NA,fun_rep)
sigmant<-rep(NA,fun_rep)
for (i in chosen_ones){
  betat[i]<-summary(rep_list[[i]])[row.names(summary(rep_list[[i]]))=="betat",][1]
  betant[i]<-summary(rep_list[[i]])[row.names(summary(rep_list[[i]]))=="betant",][1]
  theta[i]<-summary(rep_list[[i]])[row.names(summary(rep_list[[i]]))=="theta",][1]
  phit[i]<-exp(summary(rep_list[[i]])[row.names(summary(rep_list[[i]]))=="logPhit",][1])
  phint[i]<-exp(summary(rep_list[[i]])[row.names(summary(rep_list[[i]]))=="logPhint",][1])
  sigmat[i]<-exp(summary(rep_list[[i]])[row.names(summary(rep_list[[i]]))=="logSigmat",][1])
  sigmant[i]<-exp(summary(rep_list[[i]])[row.names(summary(rep_list[[i]]))=="logSigmant",][1])
}

hist_plot_high<-data.frame(value=c(betat,betant,logitpi(theta),phit,phint,
                                   sigmat,sigmant),
                           parameter=rep(c("betat","betant","theta","phit",
                                           "phint","sigmat","sigmant"),each=fun_rep),
                           true=rep(c(par$betat,par$betant,logitpi(par$theta),
                                      exp(par$logPhit),exp(par$logPhint),
                                      exp(par$logSigmat),exp(par$logSigmant))
                                    ,each=fun_rep))
hist_plot_high$run<-rep("high",length(hist_plot_high$value))
hist_plot_high$parameter<-factor(hist_plot_high$parameter,
                                 label=c(expression(beta [nt]),
                                         expression(beta [t]),
                                         expression(phi [nt]),
                                         expression(phi [t]),
                                         expression(sigma [nt]),
                                         expression(sigma [t]),
                                         expression(p [nt])))

library(ggplot2)

par_hist_high<-ggplot(hist_plot_high)+
  geom_histogram(aes(x=value),col="black",fill="grey")+
  geom_vline(aes(xintercept = true),col="red")+
  theme(strip.background = element_blank(),legend.position = "none")+
  facet_wrap(~parameter,scales="free",labeller=label_parsed)+
  xlab("Estimated Value")+ylab("Frequency")+
  theme_bw()

#Look if the predicted ldat and ldant are good

diff_ldat<-matrix(nrow=150,ncol=nrep)
diff_ldant<-matrix(nrow=150,ncol=nrep)
for (i in chosen_ones) {
  diff_ldat[,i]<-(sim_save[[i]]$ldat-rep_list[[i]]$value[which(names(rep_list[[i]]$value)=="ldat")])/sim_save[[i]]$ldat
  diff_ldant[,i]<-(sim_save[[i]]$ldant-rep_list[[i]]$value[which(names(rep_list[[i]]$value)=="ldant")])/sim_save[[i]]$ldant
}

hist_high<-ggplot()+
  geom_histogram(aes(x=as.vector(diff_ldat*100)),
                 col="red",fill="red",alpha=0.3,binwidth = 10)+
  geom_histogram(aes(x=as.vector(diff_ldant*100)),
                 col="blue",fill="blue",alpha=0.3,binwidth = 10)+
  coord_cartesian(xlim=c(-150,110))+
  xlab("Percent Difference Between Simulated and Estimated Indices")+
  ylab("Count")+
  theme_bw()
ggsave(filename="hist_high.png",plot=hist_high,width=6,height=4)

perc_ldat_rem_high<-1-((length(which(!(is.na(diff_ldat)))) - length(which(diff_ldat < (-1.5))))/length(which(!(is.na(diff_ldat)))))
perc_ldant_rem_high<-1-((length(which(!(is.na(diff_ldant)))) - length(which(diff_ldant < (-1.5))))/length(which(!(is.na(diff_ldant)))))

hist_plot_both<-rbind(hist_plot_low,hist_plot_high)

par_hist_both<-ggplot(hist_plot_both)+
  geom_histogram(aes(x=value,col=run,fill=run),alpha=0.3)+
  geom_vline(aes(xintercept = true,col=run))+
  theme(strip.background = element_blank(),legend.position = "none")+
  scale_color_manual(name="Setting",values=c("red","blue"),labels=c("1","2"))+
  scale_fill_manual(name="Setting",values=c("red","blue"),labels=c("1","2"))+
  facet_wrap(~parameter,scales="free",labeller=label_parsed)+
  xlab("Estimated Value")+ylab("Frequency")+
  theme_bw()
ggsave(filename="par_hist_both.png",plot=par_hist_both,width=8,height=6)

