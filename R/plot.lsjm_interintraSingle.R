#' @rdname plot.lsjm
#' @importFrom graphics plot par
#' @importFrom dplyr left_join
#' @importFrom ggplot2 ggplot aes geom_pointrange geom_point scale_x_continuous scale_y_continuous theme element_blank element_line element_text ggtitle coord_cartesian geom_line geom_ribbon facet_wrap scale_color_manual guide_legend guides scale_fill_manual geom_step scale_linetype_manual
#' @importFrom survminer surv_fit ggsurvplot
#' @importFrom survival Surv
#' @export

plot.lsjm_interintraSingle <- function(x, which = 'long.fit', Objectpredict = NULL, break.times = NULL, ID.ind = NULL, xlim = NULL, ylim = NULL, ...){


  Objectlsjm <- x
  Objectlsmm <- Objectlsjm$control$Objectlsmm
  #if(is.null(Objectranef)){
  #  Objectranef <- ranef(Objectlsmm)
  #}
  graph <- NULL

  oldpar <- par(no.readonly = TRUE) # code line i
  on.exit(par(oldpar)) # code line i + 1

  ObjectpredictY <- Objectpredict$predictY

  if(which == 'long.fit'){
    formFixed <- Objectlsmm$control$formFixed
    timeVar <- Objectlsmm$control$timeVar
    data.long <- Objectlsmm$control$data.long
    value.var <- as.character(formFixed[[2]])
    pred.CV <- ObjectpredictY$predY
    if(is.null(break.times)){
      timeInterv <- range(data.long[,timeVar])
      break.times <- quantile(timeInterv,prob=seq(0,1,length.out=10))
    }
    data.long$window <- cut(data.long[,timeVar], break.times, include.lowest = T)
    mean.obs <- by(data.long[,value.var], data.long$window, mean)
    sd.obs <- by(data.long[,value.var], data.long$window, sd)
    length.obs <- by(data.long[,value.var], data.long$window, length)
    IC.inf <- mean.obs - 1.96*sd.obs/sqrt(length.obs)
    IC.sup <- mean.obs + 1.96*sd.obs/sqrt(length.obs)
    window.pred <- cut(ObjectpredictY$time, break.times, include.lowest = T)
    ObjectpredictY$time.new.pred <- ObjectpredictY$time
    data.long$time.new.pred <- data.long[,timeVar]
    prediction <- left_join(ObjectpredictY[,c("id","predY", "time.new.pred")], data.long[,c("id", "window", "time.new.pred")])
    mean.pred <- by(prediction$predY, prediction$window, mean)
    obstime.mean <- by(data.long[,timeVar], data.long$window, mean)
    df <- cbind(obstime.mean, mean.obs, IC.sup, IC.inf, mean.pred)
    df <- as.data.frame(df)
    k <- ggplot(df,  aes(obstime.mean, mean.obs, ymin = IC.sup, ymax = IC.inf))
    graph.fit.long <- k +  geom_pointrange( aes(ymin = IC.sup, ymax = IC.inf), shape =1) +
      geom_point(aes(obstime.mean, mean.pred), size = 3, shape = 17) +
      scale_x_continuous(name = "Time") +
      scale_y_continuous(name = "Current Value") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     panel.background = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text=element_text(size=12),
                     axis.title=element_text(size=14,face="bold"))+
      ggtitle("Longitudinal goodness-of-fit")
    graph <- list(long.fit = graph.fit.long)
  }

  if(which == 'traj.ind'){
    if(is.null(ID.ind)){
      stop("You have to design some individual ID to plot the the individual trajectories.")
    }
    ID.ind <- as.vector(ID.ind)
    pred.CV <- as.data.frame(ObjectpredictY)
    data.long <- Objectlsmm$control$data.long
    formFixed <- Objectlsmm$control$formFixed
    value.var <- as.character(formFixed[[2]])
    graph.traj.ind <- c()
    for(ind in ID.ind){
      pred.CV.id <- pred.CV[which(pred.CV$id == ind),]
      data.idselect <- cbind(data.long$id[which(data.long$id == ind)],data.long[which(data.long$id == ind), Objectlsmm$control$timeVar],data.long[which(data.long$id == ind), value.var])
      data.idselect <- as.data.frame(data.idselect)
      colnames(data.idselect) <- c("id","time", "y")
      #pred.CV.id$y <- data.long[which(data.long$id == ind), value.var]
      pred.CV.id$CI.sup <- pred.CV.id$predY + 1.96*sqrt(pred.CV.id$predSD_inter**2 + pred.CV.id$predSD_intra**2)
      pred.CV.id$CI.inf <- pred.CV.id$predY - 1.96*sqrt(pred.CV.id$predSD_inter**2 + pred.CV.id$predSD_intra**2)

      traj_ind <- ggplot() +

        geom_line(pred.CV.id, mapping = aes(x=time, y=predY, group = id, color = 'Predicted'))+
        geom_line(pred.CV.id, mapping = aes(x=time, y=CI.sup, group = id,  color = 'Predicted'),linetype = 2)+
        geom_line(pred.CV.id, mapping = aes(x=time, y=CI.inf, group = id,  color = 'Predicted'),linetype = 2)+

        geom_ribbon( pred.CV.id,mapping=
                                aes(x=time,ymin=CI.inf,ymax=CI.sup), fill="#998ec3", alpha=0.3)+

        geom_point(data.idselect, mapping = aes(x=time, y=y, group = id,color = "Observed"),shape =17)+
        xlab("Time") + ylab("Y") +

        facet_wrap(~id, ncol = 3)+
        scale_color_manual(name='',
                                    breaks=c('Predicted', 'Observed'),
                                    values=c('Predicted'='#998ec3', 'Observed'='#000000'),
                                    guide = guide_legend(override.aes = list(
                                      linetype = c(rep("solid", 1), "blank"),
                                      shape = c(NA,  17))))+
        theme(
          panel.background = element_blank(),
          legend.position = "bottom",
          legend.box = "vertical",
          axis.title.x = element_text(color = "black", size = 13),
          axis.title.y = element_text(color = "black", size = 13),
          panel.grid = element_blank(),
          #legend.key = element_blank(),
          legend.text = element_text(color = "black", size = 13),
          axis.line = element_line(color = "black",
                                   linetype = "solid"),
          axis.text = element_text(size = 13, color = "black")
        )+coord_cartesian(xlim = xlim,ylim = ylim, expand = TRUE)

      graph[[paste("traj.ind",ind, sep = "_")]] <- traj_ind
    }

  }


  if(which == 'survival.fit'){
    data.long <- Objectlsmm$control$data.long
    data.id <- data.long[!duplicated(data.long$id),]
    data.id$e1.new <- data.id[,all.vars(Objectlsjm$control$deltas[["delta1"]])]
    C1.sort <- data.id[order(data.id[,all.vars(Objectlsjm$control$Time[["Time_T"]])]),]
    Cum.pred1 <- apply(Objectpredict$predictCum_01, 2, mean)
    Cum.pred1 <- cbind(Cum.pred1, unique(sort(data.id[,all.vars(Objectlsjm$control$Time[["Time_T"]])])))
    Cum.pred1 <- as.data.frame(Cum.pred1)
    Cum.pred1.sort <- Cum.pred1[order(Cum.pred1[,2]),]
    colnames(Cum.pred1.sort) <- c("pred","timeFormSurv")
    timeFormSurv <- Cum.pred1.sort$timeFormSurv
    pred <- Cum.pred1.sort$pred
    #arrange(Cum.pred1, V2)
    Time_Tsort <- Objectlsjm$control$Time[["Time_T"]]
    delta1sort <- Objectlsjm$control$deltas[["delta1"]]
    C1.sort$Time_Tsort <- C1.sort[all.vars(Time_Tsort)][,1]
    C1.sort$delta1sort <- C1.sort[all.vars(delta1sort)][,1]
    Surv.fit1 <- surv_fit(Surv(Time_Tsort, delta1sort) ~ 1, data = C1.sort)
    surv_plot <- ggsurvplot(Surv.fit1, data = C1.sort, fun = "cumhaz",
                                       conf.int = TRUE, legend.title = "",
                                       legend.labs = c("Survival Curve"),
                                       xlab = "Time", palette = "#B2BABB", ylim = ylim)
    surv_plot <- surv_plot$plot
    color_mapping <- c("#B2BABB","#E74C3C")
    graph.surv.1<-surv_plot +
      geom_step(aes(timeFormSurv, pred, color = "Nelson-Aalen"),
                         data = Cum.pred1.sort,
                         linetype = "3313",
                         size = 1) +
      geom_step(aes(timeFormSurv, pred, color = "Prediction"),
                         data = Cum.pred1.sort,
                         linetype = "3313",
                         size = 1)+
      scale_color_manual(name = "",
                                  values = setNames(color_mapping, c("Nelson-Aalen", "Prediction"))) +
      guides(color = guide_legend(title = "", override.aes = list(linetype = "solid", size = 2)))+
      theme(
        legend.key.size = unit(3, "lines"),  # Ajuster la taille de la clé dans la légende
        legend.text = element_text(size = 10)  # Ajuster la taille du texte dans la légende
      ) +
      theme(legend.position = c(0.1, 0.8))+
      ggtitle("1st event")
    graph <- list(graph.surv.1 = graph.surv.1)


    #print(graph.surv.1)
    #print(graph.surv.2)

    #graph <- list(graph.surv.1 = graph.surv.1, graph.surv.2 = graph.surv.2)

  }


  graph

}
