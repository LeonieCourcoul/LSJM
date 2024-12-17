#' @rdname plot
#' @import ggplot2
#' @import survminer
#' @export
#'

plot.lsjm_covDepSingle <- function(Objectlsjm, which = 'long.fit', Objectpredict, break.times = NULL, ID.ind = NULL, xlim = NULL, ylim = NULL){


  Objectlsmm <- Objectlsjm$control$Objectlsmm
  if(is.null(Objectpredict)){
    stop("Not implemented")
  }

  graph <- NULL

  oldpar <- graphics::par(no.readonly = TRUE) # code line i
  on.exit(graphics::par(oldpar)) # code line i + 1

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
    prediction <- cbind(pred.CV, data.long$window)
    mean.pred <- by(prediction[,1], prediction[,ncol(prediction)], mean)
    obstime.mean <- by(data.long[,timeVar], data.long$window, mean)
    df <- cbind(obstime.mean, mean.obs, IC.sup, IC.inf, mean.pred)
    df <- as.data.frame(df)
    k <- ggplot2::ggplot(df,  ggplot2::aes(obstime.mean, mean.obs, ymin = IC.sup, ymax = IC.inf))
    graph.fit.long <- k +  ggplot2::geom_pointrange( ggplot2::aes(ymin = IC.sup, ymax = IC.inf), shape =1)+
      ggplot2::geom_point(ggplot2::aes(obstime.mean, mean.pred), size = 3, shape = 17) +
      ggplot2::scale_x_continuous(name = "Time") +
      ggplot2::scale_y_continuous(name = "Current Value") +
      ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                     panel.background = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"),
                     axis.text=ggplot2::element_text(size=15),
                     axis.title=ggplot2::element_text(size=18),
                     plot.title = ggplot2::element_text(size = 20, face = "bold"))+
      ggplot2::ggtitle("Longitudinal goodness-of-fit")+ggplot2::coord_cartesian(xlim = xlim,ylim = ylim, expand = TRUE)
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
      pred.CV.id$y <- data.long[which(data.long$id == ind), value.var]
      pred.CV.id$CI.sup <- pred.CV.id$predY + 1.96*pred.CV.id$predSD
      pred.CV.id$CI.inf <- pred.CV.id$predY - 1.96*pred.CV.id$predSD

      #browser()
      traj_ind <- ggplot2::ggplot() +

        ggplot2::geom_line(pred.CV.id, mapping = aes(x=time, y=predY, group = id, color = 'Predicted'))+
        ggplot2::geom_line(pred.CV.id, mapping = aes(x=time, y=CI.sup, group = id,  color = 'Predicted'),linetype = 2)+
        ggplot2::geom_line(pred.CV.id, mapping = aes(x=time, y=CI.inf, group = id,  color = 'Predicted'),linetype = 2)+

        ggplot2::geom_ribbon( pred.CV.id,mapping=
                                aes(x=time,ymin=CI.inf,ymax=CI.sup), fill="#998ec3", alpha=0.3,linetype = 3)+

        ggplot2::geom_point(pred.CV.id, mapping = aes(x=time, y=y, group = id,color = "Observed"),shape =17)+
        xlab("Time") + ylab("Y") +

        ggplot2::facet_wrap(~id, ncol = 3)+
        ggplot2::scale_color_manual(name='',
                                    breaks=c('Predicted', 'Observed'),
                                    values=c('Predicted'='#998ec3', 'Observed'='#000000'),
                                    guide = guide_legend(override.aes = list(
                                      linetype = c(rep("solid", 1), "blank"),
                                      shape = c(NA,  17))))+
        ggplot2::theme(
          panel.background = element_blank(),
          legend.position = "bottom",
          legend.box = "vertical",
          #axis.title.x = element_text(color = "black", size = 10),
          #axis.title.y = element_text(color = "black", size = 10),
          panel.grid = element_blank(),
          #legend.key = element_blank(),
          axis.text=ggplot2::element_text(size=15),
          axis.title=ggplot2::element_text(size=18),
          plot.title = ggplot2::element_text(size = 20, face = "bold"),
          legend.text = element_text(color = "black", size = 14),
          axis.line = element_line(color = "black",
                                   linetype = "solid")
          #axis.text = element_text(size = 10, color = "black")
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
    Surv.fit1 <- survminer::surv_fit(Surv(Time_Tsort, delta1sort) ~ 1, data = C1.sort)
    surv_plot <- survminer::ggsurvplot(Surv.fit1, data = C1.sort, fun = "cumhaz",
                                       conf.int = TRUE, legend.title = "",
                                       legend.labs = c("Survival Curve"),
                                       xlab = "Time", palette = "#B2BABB", ylim = ylim)
    surv_plot <- surv_plot$plot
    color_mapping <- c("#B2BABB","#E74C3C")
    graph.surv.1<-surv_plot +
      ggplot2::geom_step(aes(timeFormSurv, pred, color = "Nelson-Aalen"),
                         data = Cum.pred1.sort,
                         linetype = "3313",
                         size = 1) +
      ggplot2::geom_step(aes(timeFormSurv, pred, color = "Prediction"),
                         data = Cum.pred1.sort,
                         linetype = "3313",
                         size = 1)+
      ggplot2::scale_color_manual(name = "",
                                  values = setNames(color_mapping, c("Nelson-Aalen", "Prediction"))) +
      ggplot2::guides(color = guide_legend(title = "", override.aes = list(linetype = "solid", size = 2)))+
      ggplot2::theme(
        legend.key.size = unit(3, "lines"),  # Ajuster la taille de la clé dans la légende
        legend.text = element_text(size = 10)  # Ajuster la taille du texte dans la légende
      ) +
      ggplot2::theme(legend.position = c(0.1, 0.8))+
      ggplot2::ggtitle("1st event")
    graph <- list(graph.surv.1 = graph.surv.1)


    #print(graph.surv.1)
    #print(graph.surv.2)

    #graph <- list(graph.surv.1 = graph.surv.1, graph.surv.2 = graph.surv.2)

  }


  graph

}
