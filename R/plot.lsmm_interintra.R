plot.lsmm_interintra <- function(Objectlsmm, which = 'long.fit', Objectranef = NULL, break.times = NULL, ID.ind = NULL){

  if(is.null(Objectranef)){
    Objectranef <- ranef(Objectlsmm)
  }

  graph <- NULL

  if(which == 'long.fit'){
    formFixed <- Objectlsmm$control$formFixed
    timeVar <- Objectlsmm$control$timeVar
    data.long <- Objectlsmm$control$data.long
    value.var <- as.character(formFixed[[2]])
    pred.CV <- re$cv.Pred[,3]
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
    oldpar <- graphics::par(no.readonly = TRUE) # code line i
    on.exit(graphics::par(oldpar)) # code line i + 1
    k <- ggplot2::ggplot(df,  ggplot2::aes(obstime.mean, mean.obs, ymin = IC.sup, ymax = IC.inf))
    graph.fit.long <- k +  ggplot2::geom_pointrange( ggplot2::aes(ymin = IC.sup, ymax = IC.inf), shape =1) +
      ggplot2::geom_point(ggplot2::aes(obstime.mean, mean.pred), size = 3, shape = 17) +
      ggplot2::scale_x_continuous(name = "Time") +
      ggplot2::scale_y_continuous(name = "Current Value") +
      ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                     panel.background = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"))+
      ggplot2::ggtitle("Longitudinal goodness-of-fit")
    graph <- graph.fit.long
  }

  if(which == 'traj.ind'){
    if(is.null(ID.ind)){
      stop("You have to design some individual ID to plot the the individual trajectories.")
    }
    ID.ind <- as.vector(ID.ind)
    pred.CV <- as.data.frame(re$cv.Pred)
    data.long <- Objectlsmm$control$data.long
    formFixed <- Objectlsmm$control$formFixed
    value.var <- as.character(formFixed[[2]])
    graph.traj.ind <- c()
    for(ind in ID.ind){
      pred.CV.id <- pred.CV[which(pred.CV$id == ind),]
      pred.CV.id$y <- data.long[which(data.long$id == ind), value.var]
      pred.CV.id$CI.sup <- pred.CV.id$CV + 1.96*sqrt(pred.CV.id$Residual_SD_inter**2 + pred.CV.id$Residual_SD_intra**2)
      pred.CV.id$CI.inf <- pred.CV.id$CV - 1.96*sqrt(pred.CV.id$Residual_SD_inter**2 + pred.CV.id$Residual_SD_intra**2)

      traj_ind <- ggplot() +

        geom_line(pred.CV.id, mapping = aes(x=time, y=CV, group = id, color = 'Predicted'))+
        geom_line(pred.CV.id, mapping = aes(x=time, y=CI.sup, group = id,  color = 'Predicted'))+
        geom_line(pred.CV.id, mapping = aes(x=time, y=CI.inf, group = id,  color = 'Predicted'))+

        geom_ribbon( pred.CV.id,mapping=
                       aes(x=time,ymin=CI.inf,ymax=CI.sup), fill="#998ec3", alpha=0.3)+

        geom_point(pred.CV.id, mapping = aes(x=time, y=y, group = id,color = "Observed"),shape =17)+
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
          axis.title.x = element_text(color = "black", size = 10),
          axis.title.y = element_text(color = "black", size = 10),
          panel.grid = element_blank(),
          #legend.key = element_blank(),
          legend.text = element_text(color = "black", size = 10),
          axis.line = element_line(color = "black",
                                   linetype = "solid"),
          axis.text = element_text(size = 10, color = "black")
        )

      graph.traj.ind <- c(graph.traj.ind, traj_ind)
    }
    graph <- graph.traj.ind

  }


  graph

}
