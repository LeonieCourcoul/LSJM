#' @import ggplot2
#' @export
#'


plot.lsmm_interintra <- function(object, which = 'long.fit', predictObject = NULL, break.times = NULL, ID.ind = NULL, ylim = NULL, xlim= NULL){

  if(is.null(predictObject)){
    stop("predictObject is missing.")
  }

  graph <- NULL

  oldpar <- graphics::par(no.readonly = TRUE) # code line i
  on.exit(graphics::par(oldpar)) # code line i + 1

  if(which == 'long.fit'){
    formFixed <- object$control$formFixed
    timeVar <- object$control$timeVar
    data.long <- object$control$data.long
    value.var <- as.character(formFixed[[2]])
    pred.CV <- predictObject$predY
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
    window.pred <- cut(predictObject$time, break.times, include.lowest = T)
    prediction <- cbind(pred.CV, window.pred)
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
    pred.CV <- as.data.frame(predictObject)
    data.long <- object$control$data.long
    formFixed <- object$control$formFixed
    value.var <- as.character(formFixed[[2]])
    graph.traj.ind <- c()
    for(ind in ID.ind){
      pred.CV.id <- pred.CV[which(pred.CV$id == ind),]
      data.idselect <- cbind(data.long$id[which(data.long$id == ind)],data.long[which(data.long$id == ind), object$control$timeVar],data.long[which(data.long$id == ind), value.var])
      data.idselect <- as.data.frame(data.idselect)
      colnames(data.idselect) <- c("id","time", "y")
      #pred.CV.id$y <- data.long[which(data.long$id == ind), value.var]
      pred.CV.id$CI.sup <- pred.CV.id$predY + 1.96*sqrt(pred.CV.id$predSD_inter**2 + pred.CV.id$predSD_intra**2)
      pred.CV.id$CI.inf <- pred.CV.id$predY - 1.96*sqrt(pred.CV.id$predSD_inter**2 + pred.CV.id$predSD_intra**2)

      traj_ind <- ggplot2::ggplot() +

        ggplot2::geom_line(pred.CV.id, mapping = aes(x=time, y=predY, group = id, color = 'Predicted'))+
        ggplot2::geom_line(pred.CV.id, mapping = aes(x=time, y=CI.sup, group = id,  color = 'Predicted'),linetype = 2)+
        ggplot2::geom_line(pred.CV.id, mapping = aes(x=time, y=CI.inf, group = id,  color = 'Predicted'),linetype = 2)+

        ggplot2::geom_ribbon( pred.CV.id,mapping=
                       aes(x=time,ymin=CI.inf,ymax=CI.sup), fill="#998ec3", alpha=0.3)+

        ggplot2::geom_point(data.idselect, mapping = aes(x=time, y=y, group = id,color = "Observed"),shape =17)+
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


  graph

}
