#' @rdname plot
#' @import ggplot2
#' @export
#'

plot.lsjm_interintraIDM <- function(Objectlsjm, which = 'long.fit', Objectpredict = NULL, break.times = NULL, ID.ind = NULL, ObjectSmoothHazard = NULL, xlim = NULL, ylim = NULL){


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
    #prediction <- cbind(pred.CV, data.long$window)
    ObjectpredictY$time.new.pred <- ObjectpredictY$time
    data.long$time.new.pred <- data.long[,timeVar]
    prediction <- dplyr::left_join(ObjectpredictY[,c("id","predY", "time.new.pred")], data.long[,c("id", "window", "time.new.pred")])
    mean.pred <- by(prediction$predY, prediction$window, mean)
    obstime.mean <- by(data.long[,timeVar], data.long$window, mean)
    df <- cbind(obstime.mean, mean.obs, IC.sup, IC.inf, mean.pred)
    df <- as.data.frame(df)
    k <- ggplot2::ggplot(df,  ggplot2::aes(obstime.mean, mean.obs, ymin = IC.sup, ymax = IC.inf))
    graph.fit.long <- k +  ggplot2::geom_pointrange( ggplot2::aes(ymin = IC.sup, ymax = IC.inf), shape =1) +
      ggplot2::geom_point(ggplot2::aes(obstime.mean, mean.pred), size = 3, shape = 17) +
      ggplot2::scale_x_continuous(name = "Time") +
      ggplot2::scale_y_continuous(name = "Current Value") +
      ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                     panel.background = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"),
                     axis.text=ggplot2::element_text(size=15),
                     axis.title=ggplot2::element_text(size=18),
                     plot.title = ggplot2::element_text(size = 20, face = "bold"))+
      ggplot2::ggtitle("Longitudinal goodness-of-fit")
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
      pred.CV.id$CI.sup <- pred.CV.id$CV + 1.96*sqrt(pred.CV.id$Residual_SD_inter**2 + pred.CV.id$Residual_SD_intra**2)
      pred.CV.id$CI.inf <- pred.CV.id$CV - 1.96*sqrt(pred.CV.id$Residual_SD_inter**2 + pred.CV.id$Residual_SD_intra**2)

      traj_ind <- ggplot2::ggplot() +

        ggplot2::geom_line(pred.CV.id, mapping = aes(x=time, y=CV, group = id, color = 'Predicted'))+
        ggplot2::geom_line(pred.CV.id, mapping = aes(x=time, y=CI.sup, group = id,  color = 'Predicted'))+
        ggplot2::geom_line(pred.CV.id, mapping = aes(x=time, y=CI.inf, group = id,  color = 'Predicted'))+

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
      print(traj_ind)
    }
    graph <- graph.traj.ind

  }
  if(which == 'survival.fit'){
    #Objectranef$grid.time.Cum
    Cum_01Smooth_est <- intensity(times = Objectpredict$grid.time.Cum, knots = ObjectSmoothHazard$knots01,
                                  number.knots = ObjectSmoothHazard$nknots01,
                                  theta = ObjectSmoothHazard$theta01^2)

    Cum_02Smooth_est <- intensity(times = Objectpredict$grid.time.Cum,knots = ObjectSmoothHazard$knots02,
                                  number.knots = ObjectSmoothHazard$nknots02,
                                  theta = ObjectSmoothHazard$theta02^2)

    Cum_12Smooth_est <- intensity(times = Objectpredict$grid.time.Cum, knots = ObjectSmoothHazard$knots12,
                                  number.knots = ObjectSmoothHazard$nknots12,
                                  theta = ObjectSmoothHazard$theta12^2)
    V <- ObjectSmoothHazard$V
    Cum_01.cum <- c()
    Cum_02.cum <- c()
    Cum_12.cum <- c()
    for(boot in 1:5000){
      #browser()
      tirage <- mvtnorm::rmvnorm(1, mean = c(ObjectSmoothHazard$theta01,ObjectSmoothHazard$theta02,ObjectSmoothHazard$theta12), sigma = V)
      Cum_01Smooth <- intensity(times = Objectpredict$grid.time.Cum, knots = ObjectSmoothHazard$knots01,
                                number.knots = ObjectSmoothHazard$nknots01,
                                theta = tirage[1:(ObjectSmoothHazard$nknots01+2)]^2)
      Cum_02Smooth <- intensity(times = Objectpredict$grid.time.Cum, knots = ObjectSmoothHazard$knots02,
                                number.knots = ObjectSmoothHazard$nknots02,
                                theta = tirage[(ObjectSmoothHazard$nknots02+3):(2*(ObjectSmoothHazard$nknots02+2))]^2)
      Cum_12Smooth <- intensity(times = Objectpredict$grid.time.Cum, knots = ObjectSmoothHazard$knots12,
                                number.knots = ObjectSmoothHazard$nknots12,
                                theta = tirage[(2*(ObjectSmoothHazard$nknots12+2)+1):(3*(ObjectSmoothHazard$nknots12+2))]^2)

      Cum_01.cum <- rbind(Cum_01.cum, Cum_01Smooth$cumulative.intensity)
      Cum_02.cum <- rbind(Cum_02.cum, Cum_02Smooth$cumulative.intensity)
      Cum_12.cum <- rbind(Cum_12.cum, Cum_12Smooth$cumulative.intensity)
    }

    Cum_01.quant <- apply(Cum_01.cum,2,quantile, probs = c(0.025,0.5,0.975))
    Cum_02.quant <- apply(Cum_02.cum,2,quantile, probs = c(0.025,0.5,0.975))
    Cum_12.quant <- apply(Cum_12.cum,2,quantile, probs = c(0.025,0.5,0.975))

    predictCum_01 <- Objectpredict$predictCum_01
    colnames(predictCum_01) <- c("ID", paste0("t", 1:(ncol(predictCum_01)-1)))
    predictCum_02 <- Objectpredict$predictCum_02
    colnames(predictCum_02) <- c("ID", paste0("t", 1:(ncol(predictCum_02)-1)))
    predictCum_12 <- Objectpredict$predictCum_12
    colnames(predictCum_12) <- c("ID", paste0("t", 1:(ncol(predictCum_12)-1)))
    predictCum_01 <- as.data.frame(predictCum_01)
    predictCum_02 <- as.data.frame(predictCum_02)
    predictCum_12 <- as.data.frame(predictCum_12)
    data.long <- Objectlsmm$control$data.long
    data.id <- data.long[!duplicated(data.long$id),]
    Time.R.var <- as.character(Objectlsjm$control$Time$Time_R[[2]])
    Time.L.var <- as.character(Objectlsjm$control$Time$Time_L[[2]])
    Time.T.var <- as.character(Objectlsjm$control$Time$Time_T[[2]])
    event1.var <- as.character(Objectlsjm$control$deltas$delta1[[2]])
    Cum_01_pred <- c()
    Cum_02_pred <- c()
    Cum_12_pred <- c()
    for(i in 2:ncol(predictCum_01)){
      id.01 <- data.id$id[which(data.id[,Time.R.var]>Objectpredict$grid.time.Cum[i])]
      id.02 <- data.id$id[which(data.id[,Time.T.var]>Objectpredict$grid.time.Cum[i] & (data.id[,event1.var] == 0 | (data.id[,event1.var] == 1 & data.id[,Time.R.var]>Objectpredict$grid.time.Cum[i])))]
      id.12 <- data.id$id[which(data.id[,Time.T.var]>Objectpredict$grid.time.Cum[i] & data.id[,Time.L.var]<Objectpredict$grid.time.Cum[i] & data.id[,event1.var] == 1)]
      Cum_01_pred <- c(Cum_01_pred, mean(predictCum_01[which(predictCum_01$ID %in% id.01),i]))
      Cum_02_pred <- c(Cum_02_pred, mean(predictCum_02[which(predictCum_02$ID %in% id.02),i]))
      Cum_12_pred <- c(Cum_12_pred, mean(predictCum_12[which(predictCum_01$ID %in% id.12),i]))
    }


    #Cum_01_pred <- apply(Objectpredict$predictCum_01,2,mean)
    #Cum_02_pred <- apply(Objectpredict$predictCum_02,2,mean)
    #Cum_12_pred <- apply(Objectpredict$predictCum_12,2,mean)
#

    data_tot <- cbind(Objectpredict$grid.time.Cum, Cum_01Smooth_est$cumulative.intensity, Cum_02Smooth_est$cumulative.intensity, Cum_12Smooth_est$cumulative.intensity,
                      Cum_01.quant[1,],Cum_01.quant[2,],Cum_01.quant[3,],
                      Cum_02.quant[1,],Cum_02.quant[2,],Cum_02.quant[3,],
                      Cum_12.quant[1,],Cum_12.quant[2,],Cum_12.quant[3,],
                      Cum_01_pred,Cum_02_pred,Cum_12_pred)

    data_tot <- as.data.frame(data_tot)

    colnames(data_tot) <- c("time", "Cum_01_est", "Cum_02_est", "Cum_12_est",
                            "Cum_01_2.5", "Cum_01_5", "Cum_01_97.5",
                            "Cum_02_2.5", "Cum_02_5", "Cum_02_97.5",
                            "Cum_12_2.5", "Cum_12_5", "Cum_12_97.5",
                            "Cum_01_pred", "Cum_02_pred", "Cum_12_pred")

    if(!is.null(xlim)){
      data_tot <- data_tot[which(data_tot$time < xilim[2]),]
    }


    survB <- ggplot() +
      ggplot2::geom_line(data = data_tot, mapping = aes(x = time, y = Cum_01_pred, color = "Cum_01_pred", linetype = "Joint Model"), linewidth = 1.5) +
      ggplot2::geom_line(data = data_tot, mapping = aes(x = time, y = Cum_02_pred, color = "Cum_02_pred", linetype = "Joint Model"), linewidth = 1.5) +
      ggplot2::geom_line(data = data_tot, mapping = aes(x = time, y = Cum_12_pred, color = "Cum_12_pred", linetype = "Joint Model"), linewidth = 1.5) +
      ggplot2::geom_ribbon(data = data_tot, mapping = aes(x = time, ymin = Cum_01_2.5, ymax = Cum_01_97.5, fill = "Cum_01_pred"), alpha = 0.2) +
      ggplot2::geom_ribbon(data = data_tot, mapping = aes(x = time, ymin = Cum_02_2.5, ymax = Cum_02_97.5, fill = "Cum_02_pred"), alpha = 0.2) +
      ggplot2::geom_ribbon(data = data_tot, mapping = aes(x = time, ymin = Cum_12_2.5, ymax = Cum_12_97.5, fill = "Cum_12_pred"), alpha = 0.2) +
      ggplot2::geom_line(data = data_tot, mapping = aes(x = time, y = Cum_01_est, color = "Cum_01_pred", linetype = "Smooth Hazard"), linewidth = 1.2) +
      ggplot2::geom_line(data = data_tot, mapping = aes(x = time, y = Cum_02_est, color = "Cum_02_pred", linetype = "Smooth Hazard"), linewidth = 1.2) +
      ggplot2::geom_line(data = data_tot, mapping = aes(x = time, y = Cum_12_est, color = "Cum_12_pred", linetype = "Smooth Hazard"), linewidth = 1.2) +
      xlab("Age (years)") +
      ylab("Cumulative intensity functions") +
      ggplot2::ggtitle("B")+
      ggplot2::scale_color_manual(values = c("Cum_01_pred" = "#f1a340", "Cum_02_pred" = "#998ec3", "Cum_12_pred" = "palegreen3"
      ),
      labels = c("Cum_01_pred" = "Transition 0 -> 1",
                 "Cum_02_pred" = "Transition 0 -> 2",
                 "Cum_12_pred" = "Transition 1 -> 2")) +
      ggplot2::scale_fill_manual(values = c("Cum_01_pred" = "#f1a340", "Cum_02_pred" = "#998ec3", "Cum_12_pred" = "palegreen3"),
                                 labels = c("Cum_01_pred" = "Transition 0 -> 1",
                                            "Cum_02_pred" = "Transition 0 -> 2",
                                            "Cum_12_pred" = "Transition 1 -> 2")) +
      ggplot2::scale_linetype_manual(values = c("Joint Model" = 3, "Smooth Hazard" = 1)) +

      ggplot2::theme(
        panel.background = element_blank(),
        legend.position = "bottom",
        legend.box = "vertical",
        axis.title.x = element_text(color = "black", size = 16),
        axis.title.y = element_text(color = "black", size = 16),
        panel.grid = element_blank(),
        legend.text = element_text(color = "black", size = 16),
        axis.line = element_line(color = "black", linetype = "solid"),
        axis.text = element_text(size = 14, color = "black")
      ) +
      ggplot2::guides(color = guide_legend(title = "", nrow = 1, byrow = TRUE,keywidth = 4),
                      fill = guide_legend(title = "",keywidth = 4),
                      linetype = guide_legend(title = "", keywidth = 4, keyheight = 1))


    print(survB)

    graph <- survB

  }


  graph

}
