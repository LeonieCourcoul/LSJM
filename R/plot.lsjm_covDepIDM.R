#' @rdname plot.lsjm
#' @importFrom graphics plot par
#' @importFrom dplyr left_join
#' @importFrom ggplot2 ggplot aes geom_pointrange geom_point scale_x_continuous scale_y_continuous theme element_blank element_line element_text ggtitle coord_cartesian geom_line geom_ribbon facet_wrap scale_color_manual guide_legend guides scale_fill_manual geom_step scale_linetype_manual
#' @importFrom survminer surv_fit ggsurvplot
#' @importFrom survival Surv
#' @importFrom SmoothHazard intensity
#' @importFrom mvtnorm rmvnorm
#' @export
plot.lsjm_covDepIDM <- function(x, which = 'long.fit', Objectpredict, break.times = NULL, ID.ind = NULL, ObjectSmoothHazard  = NULL, xlim = NULL, ylim = NULL, ...){


  Objectlsjm <- x
  Objectlsmm <- Objectlsjm$control$Objectlsmm
  if(is.null(Objectpredict)){
    stop("Not implemented")
  }

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
                     axis.text=element_text(size=15),
                     axis.title=element_text(size=18),
                     plot.title = element_text(size = 20, face = "bold"))+
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
      pred.CV.id$y <- data.long[which(data.long$id == ind), value.var]
      pred.CV.id$CI.sup <- pred.CV.id$predY + 1.96*pred.CV.id$predSD
      pred.CV.id$CI.inf <- pred.CV.id$predY - 1.96*pred.CV.id$predSD

      #browser()
      traj_ind <- ggplot() +

        geom_line(pred.CV.id, mapping = aes(x=time, y=predY, group = id, color = 'Predicted'))+
        geom_line(pred.CV.id, mapping = aes(x=time, y=CI.sup, group = id,  color = 'Predicted'),linetype = 2)+
        geom_line(pred.CV.id, mapping = aes(x=time, y=CI.inf, group = id,  color = 'Predicted'),linetype = 2)+

        geom_ribbon( pred.CV.id,mapping=
                                aes(x=time,ymin=CI.inf,ymax=CI.sup), fill="#998ec3", alpha=0.3,linetype = 3)+

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
          #axis.title.x = element_text(color = "black", size = 10),
          #axis.title.y = element_text(color = "black", size = 10),
          panel.grid = element_blank(),
          #legend.key = element_blank(),
          axis.text=element_text(size=15),
          axis.title=element_text(size=18),
          plot.title = element_text(size = 20, face = "bold"),
          legend.text = element_text(color = "black", size = 14),
          axis.line = element_line(color = "black",
                                   linetype = "solid")
          #axis.text = element_text(size = 10, color = "black")
        )+coord_cartesian(xlim = xlim,ylim = ylim, expand = TRUE)

      graph[[paste("traj.ind",ind, sep = "_")]] <- traj_ind
    }

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
      tirage <- rmvnorm(1, mean = c(ObjectSmoothHazard$theta01,ObjectSmoothHazard$theta02,ObjectSmoothHazard$theta12), sigma = V)
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

    Cum_01_pred <- apply(Objectpredict$predictCum_01,2,mean)
    Cum_02_pred <- apply(Objectpredict$predictCum_02,2,mean)
    Cum_12_pred <- apply(Objectpredict$predictCum_12,2,mean)

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
      geom_line(data = data_tot, mapping = aes(x = time, y = Cum_01_pred, color = "Cum_01_pred", linetype = "Joint Model"), linewidth = 1.5) +
      geom_line(data = data_tot, mapping = aes(x = time, y = Cum_02_pred, color = "Cum_02_pred", linetype = "Joint Model"), linewidth = 1.5) +
      geom_line(data = data_tot, mapping = aes(x = time, y = Cum_12_pred, color = "Cum_12_pred", linetype = "Joint Model"), linewidth = 1.5) +
      geom_ribbon(data = data_tot, mapping = aes(x = time, ymin = Cum_01_2.5, ymax = Cum_01_97.5, fill = "Cum_01_pred"), alpha = 0.2) +
      geom_ribbon(data = data_tot, mapping = aes(x = time, ymin = Cum_02_2.5, ymax = Cum_02_97.5, fill = "Cum_02_pred"), alpha = 0.2) +
      geom_ribbon(data = data_tot, mapping = aes(x = time, ymin = Cum_12_2.5, ymax = Cum_12_97.5, fill = "Cum_12_pred"), alpha = 0.2) +
      geom_line(data = data_tot, mapping = aes(x = time, y = Cum_01_est, color = "Cum_01_pred", linetype = "Smooth Hazard"), linewidth = 1.2) +
      geom_line(data = data_tot, mapping = aes(x = time, y = Cum_02_est, color = "Cum_02_pred", linetype = "Smooth Hazard"), linewidth = 1.2) +
      geom_line(data = data_tot, mapping = aes(x = time, y = Cum_12_est, color = "Cum_12_pred", linetype = "Smooth Hazard"), linewidth = 1.2) +
      xlab("Age (years)") +
      ylab("Cumulative intensity functions") +
      ggtitle("B")+
      scale_color_manual(values = c("Cum_01_pred" = "#f1a340", "Cum_02_pred" = "#998ec3", "Cum_12_pred" = "palegreen3"
      ),
      labels = c("Cum_01_pred" = "Transition 0 -> 1",
                 "Cum_02_pred" = "Transition 0 -> 2",
                 "Cum_12_pred" = "Transition 1 -> 2")) +
      scale_fill_manual(values = c("Cum_01_pred" = "#f1a340", "Cum_02_pred" = "#998ec3", "Cum_12_pred" = "palegreen3"),
                        labels = c("Cum_01_pred" = "Transition 0 -> 1",
                                   "Cum_02_pred" = "Transition 0 -> 2",
                                   "Cum_12_pred" = "Transition 1 -> 2")) +
      scale_linetype_manual(values = c("Joint Model" = 3, "Smooth Hazard" = 1)) +

      theme(
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
      guides(color = guide_legend(title = "", nrow = 1, byrow = TRUE,keywidth = 4),
             fill = guide_legend(title = "",keywidth = 4),
             linetype = guide_legend(title = "", keywidth = 4, keyheight = 1))


    print(survB)

    graph <- survB

  }


  graph

}
