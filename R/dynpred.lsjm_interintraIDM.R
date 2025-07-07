#' @rdname dynpred
#' @export

dynpred.lsjm_interintraIDM <- function(Objectlsjm,newdata,  s, horizon, IC = 95, nb.draws = 1000){

  event <- 1
  if(!is.null(IC) && (IC<=0 || IC>=100)) stop("IC must be between 0 and 100")
  if(!is.null(IC) && (is.null(nb.draws) || nb.draws <=0)) stop("draw must be higher 1")
  if(Objectlsjm$result_step1$istop != 1|| (!is.null(Objectlsjm$result_step2) && Objectlsjm$result_step2$istop !=1)){
    stop("The model didn't reach convergence.")
  }

  id <- as.integer(newdata[all.vars(Objectlsjm$control$Objectlsmm$control$formGroup)][,1])
  newdata$id <- id
  times <- horizon
  window <- times-s
  table.pred <- c()

  for(i in unique(id)){
    newdata.id <- newdata[which(newdata$id == i),]
    newdata.id <- as.data.frame(newdata.id)
    data.long.until.time.s <-subset(newdata.id, get(Objectlsjm$control$Objectlsmm$control$timeVar)<=s)
    pred.boot <- c()
    pred.ponct <- c()
    for(t in window){
      pred.ponct <- c(pred.ponct,predyn_ponct_lsjm_interintraIDM(Objectlsjm,  data.long.until.time.s, s, t, event))
      if(!is.null(IC)){
        pred.boot <- cbind(pred.boot, predyn_boot_lsjm_interintraIDM(Objectlsjm, data.long.until.time.s, s, t, event, nb.draws) )
      }
    }
    if(!is.null(IC)){
      table.pred.id <- cbind(i, times, pred.ponct, apply(pred.boot,2, function(x) quantile(x, 0.50)),
                             apply(pred.boot,2, function(x) quantile(x, ((100-IC)/2)/100)),
                             apply(pred.boot,2, function(x) quantile(x, 1-((100-IC)/2)/100)),
                             apply(pred.boot,2, sd)
      )
      table.pred.id <- as.data.frame(table.pred.id)
      colnames(table.pred.id) <- c("ID","Time","Prediction","Median","ICinf","ICsup", "Empirical SD")
    }
    else{
      table.pred.id <- cbind(i, times, pred.ponct
      )
      table.pred.id <- as.data.frame(table.pred.id)
      colnames(table.pred.id) <- c("ID","Time","Prediction")
    }

    table.pred <- rbind(table.pred, table.pred.id)

    ### Graph
    if(!is.null(IC)){
      oldpar <- graphics::par(no.readonly = TRUE) # code line i
      on.exit(graphics::par(oldpar)) # code line i + 1
      #browser()
      x.axe <- c(0,data.long.until.time.s[,Objectlsjm$control$Objectlsmm$control$timeVar],times)
      #print(x.axe)
      y.axe <- c(NA,data.long.until.time.s[,all.vars(Objectlsjm$control$Objectlsmm$control$formFixed)[1]], rep(NA,length(window)))
      #print(y.axe)
      y.axe2 <- c(NA,rep(NA,length(data.long.until.time.s[,Objectlsjm$control$Objectlsmm$control$timeVar])),table.pred.id$ICinf)
      #print(y.axe2)
      y.axe3 <- c(NA,rep(NA,length(data.long.until.time.s[,Objectlsjm$control$Objectlsmm$control$timeVar])),table.pred.id$Median)
      #print(y.axe3)
      y.axe4 <- c(NA,rep(NA,length(data.long.until.time.s[,Objectlsjm$control$Objectlsmm$control$timeVar])),table.pred.id$ICsup)
      y.axe5 <- c(NA,rep(NA,length(data.long.until.time.s[,Objectlsjm$control$Objectlsmm$control$timeVar])),table.pred.id$Prediction)
      graphics::plot(x = x.axe, y = y.axe,xlim = c(0,max(s+window)),
                     xlab = "Time", ylab = "Marker", cex.lab = 1, col = "black",
                     main = "Prediction of event", pch = 20, cex = 1, font = 1, font.lab = 1, cex.lab = 1, cex.main = 1)
      graphics::par(new = TRUE, font = 1, cex.lab = 1)
      graphics::plot(x.axe, y.axe3, axes = FALSE,col = "black", type = "l", ylim = c(0.000001, max(table.pred.id$ICsup, na.rm = T)),ylab = "", xlab = "", lwd =2, font.lab = 1, cex.lab = 1 )
      graphics::lines(x.axe, y.axe5, col = "red", lty=1, lwd = 2)
      graphics::lines(x.axe, y.axe2, col = "black", lty=2, lwd = 2)
      graphics::lines(x.axe, y.axe4, col = "black", lty=2, lwd = 2)
      graphics::axis(side= 4, cex = 2)
      graphics::abline(v = s, lty = 3)
      graphics::mtext("Probability of event", cex = 1, side = 4, line = 3, font.lab = 1)

    }

  }
  table.pred <- as.data.frame(table.pred)
  if(!is.null(IC)){
    colnames(table.pred) <- c("ID","Time","Prediction","Median","ICinf","ICsup", "Empirical SD")
  }
  else{
    colnames(table.pred) <- c("ID","Time","Prediction")
  }

  table.pred

}
