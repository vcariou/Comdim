#' @title PathComdim: Finding common dimensions in multitable data (Xk, k=1...K) related together with a path diagram
#' @description:
#' Path-ComDim aims at finding common dimensions from multibloc datasets, where the datasets at hand are assumed to have a
#' specific pattern of directed relations among them reflecting, for instance, a chain of influence.
#' The aim of Path-ComDim is to analyze these datasets taking into account the structural connections among the
#' multiblock data (Xk, k=1...K) related together with such a path diagram
#' @usage PathComdim(X,group,delta,ndim,scale,threshold,init)
#'
#' @param X :         concatenated matrices (or dataframes) (Xk) with nrow. All (Xk)s share the same rows
#' @param group :     vector containing the number of variables for each table of X with names of the tables
#' @param delta :     matrix of size K x K depicting the directed links from Xj to Xk (j,k=1, ..., K)
#' @param ndim :      number of common dimensions
#' @param scale :    scaling of variables with standardization (default FALSE)

#' @param threshold : if the difference of fit<threshold then break the iterative loop (default 1E-10)
#' @param init :  way to initialize the estimation loop, either 0 (random initialization), or 1 (average initialization), (default 0) \cr
#'
#' @return \item{group}{ input parameter group }
#' @return \item{mean}{mean values associated to the various variables}
#' @return \item{Xscale}{scaling values associated to the various variables}
#' @return \item{t}{common components associated to explanation (nrow x ndim)}
#' @return \item{u}{common components associated to prediction (nrow x ndim)}
#' @return \item{Tk}{block  components associated to explanation and related to t}
#' @return \item{Ul}{block components associated to prediction and related to u}
#' @return \item{saliences}{saliences associated to each directed link and for each dimension. It indicates the strenght of the link for the current dimension}
#' @return \item{weights}{weights to determinate the common component t at each current dimension}
#' @return \item{Wstar}{weights to determinate the common component t, regardless the deflation stage}
#' @return call : call of the method
#' @seealso ....
#'
#' @examples
#' data(RedLosses)
#' res <- PathComdim(data,Group,delta,20)
#'
#' @author  {Veronique Cariou <veronique.cariou@oniris-nantes.fr>}
#' @references {Cariou, V., Qannari, E. M., Rutledge, D. N., & Vigneau, E. (2018). ComDim: From multiblock data analysis to path modeling. Food Quality and Preference, 67, 27-34.}
#'
PathComdim=function(X,group,delta,ndim, scale=FALSE, threshold=1E-10, init=0){
  # ---------------------------------------------------------------------------
  # 0. Preliminary tests
  # ---------------------------------------------------------------------------
  if (any(is.na(X)))
    stop("No NA values are allowed")
  if (is.null(rownames(X)))
    rownames(X) <- 1:dim(X)[1]
  if (is.null(colnames(X)))
    colnames(X) <- paste("X", 1:dim(X)[2], sep=".")
  if (is.null(names(group)))
    names(group) <- paste("Tab", 1:length(group), sep=".")
  if (sum(group) != ncol(X))
    stop("Non convenient blocks parameter")
  if (init !=0 && init!=1 )
    stop("Non convenient init parameter for common components initialization")


  # ---------------------------------------------------------------------------
  # 1. Initialization
  # ---------------------------------------------------------------------------

  res=list()
  ntab=length(group)                      # Number of blocks
  nind=nrow(X)                            # # of individuals
  Lambda=array(0,dim=c(ntab,ntab,ndim));  # array of saliences
  W=array(0,dim=c(nind,nind,ntab));       # association matrices
  MatU=matrix(0,nrow=nind,ncol=ndim)      # common components associated with Xl
  MatT=T.new=matrix(0,nrow=nind,ncol=ndim)      # common components associated with Xk
  Weights= Px= matrix(0,nrow=ncol(X),ncol=ndim)
  Tk=array(0, dim=c(nind,ntab,ndim))      # array of block components associated with Xk from MatT
  Ul=array(0, dim=c(nind,ntab,ndim))      # array of block components associated with Xl from MatU

  J=rep(1:ntab , times =  group )         # J vector indicates which block each variable belongs to
  list.k<-which(apply(delta,1,sum)>0)     # Ingoing  arrows
  list.l<-which(apply(delta,2,sum)>0)     # Outgoing arrows
  Tk.weight <- vector("list",length=ntab)

  # ---------------------------------------------------------------------------
  # 2. Scaling - by default scaling put each table to unit norm
  # ---------------------------------------------------------------------------
  X0 <- X
  Xj=list()
  #X contains centered and optionnally scaled data tables and automatically unit total variance
  Xmeans <- colMeans(X)
  X <- X - rep(Xmeans, each = nind)                  # default centering of the dataset table
  if (scale) {
    scale <- sqrt(colSums((X - rep(colMeans(X), each = nind))^2) /
                    (nind-1))
    temp      <- abs(scale) < 1e-14
    if (any(temp)) {
      warning("Scaling with (near) zero standard deviation leads no standardization for some variables")
      scale[temp] <- 1
    }
    X <- X / rep(scale, each = nind)
  }else{
    scale <- rep(1,times=length(J))
  }
  # set each table to unit variance
  for(j in 1:ntab){
    ivar<-which(J==j)
    sc=sqrt(inertie(X[,ivar]))
    X[,ivar]=X[,ivar]/sc
    scale[ivar] <- scale[ivar]*sc
  }
  Xscale <- X

  # Computation of association matrices
  for(j in 1:ntab) {
    Xj[[j]]=as.matrix(X[,J==j])           # contains matrix Xj (list)
    W[,,j]=Xj[[j]]%*%t(Xj[[j]])
  }
  n.slices <- sum(delta)

  W_kl <- vector("list",n.slices)
  h <-1
  delta.vec<-matrix(0,ncol=2,nrow=n.slices)
  for (k in list.k) { #predictor
    for (l in list.l) { #to predict
      if ((delta[k,l]!=0)){
        W_kl[[h]]<- W[,,k] %*%W[,,l]
        delta.vec[h,]<-c(k,l)
        h<-h+1
      } #  if
    } # for k
  } # for l

  # ----------------------------------------------------------------------------
  # 3. Main loop on the dimensions
  # ----------------------------------------------------------------------------

  for (dimension in 1:ndim)  {            # determination of the common vs block components
    if (init==0) {
      #random initial t
      t=rnorm(nind)
      t=t1=normv(t)                            # set the global component t to norm 1

      #random initial u
      u=rnorm(nind)
      u=normv(u)                            # set the global component t to norm 1
    } else {
      #computation of an initial mean component t
      t=apply(matrix(unlist(lapply(W_kl,function(x) {apply(x,1,mean)})),nrow=nind,ncol=n.slices),1,mean) #for t
      t=t1=normv(t)                            # set the global component t to norm 1

      #computation of an initial mean component u
      u=apply(matrix(unlist(lapply(W_kl,function(x) {apply(x,2,mean)})),nrow=nind,ncol=n.slices),1,mean)
      u=normv(u)                            # set the global component t to norm 1
    }

    deltacrit=1
    crit=0

    # ----------------------------------------------------------------------------
    # 3.1. Iteration sub-loop for the current dimension
    # ----------------------------------------------------------------------------

    while(deltacrit>threshold) {
      # Updating u
      U=matrix(unlist(lapply(W_kl,function(w_kl){crossprod(w_kl,t)})),nrow=nind)
      u=U%*%t(U)%*%u
      u=normv(u)                          # u updated

      # Updating t
      T         <- matrix(unlist(lapply(W_kl,function(w_kl){w_kl%*%u})),nrow=nind)
      tk.superweight <- t(T)%*%t
      t         <- T%*%tk.superweight
      t.in      <- sqrt(sum(t^2))
      t=t/t.in                       # t updated
      tk.superweight <- tk.superweight/t.in


      # lambda_cur : saliences
      lambda_cur <- unlist(lapply(W_kl,function(w_kl){crossprod(t,w_kl%*%u)}))
      critnew=sum(lambda_cur^2)
      deltacrit=critnew-crit
      crit=critnew
    } # end while

    # ----------------------------------------------------------------------------
    # 3.2. Storing the results corresponding to the current dimension
    # ----------------------------------------------------------------------------

    MatT[,dimension]=t
    MatU[,dimension]=u

    kl<-1
    for (k in list.k) {
      for (l in list.l) {
        if (delta[k,l]!=0) {
          Lambda[k,l,dimension] <- lambda_cur[kl]
          kl <- kl+1
        }
      }
    }

    # Xk components
    for (k in 1:ntab) {
      Tk.weight[[k]]  <- rep(0,times=ncol(Xj[[k]]))
      for (l in list.l) {
        if (delta[k,l] != 0) {
          h <- which(delta.vec[,1]==k & delta.vec[,2]==l)
          Tk.weight[[k]]  <- Tk.weight[[k]]  + tk.superweight[h]*(t(Xj[[k]])%*%W[,,l]%*%u)
        }
      }
      T.new[,dimension] <-T.new[,dimension]+ Xj[[k]]%*%Tk.weight[[k]]
      Tk[,k,dimension] <- W[,,k]%*%t
    }
    Weights[,dimension]       <- unlist(sapply(1:ntab,function(k){Tk.weight[[k]]}))
    #   W[,h]                     <- W[,h] / sum(W[,h]^2)  # normalization optional and only required if global weights are set to length one
    Px[,dimension]      <- unlist(lapply(Xj,function(Xj.k) {crossprod(Xj.k,t)})) #mod1 see Tenenhaus P128


    # Xl components
    for (l in list.l)
      Ul[,l,dimension]=W[,,l]%*%u



    # ----------------------------------------------------------------------------
    # 3.3. Deflation stage
    # ----------------------------------------------------------------------------

    for(j in 1:ntab) {
      Xj[[j]]=deflation(as.matrix(X[,J==j]),t)
      W[,,j]=Xj[[j]]%*%t(Xj[[j]]);
      X[,J==j]=Xj[[j]] #updating  Xk
    }
    h <-1
    for (k in list.k) { #predictor
      for (l in list.l) { #to predict
        if (delta[k,l]!=0) {
          W_kl[[h]]<- W[,,k] %*%W[,,l]
          h<-h+1
        } #  if
      } # for k
    } # for l
  } # End main loop

  # ----------------------------------------------------------------------------
  # 4. Determination of the regression model for pure endogeneous blocks
  # ----------------------------------------------------------------------------

  # Yendo <- which(apply(delta,1,sum)==0)
  # if (is.null(Yendo)) {
  #   print("No endogeneous block")
  # }
  # else {
  #   res$Yendo <- Yendo
  #   res$Ybeta <- matrix(0,nrow=ndim,ncol=sum(group[Yendo]))
  #   Y <- as.matrix(X0[,which(J%in% Yendo)] - rep(Xmeans[which(J%in% Yendo)], each = nind))
  #   model.Yend <- lm( Y~ T.new-1)
  #   res$Ybeta <- matrix(model.Yend$coefficients,nrow=ndim)
  #   #       print((res$Ybeta))
  #   #       print(colnames(X)[J%in% Yendo])
  #   dimnames(res$Ybeta)[[1]]  <- 1:ndim
  #   dimnames(res$Ybeta)[[2]]  <- colnames(X)[J%in% Yendo]
  # }

  # ----------------------------------------------------------------------------
  # 5. Storage of the results
  # ----------------------------------------------------------------------------
#  res$T.new <- T.new

  res$group = group
  res$scale <- scale
  res$mean  <- Xmeans
#  res$Xscale <- Xscale
  res$Tk=Tk
  res$Ul=Ul
  res$t=MatT[,1:ndim]  #common components
  res$u=MatU[,1:ndim]
  res$saliences=Lambda[,,1:ndim]
  res$weights = Weights
  res$Wstar  <- Weights %*% solve(crossprod(Px,Weights),tol=1e-150) #cf p135 #modif vca 21/05/2019 #to modify whether singular matrix or not
  # ----------------------------------------------------------------------------
  # 6. Class PathComDim
  # ----------------------------------------------------------------------------

  res$call   <- match.call()
  ### rajouter les options par defaut
  class(res) <- c("PathComdim")


  return(res)
}


####################################################################################################
# computing the total variance of a dataset
inertie <-function(tab) {
  tab<- scale(tab, scale=FALSE)
  tab<-as.matrix(tab)
  # n=dim(tab)[1]
  V<-t(tab)%*%tab
  sum(diag(V))
}

#normalization of a vector
normv=function (x)
{
  normx = sqrt(t(x)%*% x)
  y = x/as.numeric(normx)
  return(y)
}

#deflation of a table with respect with a common component d
deflation=function (X, d)  {
  Y = X - d %*% t(d) %*% X
  return(Y)
}
