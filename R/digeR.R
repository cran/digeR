digeR <-
function()
{

  expData<- matrix()        # expression data (transposed,col=spot no, row =sample)
  cord<- matrix()           # spot coordinate data
  groupClass<- vector()     # group names
  group<- vector()          # group labels (factor)
  dataset<- list()          # list contains pair-wise expression data (only for >=3 groups)
  dataall<- matrix()        # expData and cord
  imagepath<- ""            # image path
  input<- list()            # list contain expData, cord, groupClass, group,dataset
  #################################################################################### 
  # function for find data file, read in and quit                                    #
  ####################################################################################
  gf.open<- function(h,...)
  {
    filepath<- gfile("Select the data file",filter = list("txt" = list(patterns = c("*.txt","*.TXT"))),type="open")
    if (filepath!="")
    {
      data.all<- read.table(filepath,sep="\t",header=TRUE)
      cord<<- data.all[,1:2]
      expData<<- t(data.all[,-c(1:2)])
      dataall<<- data.all
      colN<- rownames(expData)
      print(colN)
      for (i in 1:length(colN))
      {
        for (j in 1:10)
        {
          if(substring(colN[i],j,j)=="_")
          {
            group[i]<- substring(colN[i],1,j-1)
            break
          }
        }
      }
      groupClass<- unique.default(group)
      print(groupClass)
      for ( i in 1:length(groupClass))
      {
        dataset[[i]]<- expData[group==groupClass[i],]
      }
      # save(expData,cord,groupClass,group,dataset,dataall,file="input.RData")
      input<<- list(expData=expData,cord=cord,groupClass=groupClass,group=group,dataset=dataset)
    }
  }
  gf.quit<- function(h,..)
  {
    dispose(digeR.menu)
  }

  ####################################################################################
  # upload function for uploading the gel image                                      #
  ####################################################################################
  upload<- function(h,...)
  {
    imagepath<<- gfile("Select a jpg image",filter = list("JPEG" = list(patterns = c("*.jpg","*,jpeg"))),type="open")
  }

  ####################################################################################
  # whole profile plot PCA and PLSR                                                  #
  ####################################################################################
  whole.plot<- function(selected,inputData)   # selected = which method to use, inpuData = data
  {
    expData<- inputData$expData
    cord<- inputData$cord
    groupClass<- inputData$groupClass
    group<- inputData$group
    dataset<- inputData$dataset

    ScorePlot<- function(h,...)
    {
      if (svalue(scaling)=="No scaling"){m<- data.frame(expData)}
      if (svalue(scaling)=="Center only")  {m<- data.frame(scale(expData,center=TRUE,scale=FALSE))}
      if (svalue(scaling)=="Center and Scaling")  {m<- data.frame(scale(expData,center=TRUE,scale=TRUE))}
      if (svalue(scaling)=="Paretol Scaling")  {m<- data.frame(sqrtscale(expData))}
    ### PCA
      if (svalue(plotType)=="PCA")
      {
        fit<- prcomp(m)
        colnames(fit$x)<- paste(colnames(fit$x),paste("(",paste(summary(fit)$importance[2,]*100,"%"),")",sep=""),sep="\n")
        if (svalue(nComp)=="Top N Component")
        {
          proceed.pca<- TRUE
          if(svalue(firstCom)==1)
          {
            gmessage("The number of component should be greater than one",icon="info",cont=TRUE,handler=function(h,...){svalue(firstCom)<-2})
            proceed.pca<- FALSE
          }
          if (proceed.pca){pairs(fit$x[,1:svalue(firstCom)],pch=16,col=c(svalue(firstColor),svalue(secondColor),svalue(thirdColor))[as.factor(group)])}
        }

        if (svalue(nComp)=="pair-wise")
        {
          pca.elps<- ellipse(cov(fit$x[,c(as.numeric(svalue(firstCom)),as.numeric(svalue(secondCom)))]),level = 0.95)
          plot(fit$x[,as.numeric(svalue(firstCom))],fit$x[,as.numeric(svalue(secondCom))],col=c(svalue(firstColor),svalue(secondColor),svalue(thirdColor))[as.factor(
            group)],pch=16,xlab=paste("Component",svalue(firstCom)),ylab=paste("Component",svalue(secondCom)),xlim=range(pca.elps[,1])*1.15,ylim=range(pca.elps[,2])*1.15)
          points(pca.elps,type="l")
          abline(h=0)
          abline(v=0)
          if (svalue(withLabel)==TRUE)
          {
            text(fit$x[,as.numeric(svalue(firstCom))]+sum(abs(range(pca.elps[,1])))*svalue(moveLabel),
            fit$x[,as.numeric(svalue(secondCom))],labels=rownames(m),col=c(svalue(firstColor),svalue(secondColor),svalue(thirdColor))[as.factor(group)],cex=0.7)
          }
        }
      }
      ### PLSR
      if (svalue(plotType)=="PLSR")
      {
        l<- c(1:length(groupClass))[as.factor(group)]
        fit<- plsr(l~as.matrix(m))
        if (svalue(nComp)=="Top N Component")
        {
          proceed.pls<- TRUE
          if(svalue(firstCom)==1)
          {
            gmessage("The number of component should be greater than one",icon="info",cont=TRUE,handler=function(h,...){svalue(firstCom)<-2})
            proceed.pls<- FALSE
          }
          if (proceed.pls)
          {
            plot(fit, plottype = "scores", comps = 1:svalue(firstCom),pch=16,col=c(svalue(firstColor),svalue(secondColor),svalue(thirdColor))[as.factor(group)])
          }
        }
        if (svalue(nComp)=="pair-wise")
        {
          pls.elps<- ellipse(cov(fit$score[,c(as.numeric(svalue(firstCom)),as.numeric(svalue(secondCom)))]),level = 0.95)
          plot(fit, plottype = "scores", comps = c(as.numeric(svalue(firstCom)),as.numeric(svalue(secondCom))),col=c(svalue(firstColor),svalue(secondColor),svalue(thirdColor))[as.factor(group)],xlim=range(pls.elps[,1])*1.15,ylim=range(pls.elps[,2])*1.15,pch=16)
          points(pls.elps,type="l")
          abline(h=0)
          abline(v=0)
          if (svalue(withLabel)==TRUE)
          {
            text(fit$score[,as.numeric(svalue(firstCom))]+sum(abs(range(pls.elps[,1])))*svalue(moveLabel),
            fit$score[,as.numeric(svalue(secondCom))],labels=rownames(m),col=c(svalue(firstColor),svalue(secondColor),svalue(thirdColor))[as.factor(group)],cex=0.7)
          }
        }
      }
    }

    # GUI component
    gw.plot<<- gwindow("digeR-Score Plot")

    # split the window into left and right part
    pg <- gpanedgroup(cont = gw.plot, horizontal=TRUE)

    # define the left widgets group
    lg <- ggroup(horizontal = FALSE, cont = pg)

    # notebook on left hand side for seting the parameters for plot and classification
    # first page for plotting
    nb <- gnotebook(cont = lg)

    # setup the layout for ploting parameters
    nb.tab1<- glayout(label=" Score Plot ",cont=nb)
  
    # first line: select plot type
    method<- c("PCA","PLSR")
    nb.tab1[1,1,anchor=c(1,0)]= "Plot Type"
    nb.tab1[1,2] <- plotType<- gdroplist(method,selected= (1:length(method))[(method==selected)],editable=FALSE)
  
    # second line-- separator
    nb.tab1[2,1:4] <- gseparator(cont=nb.tab1)
  
    # fourth line, checkbox to comfirm if want to plot top N component
    # nCom: if want to plot multiple components, topN: top n components
  
    nb.tab1[4,1:3]<- nComp<- gradio(c("Top N Component","pair-wise"),horizontal= TRUE,cont=nb.tab1,handler=function(h,...)
      {
        if(svalue(nComp)=="Top N Component") {enabled(secondCom)<- FALSE;enabled(withLabel)<- FALSE}
        if(svalue(nComp)=="pair-wise") {enabled(secondCom)<- TRUE;enabled(withLabel)<- TRUE}
      })
  
    nb.tab1[5,1]<- glabel("Component 1")
    nb.tab1[5,2]<- firstCom<- gdroplist(1:10,editable=TRUE,cont=nb.tab1)
    nb.tab1[5,3]<- glabel("Component 2")
    nb.tab1[5,4]<- secondCom<- gdroplist(1:10,editable=TRUE,cont=nb.tab1)
    enabled(secondCom)<- FALSE
  
    # select the color for two components
    # define a vector with color names
    colNames<- c("black","red","blue","green","yellow","grey","purple","brown","pink")
    nb.tab1[6,1]<- glabel(groupClass[1])
    nb.tab1[6,2]<- firstColor<- gdroplist(colNames)
    nb.tab1[6,3]<- glabel(groupClass[2])
    nb.tab1[6,4]<- secondColor<- gdroplist(colNames)
    nb.tab1[7,1]<- g3label<- glabel(groupClass[3])
    nb.tab1[7,2]<- thirdColor<- gdroplist(colNames)
    enabled(g3label)<- FALSE
    enabled(thirdColor)<- FALSE
    if (length(groupClass)==3)
    {
      enabled(g3label)<- TRUE
      enabled(thirdColor)<- TRUE
    }
  
    nb.tab1[8,1] <- withLabel<- gcheckbox("With label",cont=nb.tab1,handler=ScorePlot)
    nb.tab1[8,2] <- moveLabel<- gslider(from=-0.1,to=0.1,by=0.01,value=0,handler=ScorePlot)
    enabled(withLabel)<- FALSE
    nb.tab1[8,3]<- glabel("scaling")
    nb.tab1[8,4]<- scaling<- gdroplist(c("No scaling","Center only","Center and Scaling","Paretol Scaling"),editable=FALSE)
    nb.tab1[10,1:4] <- gseparator(cont=nb.tab1)
    nb.tab1[11,3]<- gbutton("Plot    ",cont=nb.tab1,handler=ScorePlot)
  }
  ####################################################################################
  # function for boostrap classification used in classification function             #
  ####################################################################################  
  boots<- function(mdat,label,the.fit,the.predict,nb)
  {
    names(label)<- 1:length(label)
    l1<- label[label==1]
    l2<- label[label==2]
    pred.res<- c()
    label.res<- c()
    for (i in 1:nb)
    {
      select1<- unique(sample(names(l1),length(l1),replace=TRUE))     # bootstrap
      select2<- unique(sample(names(l2),length(l2),replace=TRUE))
      oob<- as.numeric(names(label)[-as.numeric(c(select1,select2))])     # out of bags
      fit<- the.fit(mdat,label,oob)
      res<- the.predict(fit,mdat,oob)
      pred.res<- c(pred.res,res)
      label.res<- c(label.res,label[oob])
    }
    res.boots<- list(pred=pred.res,label=label.res)
    return(res.boots)
  }   
  
  ####################################################################################
  # function for classification                                                      #
  ####################################################################################
  whole.class<- function(selected,inputData)
  {
    expData<- inputData$expData
    cord<- inputData$cord
    groupClass<- inputData$groupClass
    group<- inputData$group
    dataset<- inputData$dataset
    label<- (1:length(groupClass))[as.factor(group)]
    # function for the classification
    classification<- function(dat,label)
    {
        if (svalue(scaleSelect)=="None")  {m<- data.frame(dat)}
        if (svalue(scaleSelect)=="Center only")  {m<- data.frame(scale(dat,center=TRUE,scale=FALSE))}
        if (svalue(scaleSelect)=="Center and Scaling")  {m<- data.frame(scale(dat,center=TRUE,scale=TRUE))}
        if (svalue(scaleSelect)=="Paretol Scaling")  {m<- data.frame(sqrtscale(dat))}
    

        if(svalue(methodSelect)=="LDA")
        {
          theta.fit<- function(x,y,subt)
          {
            m.train<- data.frame(x,class=y)
            lda(class~.,data= m.train,method=svalue(lda.method),subset=-subt)
          }
          theta.predict<- function(fit,x,subt)
          {
            m.train2<- data.frame(x,class=1)
            predict(fit,newdata=m.train2)$post[subt,2]
          }
        }
        if(svalue(methodSelect)=="PCR")
        {
          theta.fit<- function(x,y,subt)
          {
            m.train<- data.frame(x,class=y)
            pcr(class~.,data=m.train,subset=-subt,ncomp=svalue(pcr.nComp))
          }
      		theta.predict<- function(fit,x,subt)
          {
            m.train2<- data.frame(x,class=1)
            predict(fit,newdata =m.train2,ncomp=svalue(pcr.nComp))[subt]
          }
        }
        if(svalue(methodSelect)=="PLSR")
        {
          theta.fit<- function(x,y,subt)
          {
            m.train<- data.frame(x,class=y)
            plsr(class~.,data=m.train,subset=-subt,ncomp=svalue(pls.nComp))
          }
      		theta.predict<- function(fit,x,subt)
          {
            m.train2<- data.frame(x,class=1)
            predict(fit,newdata =m.train2,ncomp=svalue(pls.nComp))[subt]
          }
        }
        if(svalue(methodSelect)=="Logistic_reg.")
        {
          theta.fit<- function(x,y,subt)
          {
            m.train<- data.frame(x,class=y-1)
            glm(class~.,data=m.train,subset=-subt,family=binomial)
          }
      		theta.predict<- function(fit,x,subt)
          {
            m.train2<- data.frame(x,class=1)
            predict(fit,newdata= m.train2,type="response")[subt]
          }
        }
        if(svalue(methodSelect)=="SVM")
        {
          theta.fit<- function(x,y,subt)
          {
            m.train<- data.frame(x,class=y)
            fit<- svm(class~.,data= m.train,subset=-subt)
          }
      		theta.predict<- function(fit,x,subt)
          {
            m.train2<- data.frame(x,class=1)
            as.numeric(predict(fit,newdata=x)[subt])
          }
        }
      ### LOO and N-fold cv
        if (svalue(classMethod)== "leave-one-out cv"|svalue(classMethod)== "N-fold cv")
        {
          if (svalue(classMethod)== "leave-one-out cv")
          {
            ngroup<- length(label)
          }
          if (svalue(classMethod)== "N-fold cv")
          {
            ngroup<- svalue(nFold)
          }
          print(ngroup)
          results<- crossvalid(m,label,theta.fit=theta.fit,theta.predict=theta.predict,ngroup=ngroup)
          return(results)
        }
     ### bootstrap  
        else if (svalue(classMethod)== "bootstrap")
        {
          res<- boots(m,label,theta.fit,theta.predict,n = svalue(nboot))  
          return(res)
        }
    }
      
    # GUI component
    gw.class<<- gwindow("digeR-Classification")
    
    # split the window into left and right part
    pg <- gpanedgroup(cont = gw.class, horizontal=TRUE)
    
    # define the left widgets group
    lg <- ggroup(horizontal = FALSE, cont = pg)
    
    # notebook on left hand side for seting the parameters for plot and classification
    # first page for plotting
    nb <- gnotebook(cont = lg)
    
    # Classification page
    nb.tab2<- glayout(label="Classification",cont=nb)
    nb.tab2[1,1,anchor=c(1,0)]= "Methods"
    method<- c("LDA","PCR","PLSR","Logistic_reg.","SVM")
    nb.tab2[1,2:4] <- methodSelect<- gdroplist(method,selected=(1:length(method))[(method==selected)],editable=FALSE,handler=
      function(h,...)
      {
        if (svalue(methodSelect)=="LDA") {svalue(arg.nb)<- 1}
        if (svalue(methodSelect)=="PCR") {svalue(arg.nb)<- 2}
        if (svalue(methodSelect)=="PLSR") {svalue(arg.nb)<- 3}
      })
    nb.tab2[2,1,anchor=c(1,0)] <- "Scaling"
    nb.tab2[2,2:4]<- scaleSelect<- gdroplist(c("None","Center only","Center and Scaling","Paretol Scaling"),editable=FALSE)
    nb.tab2[3:9,1:4]<- arg.fr<- gframe("Arguments")
    arg.nb<- gnotebook(cont= arg.fr,closebuttons=FALSE)
    # LDA arguments
    lda.arg<- glayout(label="LDA",cont=arg.nb)
    lda.arg[1,1]<- glabel("Methods")
    lda.arg[1,2:3]<-  lda.method<-gdroplist(c("moment","mle","mve","t"))
    # PCR arguments
    pcr.arg<- glayout(label="PCR",cont=arg.nb)
    pcr.arg[1,1]<- glabel("ncomp")
    pcr.arg[1,2:4]<- pcr.nComp<- gspinbutton(from =1,to = length(group)-1,by=1,value= 3,editable=FALSE)
    # PLS arguments
    pls.arg<- glayout(label="PLSR",cont=arg.nb)
    pls.arg[1,1]<- glabel("ncomp")
    pls.arg[1,2:4]<- pls.nComp<- gspinbutton(from =1,to = length(group)-1,by=1,value= 3,editable=FALSE)
    # N-fold cv
    foldcv.arg<- glayout(label="N-fold CV",cont=arg.nb)
    foldcv.arg[1,1]<- nFold<- gspinbutton(from =1,to = length(group),by=1,value= 5,editable=FALSE)
    foldcv.arg[1,2]<- glabel("Fold cross validation")
    
    # bootstrap
    btrap.arg<- glayout(label="bootstrap",cont=arg.nb)
    btrap.arg[1,1]<- glabel("nboot = ")
    btrap.arg[1,2]<- nboot<- gedit(text= "10",coerce.with= as.numeric,handler=function(h,...)
      {
        if (is.na(svalue(nboot)))
        {
          gmessage("Please input integers",icon ="info")
          svalue(nboot)<- ""
        }
      })
    
    svalue(arg.nb)<- 1
    fpath<- c()
    featureList<- c()
    feature.upload<- list()
    # Using the features
    nb.tab2[11,1]<- useFeature<- gcheckbox("Use selected feature",handler= function(h,...)
    {
      if (svalue(useFeature)){enabled(loadFeature.gb)<- TRUE;enabled(feature.gdrop)<- TRUE}
      if (!svalue(useFeature)){enabled(loadFeature.gb)<- FALSE;enabled(feature.gdrop)<- FALSE}
    })
    nb.tab2[11,2]<- loadFeature.gb<- gbutton("load features",handler=function(h,...)
    {
      fpath<<- gfile("Select the RData file containing features",filter = list(".RData" = list(patterns = c("*.RData"))),type="open")
      if(length(fpath)!=0)
      {
        load(fpath)
        feature.upload<<- feature
        featureList<- rep("",length(feature.upload))
        for (j in 1:length(feature.upload))
        {
          featureList[j]<- names(feature.upload)[j]
        }
        feature.gdrop[]<- featureList 
      }
    })
    
    nb.tab2[11,3:4]<- feature.gdrop<- gdroplist(featureList)
    enabled(loadFeature.gb)<- FALSE
    enabled(feature.gdrop)<- FALSE
    
    nb.tab2[12,1]<- classMethod<- gdroplist(c("leave-one-out cv","N-fold cv","bootstrap"),cont=nb.tab2,handler=
      function(h,..)
      {
        if (svalue(classMethod,index=TRUE)==2) {svalue(arg.nb)<- 4}
        if (svalue(classMethod,index=TRUE)==3) {svalue(arg.nb)<- 5}
      })
    
    classResult<- list()
    length(classResult)<- 10
    
    nb.tab2[12,2]<- classbutton<- gbutton("Run Classfication",cont=nb.tab2,handler=
    function(h,...)
    {
      proceed<- TRUE
      if (length(svalue(res.pred,index=TRUE))==0)
      {
        gmessage("Please select a result item to store the results",icon="info")
        proceed<- FALSE
      }
      if (length(svalue(res.pred,index=TRUE))>1)
      {
        gmessage("Please only select one result item",icon="info")
        svalue(res.pred,index=TRUE)<- svalue(res.pred,index=TRUE)[1]
        proceed<- FALSE
      }
    
      if (proceed)
      {
        if (svalue(classMethod)=="leave-one-out cv")
        {
          cv.method<- svalue(classMethod)
        }
        if (svalue(classMethod)=="N-fold cv")
        {
          cv.method<- paste(svalue(nFold),"fold cv",sep=" ")
        }
        if (svalue(classMethod)=="bootstrap")
        {
          cv.method<- paste(svalue(nboot),"times bootstrap",sep=" ")
        }
        
        i<- svalue(res.pred,index=TRUE)
        res.pred[i]<- paste(svalue(methodSelect),cv.method)
        print("Classification...")
        svalue(classbutton)<- "Processing.."
        enabled(nb)<- FALSE
        if (svalue(useFeature)==TRUE)
        {
          classResult[[i]]<<- classification(expData[,feature.upload[[svalue(feature.gdrop,index=TRUE)]]],label)
          names(classResult)[i]<<- svalue(methodSelect)
        }
        else if(svalue(useFeature)==FALSE)
        {
          classResult[[i]]<<- classification(expData,label)
          names(classResult)[i]<<- paste(svalue(methodSelect),cv.method)
        }  
        print("Finished")
        svalue(classbutton)<- "Run Classfication"
        enabled(nb)<- TRUE
      }
    })
    nb.tab2[12,3]<- gbutton("Save",handler=function(h,...)
    {
      selectPath<- gfile("Save as .RData",type="save")
      RDataPath<- paste(selectPath,".RData",sep="")
      save(classResult,file=RDataPath)
    })
    nb.tab2[14,1:4] <- gseparator(cont=nb.tab2)
    nb.tab2[1,5]<- glabel("Prediction results")
    pred<-  paste("results_",1:10)
    nb.tab2[2:10,5:7]<- res.pred<- gtable(pred,cont=nb.tab2,multiple=TRUE)
    nb.tab2[15,1]<- glabel("Legend")
    nb.tab2[15,2]<- gleg<- gdroplist(c("bottomright","bottomleft","topright","topleft","left","right","top","bottom"))
    nb.tab2[15,3]<- gbutton("ROC curve",cont=nb.tab2,handler=function(h,...)
    {
      proceed<- TRUE
      if (length(svalue(res.pred,index=TRUE))==0)
      {
        gmessage("Please run classification first",icon="info")
        proceed<- FALSE
      }
      if (proceed)
      {
        for (j in svalue(res.pred,index=TRUE))
        {
          if(is.null(classResult[[j]]))
          {
            gmessage("Can not find the selected result",icon="info")
            proceed<- FALSE
            break
          }
        }
      }
    
      if (proceed)
      {
        select.roc<- svalue(res.pred,index=TRUE)
        prd<- list()
        lb<- list()
        a<- 1
        for (i in select.roc)
        {
          prd[[a]]<- classResult[[i]]$pred
          lb[[a]]<- classResult[[i]]$label
          a<- a+1
        }
        roc(prd,lb,svalue(res.pred),svalue(gleg))
      }
    })
  }


  ####################################################################################
  # Correlation analysis                                                             #
  ####################################################################################
  gwtCor<- function(inputData,imagePath)
  {
    expData<- inputData$expData
    cord<- inputData$cord
    groupClass<- inputData$groupClass
    group<- inputData$group
    dataset<- inputData$dataset
    # the correlation plot
    plotCor <- function(...)
    {
      i<- svalue(gr,index=TRUE)
      if (i ==1){cormat<- cor(expData, method = c("pearson", "kendall", "spearman")[svalue(corrtype.gdrop,index=TRUE)])}
      if (i >1){cormat<- cor(dataset[[i-1]], method = c("pearson", "kendall", "spearman")[svalue(corrtype.gdrop,index=TRUE)])}
      if (svalue(useFeature)){j<- as.numeric(svalue(spot.gdrop))}
      if (!svalue(useFeature)){j<- as.numeric(svalue(glist))}
      par(mar=c(3.5,1,3.5,1))
      plot(cord[,1],-cord[,2],cex= 3*abs(cormat[j,])*(abs(cormat[j,])>svalue(gsl)),
        col=(cormat[j,]>0)+1,xlab="PI",ylab="MASS",xaxt="n",yaxt="n")
      points(cord[j,1],-cord[j,2],pch=3,cex=3,col="blue")
  
      shown<- (1:nrow(cord))[abs(cormat[j,])>svalue(gsl)]
      if (svalue(gch))
      {
        text(cord[shown,1],-cord[shown,2],labels= shown,cex=0.6,col= 2-(cormat[j,shown]>0))
      }
      if (svalue(g.showNo))
      {
       svalue(gn)<- shown
      }
   }
   
    # GUI component
    gw.correlation <<- gwindow("digeR-Spots Correlation")
    BigGroup<- ggroup(cont= gw.correlation)
    leftpane <- ggroup(horizontal=FALSE, cont=BigGroup)
    
    # The first frame for choose the dataset
    tmp <- gframe("Dataset", container=leftpane, expand= TRUE)
    gr <- gradio(c("All",groupClass), horizontal= FALSE, cont=tmp, handler= plotCor)
    
    # The second frame for choose the spot list
    tmp <- gframe("Spot List", container=leftpane, horizontal=FALSE, expand=TRUE)
    glist <- gspinbutton(from=1,to= ncol(expData),by=1,cont=tmp,editable=TRUE)
    featureList<- c()
    feature.upload<- list()
    useFeature<- gcheckbox("Selected feature",cont=tmp,handler= function(h,...)
    {
      if (svalue(useFeature)){enabled(loadFeature.gb)<- TRUE;enabled(feature.gdrop)<- TRUE;enabled(spot.gdrop)<-TRUE;enabled(glist)<- FALSE}
      if (!svalue(useFeature)){enabled(loadFeature.gb)<- FALSE;enabled(feature.gdrop)<- FALSE;enabled(spot.gdrop)<-FALSE;enabled(glist)<- TRUE}
    })
    loadFeature.gb<- gbutton("load features",cont=tmp,handler=function(h,...)
    {
      fpath<- gfile("Select the RData file containing features",filter = list(".RData" = list(patterns = c("*.RData"))),type="open")
      if(fpath!="")
      {
        load(fpath)
        feature.upload<<- feature
        if (length(feature.upload)==0)
        {
          gmessage("Can not find selected features",icon="info")
        }
        featureList<- rep("",length(feature.upload))
        for (j in 1:length(feature.upload))
        {
          featureList[j]<- names(feature.upload)[j]
        }
        feature.gdrop[]<- featureList 
      }
    })
    
    feature.gdrop<- gdroplist(featureList,cont=tmp,handler=function(h,...)
    {
      spot.gdrop[]<- feature.upload[[svalue(feature.gdrop,index=TRUE)]]
    })
    spot.gdrop<- gdroplist(" ",cont=tmp)
    enabled(loadFeature.gb)<- FALSE
    enabled(feature.gdrop)<- FALSE
    enabled(spot.gdrop)<- FALSE
    
    corrtype.gdrop<- gdroplist(c("pearson", "kendall", "spearman"),cont=tmp,selected = 1) 
    gb<- gbutton("Show the correlation",cont=tmp,handler=plotCor)
    gn<- gtext("",cont=tmp,editable=TRUE,horizontal=TRUE)
    
    glabel<- glabel("Correlation Coefficient",con= leftpane)
    gsl <- gslider(from = 0.01, to= 0.99, by=0.01, value = 0, cont = leftpane, handler = plotCor)
    gch<- gcheckbox("Show Spots ID", container= leftpane,handler= plotCor)
    g.showNo<- gcheckbox("Show Number",cont= leftpane,handler= plotCor)
    add(BigGroup, ggraphics())
    if (imagePath!="")
    {
      gimage(imagePath,cont=BigGroup)
    }
  }
  ####################################################################################
  # feature selection                                                                #   
  ####################################################################################
  fs<- function(selected,inputData)
  {
    expData<- inputData$expData
    cord<- inputData$cord
    groupClass<- inputData$groupClass
    group<- inputData$group
    dataset<- inputData$dataset
    l<- (1:length(groupClass))[as.factor(group)]

    # coefficient based method
    coef.ft<- function(dat,l,n,method)
    {
      # scaling
      if (svalue(scaleSelect)=="None"){m<- data.frame(dat)}
      if (svalue(scaleSelect)=="Center only")  {m<- data.frame(scale(dat,center=TRUE,scale=FALSE))}
      if (svalue(scaleSelect)=="Center and Scaling")  {m<- data.frame(scale(dat,center=TRUE,scale=TRUE))}
      if (svalue(scaleSelect)=="Paretol Scaling")  {m<- data.frame(sqrtscale(dat))}
            
      m.train<- data.frame(m,class=l)
      if (method=="lda")
      {
        lda.fit<- lda(class~.,data=m.train,method=svalue(lda.method))
        fstore<- order(-abs(lda.fit$scaling))[1:n]
        value<- abs(lda.fit$scaling)[order(-abs(lda.fit$scaling))[1:n]]
      }
      if (method=="pls")
      {
        pls.fit<- plsr(class~.,data=m.train,ncomp=svalue(pls.nComp))
        fstore<- order(-abs(pls.fit$coeff[,,ncol(pls.fit$scores)]))[1:n]
        value<-  abs(pls.fit$coeff[,,ncol(pls.fit$scores)])[order(-abs(pls.fit$coeff[,,ncol(pls.fit$scores)]))[1:n]]
      }
      names(value)<- fstore
      bplot<- barplot(value,col="green",ylim=c(0,max(value)*1.03))
      text(bplot,value*1.02,labels=fstore)
    	return(fstore)
    }
    
    ###############################################
    # RandomForest
    randf<- function(h,...)
    {
      # scaling
      if (svalue(scaleSelect)=="None"){m<- data.frame(expData)}
      if (svalue(scaleSelect)=="Center only")  {m<- data.frame(scale(expData,center=TRUE,scale=FALSE))}
      if (svalue(scaleSelect)=="Center and Scaling")  {m<- data.frame(scale(expData,center=TRUE,scale=TRUE))}
      if (svalue(scaleSelect)=="Paretol Scaling")  {m<- data.frame(sqrtscale(expData))}
      l<- (1:length(groupClass))[as.factor(group)]
      
      m.train<- data.frame(m,class=as.factor(l))
      print(svalue(rf.ntree))
      if (svalue(mtry.change)==FALSE)
      {
        m.rf<- randomForest(class ~.,data = m.train,ntree =as.numeric(svalue(rf.ntree)))
      }
      if (svalue(mtry.change)==TRUE)
      {
        m.rf<- randomForest(class ~.,data = m.train,ntree =as.numeric(svalue(rf.ntree)),mtry=as.numeric(svalue(rf.mtry)))
      }
      print(m.rf)
      par(ask=TRUE)
      plot(m.rf,main="Error rate")
      imp<- importance(m.rf)[,1]
      names(imp)<- substring(names(imp),2,6)
      bplot.rf<- barplot(sort(imp,decreasing=TRUE)[1:svalue(rf.nimp)], main = "Variable importance",col="green",ylim=c(0,max(imp)*1.03))
      text(bplot.rf,sort(imp,decreasing=TRUE)[1:svalue(rf.nimp)]*1.02,labels=names(sort(imp,decreasing=TRUE)[1:svalue(rf.nimp)]))
      feature.rf<- order(-abs(imp))[1:svalue(rf.nimp)]
      par(ask=FALSE)
      return(feature.rf)
    }
    # Adaboost                                
    ada<- function(h,...)
    {
      # scaling
      if (svalue(scaleSelect)=="None"){m<- data.frame(expData)}
      if (svalue(scaleSelect)=="Center only")  {m<- data.frame(scale(expData,center=TRUE,scale=FALSE))}
      if (svalue(scaleSelect)=="Center and Scaling")  {m<- data.frame(scale(expData,center=TRUE,scale=TRUE))}
      if (svalue(scaleSelect)=="Paretol Scaling")  {m<- data.frame(sqrtscale(expData))}
      l<- (1:length(groupClass))[as.factor(group)]
      colnames(m)<- substring(colnames(m),2,6)

      l<- as.factor(l)
      m.train<- data.frame(m,class= l)
      fit<- adaboost.M1(class~., data=m.train, boos=TRUE, mfinal=svalue(ada.mf),maxdepth=5)
      imp<- fit$importance
      bplot.ada<- barplot(sort(imp,decreasing=TRUE)[1:svalue(ada.nimp)], main = "Variable importance",col="green",ylim=c(0,max(imp)*1.03))
      text(bplot.ada,sort(imp,decreasing=TRUE)[1:svalue(ada.nimp)]*1.03,labels=names(sort(imp,decreasing=TRUE)[1:svalue(ada.nimp)]))
      feature.ada<- order(-abs(imp))[1:svalue(ada.nimp)]
      return(feature.ada)
    }
    
    
    ###################################################################
    # GUI
    # Feature selection including prefilter(ttest,lda-coef,pls-coef),method(stept,randomForest,Adaboost)
    # Main window
    gw.fs<<- gwindow("digeR-Feature Selection")
    
    # split the window into left and right part
    pg <- gpanedgroup(cont = gw.fs, horizontal=TRUE)
    
    # define the left widgets group
    lg <- ggroup(horizontal = FALSE, cont = pg)
    
    # notebook on left hand side for seting the parameters for plot and classification
    # first page for plotting
    nb <- gnotebook(cont = lg)
    ###########################################################
    
    nb.tab2<- glayout(label="Feature Selection",cont=nb)
    # select method
    nb.tab2[1,1,anchor=c(1,0)]<- "Methods"
    method<- c("LDA","PLSR","RandomForest","Adaboost")
    nb.tab2[1,2:4] <- methodSelect<- gdroplist(method,selected=(1:length(method))[(method==selected)],editable=FALSE,handler=
      function(h,...)
      {
        if (svalue(methodSelect)=="LDA") {svalue(arg.nb)<- 1}
        if (svalue(methodSelect)=="PLSR") {svalue(arg.nb)<- 2}
        if (svalue(methodSelect)=="RandomForest") {svalue(arg.nb)<- 3}
        if (svalue(methodSelect)=="Adaboost") {svalue(arg.nb)<- 4}
      })
    
    nb.tab2[2,1,anchor=c(1,0)] <- "Scaling"
    nb.tab2[2,2:4]<- scaleSelect<- gdroplist(c("None","Center only","Center and Scaling","Paretol Scaling"),editable=FALSE)
    nb.tab2[3,1:4]<- gseparator()
    
    # Arguments
    nb.tab2[5:12,1:4]<- arg.fr<- gframe("Arguments")
    arg.nb<- gnotebook(cont= arg.fr,closebuttons=FALSE)
    
    # LDA arguments
    lda.arg<- glayout(label="LDA",cont=arg.nb)
    lda.arg[1,1]<- glabel("Methods")
    lda.arg[1,2:3]<-  lda.method<-gdroplist(c("moment","mle","mve","t"))
    lda.arg[2,1]<- glabel("Top")
    lda.arg[2,2:4]<- lda.n<- gdroplist(c("10","20","30"),selected=2,editable=TRUE)
    
    # PLSR arguments
    pls.arg<- glayout(label="PLSR",cont=arg.nb)
    pls.arg[1,1]<- glabel("ncomp")
    pls.arg[1,2:4]<- pls.nComp<- gdroplist(c(1:(nrow(expData)-1)))
    pls.arg[2,1]<- glabel("Top")
    pls.arg[2,2:4]<- pls.n<- gdroplist(c("10","20","30"),selected=2,editable=TRUE)
    
    # RandomForest arguments
    rf.arg<- glayout(label="RandomForest",cont=arg.nb)
    rf.arg[1,1]<- glabel("ntree")
    rf.arg[1,2]<- rf.ntree<- gdroplist(c("100","200","300"),editable=TRUE)
    rf.arg[2,1]<- glabel("Top")
    rf.arg[2,2]<- rf.nimp<- gdroplist(c("10","20","30"),selected=2,editable=TRUE)
    rf.arg[3,1]<- mtry.change<- gcheckbox("mtry",handler= function(h,..){
      if (svalue(mtry.change)==TRUE){enabled(rf.mtry)<- TRUE}
      if (svalue(mtry.change)==FALSE){enabled(rf.mtry)<- FALSE}
      })
    rf.arg[3,2]<- rf.mtry<- gdroplist(c("5","10","15"),selected=1,editable=TRUE) 
    enabled(rf.mtry)<- FALSE
    
    # Adaboost arguments
    ada.arg<- glayout(label="Adaboost",cont=arg.nb)
    ada.arg[1,1]<- glabel("mfinal")
    ada.arg[1,2]<- ada.mf<- gdroplist(c("100","200","300"),editable=TRUE)
    ada.arg[2,1]<- glabel("Top")
    ada.arg[2,2]<- ada.nimp<- gdroplist(c("10","20","30"),selected=2,editable=TRUE)
    svalue(arg.nb)<- svalue(methodSelect,index=TRUE)
    
    feature<- list()
    nb.tab2[13,1:5] <- gseparator(cont=nb.tab2)
    nb.tab2[14,3]<- fs.button<- gbutton("Run feature selection",cont=nb.tab2,handler= function(h,...)
    {
      proceed<- TRUE
      if (length(svalue(res.pred))==0)
      {
        gmessage("Please select a result item to store the features",icon="info")
        proceed<- FALSE
      }
      if (proceed)
      {
        i<- svalue(res.pred,index=TRUE)
        res.pred[i]<- paste(svalue(methodSelect))
        svalue(fs.button)<- "Processing.."
        enabled(nb)<- FALSE
        print("Feature Selection...")
        if (svalue(methodSelect)=="LDA") {feature[[i]]<<- coef.ft(expData,l,svalue(lda.n),"lda");names(feature)[i]<<- "LDA"}
        if (svalue(methodSelect)=="PLSR") {feature[[i]]<<- coef.ft(expData,l,svalue(pls.n),"pls");names(feature)[i]<<- "PLSR"}
        if (svalue(methodSelect)=="RandomForest") {feature[[i]]<<- randf(h,...);names(feature)[i]<<- "RandomForest"}
        if (svalue(methodSelect)=="Adaboost") {feature[[i]]<<- ada(h,...);names(feature)[i]<<- "Adaboost"}
        print("Finished")
        svalue(fs.button)<- "Run feature selection"
        enabled(nb)<- TRUE
      }
    })
      
    nb.tab2[1,5]<- glabel("Selected features")
    pred<-  paste("results_",1:10)
    nb.tab2[2:12,5]<- res.pred<- gtable(pred,cont=nb.tab2,multiple=TRUE)
    nb.tab2[14,5]<- gbutton("Save features",handler=function(h,...)
    {
      proceed.output<- TRUE
      if (length(feature)==0)
      {
        gmessage("Can not find feature results",icon="info")
        proceed.output<- FALSE
      }
      if (proceed.output)
      {
        featurePath<- gfile("Save as .RData",type="save")
        RDataPath<- paste(featurePath,".RData",sep="")
        save(feature,file=RDataPath)
      }
    })
    
    # set up the intial situation
    if (svalue(methodSelect)=="LDA") {svalue(arg.nb)<- 1}
    if (svalue(methodSelect)=="PLSR") {svalue(arg.nb)<- 2}
    if (svalue(methodSelect)=="RandomForest") {svalue(arg.nb)<- 3}
    if (svalue(methodSelect)=="Adaboost") {svalue(arg.nb)<- 4}
  }
  ####################################################################################
  # function to generate ROC curves                                                  #
  ####################################################################################
  
  roc<- function(prd,lb,a,pos)           # prd=prediction, lb= label, a= method used for prediction, pos= legend postion
  {
    n<- length(prd)
   	AUC<- rep(0,n)
   	perf<- list()
    for (i in 1:n)
    {
  	 perf[[i]]<- performance(prediction(prd[[i]],lb[[i]]),"sens","fpr")
  	 AUC[i]<-  trapz(attributes(perf[[i]])$x.v[[1]],attributes(perf[[i]])$y.v[[1]])
    }
    AUC<- round(AUC,3)
    # do the plot
  	plot(c(0, 1), c(0, 1),type="n",xlab = "False Postive Rate", sub = "1-Specificity", ylab = "Sensitivity")
    for (i in 1:n)
    {
      plot(perf[[i]],lty=i,col=i,add=TRUE)
    }
  	abline(0,1,col="lightgrey",lty=3)
  	legend(pos, paste(a," (AUC =",AUC,")"),lty=1:n,col=1:n,bg = "white",cex=0.75)
  }

  ####################################################################################
  # function for doing cross validation, bootstrap                                   #
  ####################################################################################
  # x= data matrix, y=label,theta.fit= fit method,theta.predict= predict method, ngroup= N-fold crossvalidation
  crossvalid<- function(x,y,theta.fit,theta.predict,...,ngroup)  
  {
    n <- length(y)
    if (ngroup < 2|ngroup>n)
    {
      stop("2=< N <= number of observations")
    }
  
    if (ngroup == n)
    {
      groups <- 1:n
      leave.out <- 1
    }
    if (ngroup < n)
    {
      leave.out <- trunc(n/ngroup)
      o <- sample(1:n)
      groups <- vector("list", ngroup)
      for (j in 1:(ngroup - 1))
      {
        jj <- (1 + (j - 1) * leave.out)
        groups[[j]] <- (o[jj:(jj + leave.out - 1)])
      }
      groups[[ngroup]] <- o[(1 + (ngroup - 1) * leave.out):n]
    }
    print(ngroup)
    print(groups)
    res<- list()
    label<- list()
    for (i in 1:length(groups))
    {
      fit<- theta.fit(x,y,groups[[i]])
      res[[i]]<- theta.predict(fit,x,groups[[i]])
      label[[i]]<- y[groups[[i]]]
    }
   list(pred = unlist(res),label = unlist(label))
  }
  
  ####################################################################################
  # GUI for calculate the power for spots and experiment                             #
  ####################################################################################
  power.call<- function(inputData)
  {
    expData<- inputData$expData
    cord<- inputData$cord
    groupClass<- inputData$groupClass
    group<- inputData$group
    dataset<- inputData$dataset
    
    power.size<- function(h,...)
    {
      m<- expData
      l<- (1:length(groupClass))[as.factor(group)]
      ss<- as.numeric(svalue(sampleSize))
      sl<- as.numeric(svalue(sigLev))
      pw<- as.numeric(svalue(pow))
      na.no<- c(1:3)[c(is.na(ss),is.na(sl),is.na(pw))]
      if(is.na(ss)){ss<- NULL}
      if(is.na(sl)){sl<- NULL}
      if(is.na(pw)){pw<- NULL}
      proceed<-TRUE
      if (length(na.no)!=1)
      {
        gmessage("Please indicate the one to be calculated and leave it blank",icon="error",cont=TRUE,
          handler=function(h,...){svalue(sampleSize)<-"";svalue(sigLev)<-"0.05";svalue(pow)<-"0.8"})
        proceed<- FALSE
      }
      if(proceed)
      {
        if (svalue(calType)=="Single spots")
        {
          i<- as.numeric(svalue(spotNo))
          if (length(groupClass)==2)
          {
            fit<- power.t.test(n= ss,delta = mean(m[l==2,i])-mean(m[l==1,i]), sd = sd(m[,i]),
            sig.level = sl, power= pw,type = "two.sample",alternative = "two.sided")
      
            svalue(sampleSize)<- round(fit$n,3)
            svalue(sigLev)<-  round(fit$sig,3)
            svalue(pow)<- round(fit$power,3)
          }
          if (length(groupClass)>2)
          {
            anova.fit<- anova(aov(m[,i]~ l))
            fit<- power.anova.test(groups = length(groupClass), n =ss, between.var =
            anova.fit$Mean[1],within.var = anova.fit$Mean[2], sig.level =sl, power = pw)
      
            svalue(sampleSize)<- round(fit$n,3)
            svalue(sigLev)<-  round(fit$sig,3)
            svalue(pow)<- round(fit$power,3)
          }
        }
        if (svalue(calType)=="Gel")
        {
          m.train<- data.frame(m,class=l)
          ldafit<- lda(class~.,data=m.train)
          pred<- predict(ldafit,newdata=m.train,dimen=1)$x
          if (length(groupClass)==2)
          {
            fit<- power.t.test(n= ss,delta = mean(pred[l==2])-mean(pred[l==1]), sd = sd(pred),
            sig.level = sl, power= pw,type = "two.sample",alternative = "two.sided")
      
            svalue(sampleSize)<- round(fit$n,3)
            svalue(sigLev)<-  round(fit$sig,3)
            svalue(pow)<- round(fit$power,3)
          }
          if (length(groupClass)>2)
          {
            anova.fit<- anova(aov(as.numeric(pred)~ l))
            fit<- power.anova.test(groups = length(groupClass), n = ss, between.var =
            anova.fit$Mean[1],within.var = anova.fit$Mean[2], sig.level = sl, power = pw)
      
            svalue(sampleSize)<- round(fit$n,3)
            svalue(sigLev)<-  round(fit$sig,3)
            svalue(pow)<- round(fit$power,3)
          }
        }
      }
    }

    # GUI  component
    gw.power<<- gwindow("digeR-Sample size")
    pg <- gpanedgroup(cont = gw.power, horizontal=TRUE)
    lg <- ggroup(horizontal = FALSE, cont = pg)
    nb <- gnotebook(cont = lg)
    nb.tab2<- glayout(label="Power Analysis",cont=nb)
    
    nb.tab2[1,1:2]<- calType<- gradio(c("Single spots","Gel"),handler=function(h,...)
      {
        if (svalue(calType)=="Gel"){enabled(spotNo)<- FALSE}
        if (svalue(calType)=="Single spots"){enabled(spotNo)<- TRUE}
      }) 
    nb.tab2[1,4]<- glabel("Spot number")
    nb.tab2[1,5]<- spotNo<- gspinbutton(from=1,to=ncol(expData),by=1,handler=function(h,...){svalue(sampleSize)<- ""})
    
    nb.tab2[2:6,1:4]<- arg.fr<- gframe("Arguments")
    arg.set<- glayout(label="parameters",cont=arg.fr)
    arg.set[1,1]<- glabel("Effect size")
    arg.set[1,2]<- effectSize<- gedit()
    enabled(effectSize)<- FALSE
    arg.set[1,3]<- glabel("Significant level")
    arg.set[1,4]<- sigLev<- gedit("0.05")
    arg.set[2,1]<- glabel("Power")
    arg.set[2,2]<- pow<- gedit("0.8")
    arg.set[2,3]<- glabel("Sample size Per group")
    arg.set[2,4]<- sampleSize<- gedit()
    
    nb.tab2[7,1:4]<- gseparator()
    nb.tab2[9,1:2]<- glabel("Effect size is estimated from data")
    nb.tab2[9,3]<- gbutton("Calculate",handler= power.size)
  }

  ####################################################################################
  # function for square root scaling                                                 #
  ####################################################################################
  sqrtscale<- function(x)
  {
    x.s<- x
    for (i in 1:ncol(x))
    {
    x.s[,i]<- x[,i]-mean(x[,i])
    x.s[,i]<- x[,i]/sqrt(sd(x[,i]))
    }
    return(x.s)
  }
  
  ####################################################################################
  # digeR GUI component                                                              # 
  ####################################################################################
  
  digeR.menu<- gwindow("digeR",cont=TRUE)
  addHandlerDestroy(digeR.menu,handler=function(h,...)
    {
      dispose(gw.plot)
      dispose(gw.class)
      dispose(gw.correlation)
      dispose(gw.fs)
      dispose(gw.power)    
    })

  gw.plot<- gwindow(visible=FALSE)
  gw.class<- gwindow(visible=FALSE)
  gw.correlation<- gwindow(visible=FALSE)
  gw.fs<- gwindow(visible=FALSE)
  gw.power<- gwindow(visible=FALSE)
  
  
  gimage(paste(file.path(.path.package(package="digeR")[1]),"/data/cover.jpg",sep=""),cont= digeR.menu)
  mblst = list(
      File=list(
        Open = list(handler=gf.open,icon="open"),
        Upload_gel_image = list(handler= upload,icon="newplot"),
        Quit = list(handler=gf.quit,icon="close")
        ),                             
      Correlation = list(
        Gel = list(handler=function(h,...)
          {if (length(input)==0){gmessage("Please upload the data first",icon="info")}
           else{gwtCor(input,imagepath)}})),
      Score_Plot = list(                     
        PCA = list(handler= function(h,..)
          {if (length(input)==0){gmessage("Please upload the data first",icon="info")}
           else{whole.plot(selected="PCA",input)}}),
        PLSR = list(handler=function(h,..)
          {if (length(input)==0){gmessage("Please upload the data first",icon="info")}
           else{whole.plot(selected="PLSR",input)}})),
      Classification = list(
        LDA = list(handler= function(h,..)
          {if (length(input)==0){gmessage("Please upload the data first",icon="info")}
           else{whole.class(selected="LDA",input)}}),
        PCR = list(handler= function(h,..)
          {if (length(input)==0){gmessage("Please upload the data first",icon="info")}
           else{whole.class(selected="PCR",input)}}),
        PLSR = list(handler= function(h,..)
          {if (length(input)==0){gmessage("Please upload the data first",icon="info")}
           else{whole.class(selected="PLSR",input)}}),
        Logistic_reg. = list(handler= function(h,..)
          {if (length(input)==0){gmessage("Please upload the data first",icon="info")}
           else{whole.class(selected="Logistic_reg.",input)}}),
        SVM = list(handler= function(h,..)
          {if (length(input)==0){gmessage("Please upload the data first",icon="info")}
           else{whole.class(selected="SVM",input)}})),
      Feature_Selection = list(                  
        LDA = list(handler=  function(h,..)
          {if (length(input)==0){gmessage("Please upload the data first",icon="info")}      
           else{fs(selected="LDA",input)}}),
        PLSR = list(handler=  function(h,..)
          {if (length(input)==0){gmessage("Please upload the data first",icon="info")}     
           else{fs(selected="PLSR",input)}}),
        RandomForest = list(handler= function(h,..)
          {if (length(input)==0){gmessage("Please upload the data first",icon="info")}     
           else{fs(selected="RandomForest",input)}}),
        Adaboost = list(handler=function(h,..)
          {if (length(input)==0){gmessage("Please upload the data first",icon="info")}        
           else{fs(selected="Adaboost",input)}})),
      Power = list(Sample_Size = list(handler=function(h,..)
          {if (length(input)==0){gmessage("Please upload the data first",icon="info")}  
          else{power.call(input)}}))
  )
  menu<- gmenu(mblst, cont = digeR.menu)
}

