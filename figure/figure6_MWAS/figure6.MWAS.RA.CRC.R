	#### pub data
		dat <- read.table("abundance_profile_strain_pub.jgi.tsv",
			sep="\t",head=T,row.names=1)

		dat[1:4,1:4]
		dim(dat)

	### pub data phenotype
		phe<-read.table("abundance_profile_strain_pub.add_phenotype.jgi.tsv",
				sep="\t",head=T,as.is=T )
		head(phe)
		dim(phe)
		library(data.table)
		dat<-dcast(phe[,1:3],lineages_strain_new~sample_id)
		dat[1:4,1:4]
		rownames(dat)=dat[,1]
		dat=dat[,-1]
		dim(dat)

		phe1=phe[,c(2,4:ncol(phe))]
		phe1=phe1[!duplicated(phe1$sample_id),]
		dim(phe1)
		head(phe)
		tem=aggregate(phe1$Diagnosis,list(study=phe1$Study,diag=phe1$Diagnosis),length)
		tem


	
		ph_sa=phe1[phe1$Study=="ZellerG_2014" & phe1$Diagnosis !="small adenoma" & phe1$Country!="Germany",]
		table(ph_sa[,c("Diagnosis","Country")])
		ph_sa
		sp_sa=t(dat[,phe1$Study=="ZellerG_2014"& phe1$Diagnosis !="small adenoma"& phe1$Country!="Germany"])
		dim(sp_sa)


	### abundance are asin sqrt transform
	occu<-function(sp1){
		asin(sqrt(sp1[,apply(sp1,2,function(x){ mean(x>0)>0.1})]/100))
	}
	sp_sa2=occu(sp_sa)
	dim(sp_sa2)

	## rename the data
	ta=unlist(lapply(strsplit(as.character(colnames(sp_sa2)),"__",fixed=T),function(x){ 
					id=unlist(lapply(strsplit(x[length(x)],"_",fixed=T),function(j){ j[length(j)]}));
					x=x[-1];
					x1=x[!grepl("unclassified",x)];
					ifelse(grepl("mgs_", x1[length(x1)]), gsub(".t$","",x1[length(x1)]),
						paste( gsub(".t$","",x1[length(x1)]),"_mgs_",id,sep=""))
					#paste(x1[length(x1)]	, 
					#unlist(lapply(strsplit(x[length(x)],"_"),function(j){j[length(j)]})),sep="_"  )
				 }))

	colnames(sp_sa2)=ta	
	head(ta)
	## deal with R name system
	write.csv(sp_sa2,"tmp.csv")
	sp_sa2=read.csv("tmp.csv",row.names=1,head=T)	
	head(colnames(sp_sa2))

	### cross validation prediction




		### load package
		library(caret) # for model-building
		library(DMwR) # for smote implementation
		library(pROC) # for AUC calculations
		library(gbm)
		library(ppcor)
		library(randomForest)
		library(e1071)
		set.seed(0)
		library(doMC)
		registerDoMC(cores = 12)

		##### split data 
		 disease=factor(factor(ph_sa$Diagnosis,labels=c("D","C")),levels=c("C","D"))
		table( disease)

		### 10 fold 5 repeat
		set.seed(0)
		index=createMultiFolds(disease,k=10,times=5)
		


		## multiple core run for each cv
		res=foreach(ind= 1:1:length(index), .combine=rbind ) %dopar% {
			trainindex=index[[ind]]  #trainset
			testindex=c(1:length(disease))[-trainindex] #testset
			table(disease[trainindex])
			table(disease[testindex])

			### glm for biomaker , adjust age, gender, bmi.
			### use train data only
			pvalue=c()
			beta=c() 
			for( i in 1:ncol(sp_sa2)){
				tmp=data.frame(Gender=factor(ph_sa[,c("Sex")]),Age=ph_sa$Age,
						disease=factor(factor(ph_sa$Diagnosis,labels=c("D","C")),levels=c("C","D")) )
				tmp$y=as.numeric(t(sp_sa2[,i]))
					#wil_p=c(wil_p,wilcox.test(tmp$y ~ tmp$disease)$p.value)
					#pcor_p=c(pcor_p,pcor(data.matrix(tmp),method="spearman")$p.value["y","disease"])
					#pcor_r=c(pcor_r,pcor(data.matrix(tmp),method="spearman")$estimate["y","disease"])
				tmp=tmp[trainindex,]  # train data only use for biomarker
				mod=glm(y ~ disease + Gender  + Age, data=tmp) 
				pvalue=c(pvalue,coef(summary(mod))["diseaseD","Pr(>|t|)"])
				beta=c(beta,coef(summary(mod))["diseaseD","Estimate"])
			}
			### fdr 
			qvalue=p.adjust(pvalue,method="fdr")
			res=data.frame(beta,pvalue,qvalue)
			rownames(res)=colnames(sp_sa2)

			### keep q<0.01 biomarker for next prediction
			px=sp_sa2[,pmatch( intersect(rownames(res),colnames(sp_sa2)),colnames(sp_sa2) )]

			tmp=data.frame( disease=factor(factor(ph_sa$Diagnosis,labels=c("D","C")),levels=c("C","D")) ,
					px[, res[pmatch(colnames(px),rownames(res)),"qvalue"] <0.01 ] )
			tmp=tmp[trainindex,]  # train data only 
		
			#### gbm model 
			set.seed(2019)
			ctrl <- trainControl(method = "cv",
                     number = 2,
                     summaryFunction = twoClassSummary,
                     classProbs = TRUE,
                     allowParallel = TRUE,
                     savePredictions = TRUE )

			ctrl$sampling <- "smote" 

			grid <- expand.grid(n.trees=c(4000),
				shrinkage=c(0.002),
				n.minobsinnode = c(30),	
				interaction.depth=c(1,3)
				)

			smote_fit01 <- train(disease ~ .,
                  	 data = tmp,
                   	method = "gbm",
                  	 verbose = FALSE,
                   	tuneGrid=grid,
				bag.fraction =0.8)

			#### prediction held out samples
			testset=px[testindex,]
			pred=predict(smote_fit01,px[testindex,],type="prob")
			
			### rf model
			set.seed(2019)
			rf=randomForest(disease ~ ., data = tmp, ntree=1000,importance=F)			
			## predcit held out samples
			predrf=predict(rf,px[testindex,],type="prob")

			### merge the prediction , target for output
			out=data.frame( cbind(pred,predrf),target=disease[testindex],run=names(index)[ind])
			out
		}
		

	write.csv(res,"crc.cv.pred.csv")
	res=read.csv("crc.cv.pred.csv")
	res$run1=unlist(lapply(strsplit(as.character(res[,"run"]),".",fixed=T),function(x){x[length(x)]}))
	auc=c()
	auc1=c()
	for(i in unique(res$run1)){
		res1=res[res$run1==i,]
		auc=rbind(auc,as.numeric(roc(res1$target,res1[,"D"],ci=T)$ci) )
		auc1=rbind(auc1,as.numeric(roc(res1$target,res1[,"D.1"],ci=T)$ci))
	}
	#### report auc for the cv 
	auc




########################################################
#### fig
######################################################
	#### pub data
		dat <- read.table("abundance_profile_strain_pub.jgi.tsv",
			sep="\t",head=T,row.names=1)

		dat[1:4,1:4]
		dim(dat)

	### pub data phenotype
		phe<-read.table("abundance_profile_strain_pub.add_phenotype.jgi.tsv",
				sep="\t",head=T,as.is=T )
		head(phe)
		dim(phe)
		library(data.table)
		dat<-dcast(phe[,1:3],lineages_strain_new~sample_id)
		dat[1:4,1:4]
		rownames(dat)=dat[,1]
		dat=dat[,-1]
		dim(dat)

		phe1=phe[,c(2,4:ncol(phe))]
		phe1=phe1[!duplicated(phe1$sample_id),]
		dim(phe1)
		head(phe)
		tem=aggregate(phe1$Diagnosis,list(study=phe1$Study,diag=phe1$Diagnosis),length)
		tem


	
		ph_sa=phe1[phe1$Study=="ZellerG_2014" & phe1$Diagnosis !="small adenoma" & phe1$Country!="Germany",]
		table(ph_sa[,c("Diagnosis","Country")])
		ph_sa
		sp_sa=t(dat[,phe1$Study=="ZellerG_2014"& phe1$Diagnosis !="small adenoma"& phe1$Country!="Germany"])
		dim(sp_sa)


	### abundance are asin sqrt transform
	occu<-function(sp1){
		asin(sqrt(sp1[,apply(sp1,2,function(x){ mean(x>0)>0.1})]/100))
	}
	sp_sa2=occu(sp_sa)
	dim(sp_sa2)
	head(colnames(sp_sa2))

		### load package
		library(caret) # for model-building
		library(DMwR) # for smote implementation
		library(pROC) # for AUC calculations
		library(gbm)
		library(ppcor)
		library(randomForest)
		library(e1071)
		set.seed(0)
		library(doMC)
		registerDoMC(cores = 12)


	library(ppcor)
	pvalue=c()
	beta=c()
	pcor_p=c()
	pcor_r=c()
	wil_p=c()

		for( i in 1:ncol(sp_sa2)){

			tmp=data.frame(Gender=factor(ph_sa[,c("Sex")]),Age=ph_sa$Age,
				disease=factor(factor(ph_sa$Diagnosis,labels=c("D","C")),levels=c("C","D")) )
			tmp$y=as.numeric(t(sp_sa2[,i]))
		
			wil_p=c(wil_p,wilcox.test(tmp$y ~ tmp$disease)$p.value)
			pcor_p=c(pcor_p,pcor(data.matrix(tmp),method="spearman")$p.value["y","disease"])
			pcor_r=c(pcor_r,pcor(data.matrix(tmp),method="spearman")$estimate["y","disease"])
			mod=glm(y ~ disease + Gender  + Age, data=tmp) 
			pvalue=c(pvalue,coef(summary(mod))["diseaseD","Pr(>|t|)"])
			beta=c(beta,coef(summary(mod))["diseaseD","Estimate"])

		}
			qvalue=p.adjust(pvalue,method="fdr")
			wil_q=p.adjust(wil_p,method="fdr")
			pcor_q=p.adjust(pcor_p,method="fdr")
			res=data.frame(beta,pcor_r,pvalue,qvalue,wil_p,pcor_p,wil_q,pcor_q)
			rownames(res)=colnames(sp_sa2)




	px=sp_sa2[,pmatch( intersect(rownames(res),colnames(sp_sa2)),colnames(sp_sa2) )]
	tmp=data.frame( disease=factor(factor(ph_sa$Diagnosis,labels=c("D","C")),levels=c("C","D")) ,
			px[, res[pmatch(colnames(px),rownames(res)),"qvalue"] <0.01 ] )

		

	set.seed(2019)
	ctrl <- trainControl(method = "repeatedcv",
                     number = 10,
                     repeats = 5,
                     summaryFunction = twoClassSummary,
                     classProbs = TRUE,
                     allowParallel = TRUE,
                     savePredictions = TRUE )

	ctrl$sampling <- "smote" 

	grid <- expand.grid(n.trees=c(4000),
		shrinkage=c(0.002),
		n.minobsinnode = c(30),	
		interaction.depth=c(1,3)
		)
	smote_fit01 <- train(disease ~ .,
                   data = tmp,
                   method = "gbm",
                   verbose = FALSE,
                   metric = "ROC",
                   tuneGrid=grid,
			bag.fraction =0.8,
                   trControl = ctrl)


	set.seed(2019)
	rf=randomForest(disease ~ ., data = tmp, ntree=4000,importance=T)
	roc(tmp$disease,rf$votes[,2])

	rfimp=rf$importance


saveRDS(smote_fit01,"crc.smote_fit01.rds")
saveRDS(rf,"crc.rf.rds")

	### importance 
	var_imp01 = varImp(smote_fit01)
	res=data.frame(res)
	res$gbmimportance= var_imp01$importance$Overall[pmatch( 
					unlist(lapply(strsplit(rownames(res),"_"),function(x){x[length(x)]})) ,
				unlist(lapply(strsplit(rownames(var_imp01$importance),"_"),function(x){x[length(x)]}))    )]

	res$rfimportance= rfimp[pmatch( 
					unlist(lapply(strsplit(rownames(res),"_"),function(x){x[length(x)]})) ,
				unlist(lapply(strsplit(rownames(rfimp),"_"),function(x){x[length(x)]}))    ),4]

	res$rfrank=rank(res$rfimportance,na.last=F)
	res$rfrank=res$rfrank/max(res$rfrank)

	res=res[order(res$gbmimportance,decreasing=T),]
	#res=res[order(res$rfimportance,decreasing=T),]
	write.csv(res,"fig5.crc.csv")

	#### select by 3 method
		## mag, ref, size of new mgs
		size=read.table("oral_mgs_representative.tsv",sep="\t",head=T,as.is=T)
		size=size[,1:6]

	res1= data.frame( res, size[pmatch( 
					unlist(lapply(strsplit(rownames(res),"_"),function(x){x[length(x)]})) ,
				unlist(lapply(strsplit(size[,1],"_"),function(x){x[length(x)]}))    ),] )

	 tax=matrix(unlist(strsplit(rownames(res1),"|",fixed=T)),ncol=8,byrow=T)
	
	 res1=data.frame(res1,tax)	

		#rk=rank(abs(res1$pcor_r),na.last=F)+
		#rank(res1$gbmimportance,na.last=F)+
		#rank(res1$rfimportance,na.last=F)

	#res1=res1[order(rk,decreasing=T),]
	
	xor=read.table("oral.mgs.order",sep="\t",as.is=T)
	xor[,1]=gsub("OTU","mgs",xor[,1])
	res1$xorder=xor[pmatch(res1$mgs_id,xor[,1],duplicates.ok=T),2]
	res1$xorder[is.na(res1$xorder)]=max(res1$xorder,na.rm=T)
		res1$xorder=seq(1,1.5*14,1.5)[as.numeric(res1$X2)]+as.numeric(res1$pcor_r)

	write.csv(res1,"fig5_manhplot.crc.csv")
	res1=read.csv("fig5_manhplot.crc.csv")
nCHR <- length(unique(res1$X2))
phy_type = c("#004cce", "#059633", "#bdf2c1", "#CB1414", "#C3C00E", 
             "#9c5999", "#ffc49c", "#7a4b03", "#95cde8", "#e59735", "#ababab",
			"steelblue", "red3",'#999999','#E69F00', '#56B4E9')
sig=0.01
 manhplot <-ggplot(res1, aes(x=xorder, y=-log10(qvalue), size=-log10(qvalue) )) +
  geom_point(aes(color=as.factor(X2),shape=as.factor(mtype))) +
  scale_color_manual(values = phy_type ) +
  scale_shape_manual(values = c(1,19) ) +
	scale_size_continuous(range = c(0.5,3))+
  #scale_x_continuous(label = axis.set$CHR, breaks = axis.set$center) +
  #scale_y_continuous(expand = c(0,0), limits = c(0,ylim)) +
  geom_hline(yintercept = -log10(sig), color = "grey40", linetype = "dashed") + 
  labs(x = NULL, y = "-log10(p)") + 
  theme_minimal() +
  theme( 
    #legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5)
  )

manhplot


ggsave("fig5.manhattanplot.crc.pdf", manhplot, width = 30, height = 15, units = "cm")

	### importance plot
	
	mtcars2 <- data.frame(type = factor(mtcars$cyl), group = factor(mtcars$gear))
	library(plyr); library(ggplot2)
	dat <- rbind(ddply(mtcars2, .(type, group), summarise, count = length(group)), c(8, 4, NA))

	res2=res1[res1$size>10,]
	res2=res2[order(res2$gbmimportance,decreasing=T),][1:50,]
	res2$tax2= unlist(lapply(strsplit(rownames(res2),"|",fixed=T),function(x){ 
				x1=x[!grepl("unclassified",x)];
				paste(x1[length(x1)]	, 
				unlist(lapply(strsplit(x[length(x)],"_"),function(j){j[length(j)]})),sep="_"  )
			 }))
	res2$size1=paste(res2$MAG,"/",res2$size,sep="")

	write.csv(res2,"fig5_barplot.crc.csv")
	res2=read.csv("fig5_barplot.crc.csv")

	p2 <- ggplot(res2, aes(x = reorder(tax2,pcor_r*(gbmimportance+1)) ,y = -log10(qvalue),fill = factor(sign(beta)))) + 
  		geom_bar(stat = "identity", position = "dodge", width = 0.8) +
		geom_point( aes(y = sqrt(gbmimportance), x = reorder(tax2,pcor_r) ),pch=15,col="red3") +
		geom_point( aes(y = rfimportance*10, x = reorder(tax2,pcor_r) ),pch=17) +
 		 geom_hline(yintercept = -log10(sig), color = "grey40", linetype = "dashed") +
 	 	ylim(0, 14) +
  		geom_text(   aes(x =  reorder(tax2,pcor_r) ,y = -log10(qvalue), label =paste(MAG,"/",size,sep="") ), 
   				 hjust = -0.5, size = 2,
   				 position = position_dodge(width = 1),
   				 inherit.aes = TRUE )+
		 #scale_fill_brewer(palette="Greens") +
			scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
		 coord_flip()+
 		 theme_minimal() +
 		 theme( 
   		 #legend.position = "none",
   		 panel.border = element_blank(),
   		 panel.grid.major.x = element_blank(),
   		 panel.grid.minor.x = element_blank(),
   		 axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5)
		)
	p2
ggsave("fig5.a.barplot.crc.pdf", p2, width = 30, height = 15, units = "cm")


