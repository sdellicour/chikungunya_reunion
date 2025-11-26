treeExtractions_DTA = function(dta_tre, mostRecentSamplingDatum)
	{
		dta_tab = data.frame(matrix(nrow=dim(dta_tre$edge)[1], ncol=5))
		colnames(dta_tab) = c("node1","node2","length","startLoc","endLoc")
		dta_tab[,c("node1","node2")] = as.numeric(dta_tre$edge)
		dta_tab[,c("length")] = as.numeric(dta_tre$edge.length)
		for (i in 1:length(dta_tre$annotations))
			{
				annotations = dta_tre$annotations[[i]]
				if (!is.null(dta_tre$annotations[[i]]))
					{
						dta_tab[i,"endLoc"] = annotations$location
					}
			}
		for (i in 1:length(dta_tre$annotations))
			{
				index = which(dta_tab[,"node2"] == dta_tab[i,"node1"])
				if (length(index) > 0)
					{
						dta_tab[i,"startLoc"] = dta_tab[index,"endLoc"]
					}	else	{
						annotations = dta_tre$root.annotation
						dta_tab[i,"startLoc"] = annotations$location
					}
			}	
		l = length(dta_tab[,1]); ll = matrix(1:l,nrow=l,ncol=l); ll[] = 0
		for (j in 1:l)
			{
				subMat = dta_tab[j,2]
				subMat = subset(dta_tab,dta_tab[,2]==subMat)
				ll[j,1] = subMat[,3]
				subMat = subMat[1,1]
				subMat1 = subset(dta_tab,dta_tab[,2]==subMat)
				for (k in 1:l)
					{
						if (nrow(subMat1) > 0)
							{
								ll[j,k+1] = subMat1[,3]
	 							subMat2 = subMat1[1,1]
	 							subMat1 = subset(dta_tab,dta_tab[,2]==subMat2)
	 						}
	 				}
			}
		endNodeL = rowSums(ll)
		dta_tab = cbind(dta_tab, endNodeL)
		startNodeL = matrix(1:l,nrow=l,ncol=1)
		startNodeL[] = 0
		for (j in 1:l)
			{
				r = dta_tab[j,1]
				s = subset(dta_tab,dta_tab[,2]==r)
				for (k in 1:l)
					{
						if (nrow(s) > 0)
							{
								startNodeL[j,1] = s[,dim(s)[2]]
	 						}
	 				}	
			}
		colnames(startNodeL) = "startNodeL"
		dta_tab = cbind(dta_tab,startNodeL)
		maxEndLIndice = which.max(dta_tab[,"endNodeL"])
		maxEndL = dta_tab[maxEndLIndice,"endNodeL"]
	 	endYear = matrix(dta_tab[,"endNodeL"]-maxEndL)
	 	endYear = matrix(mostRecentSamplingDatum+(endYear[,1]))
		startYear = matrix(dta_tab[,"startNodeL"]-maxEndL)
	 	startYear = matrix(mostRecentSamplingDatum+(startYear[,1]))
	 	colnames(startYear) = "startYear"; colnames(endYear) = "endYear"
	 	dta_tab = cbind(dta_tab,startYear,endYear)
		dta_tab = dta_tab[order(dta_tab[,"startYear"],decreasing=F),]
		dta_tab1 = dta_tab[1,]; dta_tab2 = dta_tab[2:dim(dta_tab)[1],]
		dta_tab2 = dta_tab2[order(dta_tab2[,"endYear"],decreasing=F),]
		dta_tab = rbind(dta_tab1, dta_tab2); dta_tab = dta_tab[,c(1,2,4,5,3,6:9)]
		tipLabels = matrix(nrow=dim(dta_tab)[1], ncol=1); colnames(tipLabels) = "tipLabel"
		for (i in 1:length(dta_tre$tip.label))
			{
				tipLabels[which(dta_tab[,"node2"]==i),"tipLabel"] = gsub("'","",dta_tre$tip.label[i])
			}
		dta_tab = cbind(dta_tab, tipLabels)
		return(dta_tab)
	}
