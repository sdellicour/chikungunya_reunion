library(adephylo)
library(colorspace)
library(diagram)
library(EpiEstim)
library(HDInterval)
library(lubridate)
library(MetBrewer)
library(nlme)
library(seraphim)
library(sf)
library(treeio)

analysis = "101025"
showingPlots = FALSE

# 1. Investigating the temporal signal and sampling representativity
# 2. Preliminary BEAST analysis to estimate the substitution rate
# 3. Preparing the skygrid and discrete phylogeographic analyses
# 4. Preparing the different predictors for the discrete-GLM analysis
# 5. Preparing climatic predictors for the skygrid/EBDS-GLM analyses
# 6. Evolution of Ne and R(t): skygrid and sampling aware analyses
# 7. Discrete phylogeographic analysis based on the municipalities
# 8. Analysing and reporting the results of the discrete-GLM analysis
# 9. Conducting complementary isolation-by-distance (IBD) analyses

tab = read.table(paste0("BEAST_CTA_analysis/Alignment_",analysis,".txt"), head=T)
mostRecentSamplingDatum = max(decimal_date(ymd(tab[,"collection_date"]))) # 2025.595
admin0 = shapefile("GADM_REU_shapefile/GADM_REU_0.shp") # island borders
admins2 = shapefile("GADM_REU_shapefile/GADM_REU_2.shp") # municipalities
locations = unique(admins2@data[,"GID_2"]) # municipality IDs
elevation = raster::mask(crop(raster("Copernicus_DEM30.tif"), admin0), admin0)
elevation[elevation[]<0] = 0; elevation[elevation[]>3000] = 3000
osm_lines = st_read("OpenStreetMap_files/OSM_25102014.pbf", layer="lines")
osm_polygons = st_read("OpenStreetMap_files/OSM_25102014.pbf", layer="multipolygons")
main_roads = subset(osm_lines, highway%in%c("trunk","primary"))
residential_areas = subset(osm_polygons, landuse=="residential")
construction_areas = subset(osm_polygons, landuse=="construction")
forest_areas = subset(osm_polygons, landuse=="forest")
elevation_cols = paste0(divergingx_hcl(14,"fall")[4:14],"BF") # 75% transparency
collection_week_cols1 = met.brewer(name="Hiroshige", n=60, type="continuous")[1:51]
collection_week_cols2 = met.brewer(name="Hiroshige", n=60, type="continuous")[1:52]
collection_month_cols = met.brewer(name="Hiroshige", n=14, type="continuous")[1:12]

# 1. Investigating the temporal signal and sampling representativity

	# 1.1. Generating a second alignment and associated metadata files without outliers identified with the first TempEst analysis

tre1 = read.tree(paste0("Temporal_signal_test1/Alignment_",analysis,".tre")); fas2 = c()
tab1 = read.csv(paste0("Temporal_signal_test1/Alignment_",analysis,".csv"), head=T, sep=";")
fas1 = scan(paste0("Temporal_signal_test1/Alignment_",analysis,".fas"), what="", sep="\n", quiet=T)
outliers = read.table(paste0("Temporal_signal_test1/Outliers_",gsub("25","2025",analysis),".txt"))
tre2 = drop.tip(tre1, outliers[,1]); tab2 = tab1[which(!tab1[,"trait"]%in%outliers[,1]),]
for (i in 1:dim(tab2)[1])
	{
		if ((length(unlist(strsplit(tab2[i,2],"-"))) == 3) && (nchar(unlist(strsplit(tab2[i,2],"-"))[3]) == 4))
			{
				txt = unlist(strsplit(tab2[i,2],"-")); tab2[i,2] = paste(txt[3],txt[2],txt[1],sep="-")
			}
	}
for (i in 1:length(fas1))
	{
		if ((grepl(">",fas1[i])) && (!gsub(">","",gsub("\\|","_",fas1[i]))%in%outliers[,1])) fas2 = c(fas2, fas1[c(i,i+1)])
	}
for (i in 1:length(fas2)) fas2[i] = gsub("\\|","_",fas2[i])
write.csv(tab2, paste0("Temporal_signal_test2/Alignment_",analysis,".csv"), row.names=F, quote=F)
write.table(tab2, paste0("Temporal_signal_test2/Alignment_",analysis,".txt"), row.names=F, quote=F, sep="\t")
write.tree(tre2, paste0("Temporal_signal_test2/Alignment_",analysis,".tre"))
write(fas2, paste0("Temporal_signal_test2/Alignment_",analysis,".fas"))

	# 1.2. Investigating the correlation between the number of genomic samples and the number of cases per week and per municipality

tab1 = read.table(paste0("Temporal_signal_test1/TempEst_regression.txt"), head=T, sep="\t")
tab2 = read.table(paste0("Temporal_signal_test2/TempEst_regression.txt"), head=T, sep="\t")
rP = cor(tab1[,"date"], tab1[,"distance"], method="pearson") # rP (Pearson) = 0.82
lr1 = lm(distance ~ date, data=tab1); R2b = 0.67 # from the 1° root-to-tips regression analysis
rP = cor(tab2[,"date"], tab2[,"distance"], method="pearson") # rP (Pearson) = 0.84
lr2 = lm(distance ~ date, data=tab2); R2b = 0.71 # from the 2° root-to-tips regression analysis
tab3 = read.csv("Cases_symptoms.csv", head=T); buffer = matrix(nrow=dim(tab3)[1], ncol=1)
colnames(buffer) = "sequences"; tab3 = cbind(tab3, buffer)
for (i in 1:dim(tab3)[1])
	{
		if (i < dim(tab3)[1]) indices = which((ymd(tab[,"collection_date"])>=dmy(tab3[i,"start"]))&(ymd(tab[,"collection_date"])<dmy(tab3[i+1,"start"])))
		if (i == dim(tab3)[1]) indices = which(ymd(tab[,"collection_date"])>=dmy(tab3[i,"start"]))
		tab3[i,"sequences"] = length(indices)
	}
rP = cor(tab3[,"cases"], tab3[,"sequences"], method="pearson") # rP (Pearson) = 0.95
lr3 = lm(sequences ~ cases, data=tab3); R2c = 0.90 # pretty good correlation between the two
tab4 = read.csv("Cases_municipality.csv", head=T)
rP = cor(tab4[,"cases"], tab4[,"sequences"], method="pearson") # rP (Pearson) = 0.86
lr4 = lm(sequences ~ cases, data=tab4); R2d = 0.75 # good correlation between the two

	# 1.3. Visualising the root-to-tip regression analyses and the associations between the number of genomic samples and cases 

pdf(paste0("Figure_S2_",analysis,"_NEW.pdf"), width=8, height=5) # dev.new(width=8, height=5)
par(mfrow=c(2,2), oma=c(0,0.6,0.5,0), mar=c(3.3,4.5,0.5,1), lwd=0.3, bty="o", col="gray30", col.axis="gray30", fg="gray30")
plot(tab3[,"cases"], tab3[,"sequences"], col=NA, axes=F, ann=F, xlim=c(-170,8500), ylim=c(-17,305))
abline(lr3, col=rgb(222,67,39,255,maxColorValue=255), lwd=0.75, lty=2)
points(tab3[,"cases"], tab3[,"sequences"], pch=16, cex=1.0, col=rgb(70,118,187,100,maxColorValue=255))
points(tab3[,"cases"], tab3[,"sequences"], pch=1, cex=1.0, lwd=0.5, col=rgb(70,118,187,255,maxColorValue=255))
ats1 = seq(-1500, 9000, 1500); ats2 = seq(-50, 350, 50)
axis(side=1, lwd=0.5, cex.axis=0.85, mgp=c(0,0.30,0), lwd.tick=0.5, col.lab="gray30", col="gray30", tck=-0.03, las=1, at=ats1)
axis(side=2, lwd=0.5, cex.axis=0.85, mgp=c(0,0.55,0), lwd.tick=0.5, col.lab="gray30", col="gray30", tck=-0.03, las=1, padj=0.4, at=ats2)
title(xlab="Number of weekly cases", mgp=c(1.5,0,0), cex.lab=1.0, col.lab="gray30")
title(ylab="Genomic samples/week     ", mgp=c(2.7,0,0), cex.lab=1.0, col.lab="gray30")
mtext(expression(bold(A)), at=-2450, line=-0.95, cex.lab=0.80, col="gray30")
mtext(expression("Pearson correlation"~"="~"0.82"), at=2160, line=-0.70, cex=0.75, col=rgb(222,67,39,255,maxColorValue=255))
mtext(expression(R^2~"="~"0.90"), at=700, line=-1.80, cex=0.75, col=rgb(222,67,39,255,maxColorValue=255))
plot(tab4[,"cases"], tab4[,"sequences"], col=NA, axes=F, ann=F, xlim=c(-170,8500), ylim=c(-17,500))
abline(lr4, col=rgb(222,67,39,255,maxColorValue=255), lwd=0.75, lty=2)
points(tab4[,"cases"], tab4[,"sequences"], pch=16, cex=1.0, col=rgb(70,118,187,100,maxColorValue=255))
points(tab4[,"cases"], tab4[,"sequences"], pch=1, cex=1.0, lwd=0.5, col=rgb(70,118,187,255,maxColorValue=255))
ats1 = seq(-1500, 9000, 1500); ats2 = seq(-100, 600, 100)
axis(side=1, lwd=0.5, cex.axis=0.85, mgp=c(0,0.30,0), lwd.tick=0.5, col.lab="gray30", col="gray30", tck=-0.03, las=1, at=ats1)
axis(side=2, lwd=0.5, cex.axis=0.85, mgp=c(0,0.55,0), lwd.tick=0.5, col.lab="gray30", col="gray30", tck=-0.03, las=1, padj=0.4, at=ats2)
title(xlab="Number of cases per municipality", mgp=c(1.5,0,0), cex.lab=1.0, col.lab="gray30")
title(ylab="Gen. samples/municipality       ", mgp=c(2.7,0,0), cex.lab=1.0, col.lab="gray30")
mtext(expression(bold(B)), at=-2450, line=-0.95, cex.lab=0.80, col="gray30")
mtext(expression("Pearson correlation"~"="~"0.84"), at=2160, line=-0.70, cex=0.75, col=rgb(222,67,39,255,maxColorValue=255))
mtext(expression(R^2~"="~"0.86"), at=700, line=-1.80, cex=0.75, col=rgb(222,67,39,255,maxColorValue=255))
plot(tab1[,"date"], tab1[,"distance"], col=NA, axes=F, ann=F, xlim=c(1961.5,2025), ylim=c(0.0045,0.022))
abline(lr1, col=rgb(222,67,39,255,maxColorValue=255), lwd=0.75, lty=2)
points(tab1[,"date"], tab1[,"distance"], pch=16, cex=1.0, col=rgb(70,118,187,100,maxColorValue=255))
points(tab1[,"date"], tab1[,"distance"], pch=1, cex=1.0, lwd=0.5, col=rgb(70,118,187,255,maxColorValue=255))
ats1 = seq(1950, 2030, 10); ats2 = seq(0, 0.025, 0.005)
axis(side=1, lwd=0.5, cex.axis=0.85, mgp=c(0,0.30,0), lwd.tick=0.5, col.lab="gray30", col="gray30", tck=-0.03, las=1, at=ats1)
axis(side=2, lwd=0.5, cex.axis=0.85, mgp=c(0,0.55,0), lwd.tick=0.5, col.lab="gray30", col="gray30", tck=-0.03, las=1, padj=0.4, at=ats2)
title(xlab="Time", mgp=c(1.5,0,0), cex.lab=1.0, col.lab="gray30")
title(ylab="Root-to-tip divergence    ", mgp=c(2.7,0,0), cex.lab=1.0, col.lab="gray30")
mtext(expression(bold(C)), at=1945, line=-0.95, cex.lab=0.80, col="gray30")
mtext(expression("Pearson correlation"~"="~"0.95"), at=1978.2, line=-0.70, cex=0.75, col=rgb(222,67,39,255,maxColorValue=255))
mtext(expression(R^2~"="~"0.67"), at=1967.5, line=-1.80, cex=0.75, col=rgb(222,67,39,255,maxColorValue=255))
plot(tab2[,"date"], tab2[,"distance"], col=NA, axes=F, ann=F, xlim=c(1961.5,2025), ylim=c(0.0045,0.022))
abline(lr2, col=rgb(222,67,39,255,maxColorValue=255), lwd=0.75, lty=2)
points(tab2[,"date"], tab2[,"distance"], pch=16, cex=1.0, col=rgb(70,118,187,100,maxColorValue=255))
points(tab2[,"date"], tab2[,"distance"], pch=1, cex=1.0, lwd=0.5, col=rgb(70,118,187,255,maxColorValue=255))
ats1 = seq(1950, 2030, 10); ats2 = seq(0, 0.025, 0.005)
axis(side=1, lwd=0.5, cex.axis=0.85, mgp=c(0,0.30,0), lwd.tick=0.5, col.lab="gray30", col="gray30", tck=-0.03, las=1, at=ats1)
axis(side=2, lwd=0.5, cex.axis=0.85, mgp=c(0,0.55,0), lwd.tick=0.5, col.lab="gray30", col="gray30", tck=-0.03, las=1, padj=0.4, at=ats2)
title(xlab="Time", mgp=c(1.5,0,0), cex.lab=1.0, col.lab="gray30")
title(ylab="Root-to-tip divergence    ", mgp=c(2.7,0,0), cex.lab=1.0, col.lab="gray30")
mtext(expression(bold(D)), at=1945, line=-0.95, cex.lab=0.80, col="gray30")
mtext(expression("Pearson correlation"~"="~"0.86"), at=1978.2, line=-0.70, cex=0.75, col=rgb(222,67,39,255,maxColorValue=255))
mtext(expression(R^2~"="~"0.71"), at=1967.5, line=-1.80, cex=0.75, col=rgb(222,67,39,255,maxColorValue=255))
dev.off()

# 2. Preliminary BEAST analysis to estimate the substitution rate

tab = read.csv(paste0("BEAST_first_analysis/Alignment_",analysis,".csv"), head=T)
log = read.table(paste0("BEAST_first_analysis/Alignment_",analysis,".log"), head=T, sep="\t")
ucld_mean = log[ceiling(0.1*dim(log)[1]):dim(log)[1],"ucld.mean"] # substitution rate
tree = readAnnotatedNexus(paste0("BEAST_first_analysis/Alignment_",analysis,".tree")) # MCC tree
rootHeight = max(nodeHeights(tree)); root_time = mostRecentSamplingDatum-rootHeight
minYear = mostRecentSamplingDatum-tree$root.annotation$`height_95%_HPD`[[2]]
maxYear = mostRecentSamplingDatum; nodes_cols = c(); locations = c()
locations = unique(tab[,"location"]); locations = locations[order(locations)]
location_cols = colorRampPalette(brewer.pal(11,"Spectral"))(length(locations))
location_cols = c("#FAA521","#4676BB","#D1E5F0","#C0C0C0","#6A3D9A","#DE4327","#BF812D")

pdf(paste0("Figure_S3A_",analysis,"_NEW.pdf"), width=8, height=3.5) # dev.new(width=8, height=4)
par(oma=c(0,0,0,0), mar=c(0.8,2.2,4.3,0), lwd=0.3, bty="o", col="gray30", col.axis="gray30", fg="gray30", lheight=0.85)
plot(tree, show.tip.label=F, show.node.label=F, edge.width=0.7, cex=0.6, align.tip.label=3, direction="downwards",
	 y.lim=c(0,rootHeight), col="gray30", edge.color="gray30")
tree_obj = get("last_plot.phylo", envir=.PlotPhyloEnv); rootBarPlotted = FALSE
gray90_transparent = rgb(229, 229, 229, 150, names=NULL, maxColorValue=255)
for (j in 1:dim(tree$edge)[1])
	{
		if ((tree$edge[j,2]%in%tree$edge[,1])&(length(tree$annotations[[j]]$`height_95%_HPD`) > 1))
			{
				y1 = mostRecentSamplingDatum-(mostRecentSamplingDatum-tree$annotations[[j]]$`height_95%_HPD`[[2]])
				y2 = mostRecentSamplingDatum-(mostRecentSamplingDatum-tree$annotations[[j]]$`height_95%_HPD`[[1]])
				lines(x=rep(tree_obj$xx[tree$edge[j,2]],2), y=c(y2,y1), lwd=3.5, lend=0, col=gray90_transparent)
			}
		if ((rootBarPlotted == FALSE)&&(!tree$edge[j,1]%in%tree$edge[,2]))
			{
				y1 = mostRecentSamplingDatum-(mostRecentSamplingDatum-tree$root.annotation$`height_95%_HPD`[[2]])
				y2 = mostRecentSamplingDatum-(mostRecentSamplingDatum-tree$root.annotation$`height_95%_HPD`[[1]])
				lines(x=rep(tree_obj$xx[tree$edge[j,1]],2), y=c(y2,y1), lwd=3.5, lend=0, col=gray90_transparent)
				rootBarPlotted = TRUE
			}
	}
for (i in 1:dim(tree$edge)[1])
	{
		if (tree$edge[i,2]%in%tree$edge[,1])
			{
				nodelabels(node=tree$edge[i,2], pch=16, cex=0.4, col="gray30")
			}
		if (!tree$edge[i,2]%in%tree$edge[,1])
			{
				tip_col = location_cols[which(locations==tab[which(tab[,"sequence_ID"]==tree$tip.label[tree$edge[i,2]]),"location"])]
				nodelabels(node=tree$edge[i,2], pch=16, cex=0.8, col=tip_col)
				nodelabels(node=tree$edge[i,2], pch=1, cex=0.8, col="gray30", lwd=0.25)
			}
	}
dates_to_plot = seq(1900,2030,10); dates_to_print = dates_to_plot; indices = c(1,2,5,7,3,4,6)
axis(lwd=0.5, at=mostRecentSamplingDatum-dates_to_plot, labels=dates_to_print, cex.axis=0.60, 
	 mgp=c(0,0.1,-0.3), lwd.tick=0.5, col.lab="gray30", col="gray30", tck=-0.013, side=2, las=1)
legend(x=130, y=60, gsub("_"," ",locations)[indices], col=location_cols[indices], text.col="gray30", pch=16, pt.cex=1.25, box.lty=0, cex=0.65, x.intersp=0.8, y.intersp=1.1)
legend(x=130, y=60, gsub("_"," ",locations)[indices], col="gray30", text.col=rgb(0,0,0,0), pch=1, pt.cex=1.25, box.lty=0, cex=0.65, pt.lwd=0.3, x.intersp=0.8, y.intersp=1.1)
dev.off() # two manual edits to conduct in Illustrator: (i) removing the clipping mask on the tree, and (ii) putting the horizontal bars behind

pdf(paste0("Figure_S3B_",analysis,"_NEW.pdf"), width=3, height=1.5) # dev.new(width=3, height=2)
par(oma=c(0,0,0,0), mar=c(1.5,2,0.7,1), lwd=0.3, bty="o", col="gray30", col.axis="gray30", fg="gray30")
plot(density(ucld_mean), axes=F, ann=F, col=NA); polygon(stats::density(ucld_mean), col="#DE43274D", border=NA); lines(density(ucld_mean), col="#DE4327", lwd=1)
ats = seq(0.0002,0.0007,0.0001); labels = c("",expression("3x10"^"-4"),expression("4x10"^"-4"),expression("5x10"^"-4"),expression("6x10"^"-4"),"")
axis(side=1, lwd=0.5, cex.axis=0.60, mgp=c(0,0,0), lwd.tick=0.5, col.lab="gray30", col="gray30", tck=-0.03, las=1, at=ats, label=labels)
axis(side=2, lwd=0.5, cex.axis=0.60, mgp=c(0,0.1,-0.3), lwd.tick=0.5, col.lab="gray30", col="gray30", tck=-0.03, las=1, padj=0.22)
dev.off()

# 3. Preparing the skygrid and discrete phylogeographic analyses

fas = scan(paste0("Alignment_",analysis,".fas"), what="", sep="\n", quiet=T)
tab1 = read.table(paste0("Alignment_",analysis,".txt"), head=T)
tab2 = read.csv(paste0("Alignment_",analysis,".csv"), head=T, sep=";")
tab = matrix(nrow=dim(tab2)[1], ncol=5); colnames(tab) = c("trait","collection_date","location","longitude","latitude")
for (i in 1:dim(tab)[1])
	{
		tab[i,"trait"] = tab1[i,"trait"]; tab[i,"collection_date"] = tab1[i,"collection_date"]
		index = which(tab2[,"trait"]==tab[i,"trait"]); J = NA
		for (j in 1:length(admins2))
			{
				pol_coords = admins2@polygons[[j]]@Polygons[[1]]@coords
				if (point.in.polygon(tab2[index,"longitude"],tab2[index,"latitude"],pol_coords[,"x"],pol_coords[,"y"]) == 1) J = j
			}
		tab[i,"location"] = admins2@data[J,"GID_2"]
		tab[i,"longitude"] = tab2[index,"longitude"]
		tab[i,"latitude"] = tab2[index,"latitude"]
	}
write.table(tab, paste0("BEAST_DTA_analysis/Alignment_",analysis,".txt"), row.names=F, quote=F, sep="\t")
write(fas, paste0("BEAST_DTA_analysis/Alignment_",analysis,".fas"))

# 4. Preparing the different predictors for the discrete-GLM analysis

	# 4.1. Population count at the municipality of origin and at the destination municipality

population_counts = matrix(nrow=dim(admins2)[1], ncol=1); correspondences = matrix(nrow=dim(admins2@data)[1], ncol=2)
colnames(population_counts) = c("population_count"); row.names(population_counts) = admins2@data[,"GID_2"]
population_data = read.csv("GLM_predictors_data/Original_data_files/Municipality_population_INSEE.csv", head=T, sep=";")
population_data = population_data[which(population_data[,"Département"]==974),c("Code.géographique","Libellé.géographique","Population.en.2021")]
admins2@data[which(admins2@data[,"NAME_2"]=="L'Entre-Deux"),"NAME_2"] = "Entre-Deux"
for (i in 1:dim(population_counts)[1])
	{
		index = which(population_data[,"Libellé.géographique"]==admins2@data[i,"NAME_2"])
		population_counts[i] = as.numeric(gsub(" ","",population_data[index,"Population.en.2021"]))
		correspondences[i,] = cbind(admins2@data[index,"GID_2"], population_data[i,"Code.géographique"])
	}
write.csv(population_counts, "GLM_predictors_data/Prepared_variables/Discrete-GLM_analyses/Population_counts.csv", quote=F)

	# 4.2. Pairwise great-circle geographic distance between population-weighted centroid points

centroids = matrix(nrow=dim(admins2@data)[1], ncol=2); row.names(centroids) = admins2@data[,"GID_2"]
population = raster("GLM_predictors_data/Original_data_files/Population_2025_100m_R2025A.tif")
for (i in 1:dim(centroids)[1])
	{
		p = Polygon(admins2@polygons[[i]]@Polygons[[1]]@coords); ps = Polygons(list(p),1); sps = SpatialPolygons(list(ps))
		r = mask(crop(population, sps), sps); tab = matrix(nrow=sum(!is.na(r[])), ncol=3); colnames(tab) = c("pop","x","y")
		indices = which(!is.na(r[]))
		for (j in 1:length(indices))
			{
				tab[j,"pop"] = r[indices[j]]; tab[j,c("x","y")] = xyFromCell(r, indices[j])
			}
		centroids[i,1] = sum(tab[,"x"]*(tab[,"pop"]/sum(tab[,"pop"])))
		centroids[i,2] = sum(tab[,"y"]*(tab[,"pop"]/sum(tab[,"pop"])))
	}
geographic_distances = matrix(nrow=dim(centroids)[1], ncol=dim(centroids)[1])
row.names(geographic_distances) = admins2@data[,"GID_2"]
colnames(geographic_distances) = admins2@data[,"GID_2"]
for (i in 1:dim(centroids)[1])
	{
		for (j in 1:dim(centroids)[1])
			{
				if (i == j)
					{
						geographic_distances[i,j] = 0
					}	else	{
						geographic_distances[i,j] = rdist.earth(cbind(centroids[i,1],centroids[i,2]), cbind(centroids[j,1],centroids[j,2]), miles=F)
					}
			}
	}
write.csv(geographic_distances, "GLM_predictors_data/Prepared_variables/Discrete-GLM_analyses/Geographic_distances.csv", quote=F)

	# 4.3. Pairwise measure of the professional (commuting) mobility flux between municipalities (symetric; for the main GLM analysis)

mobility_metric = matrix(nrow=dim(admins2@data)[1], ncol=dim(admins2@data)[1])
row.names(mobility_metric) = admins2@data[,"GID_2"]; colnames(mobility_metric) = admins2@data[,"GID_2"]
mobility_data = read.csv("GLM_predictors_data/Original_data_files/Mobility_data_INSEE_for_2021.csv", head=T, sep=";")
for (i in 1:dim(admins2@data)[1])
	{
		for (j in 1:dim(admins2@data)[1])
			{
				index1 = which((mobility_data[,"CODGEO"]==correspondences[which(correspondences[,1]==admins2@data[i,"GID_2"]),2])
							  &(mobility_data[,"DCLT"]==correspondences[which(correspondences[,1]==admins2@data[j,"GID_2"]),2]))
				index2 = which((mobility_data[,"CODGEO"]==correspondences[which(correspondences[,1]==admins2@data[j,"GID_2"]),2])
							  &(mobility_data[,"DCLT"]==correspondences[which(correspondences[,1]==admins2@data[i,"GID_2"]),2]))
				if ((length(index1) != 0)&(length(index2) != 0))
					{
						mobility_metric[i,j] = (mobility_data[index1,"NBFLUX_C21_ACTOCC15P"]+mobility_data[index2,"NBFLUX_C21_ACTOCC15P"])/2
					}
				if ((length(index1) != 0)&(length(index2) == 0) )mobility_metric[i,j] = mobility_data[index1,"NBFLUX_C21_ACTOCC15P"]
				if ((length(index1) == 0)&(length(index2) != 0)) mobility_metric[i,j] = mobility_data[index2,"NBFLUX_C21_ACTOCC15P"]
				if ((length(index1) == 0)&(length(index2) == 0)) mobility_metric[i,j] = 0.1
			}
	}
write.csv(mobility_metric, "GLM_predictors_data/Prepared_variables/Discrete-GLM_analyses/Commuting_metric.csv", quote=F)

# 5. Preparing climatic predictors for the skygrid/EBDS-GLM analyses

	# 5.1. Preparation of monthly temperature and precipitation data (source: https://meteo.data.gouv.fr/datasets/6569b3d7d193b4daf2b43edc)

selected_months = c("202406","202407","202408","202409","202410","202411","202412","202501","202502","202503","202504","202505","202506","202507","202508")
climatic_variables = matrix(nrow=length(selected_months), ncol=5); colnames(climatic_variables) = c("year","month","date","temperature","precipitation")
climatic_data = read.csv("GLM_predictors_data/Original_data_files/Monthly_climate_data_2024-25.csv", head=T, sep=";")
climatic_data = climatic_data[which(climatic_data[,"AAAAMM"]%in%selected_months),]; pol_coords = admin0@polygons[[1]]@Polygons[[1]]@coords
climatic_data = climatic_data[which(point.in.polygon(climatic_data[,"LON"],climatic_data[,"LAT"],pol_coords[,"x"],pol_coords[,"y"])==1),]
lon_lat = paste0(climatic_data[,"LON"],"_",climatic_data[,"LAT"]); indices = which((is.na(climatic_data[,"TMM"]))|(is.na(climatic_data[,"RR"])))
lon_lat_to_discard = unique(lon_lat[indices]); climatic_data = climatic_data[which(!lon_lat%in%lon_lat_to_discard),]
for (i in 1:length(selected_months))
	{
		if (selected_months[i] != "202412")
			{
				climatic_variables[i,"date"] = mean(decimal_date(ym(c(selected_months[i],as.character(as.numeric(selected_months[i])+1)))))
			}	else	{
				climatic_variables[i,"date"] = mean(decimal_date(ym(c(selected_months[i],"202501"))))
			}
		climatic_variables[i,"year"] = substr(selected_months[i], 1, 4); climatic_variables[i,"month"] = substr(selected_months[i], 5, 6)
		climatic_variables[i,"temperature"] = mean(climatic_data[which(climatic_data[,"AAAAMM"]==selected_months[i]),"TMM"]) # °C
		climatic_variables[i,"precipitation"] = mean(climatic_data[which(climatic_data[,"AAAAMM"]==selected_months[i]),"RR"]) # mm
	}
write.csv(climatic_variables, "GLM_predictors_data/Prepared_variables/EBDS_skygrid_analyses/Climatic_variables_monthly.csv", row.names=F, quote=F)

	# 5.2. Preparation of weekly temperature and precipitation data (source: https://cds.climate.copernicus.eu/; retrieved with a Python script)

temperatures = brick("GLM_predictors_data/Original_data_files/ERA5-Land_daily_temperatures.nc")
precipitations = brick("GLM_predictors_data/Original_data_files/ERA5-Land_daily_precipitations.nc")
dates = as.POSIXct("1970-01-01 00:00:00", format="%Y-%m-%d %H:%M:%OS") + seconds(data@z$valid_time)

		# 5.2.1. Prepration of the climate data for the skygrid-GLM analyses

climatic_variables_xml = rep("\t\t\tvalue=\"",2); epoch_width = 0.01960784
climatic_variables = matrix(nrow=60, ncol=3); colnames(climatic_variables) = c("date","temperature","precipitation")
for (i in 1:60)
	{
		date1 = mostRecentSamplingDatum-(i*epoch_width); date2 = date1+epoch_width
		climatic_variables[i,"date"] = (date1+date2)/2
		indices = which((decimal_date(dates)>date1)&(decimal_date(dates)<=date2))
		vS = rep(NA, length(indices))
		for (j in 1:length(vS))
			{
				vS[j] = mean(temperatures[[indices[j]]][], na.rm=T)-273.15 # Kelvins to °C
			}
		climatic_variables[i,"temperature"] = mean(vS)
		climatic_variables_xml[1] = paste0(climatic_variables_xml[1],log(mean(vS))," ")
		for (j in 1:length(vS))
			{
				vS[j] = mean(precipitations[[indices[j]]][], na.rm=T)*1000 # m to mm
			}
		climatic_variables[i,"precipitation"] = sum(vS)
		climatic_variables_xml[2] = paste0(climatic_variables_xml[2],log(sum(vS)*4)," ")
	}
write.csv(climatic_variables, "GLM_predictors_data/Prepared_variables/EBDS_skygrid_analyses/Climatic_variables_skygrid.csv", row.names=F, quote=F)
write(climatic_variables_xml, "GLM_predictors_data/Prepared_variables/EBDS_skygrid_analyses/Climatic_variables_skygrid.txt")
climatic_variables_monthly = read.csv("GLM_predictors_data/Prepared_variables/EBDS_skygrid_analyses/Climatic_variables_monthly.csv", head=T)
# plot(climatic_variables_monthly[,c("date","temperature")], type="l"); lines(climatic_variables[,c("date","temperature")], col="red")
# plot(climatic_variables_monthly[,c("date","precipitation")], type="l"); lines(climatic_variables[,"date"],climatic_variables[,"precipitation"]*4, col="red")

		# 5.2.2. Prepration of the climate data for the EBDS-GLM analyses

epoch_width = 1/52; climatic_variables = matrix(nrow=60, ncol=3)
colnames(climatic_variables) = c("date","temperature","precipitation")
for (i in 1:60)
	{
		date1 = mostRecentSamplingDatum-(i*epoch_width); date2 = date1+epoch_width
		climatic_variables[i,"date"] = (date1+date2)/2
		indices = which((decimal_date(dates)>date1)&(decimal_date(dates)<=date2))
		vS = rep(NA, length(indices))
		for (j in 1:length(vS))
			{
				vS[j] = mean(temperatures[[indices[j]]][], na.rm=T)-273.15 # Kelvins to °C
			}
		climatic_variables[i,"temperature"] = mean(vS)
		for (j in 1:length(vS))
			{
				vS[j] = mean(precipitations[[indices[j]]][], na.rm=T)*1000 # m to mm
			}
		climatic_variables[i,"precipitation"] = sum(vS)
	}
write.csv(climatic_variables, "GLM_predictors_data/Prepared_variables/EBDS_skygrid_analyses/Climatic_variables_EBDS.csv", row.names=F, quote=F)
climatic_variables_monthly = read.csv("GLM_predictors_data/Prepared_variables/EBDS_skygrid_analyses/Climatic_variables_monthly.csv", head=T)
# plot(climatic_variables_monthly[,c("date","temperature")], type="l"); lines(climatic_variables[,c("date","temperature")], col="red")
# plot(climatic_variables_monthly[,c("date","precipitation")], type="l"); lines(climatic_variables[,"date"],climatic_variables[,"precipitation"]*4, col="red")

# 6. Evolution of Ne and R(t): skygrid and sampling aware analyses

		# For the generation (or serial) time distribution, the mean was drawn from a uniform distribution ranging from 8 to 23 days and the standard
		# deviation from a uniform distribution ranging from 4 to 8 days, reflecting a range of values reported and considered in the literature:
		#     - Cauchemez et al. (2014, Eurosurveillance): 23 days (SD ~6 days)
		#     - Riou et al. (2017, Epidemics): 1.5-2.7 weeks (≈10.5–18.9 days)
		#     - Meyer et al. (2023, Epidemics): 13.8 for Italy, 12.2 for Cambodia, and 10.3 days for Bangladesh
		#     - Meyer et al. (2025, Sci. Adv.): median of 10 days, as well as 9 and 12 days for the warmest and collest populations, respectively

	# 6.1. Estimation of R0 based on the exponential growth rate (Grassly & Fraser, 2008)

n = 1000; mean_range = c(9, 23); sd_range = c(4, 8) # n = number of iterations, and ranges are in days
log = read.table(paste0("BEAST_DTA_analysis/Without_DTA_model/Alignment_101025_exp.log"), head=T)
rS1 = log[(((dim(log)[1]-1)/10)+2):dim(log)[1],"exponential.growthRate"]/365.25 # exponential growth rate (per day)
rS2 = sample(rS1, 1000, replace=F); R0s_list = list(); Rfs_list = list()
R0_lower_95ci = rep(NA, length(rS2)); R0_upper_95ci = rep(NA, length(rS2))
Rf_lower_95ci = rep(NA, length(rS2)); Rf_upper_95ci = rep(NA, length(rS2))
pC = 0; pI = 0.66 # as of early November 2025 (see SPF report of 03/11/25)
pC = 0; pI = (0.66-0.20)/(1.00-0.20) # correction for the study revision
for (i in 1:length(rS2))
	{
		R0s = rep(NA, n); Rfs = rep(NA, n); r = rS2[i]
		for (j in 1:n)
			{
				mean_si_i = runif(1, mean_range[1], mean_range[2]); sd_si_i = runif(1, sd_range[1], sd_range[2])
				a = (mean_si_i^2)/(sd_si_i^2); b = mean_si_i/(sd_si_i^2); R0s[j] = (1+(r/b))^a # https://beast.community/estimating_R0
				Rfs[j] = R0s[j]*(1-pC)*(1-pI) # as summarised in Fontanet & Cauchemez (2020; with pC set to 0)
			}
		R0s_list[[i]] = R0s; R0_lower_95ci[i] = quantile(R0s, prob=0.025); R0_upper_95ci[i] = quantile(R0s, prob=0.975)
		Rfs_list[[i]] = Rfs; Rf_lower_95ci[i] = quantile(Rfs, prob=0.025); Rf_upper_95ci[i] = quantile(Rfs, prob=0.975)
	}
R0_lower_95ci_median = median(R0_lower_95ci); R0_lower_95ci_95hpd = hdi(R0_lower_95ci)[1:2]
R0_upper_95ci_median = median(R0_upper_95ci); R0_upper_95ci_95hpd = hdi(R0_upper_95ci)[1:2]
	# --> R0 estimates range from 1.39 (95% HPD = [1.36, 1.41]) to 2.27 (95% HPD = [2.18, 2.38])
Rf_lower_95ci_median = median(Rf_lower_95ci); Rf_lower_95ci_95hpd = hdi(Rf_lower_95ci)[1:2]
Rf_upper_95ci_median = median(Rf_upper_95ci); Rf_upper_95ci_95hpd = hdi(Rf_upper_95ci)[1:2]
	# --> Rf estimates range from 0.47 (95% HPD = [0.46, 0.48]) to 0.77 (95% HPD = [0.74, 0.81]) if pI = 0.66
	# --> Rf estimates range from 0.59 (95% HPD = [0.58, 0.60]) to 0.97 (95% HPD = [0.93, 1.01]) if pI = 0.46/0.80

dev.new(width=8/2.8, height=8/3); par(oma=c(0,0,0,0), mar=c(2.8,3.0,1.0,1.0), lwd=0.3, bty="o", col="gray30", col.axis="gray30", fg="gray30")
plot(density(R0s_list[[1]]), col=NA, xlim=c(0.5,3.5), ylim=c(0,1.65), axes=F, ann=F)
for (i in 1:length(R0s_list)) lines(density(R0s_list[[i]]), lwd=0.1, col="gray70")
axis(side=1, lwd=0.5, cex.axis=0.7, mgp=c(0,0.17,0), lwd.tick=0.5, col.lab="gray30", col="gray30", tck=-0.03, las=1)
axis(side=2, lwd=0.5, cex.axis=0.7, mgp=c(0,0.4,-0.1), lwd.tick=0.5, col.lab="gray30", col="gray30", tck=-0.03, las=1, padj=0.4, at=seq(0,2,0.5))
mtext("Density       ", side=2, col="gray30", cex=0.8, line=1.7, las=3)
mtext("R0", side=1, col="gray30", cex=0.8, line=1.2)

	# 6.2. Setting the parameters for an episodic birth-death-sampling (EBDS) analysis

		# Estimation of the transmission period used to anchor the death rate: 8-15 days = 5-8 days of human infectiousness (Cauchemez et al. 2014;
		# Eurosurveillance, SI data) + 3-7 days to take into account the mosquito extrinsic incubation period (EIP; Zhao et al. 2025, Biosaf. Health).
		# Note on the distinction between the death (mu) and sampling (psi) rates: mu and psi are distinct hazards running in parallel, with mu removing
		# infectious lineages without observation (recovery, death, isolation), whereas psi governs observation events that create tips in the tree
		# (and, under the usual “sampling with removal” setup, also terminate that lineage). With two competing exponential clocks, the chance a lineage
		# ends by being sampled rather than silently terminating is psi relative to the total exit pressure mu+psi, hence rho being equal to psi/(mu+psi).

root_age = 1; num_taxa = 3114 # root age estimate and number of taxa
min_duration_days = 8 # lower bound of infectious duration (days)
max_duration_days = 15 # upper bound of infectious duration (days)
sampling_fraction_lower = 0.015 # lower bound of sampling fraction at current time 
sampling_fraction_higher = 0.015 # upper bound of sampling fraction at current time
birth_rate_logmean = log((log(num_taxa+2)-log(2))/root_age) # birth rate (lambda; Magee et al. 2020, PLoS Comp. Biol.)
birth_rate_logsd = 2*0.587405 # to get a weakly informative prior, spanning roughly 0.1 times to 10 times around the center:
							  # 		Pr((lambda_mean/10) <= lambda < (10*lambda_mean))
compute_death_rate_log_params = function(min_duration_days, max_duration_days)
	{
		min_duration_years = min_duration_days/365.25; max_death_rate = 1/min_duration_years
		max_duration_years = max_duration_days/365.25; min_death_rate = 1/max_duration_years
		log_mean = (log(min_death_rate)+log(max_death_rate))/2
		log_sd = (log(max_death_rate)-log(min_death_rate))/(2*1.96)
		 	# to fit a lognormal prior whose central 95% covers the lower and upper bounds
		return(list(death_meanlog=log_mean, death_sdlog=log_sd))
	}
death_rate_log_params = compute_death_rate_log_params(min_duration_days, max_duration_days)
death_rate_logmean = as.numeric(death_rate_log_params[1]) # death_rate (mu) = 1/infectious_duration
death_rate_logsd = as.numeric(death_rate_log_params[2])
compute_log_sampling_rate_params = function(death_rate, lower_bound, upper_bound) # sampling rate (psi)
	{
		sampling_rate_lower = death_rate*(lower_bound/(1-lower_bound)) # because sampling proportion = psi/(mu+psi)
		sampling_rate_upper = death_rate*(upper_bound/(1-upper_bound)) # because sampling proportion = psi/(mu+psi)
		log_mean = (log(sampling_rate_lower)+log(sampling_rate_upper))/2; log_sd = 0.587405
		return(list(sampling_meanlog=log_mean, sampling_sdlog=log_sd))
	}
death_rate = exp(death_rate_logmean)
sampling_rate_log_params = compute_log_sampling_rate_params(death_rate, sampling_fraction_lower, sampling_fraction_higher)
sampling_rate_logmean = as.numeric(sampling_rate_log_params[1])
sampling_rate_logsd = as.numeric(sampling_rate_log_params[2])
print(c(birth_rate_logmean, birth_rate_logsd)) # 1.994858 and 1.17481
print(c(death_rate_logmean, death_rate_logsd)) # 3.506836 and 0.1603594
print(c(sampling_rate_logmean, sampling_rate_logsd)) # -0.6777553 and 0.587405
	# To check the compound parameters:
lambda = rlnorm(100000, meanlog=birth_rate_logmean, sdlog=birth_rate_logsd) # birth_rate
mu = rlnorm(100000, meanlog=death_rate_logmean, sdlog=death_rate_logsd) # death_rate
psi = rlnorm(10000,  meanlog=sampling_rate_logmean, sdlog=sampling_rate_logsd) # sampling_rate
Re = lambda/(mu+psi); del = mu+psi; rho = psi/(mu+psi)
if (showingPlots)
	{
		dev.new(width=6, height=4); par(mfrow=c(1,3))
		hist(Re, freq=F, breaks=1000000, xlim=c(0,10), main="Re")
		hist(del, freq=F, breaks=10000, xlim=c(0,100), main="Become uninfectious rate")
		hist(rho, freq=F, breaks=10000, xlim=c(0,1), ylim=c(0,20), main="Sampling proportion")
	}
quantile(Re, probs=c(0.025,0.5,0.975))
quantile(del, probs=c(0.025,0.5,0.975))
quantile(rho, probs=c(0.025,0.5,0.975))

	# 6.3. Visualisation of the evolution through time of the number of cases, Ne, and R(t)

climatic_variables = read.csv("GLM_predictors_data/Prepared_variables/EBDS_skygrid_analyses/Climatic_variables_monthly.csv", head=T)
tab = read.table(paste0("BEAST_DTA_analysis/Alignment_",analysis,".txt"), head=T, sep="\t")
tab_weeks = interval(min(ymd(tab[,"collection_date"])),ymd(tab[,"collection_date"]))%/%weeks(1)+1
skg = read.csv(paste0("BEAST_DTA_analysis/Without_DTA_model/Alignment_",analysis,"_skg.csv"), head=T)[1:51,]
skg = skg[,c("time","median","lower","upper")]; colnames(skg) = c("time","median","95pHDP_lower","95pHDP_upper")
skg[,"median"] = log(skg[,"median"]+1); timeSlice = skg[1,"time"]-skg[2,"time"]; skg_weeks = 51:1
skg[,"95pHDP_lower"] = log(skg[,"95pHDP_lower"]+1); skg[,"95pHDP_upper"] = log(skg[,"95pHDP_upper"]+1)
skg_sa = read.csv(paste0("BEAST_DTA_analysis/Without_DTA_model/Alignment_",analysis,"_sa.csv"), head=T)[1:51,]
skg_sa = skg_sa[,c("time","median","lower","upper")]; colnames(skg_sa) = c("time","median","95pHDP_lower","95pHDP_upper")
skg_sa[,"median"] = log(skg_sa[,"median"]+1); timeSlice = skg_sa[1,"time"]-skg_sa[2,"time"]; skg_weeks = 51:1
skg_sa[,"95pHDP_lower"] = log(skg_sa[,"95pHDP_lower"]+1); skg_sa[,"95pHDP_upper"] = log(skg_sa[,"95pHDP_upper"]+1)
ebdsRt = read.table(paste0("BEAST_DTA_analysis/Without_DTA_model/Alignment_",analysis,"_ebds.log"), head=T, sep="\t")
cutOff = 1; nberOfPoints = sum(grepl("effectiveReproductiveNumber",colnames(ebdsRt)))
Rts = matrix(nrow=nberOfPoints, ncol=4); colnames(Rts) = c("time","median","95pHDP_lower","95pHDP_upper")
minYear = mostRecentSamplingDatum-1; maxYear = mostRecentSamplingDatum; interval = 1/nberOfPoints
Rt_times = seq(minYear+(interval/2), maxYear-(interval/2), interval); Rts[,"time"] = rev(Rt_times)
for (i in 1:dim(Rts)[1])
	{
		vS = ebdsRt[,paste0("effectiveReproductiveNumber",i)]; hpd95 = HDInterval::hdi(vS)
		Rts[i,"median"] = median(vS); Rts[i,"95pHDP_lower"] = hpd95[1]; Rts[i,"95pHDP_upper"] = hpd95[2]
	}
Rts = Rts[dim(Rts)[1]:1,]

pdf(paste0("Figure_1B_",analysis,"_NEW.pdf"), width=8/2, height=8/2) # dev.new(width=8/2, height=8/2)
par(mfrow=c(3,1), oma=c(0,0,0,0), mar=c(2.0,3.5,0.1,0.5), lwd=0.3, bty="o", col="gray30", col.axis="gray30", fg="gray30")
hist(decimal_date(ymd(tab[,"collection_date"])), breaks=52, col="gray70", border=NA, axes=F, ann=F, ylim=c(0,380), xlim=c(2024.62,2025.59))
hist(decimal_date(ymd(tab[,"collection_date"])), breaks=52, col=paste0(collection_week_cols1,"BF"), add=T)
dates = c("2024-08-07","2024-09-01","2024-10-01","2024-11-01","2024-12-01","2025-01-01","2025-02-01","2025-03-01","2025-04-01","2025-05-01","2025-06-01","2025-07-01","2025-08-08")
ats1 = decimal_date(ymd(dates)); ats2 = ats1[2:(length(ats1)-1)]; labels1 = gsub("-","\\/",dates[2:(length(dates)-1)])
axis(side=1, lwd=0.5, cex.axis=0.7, mgp=c(0,0.17,0), lwd.tick=0, col.lab="gray30", col="gray30", tck=-0.04, las=1, at=ats1, label=rep("",length(ats1)))
axis(side=1, lwd=0, cex.axis=0.7, mgp=c(0,0.17,0), lwd.tick=0.5, col.lab="gray30", col="gray30", tck=-0.04, las=1, at=ats2, label=labels1)
axis(side=2, lwd=0.5, cex.axis=0.7, mgp=c(0,0.4,-0.1), lwd.tick=0.5, col.lab="gray30", col="gray30", tck=-0.04, las=1, padj=0.4, at=seq(0,300,100))
mtext("Genomic samples", side=2, col="gray30", cex=0.6, line=1.7, las=3)
plot(skg_sa[,"time"], skg_sa[,"median"], lwd=0.7, type="l", cex.axis=0.8, cex.lab=0.8, col="gray30", axes=F, xlab=NA, ylab=NA, ylim=c(0,7), xlim=c(2024.62,2025.59))
xx_l = c(skg_sa[,c("time")],rev(skg_sa[,c("time")])); yy_l = c(skg_sa[,"95pHDP_lower"],rev(skg_sa[,"95pHDP_upper"]))
getOption("scipen"); opt = options("scipen"=20); polygon(xx_l,yy_l,col=rgb(187/255,187/255,187/255,0.25),border=0)
for (i in 1:length(skg_sa[,"time"]))
	{
		colour = paste0(collection_week_cols1[skg_weeks[i]],"BF")
		x1 = skg_sa[i,"time"]-(timeSlice/2); x2 = skg_sa[i,"time"]+(timeSlice/2)
		y1 = skg_sa[i,"95pHDP_lower"]-1; y2 = skg_sa[i,"95pHDP_upper"]+1
		polygon(c(x1,x2,x2,x1), c(y1,y1,y2,y2), col="gray70", border=NA)
		polygon(c(x1,x2,x2,x1), c(y1,y1,y2,y2), col=colour, border=NA)
	}
getOption("scipen"); opt = options("scipen"=20); polygon(xx_l,yy_l,col=NA,border="gray30")
lines(skg_sa[,"time"], skg_sa[,"median"], lwd=0.3, type="l", cex.axis=0.8, cex.lab=0.8, col="gray30")
lines(skg[,"time"], skg[,"median"], lwd=0.3, type="l", cex.axis=0.8, cex.lab=0.8, col="gray30", lty=2)
lines(skg[,"time"], skg[,"95pHDP_lower"], lwd=0.1, type="l", cex.axis=0.8, cex.lab=0.8, col="gray30", lty=2)
lines(skg[,"time"], skg[,"95pHDP_upper"], lwd=0.1, type="l", cex.axis=0.8, cex.lab=0.8, col="gray30", lty=2)
labels2 = c(1,5,10,20,50,100,200,500,2000); ats3 = c(log(labels2))
axis(side=1, lwd=0.5, cex.axis=0.7, mgp=c(0,0.17,0), lwd.tick=0, col.lab="gray30", col="gray30", tck=-0.04, las=1, at=ats1, label=rep("",length(ats1)))
axis(side=1, lwd=0, cex.axis=0.7, mgp=c(0,0.17,0), lwd.tick=0.5, col.lab="gray30", col="gray30", tck=-0.04, las=1, at=ats2, label=labels1)
axis(side=2, lwd=0.5, cex.axis=0.7, mgp=c(0,0.4,-0.1), lwd.tick=0.5, col.lab="gray30", col="gray30", tck=-0.04, las=1, padj=0.4, at=ats3, label=labels2)
mtext("Effective population size", side=2, col="gray30", cex=0.55, line=1.7, las=3)
plot(Rts[,"time"], Rts[,"median"], col=NA, axes=F, xlab=NA, ylab=NA, ylim=c(0,3.5), xlim=c(2024.62,2025.59)); timeSlice = 1/52
for (i in 1:dim(Rts)[1])
	{
		colour = paste0(collection_week_cols2[i],"BF")
		x1 = Rts[i,"time"]-(timeSlice/2); x2 = Rts[i,"time"]+(timeSlice/2)
		y1 = Rts[i,"95pHDP_lower"]; y2 = Rts[i,"95pHDP_upper"]
		polygon(c(x1,x2,x2,x1), c(y1,y1,y2,y2), col="gray70", border=NA)
		polygon(c(x1,x2,x2,x1), c(y1,y1,y2,y2), col=colour, border="gray30")
		x1 = Rts[i,"time"]-(timeSlice/2); x2 = Rts[i,"time"]+(timeSlice/2)
		y1 = Rts[i,"median"]; y2 = Rts[i,"median"]
		polygon(c(x1,x2,x2,x1), c(y1,y1,y2,y2), col="gray30", border=1)
	}
axis(side=1, lwd=0.5, cex.axis=0.7, mgp=c(0,0.17,0), lwd.tick=0, col.lab="gray30", col="gray30", tck=-0.04, las=1, at=ats1, label=rep("",length(ats1)))
axis(side=1, lwd=0, cex.axis=0.7, mgp=c(0,0.17,0), lwd.tick=0.5, col.lab="gray30", col="gray30", tck=-0.04, las=1, at=ats2, label=labels1)
axis(side=2, lwd=0.5, cex.axis=0.7, mgp=c(0,0.4,-0.1), lwd.tick=0.5, col.lab="gray30", col="gray30", tck=-0.04, las=1, padj=0.4)
mtext(expression("Reproduction number R"[t]), side=2, col="gray30", cex=0.55, line=1.7, las=3)
dev.off()

	# 6.4. Visualisation and summary of the skygrid- and EBDS-GLM analyses with covariates

analyses = rep(NA, 6); burnIns = rep(NA, 6); skygrids = list()
analyses[1] = "Alignment_101025_sa_covs1/Alignment_101025_sa_cov3_lag0a"; burnIns[1] = 100000000 # univariate analysis that only includes the log-transformed temperature variable
analyses[2] = "Alignment_101025_sa_covs1/Alignment_101025_sa_cov4_lag0a"; burnIns[2] = 100000000 # univariate analysis that only includes the log-transformed precipitation variable
analyses[3] = "Alignment_101025_sa_covs1/Alignment_101025_sa_cov2_lag0a"; burnIns[3] = 200000000 # multivariate analysis while assuming that the effective population size is associated 
											   						# with one of the two covariates, using indicator variables to select between the two covariates (not used anymore)
analyses[3] = "Alignment_101025_sa_covs2/Alignment_101025_sa_cov1_lag0"; burnIns[3] = 0 # (combined) multivariate analysis without assuming that the effective population size is associated 
											   						# with one of the two covariates, using indicator variables to select between the two covariates
analyses[4] = "Alignment_101025_sa_covs1/Alignment_101025_sa_cov3_lag1a"; burnIns[4] = 170000000 # same as "Alignment_101025_sa_cov3_lag0" but condering a 1-month lag with the covariates
analyses[5] = "Alignment_101025_sa_covs1/Alignment_101025_sa_cov4_lag1a"; burnIns[5] = 100000000 # same as "Alignment_101025_sa_cov4_lag0" but condering a 1-month lag with the covariates
analyses[6] = "Alignment_101025_sa_covs1/Alignment_101025_sa_cov2_lag1a"; burnIns[6] = 100000000 # same as "Alignment_101025_sa_cov2_lag0" but condering a 1-month lag with the covariates
analyses[6] = "Alignment_101025_sa_covs2/Alignment_101025_sa_cov1_lag1"; burnIns[6] = 0 # (combined) same as "Alignment_101025_sa_cov1_lag0" but condering a 1-month lag with the covariates
epoch_widths = c(0.09863014,0.08219178,0.08493151,0.08219178,0.08493151,0.07671233,0.08493151,0.08469945,0.08196721,0.08469945,0.08196721)
for (i in 1:length(analyses))
	{
		tab = data.frame(matrix(nrow=12, ncol=6)); colnames(tab) = c("year","month","date","median_Ne","lower_Ne","upper_Ne")
		log = read.table(paste0("BEAST_DTA_analysis/Without_DTA_model/Alignment_101025_sa_covs/",analyses[i],".log"), head=T)
		if (burnIns[i] != 0) log = log[(which(log[,"state"]==burnIns[i])+1):dim(log)[1],]
		for (j in 1:dim(tab)[1])
			{
				tab[j,"date"] = mostRecentSamplingDatum-sum(epoch_widths[1:j])+(15/365)
				if (j == dim(tab)[1]) tab[j,"date"] = decimal_date(ymd("2024-08-15"))
				date = date_decimal(tab[j,"date"]); tab[j,"year"] = format(date,"%Y")
				tab[j,"month"] = as.numeric(format(date,"%m"))
				vS = log[,paste0("effPop",j)]
				tab[j,"median_Ne"] = log(median(vS)+1); hpd95 = log(HDInterval::hdi(vS)+1)
				tab[j,"lower_Ne"] = hpd95["lower"]; tab[j,"upper_Ne"] = hpd95["upper"]
			}
		tab = tab[order(as.numeric(tab[,"date"])),]; skygrids[[i]] = tab
	}

analyses = rep(NA, 4); burnIns = rep(NA, 4); EBDSs = list()
analyses[1] = "Alignment_101025_ebds_GP_a_lag0"; burnIns[1] = 0 # analysis with a pure log-linear regression baseline
analyses[2] = "Alignment_101025_ebds_GP_b_lag0"; burnIns[2] = 0 # analysis with a a proper squared-exponential kernel
analyses[3] = "Alignment_101025_ebds_GP_a_lag1"; burnIns[3] = 0 # same as the 1° but condering a 1-month lag with the covariates
analyses[4] = "Alignment_101025_ebds_GP_b_lag1"; burnIns[4] = 0 # same as the 2° but condering a 1-month lag with the covariates
epoch_widths = c(0.09863014,0.08219178,0.08493151,0.08219178,0.08493151,0.07671233,0.08493151,0.08469945,0.08196721,0.08469945,0.08196721)
for (i in 1:length(analyses))
	{
		tab = data.frame(matrix(nrow=12, ncol=6)); colnames(tab) = c("year","month","date","median_Rt","lower_Rt","upper_Rt")
		log = read.table(paste0("BEAST_DTA_analysis/Without_DTA_model/Alignment_101025_ebds_GP/",analyses[i],".log"), head=T)
		if (burnIns[i] != 0) log = log[(which(log[,"state"]==burnIns[i])+1):dim(log)[1],]
		for (j in 1:dim(tab)[1])
			{
				tab[j,"date"] = mostRecentSamplingDatum-sum(epoch_widths[1:j])+(15/365)
				if (j == dim(tab)[1]) tab[j,"date"] = decimal_date(ymd("2024-08-15"))
				date = date_decimal(tab[j,"date"]); tab[j,"year"] = format(date,"%Y")
				tab[j,"month"] = as.numeric(format(date,"%m"))
				vS = log[,paste0("effectiveReproductiveNumber",j)]
				tab[j,"median_Rt"] = median(vS); hpd95 = HDInterval::hdi(vS)
				tab[j,"lower_Rt"] = hpd95["lower"]; tab[j,"upper_Rt"] = hpd95["upper"]
			}
		tab = tab[order(as.numeric(tab[,"date"])),]; EBDSs[[i]] = tab
	}

dev.new(width=8, height=4) # pdf(paste0("Figure_SX_",analysis,"_NEW.pdf"), width=8, height=4)
par(mfrow=c(2,3), oma=c(0.2,0.4,1.0,0.5), mar=c(2.5,3.5,0.3,0.5), lwd=0.3, bty="o", col="gray30", col.axis="gray30", fg="gray30")
labels = c("A","B","C","D","E","F"); months = c("J","J","M","A","M","F","J","D","N","O","S","A")
for (i in 1:length(skygrids))
	{
		plot(skygrids[[i]][,"date"], skygrids[[i]][,"median_Ne"], col=NA, axes=F, xlab=NA, ylab=NA, xlim=c(2024.56,2025.62), ylim=c(0,max(skygrids[[i]][,"upper_Ne"])))
		timeSlice = 1/13
		for (j in dim(skygrids[[i]])[1]:1)
			{
				colour = paste0(collection_month_cols[j],"BF")
				x1 = skygrids[[i]][j,"date"]-(timeSlice/2); x2 = skygrids[[i]][j,"date"]+(timeSlice/2)
				y1 = skygrids[[i]][j,"lower_Ne"]; y2 = skygrids[[i]][j,"upper_Ne"]
				polygon(c(x1,x2,x2,x1), c(y1,y1,y2,y2), col="gray70", border=NA)
				polygon(c(x1,x2,x2,x1), c(y1,y1,y2,y2), col=colour, border="gray30")
				x1 = skygrids[[i]][j,"date"]-(timeSlice/2); x2 = skygrids[[i]][j,"date"]+(timeSlice/2)
				y1 = skygrids[[i]][j,"median_Ne"]; y2 = skygrids[[i]][j,"median_Ne"]
				polygon(c(x1,x2,x2,x1), c(y1,y1,y2,y2), col="gray30", border=1)
			}
		axis(side=1, lwd=0.5, cex.axis=0.85, mgp=c(0,0.40,0), lwd.tick=0.5, col.lab="gray30", col="gray30", tck=-0.03, las=1, at=skygrids[[i]][,"date"], label=rev(months))
		axis(side=2, lwd=0.5, cex.axis=0.85, mgp=c(0,0.65,0), lwd.tick=0.5, col.lab="gray30", col="gray30", tck=-0.03, las=1, padj=0.4)
		title(ylab="Effective population size", mgp=c(1.7,0,0), cex.lab=1.1, col.lab="gray30")
		mtext(expression(bold(labels[i])), at=-2450, line=-0.95, cex.lab=0.80, col="gray30")
	}

climatic_variables = read.csv("GLM_predictors_data/Prepared_variables/EBDS_skygrid_analyses/Climatic_variables_monthly.csv", head=T)
climatic_variable_names = c("temperature","precipitation"); climatic_variable_titles = c("Temperature (°C)","Precipitation (mm)")
skygrid_tabs = list()
for (i in 1:8)
	{
		tab = data.frame(matrix(nrow=12, ncol=6))
		if (i == 1) { skygrid = skygrids[[1]]; h = 1; lag = 0 } # skygrid with the temperature covariate, 0-month lag
		if (i == 2) { skygrid = skygrids[[2]]; h = 2; lag = 0 } # skygrid with the precipitation covariate, 0-month lag
		if (i == 3) { skygrid = skygrids[[3]]; h = 1; lag = 0 } # skygrid with both covariates and indicators, 0-month lag
		if (i == 4) { skygrid = skygrids[[3]]; h = 2; lag = 0 } # skygrid with both covariates and indicators, 0-month lag
		if (i == 5) { skygrid = skygrids[[4]]; h = 1; lag = 1 } # skygrid with the temperature covariate, 1-month lag
		if (i == 6) { skygrid = skygrids[[5]]; h = 2; lag = 1 } # skygrid with the precipitation covariate, 1-month lag
		if (i == 7) { skygrid = skygrids[[6]]; h = 1; lag = 1 } # skygrid with both covariates and indicators, 1-month lag
		if (i == 8) { skygrid = skygrids[[6]]; h = 2; lag = 1 } # skygrid with both covariates and indicators, 1-month lag
		colnames(tab) = c("year","month",climatic_variable_names[h],"median_Ne","lower_Ne","upper_Ne")
		tab[,"year"] = as.numeric(skygrid[,"year"]); tab[,"month"] = skygrid[,"month"]; indices = rep(NA, 12)
		for (j in 1:dim(tab)[1])
			{
				indices[j] = which((climatic_variables[,"year"]==tab[j,"year"])&(climatic_variables[,"month"]==tab[j,"month"]))-lag
			}
		tab[,climatic_variable_names[h]] = climatic_variables[indices,climatic_variable_names[h]]; tab[,c("median_Ne")] = skygrid[,c("median_Ne")]
		tab[,c("lower_Ne")] = skygrid[,c("lower_Ne")]; tab[,c("upper_Ne")] = skygrid[,c("upper_Ne")]; skygrid_tabs[[i]] = tab
	}
EBDS_tabs = list()
for (i in 1:8)
	{
		tab = data.frame(matrix(nrow=12, ncol=6))
		if (i == 1) { EBDS = EBDSs[[1]]; h = 1; lag = 0 } # analysis with a pure log-linear regression baseline, 0-month lag
		if (i == 2) { EBDS = EBDSs[[1]]; h = 2; lag = 0 } # analysis with a pure log-linear regression baseline, 0-month lag
		if (i == 3) { EBDS = EBDSs[[2]]; h = 1; lag = 0 } # analysis with a proper squared-exponential kernel, 0-month lag
		if (i == 4) { EBDS = EBDSs[[2]]; h = 2; lag = 0 } # analysis with a proper squared-exponential kernel, 0-month lag
		if (i == 5) { EBDS = EBDSs[[3]]; h = 1; lag = 1 } # analysis with a pure log-linear regression baseline, 1-month lag
		if (i == 6) { EBDS = EBDSs[[3]]; h = 2; lag = 1 } # analysis with a pure log-linear regression baseline, 1-month lag
		if (i == 7) { EBDS = EBDSs[[4]]; h = 1; lag = 1 } # analysis with a proper squared-exponential kernel, 1-month lag
		if (i == 8) { EBDS = EBDSs[[4]]; h = 2; lag = 1 } # analysis with a proper squared-exponential kernel, 1-month lag
		colnames(tab) = c("year","month",climatic_variable_names[h],"median_Rt","lower_Rt","upper_Rt")
		tab[,"year"] = as.numeric(EBDS[,"year"]); tab[,"month"] = EBDS[,"month"]; indices = rep(NA, 12)
		for (j in 1:dim(tab)[1])
			{
				indices[j] = which((climatic_variables[,"year"]==tab[j,"year"])&(climatic_variables[,"month"]==tab[j,"month"]))-lag
			}
		tab[,climatic_variable_names[h]] = climatic_variables[indices,climatic_variable_names[h]]; tab[,c("median_Rt")] = EBDS[,c("median_Rt")]
		tab[,c("lower_Rt")] = EBDS[,c("lower_Rt")]; tab[,c("upper_Rt")] = EBDS[,c("upper_Rt")]; EBDS_tabs[[i]] = tab
	}

dev.new(width=8, height=7.1) # pdf(paste0("Figure_2_",analysis,"_NEW.pdf"), width=8, height=7.1)
par(mfrow=c(4,4), oma=c(0.2,0.3,1.0,0.5), mar=c(3.2,3.2,0.3,0.7), lwd=0.5, bty="o", col="gray30", col.axis="gray30", fg="gray30")
for (i in 1:length(skygrid_tabs))
	{
		if ((i %% 2) == 0)
			{
				h = 2
			}	else	{
				h = 1
			}
		maxY = max(skygrid_tabs[[i]][,"upper_Ne"])+(0.05*max(skygrid_tabs[[i]][,"upper_Ne"]))
		minX = min(skygrid_tabs[[i]][,climatic_variable_names[h]])-(0.05*((max(skygrid_tabs[[i]][,climatic_variable_names[h]])-min(skygrid_tabs[[i]][,climatic_variable_names[h]]))))
		maxX = max(skygrid_tabs[[i]][,climatic_variable_names[h]])+(0.05*((max(skygrid_tabs[[i]][,climatic_variable_names[h]])-min(skygrid_tabs[[i]][,climatic_variable_names[h]]))))
		plot(skygrid_tabs[[i]][,climatic_variable_names[h]], skygrid_tabs[[i]][,"median_Ne"], col="gray30", type="l", lwd=0.5, axes=F, xlab=NA, ylab=NA, xlim=c(minX,maxX), ylim=c(0,maxY))
		for (j in 1:dim(tab)[1])
			{
				plotrix::plotCI(skygrid_tabs[[i]][j,climatic_variable_names[h]], skygrid_tabs[[i]][j,"median_Ne"],
							    li=skygrid_tabs[[i]][j,"lower_Ne"], ui=skygrid_tabs[[i]][j,"upper_Ne"], add=T, lwd=0.5, col="gray30", pt.bg=NA, pch=16, cex=0.2)
				points(skygrid_tabs[[i]][j,climatic_variable_names[h]], skygrid_tabs[[i]][j,"median_Ne"], pch=16, cex=1.3, col="white")
				points(skygrid_tabs[[i]][j,climatic_variable_names[h]], skygrid_tabs[[i]][j,"median_Ne"], pch=16, cex=1.3, col=paste0(collection_month_cols[j],"BF"))
				points(skygrid_tabs[[i]][j,climatic_variable_names[h]], skygrid_tabs[[i]][j,"median_Ne"], pch=1, cex=1.3, col="gray30", lwd=0.5)
			}
		box(lwd=0.5, col="gray30")
		axis(side=1, lwd=0, cex.axis=0.85, mgp=c(0,0.35,0), lwd.tick=0.5, col.lab="gray30", col="gray30", tck=-0.04, las=1)
		axis(side=2, lwd=0, cex.axis=0.85, mgp=c(0,0.65,0), lwd.tick=0.5, col.lab="gray30", col="gray30", tck=-0.04, las=1, padj=0.4)
		title(xlab=climatic_variable_titles[h], mgp=c(1.5,0,0), cex.lab=0.9, col.lab="gray30")
		title(ylab="Effective population size (log)", mgp=c(1.7,0,0), cex.lab=0.9, col.lab="gray30")
	}
for (i in 1:length(EBDS_tabs))
	{
		if ((i %% 2) == 0)
			{
				h = 2
			}	else	{
				h = 1
			}
		maxY = max(EBDS_tabs[[i]][,"upper_Rt"])+(0.05*max(EBDS_tabs[[i]][,"upper_Rt"]))
		minX = min(EBDS_tabs[[i]][,climatic_variable_names[h]])-(0.05*((max(EBDS_tabs[[i]][,climatic_variable_names[h]])-min(EBDS_tabs[[i]][,climatic_variable_names[h]]))))
		maxX = max(EBDS_tabs[[i]][,climatic_variable_names[h]])+(0.05*((max(EBDS_tabs[[i]][,climatic_variable_names[h]])-min(EBDS_tabs[[i]][,climatic_variable_names[h]]))))
		plot(EBDS_tabs[[i]][,climatic_variable_names[h]], EBDS_tabs[[i]][,"median_Rt"], col="gray30", type="l", lwd=0.5, axes=F, xlab=NA, ylab=NA, xlim=c(minX,maxX), ylim=c(0,maxY))
		for (j in 1:dim(tab)[1])
			{
				plotrix::plotCI(EBDS_tabs[[i]][j,climatic_variable_names[h]], EBDS_tabs[[i]][j,"median_Rt"],
							    li=EBDS_tabs[[i]][j,"lower_Rt"], ui=EBDS_tabs[[i]][j,"upper_Rt"], add=T, lwd=0.5, col="gray30", pt.bg=NA, pch=16, cex=0.2)
				points(EBDS_tabs[[i]][j,climatic_variable_names[h]], EBDS_tabs[[i]][j,"median_Rt"], pch=16, cex=1.3, col="white")
				points(EBDS_tabs[[i]][j,climatic_variable_names[h]], EBDS_tabs[[i]][j,"median_Rt"], pch=16, cex=1.3, col=paste0(collection_month_cols[j],"BF"))
				points(EBDS_tabs[[i]][j,climatic_variable_names[h]], EBDS_tabs[[i]][j,"median_Rt"], pch=1, cex=1.3, col="gray30", lwd=0.5)
			}
		box(lwd=0.5, col="gray30")
		axis(side=1, lwd=0, cex.axis=0.85, mgp=c(0,0.35,0), lwd.tick=0.5, col.lab="gray30", col="gray30", tck=-0.04, las=1)
		axis(side=2, lwd=0, cex.axis=0.85, mgp=c(0,0.65,0), lwd.tick=0.5, col.lab="gray30", col="gray30", tck=-0.04, las=1, padj=0.4)
		title(xlab=climatic_variable_titles[h], mgp=c(1.5,0,0), cex.lab=0.9, col.lab="gray30")
		title(ylab="Effective reproduction number", mgp=c(1.7,0,0), cex.lab=0.9, col.lab="gray30")	
	}

# 7. Discrete phylogeographic analysis based on the municipalities

nberOfTreesToSample = 100; nberOfExtractionFiles = nberOfTreesToSample

	# 7.1. Extracting the spatio-temporal information embedded in posterior trees

		# 7.1.1. Extracting the lineage transition events for the visualisations

trees = scan(paste0("BEAST_DTA_analysis/With_empirical_trees/Alignment_",analysis,"_emp1.trees"), what="", sep="\n", quiet=T, blank.lines.skip=F)
burnIn = 251; index1 = which(trees=="\t\t;")[length(which(trees=="\t\t;"))]; index2 = index1 + burnIn + 1
indices3 = which(grepl("tree STATE",trees)); index3 = indices3[length(indices3)]
interval = floor((index3-(index1+burnIn))/nberOfTreesToSample)
indices = seq(index3-((nberOfTreesToSample-1)*interval),index3,interval)
selected_trees = c(trees[c(1:index1,indices)],"End;")
write(selected_trees, paste0("BEAST_DTA_analysis/With_empirical_trees/Alignment_",analysis,"_",nberOfExtractionFiles,"_1.trees"))

source("treeExtractions_DTA.r")
trees = readAnnotatedNexus(paste0("BEAST_DTA_analysis/With_empirical_trees/Alignment_",analysis,"_",nberOfExtractionFiles,"_1.trees"))
for (i in 1:nberOfExtractionFiles)
	{
		if (length(trees) == 1) dta_tab = treeExtractions_DTA(trees, mostRecentSamplingDatum)	
		if (length(trees) > 1) dta_tab = treeExtractions_DTA(trees[[i]], mostRecentSamplingDatum)	
		write.csv(dta_tab, paste0("BEAST_DTA_analysis/With_empirical_trees/Alignment_",analysis,"_",nberOfTreesToSample,"_ext/TreeExtractions_",i,".csv"), row.names=F, quote=F)
	}
samplingCoordinates = read.csv(paste0("Alignment_",analysis,".csv"), head=T, sep=";")
for (i in 1:nberOfExtractionFiles)
	{
		tab1 = read.csv(paste0("BEAST_DTA_analysis/With_empirical_trees/Alignment_",analysis,"_",nberOfTreesToSample,"_ext/TreeExtractions_",i,".csv"), head=T)
		if (!"tipLabel"%in%colnames(tab1))
			{
				tipLabels = matrix(nrow=dim(tab1)[1], ncol=1); colnames(tipLabels) = "tipLabel"
				for (j in 1:length(trees[[i]]$tip.label))
					{
						tipLabels[which(tab1[,"node2"]==j),"tipLabel"] = gsub("'","",trees[[i]]$tip.label[j])
					}
				tab2 = cbind(tab1, tipLabels)
			}	else	{
				tab2 = tab1
			}
		if (!"startLon"%in%colnames(tab2))
			{
				endLonLat = matrix(nrow=dim(tab2)[1], ncol=2); colnames(endLonLat) = c("endLon","endLat")
				for (j in 1:dim(tab2)[1])
					{
						if (!tab2[j,"node2"]%in%tab2[,"node1"])
							{
								index = which(samplingCoordinates[,"trait"]==tab2[j,"tipLabel"])
								if (length(index) == 1)
									{
										endLonLat[j,"endLon"] = samplingCoordinates[index,"longitude"]
										endLonLat[j,"endLat"] = samplingCoordinates[index,"latitude"]
									}	else	{
										print(c(i,j))
									}
							}
					}
				tab3 = cbind(tab2, endLonLat)
				write.csv(tab3, paste0("BEAST_DTA_analysis/With_empirical_trees/Alignment_",analysis,"_",nberOfTreesToSample,"_ext/TreeExtractions_",i,".csv"), row.names=F, quote=F)
			}
	}

		# 7.1.2. Extracting the Markov jumps for the Shannon entropy estimates

trees = treeio::read.beast(paste0("BEAST_DTA_analysis/With_empirical_trees/Alignment_",analysis,"_emp2.trees")); tabs = list(); ii = 0
for (i in 251:length(trees)) # to discard the 250 first sampled trees as burn-in
	{
		tip_depths = node.depth.edgelength(as.phylo(trees[[i]])); max_depth = max(tip_depths)
		tab = as.data.frame(as_tibble(trees[[i]]))
		tab$height = max_depth - tip_depths
		history = tab[,c("node","history")]
		tab = tab[-which(tab[,"parent"]==tab[,"node"]),]
		colnames(tab) = c("node1","node2","length","tipLabel","history","endLoc","height2")
		tmp = matrix(nrow=dim(tab)[1], ncol=2)
		colnames(tmp) = c("height1","startLoc"); tab = cbind(tab,tmp)
		tab = tab[,c("node1","node2","height1","height2","length","startLoc","endLoc","tipLabel","history")]
		for (j in 1:dim(tab)[1])	
			{
				index = which(tab[,"node2"]==tab[j,"node1"])
				if (length(index) == 1)
					{
						tab[j,"height1"] = tab[index,"height2"]; tab[j,"startLoc"] = tab[index,"endLoc"]
					}
			}
		tmp = c(); n = 0
		for (j in 1:dim(tab)[1])
			{
				if (is.na(tab$history[[j]])[1])
					{
						tmp = rbind(tmp, tab[j,])
					}
				if (!is.na(tab$history[[j]])[1])
					{
						MJs = c()
						for (k in 1:length(tab$history[[j]]))
							{
								if (grepl("\\{",tab$history[[j]][k]))
									{
										MJs = rbind(MJs,cbind(tab$history[[j]][k+1],gsub("\\}","",tab$history[[j]][k+2]),gsub("\\{","",tab$history[[j]][k])))
									}
							}
						sub_branches = c()
						for (k in 1:dim(MJs)[1])
							{
								sub_branch = tab[j,]; n = n+1
								if (k == 1)
									{
										sub_branch[,"node2"] = paste0("i",n)
										sub_branch[,"height2"] = MJs[k,3]
										sub_branch[,"endLoc"] = MJs[k,2]	
									}
								if (k > 1)
									{
										sub_branch[,"node1"] = sub_branches[k-1,"node2"]
										sub_branch[,"height1"] = as.numeric(sub_branches[k-1,"height2"])
										sub_branch[,"startLoc"] = sub_branches[k-1,"endLoc"]
										sub_branch[,"node2"] = paste0("i",n)
										sub_branch[,"height2"] = MJs[k,3]
										sub_branch[,"endLoc"] = MJs[k,2]
									}
								sub_branch[,"length"] = sub_branch[,"height1"]-as.numeric(sub_branch[,"height2"])
								sub_branches = rbind(sub_branches, sub_branch)
							}
						sub_branch = tab[j,]
						sub_branch[,"node1"] = sub_branches[dim(MJs)[1],"node2"]
						sub_branch[,"height1"] = as.numeric(sub_branches[dim(MJs)[1],"height2"])
						sub_branch[,"startLoc"] = sub_branches[dim(MJs)[1],"endLoc"]
						sub_branch[,"length"] = sub_branch[,"height1"]-as.numeric(sub_branch[1,"height2"])
						sub_branches = rbind(sub_branches, sub_branch)
						tmp = rbind(tmp, sub_branches)
					}
			}
		ii = ii+1; tabs[[ii]] = tmp
	}
saveRDS(tabs, paste0("BEAST_DTA_analysis/With_empirical_trees/Alignment_",analysis,"_emp3.rds"))

		# 7.1.3. Averaging the Markov jumps for the overall and four time periods

if (!file.exists(paste0("BEAST_DTA_analysis/With_empirical_trees/Alignment_",analysis,"_",nberOfExtractionFiles,"_3.csv")))
	{
		tab = read.table(paste0("BEAST_DTA_analysis/With_empirical_trees/Alignment_",analysis,"_emp3.txt"), head=T)
		MJs = matrix(nrow=dim(admins2@data)[1], ncol=dim(admins2@data)[1]); states = unique(tab[,"state"])
		for (i in 1:dim(admins2@data)[1])
			{
				for (j in 1:dim(admins2@data)[1])
					{
						sub = tab[which((tab[,"from"]==admins2@data[i,"GID_2"])&(tab[,"to"]==admins2@data[j,"GID_2"])),]
						MJs_post = c(); # print(c(i,j))
						for (k in 1:length(states))
							{
								MJ = which(sub[,"state"]==states[k])
								MJs_post = c(MJs_post, length(MJ))
							}
						MJs[i,j] = mean(MJs_post)
					}
			}
		row.names(MJs) = admins2@data[,"GID_2"]; colnames(MJs) = admins2@data[,"GID_2"]
		write.table(MJs, paste0("BEAST_DTA_analysis/With_empirical_trees/Alignment_",analysis,"_emp3.csv"), quote=F, sep=",")
	}
MJs = read.csv(paste0("BEAST_DTA_analysis/With_empirical_trees/Alignment_",analysis,"_",nberOfExtractionFiles,"_3.csv"), header=T)
log = read.table(paste0("BEAST_DTA_analysis/With_empirical_trees/Alignment_",analysis,"_emp1.log"), header=T)
log = log[(which(log[,"state"]==25100000)+1):dim(log)[1],]
BFs = matrix(nrow=dim(admins2@data)[1], ncol=dim(admins2@data)[1])
row.names(BFs) = admins2@data[,"GID_2"]; colnames(BFs) = admins2@data[,"GID_2"]
for (i in 1:dim(admins2@data)[1])
	{
		for (j in 1:dim(admins2@data)[1])
			{
				if (i != j)
					{
						colName = paste0("location.indicators.",admins2@data[i,"GID_2"],".",admins2@data[j,"GID_2"])
						index1 = which(colnames(log)==colName)
						p = sum(log[,index1]==1)/dim(log)[1]
						K = dim(admins2@data)[1]
						q = (log(2)+K-1)/(K*(K-1))
						BFs[i,j] = (p/(1-p))/(q/(1-q))
					}
			}
	}
for (i in 1:dim(admins2@data)[1])
	{
		for (j in 1:dim(admins2@data)[1])
			{
				if ((MJs[i,j] > 1)&(BFs[i,j] < 3)) cat(admins2@data[i,"GID_2"]," - ",admins2@data[j,"GID_2"],": MJ = ",MJs[i,j],"\n",sep="") # none
			}
	}
MJs[BFs<3] = NA; BFs = round(BFs, 1)
cutOffs = c(decimal_date(ymd(c("2025-01-01","2025-03-01","2025-05-01","2025-09-01"))))
for (h in 1:length(cutOffs))
	{
		tab = read.table(paste0("BEAST_DTA_analysis/With_empirical_trees/Alignment_",analysis,"_emp3.txt"), head=T)
		tab[,"time"] = mostRecentSamplingDatum-tab[,"time"]
		if (h == 1) tab = tab[which(tab[,"time"]<cutOffs[h]),]
		if (h > 1) tab = tab[which((tab[,"time"]>=cutOffs[h-1])&(tab[,"time"]<cutOffs[h])),]
		MJs = matrix(nrow=dim(admins2@data)[1], ncol=dim(admins2@data)[1]); states = unique(tab[,"state"])
		for (i in 1:dim(admins2@data)[1])
			{
				for (j in 1:dim(admins2@data)[1])
					{
						sub = tab[which((tab[,"from"]==admins2@data[i,"GID_2"])&(tab[,"to"]==admins2@data[j,"GID_2"])),]
						MJs_post = c(); # print(c(i,j))
						for (k in 1:length(states))
							{
								MJ = which(sub[,"state"]==states[k])
								MJs_post = c(MJs_post, length(MJ))
							}
						MJs[i,j] = mean(MJs_post)
					}
			}
		row.names(MJs) = admins2@data[,"GID_2"]; colnames(MJs) = admins2@data[,"GID_2"]
		write.table(MJs, paste0("BEAST_DTA_analysis/With_empirical_trees/Alignment_",analysis,"_MJt",h,".csv"), quote=F, sep=",")
	}

	# 7.2. Computing a normalised Shannon entropy through time for each municipality

		# 7.2.1. First approach based on lineage transition events transition events

entropies = list()
selected_periods = c("2024-08","2024-09","2024-10","2024-11","2024-12","2025-01","2025-02","2025-03","2025-04","2025-05","2025-06","2025-07","2025-08")
selected_periods = c("2024-08"); last_period = "2025-09" # to select the overall period
for (i in 1:dim(admins2@data)[1])
	{
		admin2 = admins2@data[i,"GID_2"]; buffer = matrix(nrow=nberOfExtractionFiles, ncol=length(selected_periods))
		for (j in 1:nberOfExtractionFiles)
			{
				tab = read.csv(paste0("BEAST_DTA_analysis/With_empirical_trees/Alignment_",analysis,"_",nberOfTreesToSample,"_ext/TreeExtractions_",j,".csv"), head=T)
				for (k in 1:length(selected_periods))
					{
						date1 = decimal_date(ymd(paste0(selected_periods[k],"-01")))
						if (k < length(selected_periods))
							{
								date2 = decimal_date(ymd(paste0(selected_periods[k+1],"-01")))
							}	else	{
								date2 = decimal_date(ymd(paste0(last_period,"-01")))
							}
						sub1 = tab[which((tab[,"endYear"]>=date1)&(tab[,"endYear"]<date2)),]
						if (dim(sub1)[1] > 0)
							{
								sub2 = sub1[which((sub1[,"startLoc"]==admin2)|(sub1[,"endLoc"]==admin2)),]
								if (dim(sub2)[1] == 0) print(c(i,j,k))
								for (l in 1:dim(sub2)[1])
									{
										if (sub2[l,"startLoc"] != sub2[l,"endLoc"]) sub2[l,"length"] = sub2[l,"length"]/2
									}
								if (dim(sub2)[1] > 0)
									{			
										t_total = sum(sub2[,"length"])
										partitions = matrix(nrow=dim(sub2)[1], ncol=1)
										roots = which(!sub2[,"node1"]%in%sub2[,"node2"])
										for (l in 1:length(roots))
											{
												partitions[roots[l],1] = l
											}
										while (sum(is.na(partitions[,1])) != 0)
											{
												for (l in 1:dim(partitions)[1])
													{
														if (is.na(partitions[l,1]))
															{
																ancestor = which(sub2[,"node2"]==sub2[l,"node1"])
																if (!is.na(partitions[ancestor,1]))
																	{
																		partitions[l,1] = partitions[ancestor,1]
																	}
															}
													}
											}
										entropy = 0; n = length(unique(partitions[,1]))
										for (l in 1:n)
											{
												t = sum(sub2[which(partitions[,1]==l),"length"]); p = t/t_total
												entropy = entropy + (p*log(p))
											}
										buffer[j,k] = -entropy/log(n)
									}
							}
					}
			}
		if (length(selected_periods) == 1)
			{
				meanV = round(mean(buffer[,1], na.rm=T),3); hpd95 = round(HDInterval::hdi(buffer[,1]),3)
				cat(admin2,": mean Shannon entropy = ",meanV,", 95% HPD = [",hpd95[1],", ",hpd95[2],"]\n",sep="")
			}
		entropies[[i]] = buffer
	}

		# 7.2.2. Second approach based on Markov jumps (like in Lemey et al. 2021, Nature)

entropies = list()
tabs = readRDS(paste0("BEAST_DTA_analysis/With_empirical_trees/Alignment_",analysis,"_emp3.rds"))
for (i in 1:dim(admins2@data)[1])
	{
		admin2 = admins2@data[i,"GID_2"]; buffer = matrix(nrow=length(tabs), ncol=1)
		for (j in 1:length(tabs))
			{
				sub = tabs[[j]][which(tabs[[j]][,"startLoc"]==admin2),]
				if (dim(sub)[1] == 0) print(c(i,j,k))
				if (dim(sub)[1] > 0)
					{			
						t_total = sum(sub[,"length"])
						partitions = matrix(nrow=dim(sub)[1], ncol=1)
						roots = which(!sub[,"node1"]%in%sub[,"node2"])
						for (l in 1:length(roots))
							{
								partitions[roots[l],1] = l
							}
						while (sum(is.na(partitions[,1])) != 0)
							{
								for (l in 1:dim(partitions)[1])
									{
										if (is.na(partitions[l,1]))
											{
												ancestor = which(sub[,"node2"]==sub[l,"node1"])
												if (!is.na(partitions[ancestor,1]))
													{
														partitions[l,1] = partitions[ancestor,1]
													}
											}
									}
							}
						entropy = 0; n = length(unique(partitions[,1]))
						for (l in 1:n)
							{
								t = sum(sub[which(partitions[,1]==l),"length"]); p = t/t_total
								entropy = entropy + (p*log(p))
							}
						buffer[j,k] = -entropy/log(n)
					}
			}
		meanV = round(mean(buffer[,1], na.rm=T),3); hpd95 = round(HDInterval::hdi(buffer[,1]),3)
		cat(admin2,": mean Shannon entropy = ",meanV,", 95% HPD = [",hpd95[1],", ",hpd95[2],"]\n",sep="")
		entropies[[i]] = buffer
	}

		# 7.2.3. Investigating the correlation between entropy estimates and population counts

population_counts = read.csv("GLM_predictors_data/Prepared_variables/Population_counts.csv", head=T)
mean_entropies = population_counts; mean_entropies[,2] = NA; colnames(mean_entropies)[2] = "mean_entropy"
for (i in 1:dim(admins2@data)[1]) mean_entropies[i,2] = mean(entropies[[i]][,1], na.rm=T)
indices = which(!is.na(mean_entropies[,2])); v = cor(population_counts[indices,2], mean_entropies[indices,2], method="spearman")
cor.test(population_counts[indices,2], mean_entropies[indices,2], method="spearman", alternative="less") # rS = 0.190
n = 10000; vS_permutations = rep(NA, n); x = population_counts[indices,2]; y = mean_entropies[indices,2]
for (i in 1:n) vS_permutations[i] = cor(x, sample(y, length(y), replace=F), method="spearman")
pValue = sum(v > vS_permutations)/n # p-value = 0.818

		# 7.2.4. Generation of the graphic with the median and 95% HPD for each municipality

tab_entropies = as.data.frame(matrix(nrow=dim(admins2@data)[1], ncol=5)); writingFiles = FALSE
tab_entropies[,1] = c(1:dim(admins2@data)[1]); tab_entropies[,2] = admins2@data[,"GID_2"]
colnames(tab_entropies) = c("municipality_index","municipality","median","95pHDP_lower","95pHDP_upper")
for (i in 1:dim(admins2@data)[1])
	{
		tab_entropies[i,"median"] = median(entropies[[i]][,1]); hpd95 = HDInterval::hdi(entropies[[i]][,1])
		tab_entropies[i,"95pHDP_lower"] = hpd95[1]; tab_entropies[i,"95pHDP_upper"] = hpd95[2]
	}
higherThan050 = sum(tab_entropies[,"95pHDP_lower"]>0.5); higherThan075 = sum(tab_entropies[,"95pHDP_lower"]>0.75) # 19 and 13/24
buffer = tab_entropies[,2:dim(tab_entropies)[2]]; buffer[,2] = round(as.numeric(buffer[,2]),4)
buffer[,3] = round(as.numeric(buffer[,3]),4); buffer[,4] = round(as.numeric(buffer[,4]),4)
if (writingFiles) write.csv(buffer, "Shannon_entropies.csv", row.names=F, quote=F)

pdf(paste0("Figure_S5_",analysis,"_NEW.pdf"), width=8, height=2.2) # dev.new(width=8, height=2.2)
par(oma=c(0,0,0,0), mar=c(3.0,3.0,0.5,0.0), lwd=0.3, bty="o", col="gray30", col.axis="gray30", fg="gray30")
plot(tab_entropies[,"municipality_index"], tab_entropies[,"median"], col=NA, axes=F, xlab=NA, ylab=NA, ylim=c(0,1))
abline(h=0.5, col="gray50", lwd=1, lty=2); abline(h=0.75, col="gray50", lwd=1, lty=2)
for (i in 1:dim(tab_entropies)[1])
	{
		colour = rgb(70,118,187,100,maxColorValue=255) # blue
		x1 = tab_entropies[i,"municipality_index"]-0.45
		x2 = tab_entropies[i,"municipality_index"]+0.45
		y1 = tab_entropies[i,"95pHDP_lower"]; y2 = tab_entropies[i,"95pHDP_upper"]
		polygon(c(x1,x2,x2,x1), c(y1,y1,y2,y2), col=colour, border="gray30")
		x1 = tab_entropies[i,"municipality_index"]-0.44
		x2 = tab_entropies[i,"municipality_index"]+0.44
		y1 = tab_entropies[i,"median"]; y2 = tab_entropies[i,"median"]
		lines(c(x1,x2), cbind(y1,y2), lwd=1, col="gray30")
		# polygon(c(x1,x2,x2,x1), c(y1,y1,y2,y2), col="gray30", border=2)
	}
ats = c(1:dim(tab_entropies)[1]); labels = admins2@data[,"NAME_2"]
axis(side=1, lwd=0.5, cex.axis=0.7, mgp=c(0,0.17,0), lwd.tick=0, col.lab="gray30", col="gray30", tck=-0.04, las=1, at=c(-2,dim(tab_entropies)[1]+2))
axis(side=1, lwd=0, cex.axis=0.7, mgp=c(0,0.50,0), lwd.tick=0.5, col.lab="gray30", col="gray30", tck=-0.04, las=2, at=ats, label=labels)
axis(side=2, lwd=0.5, cex.axis=0.7, mgp=c(0,0.4,-0.1), lwd.tick=0.5, col.lab="gray30", col="gray30", tck=-0.04, las=1, padj=0.4)
mtext("Normalised Shannon entropy", side=2, col="gray30", cex=0.7, line=1.7, las=3)
dev.off()

	# 7.3. Visualisation of the time-scaled phylogenetic inference (based on the MCC tree)

tree = readAnnotatedNexus(paste0("BEAST_DTA_analysis/With_empirical_trees/Alignment_",analysis,"_",nberOfExtractionFiles,".tree"))
tab = read.table(paste0("Alignment_",analysis,".txt"), head=T)
rootHeight = max(nodeHeights(tree)); root_time = mostRecentSamplingDatum-rootHeight
node_dates = root_time+nodeHeights(tree)[,2]
node_weeks = interval(min(ymd(tab[,"collection_date"])),date_decimal(node_dates))%/%weeks(1)+1
node_weeks[node_weeks<1] = 1; nodes_cols = paste0(collection_week_cols1[node_weeks],"BF")

pdf(paste0("Figure_S4_",analysis,"_NEW.pdf"), width=8, height=8) # dev.new(width=8, height=8)
par(oma=c(0,0,0,0), mar=c(0.0,0.0,0.0,0), lwd=0.3, bty="o", col="gray30", col.axis="gray30", fg="gray30", lheight=0.85)
plot(tree, type="fan", show.tip.label=F, show.node.label=F, edge.width=0.2, cex=0.6, align.tip.label=3, col="gray30", edge.color="gray30")
for (i in 1:dim(tree$edge)[1])
	{
		if (!tree$edge[i,2]%in%tree$edge[,1])
			{
				nodelabels(node=tree$edge[i,2], pch=16, cex=0.6, col="white")
				nodelabels(node=tree$edge[i,2], pch=16, cex=0.6, col=nodes_cols[i])
				nodelabels(node=tree$edge[i,2], pch=1, cex=0.6, col="gray30", lwd=0.3)
			}
	}
add.scale.bar(x=0.37, y=-0.9, length=NULL, ask=F, lwd=0.5 , lcol="gray30", cex=0.7)
dev.off()

	# 7.4. Visualising the dispersal history of viral lineages among municipalities

nberOfExtractionFiles = 100; tMRCAs = rep(NA, nberOfExtractionFiles)
for (i in 1:nberOfExtractionFiles)
	{
		dta_tab = read.csv(paste0("BEAST_DTA_analysis/With_empirical_trees/Alignment_",analysis,"_",nberOfExtractionFiles,"_ext/TreeExtractions_",i,".csv"), head=T)
		tMRCAs[i] = min(dta_tab[,"startYear"])
	}
minYear = as.numeric(HDInterval::hdi(tMRCAs)[1]); maxYear = mostRecentSamplingDatum
cutOffs = c(decimal_date(ymd(c("2025-01-01","2025-03-01","2025-05-01","2025-09-01"))))
titles1 = c("Aug. 2024 -","Jan. 2025 -","Mar. 2025 -","May 2025 -")
titles2 = c("Dec. 2024","Feb. 2025","Apr. 2025","Aug 2025")
matrices_list = list(); matrix_means = list()
for (h in 1:length(cutOffs))
	{
		matrices = list()
		for (i in 1:nberOfExtractionFiles)
			{
				mat = matrix(0, nrow=length(locations), ncol=length(locations))
				row.names(mat) = locations; colnames(mat) = locations
				dta_tab = read.csv(paste0("BEAST_DTA_analysis/With_empirical_trees/Alignment_",analysis,"_",nberOfExtractionFiles,"_ext/TreeExtractions_",i,".csv"), head=T)
				if (h == 1) dta_tab = dta_tab[which(dta_tab[,"endYear"]<cutOffs[h]),]
				if (h > 1) dta_tab = dta_tab[which((dta_tab[,"endYear"]>=cutOffs[h-1])&(dta_tab[,"endYear"]<cutOffs[h])),]
				for (j in 1:dim(dta_tab)[1])
					{
						index1 = which(locations==dta_tab[j,"startLoc"])
						index2 = which(locations==dta_tab[j,"endLoc"])
						mat[index1,index2] = mat[index1,index2]+1
					}
				matrices[[i]] = mat
			}
		matrices_list[[h]] = matrices
	}
saveRDS(matrices_list, paste0("BEAST_DTA_analysis/With_empirical_trees/Alignment_",analysis,"_",nberOfExtractionFiles,".rds"))
matrices_list = readRDS(paste0("BEAST_DTA_analysis/With_empirical_trees/Alignment_",analysis,"_",nberOfExtractionFiles,".rds"))
minVals1 = 9999; minVals2 = 9999; maxVals1 = -9999; maxVals2 = -9999
for (h in 1:length(cutOffs))
	{
		matrix_mean = matrix(0, nrow=length(locations), ncol=length(locations))
		for (i in 1:nberOfExtractionFiles) matrix_mean = matrix_mean+matrices_list[[h]][[i]]
		matrix_mean = matrix_mean/nberOfExtractionFiles
		matrix_means[[h]] = matrix_mean; mat = matrix_mean; diag(mat) = NA
		if (minVals1 > min(diag(matrix_mean))) minVals1 = min(diag(matrix_mean))
		if (maxVals1 < max(diag(matrix_mean))) maxVals1 = max(diag(matrix_mean))
		if (minVals2 > min(mat, na.rm=T)) minVals2 = min(mat, na.rm=T)
		if (maxVals2 < max(mat, na.rm=T)) maxVals2 = max(mat, na.rm=T)
	}
centroids = coordinates(admins2); row.names(centroids) = admins2@data[,"GID_2"] # not used
tab = read.table(paste0("BEAST_DTA_analysis/Alignment_",analysis,".txt"), head=T, sep="\t")
for (i in 1:length(admins2))
	{
		pts = tab[which(tab[,"location"]==admins2@data[i,"GID_2"]),c("longitude","latitude")]
		centroids[i,1] = mean(pts[,1]); centroids[i,2] = mean(pts[,2])
	}
main_cities = c("Saint-\nDenis","Saint-\nPaul","Saint-Pierre","Le Tampon","Saint-\nAndré",
				"Saint-Louis","Le Port","Saint-Benoît","Saint-Joseph","Sainte-Marie")
main_city_coordinates = rbind(cbind(55.4485,-20.8785),cbind(55.2792,-21.0093),
	cbind(55.4777,-21.3409),cbind(55.5153,-21.2709),cbind(55.6509,-20.9604),cbind(55.4079,-21.2866),
	cbind(55.2903,-20.9357),cbind(55.7130,-21.0335),cbind(55.6190,-21.3788),cbind(55.5490,-20.8974))
city_name_coordinates = rbind(cbind(55.4185,-20.9145),cbind(55.2492,-21.0033),
	cbind(55.4277,-21.3689),cbind(55.5153,-21.2709),cbind(55.6179,-20.9854),cbind(55.3499,-21.3176),
	cbind(55.2473,-20.9507),cbind(55.7690,-21.0375),cbind(55.5900,-21.4130),cbind(55.5790,-20.8934))

pdf(paste0("Figure_3_",analysis,"_NEW.pdf"), width=8, height=7.2); using_MJs = FALSE # dev.new(width=8, height=7.2)
par(mfrow=c(2,2), oma=c(0,0,0,0), mar=c(0.0,0.0,0.0,0), lwd=0.3, bty="o", col="gray30", col.axis="gray30", fg="gray30", lheight=0.85)
for (h in 1:length(cutOffs))
	{
		multiplier1 = 500; multiplier2 = 2; multiplier3 = 0.1
		plot(admin0, col=NA, border=NA); plot(elevation, add=T, legend=F, col=elevation_cols)
		plot(residential_areas$geometry, add=T, border=NA, col=rgb(77,77,77,180,maxColorValue=255)) # previously alpha = 120
		if (h == 1) plot(admins2, add=T, col=NA, border="white", lwd=0.3)
		plot(main_roads$geometry, add=T, lwd=0.7, col=rgb(222,67,39,255,maxColorValue=255))
		plot(admin0, col=NA, border="gray30", lwd=0.3, add=T)
		mtext(titles1[h], at=55.265, line=-1.6, cex=0.7, col="gray30")
		mtext(titles2[h], at=55.257, line=-2.4, cex=0.7, col="gray30")
		mat = matrix_means[[h]]
		points(centroids, cex=sqrt((multiplier1*((diag(mat)-minVals1)/(maxVals1-minVals1)))/pi), pch=16, col="#DE432780")
		if (using_MJs == TRUE)
			{
				mat = read.csv(paste0("BEAST_DTA_analysis/With_empirical_trees/Alignment_",analysis,"_MJt",h,".csv"))
			}
		for (i in 1:dim(centroids)[1])
			{
				for (j in 1:dim(centroids)[1])
					{
						if ((i!=j)&(mat[i,j]>=0.95)&(mat[i,j]<30))
							{
								LWD = (((mat[i,j]-minVals2)/(maxVals2-minVals2))*multiplier2)+0.2; arrow = 0
								curvedarrow(centroids[i,], centroids[j,], arr.length=arrow*1.3, arr.width=arrow, lwd=LWD, lty=1,
											lcol="black", arr.col=NA, arr.pos=0.5, curve=0.15, dr=NA, endhead=F, arr.type="none")
							}
						if ((i!=j)&(mat[i,j]>=0.95)&(mat[i,j]>=30))
							{
								LWD = (((mat[i,j]-minVals2)/(maxVals2-minVals2))*multiplier2)+0.2; arrow = (multiplier3*(mat[i,j]/maxVals2))+0.02
								curvedarrow(centroids[i,], centroids[j,], arr.length=arrow*1.3, arr.width=arrow, lwd=LWD, lty=1,
											lcol="black", arr.col="black", arr.pos=0.5, curve=0.15, dr=NA, endhead=F, arr.type="triangle")
							}
					}
			}
		if (h == 1)
			{
				buffer = city_name_coordinates; buffer[,2] = buffer[,2]+0.015
				points(main_city_coordinates, pch=16, cex=0.8, col="black")
				text(buffer, paste0(main_cities,""), cex=0.7, col="black", adj=c(0.5,0.5))
			}
		if (h == 2)
			{
				points(rbind(cbind(55.242,-21.34),cbind(55.319,-21.365),cbind(55.358,-21.376)), cex=sqrt((multiplier1*((c(500,100,20)-minVals1)/(maxVals1-minVals1)))/pi), pch=16, col="#DE432780")				
				text(rbind(cbind(55.242,-21.34),cbind(55.319,-21.365)), labels=c(500,100), cex=0.8, col="gray30", adj=c(0.5,0.5))
				text(cbind(55.428,-21.376), labels=c("20 transitions"), cex=0.8, col="gray30", adj=c(0.5,0.5))
				labels = rev(c("0-300m","300-600m","600-900m","900-1200m","1200-1500m","1500-1800m","1800-2100m","2100-2400m","2400-2700m","2700-3050m"))
				legend(x=55.725, y=-20.867, labels, text.col="gray30", pch=15, pt.cex=1.3, col=rev(elevation_cols), box.lty=0, cex=0.65, x.intersp=0.75, y.intersp=0.80)
				legend(x=55.737, y=-21.013, c("",""), text.col=NA, pch=16, pt.cex=1.3, col=c(elevation_cols[1],NA), box.lty=0, cex=0.65, x.intersp=0.75, y.intersp=0.75)
				legend(x=55.737, y=-21.013, c("Residential","areas"), text.col="gray30", pch=16, pt.cex=1.3, col=c(rgb(77,77,77,180,maxColorValue=255),NA),
					   box.lty=0, cex=0.65, x.intersp=0.75, y.intersp=0.8)
			}
		if (h == 3)
			{
				LWD = (((10-minVals2)/(maxVals2-minVals2))*multiplier2)+0.2; arrow = (multiplier3*(10/maxVals2))+0.02
				curvedarrow(cbind(55.225,-21.330), cbind(55.300,-21.330), arr.length=arrow*1.3, arr.width=arrow, lwd=LWD, lty=1,
							lcol="gray30", arr.col="gray30", arr.pos=0.5, curve=0, dr=NA, endhead=F, arr.type="triangle")
				LWD = (((30-minVals2)/(maxVals2-minVals2))*multiplier2)+0.2; arrow = (multiplier3*(30/maxVals2))+0.02
				curvedarrow(cbind(55.225,-21.350), cbind(55.300,-21.350), arr.length=arrow*1.3, arr.width=arrow, lwd=LWD, lty=1,
							lcol="gray30", arr.col="gray30", arr.pos=0.5, curve=0, dr=NA, endhead=F, arr.type="triangle")
				LWD = (((90-minVals2)/(maxVals2-minVals2))*multiplier2)+0.2; arrow = (multiplier3*(90/maxVals2))+0.02
				curvedarrow(cbind(55.225,-21.370), cbind(55.300,-21.370), arr.length=arrow*1.3, arr.width=arrow, lwd=LWD, lty=1,
							lcol="gray30", arr.col="gray30", arr.pos=0.5, curve=0, dr=NA, endhead=F, arr.type="triangle")
				text(rbind(cbind(55.318,-21.330),cbind(55.318,-21.350)), labels=c(10,30), cex=0.8, col="gray30", adj=c(0.5,0.5))
				text(cbind(55.363,-21.370), labels=c("90 transitions"), cex=0.8, col="gray30", adj=c(0.5,0.5))
			}
		if (h == 4)
			{
				points(rbind(cbind(55.215,-21.360),cbind(55.350,-21.360)), pch=16, cex=0.8, col="gray30")
				LWD = (((35-minVals2)/(maxVals2-minVals2))*multiplier2)+0.2; arrow = (multiplier3*(90/maxVals2))+0.02
				curvedarrow(cbind(55.215,-21.360), cbind(55.350,-21.360), arr.length=arrow*1.3, arr.width=arrow, lwd=LWD, lty=1,
							lcol="gray30", arr.col="gray30", arr.pos=0.5, curve=0.15, dr=NA, endhead=F, arr.type="triangle")
				curvedarrow(cbind(55.350,-21.360), cbind(55.215,-21.360), arr.length=arrow*1.3, arr.width=arrow, lwd=LWD, lty=1,
							lcol="gray30", arr.col="gray30", arr.pos=0.5, curve=0.15, dr=NA, endhead=F, arr.type="triangle")
				text(cbind(55.250,-21.306), labels=c("Dispersal\ndirection:"), cex=0.8, col="gray30", adj=c(0.5,0.5))
			}
	}
dev.off()

pdf(paste0("Figure_S1_",analysis,"_NEW.pdf"), width=8, height=7.2) # dev.new(width=8, height=7.2)
par(mfrow=c(2,2), oma=c(0,0,0,0), mar=c(0.0,0.0,0.0,0), lwd=0.3, bty="o", col="gray30", col.axis="gray30", fg="gray30", lheight=0.85)
months = interval(min(ymd(tab[,"collection_date"])),ymd(tab[,"collection_date"]))%/%months(1)+1
dates = decimal_date(ymd(tab[,"collection_date"])); minYear = min(decimal_date(ymd("2024-08-01")))
maxYear = max(dates); endYears_colours = collection_month_cols[(((dates-minYear)/(maxYear-minYear))*100)+1]
cexs = c(1.0, 0.8, 0.6, 0.7)
for (h in 1:length(cutOffs))
	{
		plot(admin0, col=NA, border=NA); plot(elevation, add=T, legend=F, col=elevation_cols)
		plot(residential_areas$geometry, add=T, border=NA, col=rgb(77,77,77,180,maxColorValue=255))
		if (h == 1) plot(admins2, add=T, col=NA, border="white", lwd=0.3)
		plot(main_roads$geometry, add=T, lwd=0.7, col=rgb(222,67,39,255,maxColorValue=255))
		plot(admin0, col=NA, border="gray30", lwd=0.3, add=T)
		mtext(titles1[h], at=55.265, line=-1.6, cex=0.7, col="gray30")
		mtext(titles2[h], at=55.257, line=-2.4, cex=0.7, col="gray30")
		if (h == 1)
			{
				buffer = city_name_coordinates; buffer[,2] = buffer[,2]+0.015
				points(main_city_coordinates, pch=16, cex=0.8, col="black")
				text(buffer, paste0(main_cities,""), cex=0.7, col="black", adj=c(0.5,0.5))
			}
		if (h == 1) indices = which(dates<cutOffs[h])
		if (h > 1) indices = which((dates<cutOffs[h])&(dates>=cutOffs[h-1]))
		sub = tab[indices,]; sub = sub[order(sub[,"collection_date"]),]
		months = interval(min(ymd(tab[,"collection_date"])),ymd(sub[,"collection_date"]))%/%months(1)+1
		for (i in dim(sub)[1]:1)
			{
				points(sub[i,c("longitude","latitude")], pch=16, cex=cexs[h], col=collection_month_cols[months[i]])
				points(sub[i,c("longitude","latitude")], pch=1, cex=cexs[h], lwd=0.4, col="black")
			}
		if (h == 2)
			{
				labels = rev(c("0-300m","300-600m","600-900m","900-1200m","1200-1500m","1500-1800m","1800-2100m","2100-2400m","2400-2700m","2700-3050m"))
				legend(x=55.725, y=-20.867, labels, text.col="gray30", pch=15, pt.cex=1.3, col=rev(elevation_cols), box.lty=0, cex=0.65, x.intersp=0.75, y.intersp=0.80)
				legend(x=55.737, y=-21.013, c("",""), text.col=NA, pch=16, pt.cex=1.3, col=c(elevation_cols[1],NA), box.lty=0, cex=0.65, x.intersp=0.75, y.intersp=0.75)
				legend(x=55.737, y=-21.013, c("Residential","areas"), text.col="gray30", pch=16, pt.cex=1.3, col=c(rgb(77,77,77,180,maxColorValue=255),NA), box.lty=0, cex=0.65, x.intersp=0.75, y.intersp=0.8)
			}
		if (h == 4)
			{
				labels = c("2024-08","2024-09","2024-10","2024-11","2024-12","2025-01","2025-02","2025-03","2025-04","2025-05","2025-06",">")
				legend(x=55.755, y=-20.870, labels, text.col="gray30", pch=15, pt.cex=1.3, col=collection_month_cols, box.lty=0, cex=0.65, x.intersp=0.75, y.intersp=0.80)
			}
	}
dev.off()

pdf(paste0("Figure_1A_",analysis,"_NEW.pdf"), width=8, height=7.2) # dev.new(width=8, height=7.2)
par(mfrow=c(2,2), oma=c(0,0,0,0), mar=c(0.0,0.0,0.0,0), lwd=0.3, bty="o", col="gray30", col.axis="gray30", fg="gray30", lheight=0.85)
months = interval(min(ymd(tab[,"collection_date"])),ymd(tab[,"collection_date"]))%/%months(1)+1
dates = decimal_date(ymd(tab[,"collection_date"])); minYear = min(decimal_date(ymd("2024-08-01")))
maxYear = max(dates); endYears_colours = collection_month_cols[(((dates-minYear)/(maxYear-minYear))*100)+1]
for (h in 1:2)
	{
		plot(admin0, col=NA, border=NA); plot(elevation, add=T, legend=F, col=elevation_cols)
		plot(residential_areas$geometry, add=T, border=NA, col=rgb(77,77,77,200,maxColorValue=255))
		plot(main_roads$geometry, add=T, lwd=0.7, col=rgb(222,67,39,255,maxColorValue=255))
		plot(admin0, col=NA, border="gray30", lwd=0.3, add=T)
		sub = tab; sub = sub[order(sub[,"collection_date"]),]
		weeks = interval(min(ymd(tab[,"collection_date"])),ymd(sub[,"collection_date"]))%/%weeks(1)+1
		for (i in dim(sub)[1]:1)
			{
				points(sub[i,c("longitude","latitude")], pch=16, cex=0.6, col=collection_week_cols1[weeks[i]])
				points(sub[i,c("longitude","latitude")], pch=1, cex=0.6, lwd=0.4, col="black")
			}
		labels = rev(c("0-300m","300-600m","600-900m","900-1200m","1200-1500m","1500-1800m","1800-2100m","2100-2400m","2400-2700m","2700-3050m"))
		legend(x=55.725, y=-20.867, labels, text.col="gray30", pch=15, pt.cex=1.3, col=rev(elevation_cols), box.lty=0, cex=0.65, x.intersp=0.75, y.intersp=0.80)
		legend(x=55.737, y=-21.013, c("",""), text.col=NA, pch=16, pt.cex=1.3, col=c(elevation_cols[1],NA), box.lty=0, cex=0.65, x.intersp=0.75, y.intersp=0.75)
		legend(x=55.737, y=-21.013, c("Residential","areas"), text.col="gray30", pch=16, pt.cex=1.3, col=c(rgb(77,77,77,200,maxColorValue=255),NA), box.lty=0, cex=0.65, x.intersp=0.75, y.intersp=0.8)
	}
dev.off()

	# 7.5. Preparing the input files to generate an animated visualisation with spread.gl

buffer = matrix(nrow=dim(centroids)[1], ncol=3); buffer[,1] = row.names(centroids)
buffer[,3] = centroids[,1]; buffer[,2] = centroids[,2]; colnames(buffer) = c("location","latitude","longitude")
write.csv(buffer, paste0("BEAST_DTA_analysis/With_empirical_trees/Alignment_",analysis,"_sgl.csv"), row.names=F, quote=F)
mcc_tre = readAnnotatedNexus(paste0("BEAST_DTA_analysis/With_empirical_trees/Alignment_",analysis,"_",nberOfExtractionFiles,".tree"))
source("treeExtractions_MCC.r"); mcc_tab = treeExtractions_MCC(mcc_tre, mostRecentSamplingDatum)
write.csv(mcc_tab, paste0("BEAST_DTA_analysis/With_empirical_trees/Alignment_",analysis,"_",nberOfExtractionFiles,".csv"), row.names=F, quote=F)
	# --> problem: some start/endLocation are associated with NA values in the MCC tree (tree STATE_107000000, i.e. tree number 35):
mcc_tab = read.csv(paste0("BEAST_DTA_analysis/With_empirical_trees/Alignment_",analysis,"_",nberOfExtractionFiles,"_ext/TreeExtractions_35.csv"), head=T)
colNames = c("id","type","length","start_time","end_time","start_name","start_lat","start_lon","end_name","end_lat","end_lon")
buffer = matrix(nrow=dim(mcc_tab)[1], ncol=11); colnames(buffer) = colNames
for (i in 1:dim(buffer)[1])
	{
		buffer[i,"id"] = i; buffer[i,"type"] = "internal"; buffer[i,"length"] = mcc_tab[i,"length"]
		if (!mcc_tab[i,"node2"]%in%mcc_tab[,"node1"]) buffer[i,"type"] = "external"
		buffer[i,"start_time"] = as.character(date_decimal(mcc_tab[i,"startYear"]))
		buffer[i,"end_time"] = as.character(date_decimal(mcc_tab[i,"endYear"]))
		buffer[i,"start_time"] = unlist(strsplit(buffer[i,"start_time"],"\\."))[1]
		buffer[i,"end_time"] = unlist(strsplit(buffer[i,"end_time"],"\\."))[1]
		buffer[i,"start_name"] = mcc_tab[i,"startLoc"]; buffer[i,"end_name"] = mcc_tab[i,"endLoc"]
		index1 = which(row.names(centroids)==mcc_tab[i,"startLoc"])
		index2 = which(row.names(centroids)==mcc_tab[i,"endLoc"])
		buffer[i,"start_lat"] = centroids[index1,2]; buffer[i,"start_lon"] = centroids[index1,1]
		buffer[i,"end_lat"] = centroids[index2,2]; buffer[i,"end_lon"] = centroids[index2,1]
	}
write.csv(buffer, paste0("BEAST_DTA_analysis/With_empirical_trees/Alignment_",analysis,"_mcc.csv"), row.names=F, quote=F)

# 8. Analysing and reporting the results of the discrete-GLM analysis

predictors = c("population_count_origin","population_count_destination","geographic_distances","pairwise_mobility_metric",
			   "temperature_origin","temperature_destination","precipitation_origin","precipitation_destination")
log = read.table(paste0("BEAST_DTA_analysis/With_empirical_trees/Alignment_",analysis,"_glm2.log"), head=T, sep="\t")
burnIn = ((dim(log)[1]-1)/10)+1; indices = (burnIn+1):dim(log)[1]; log = log[indices,]
glm_results = matrix(nrow=length(predictors), ncol=3); row.names(glm_results) = predictors
colnames(glm_results) = c("inclusion_probability","BF","coefficient_if_included")
q = 1-(0.5^(1/length(predictors))) # corresponds to a probability of zero predictor being included equal to 0.5:
				 # if K = number of predictors, P(zero predictor included) = 0.5 = (1-q)^K --> q = 1-(0.5^(1/K))
for (i in 1:length(predictors))
	{
		glm_results[i,"inclusion_probability"] = round(mean(log[,paste0("location.coefIndicators",i)]),3)
		p = as.numeric(glm_results[i,"inclusion_probability"]); glm_results[i,"BF"] = round((p/(1-p))/(q/(1-q)),1)
		median = round(median(log[,paste0("location.coefficientsTimesIndicators",i)]),3)
		hpd95 = round(hdi(log[,paste0("location.coefficientsTimesIndicators",i)])[1:2],3)
		glm_results[i,"coefficient_if_included"] = paste0(median," [",hpd95[1],", ",hpd95[2],"]")
		if (as.numeric(glm_results[i,"BF"]) > 999) glm_results[i,"BF"] = ">999"
		
	}
write.csv(glm_results, paste0("BEAST_DTA_analysis/With_empirical_trees/Alignment_",analysis,"_glm.csv"), quote=F)

# 9. Conducting complementary isolation-by-distance (IBD) analyses

	# 9.1. Approach based on the code within the "spreadStatistics" function (too slow)

nberOfExtractionFiles = 100; rP2s = rep(NA, nberOfExtractionFiles)
for (i in 1:nberOfExtractionFiles)
	{
		tab = read.csv(paste0("BEAST_DTA_analysis/With_empirical_trees/Alignment_",analysis,"_",nberOfExtractionFiles,"_ext/TreeExtractions_",i,".csv"), head=T)
		tipNodeIndices = which(!tab[,"node2"]%in%tab[,"node1"])
		distTree = matrix(nrow=length(tipNodeIndices), ncol=length(tipNodeIndices))
		for (k in 2:dim(distTree)[1]) # extracted from the "spreadStatistics" function of "seraphim"
			{
				for (l in 1:(k-1))
					{
						index1 = tipNodeIndices[k]
						index2 = tipNodeIndices[l]
						indices1 = index1; root = FALSE
						while (root == FALSE)
							{	
								if (tab[indices1[length(indices1)],"node1"]%in%tab[,"node2"])
									{
										indices1 = c(indices1, which(tab[,"node2"]==tab[indices1[length(indices1)],"node1"]))
									}	else	{
										root = TRUE
									}
							}
						indices2 = index2; root = FALSE
						while (root == FALSE)
							{	
								if (tab[indices2[length(indices2)],"node1"]%in%tab[,"node2"])
									{
										indices2 = c(indices2, which(tab[,"node2"]==tab[indices2[length(indices2)],"node1"]))
									}	else	{
										root = TRUE
									}
							}
						indices3 = indices1[which(indices1%in%indices2)]; patristic_dis = NULL
						if (length(indices3) == 0)
							{
								patristic_dis = sum(tab[c(indices1,indices2),"length"])
							}	else	{
								patristic_dis = sum(tab[c(indices1[which(!indices1%in%indices3)],indices2[which(!indices2%in%indices3)]),"length"])
							}
						distTree[k,l] = patristic_dis; distTree[l,k] = patristic_dis
					}
			}
		distsGeo = matrix(nrow=dim(distTree)[1], ncol=dim(distTree)[2])
		for (k in 2:dim(distsGeo)[1]) # extracted from the "spreadStatistics" function of "seraphim"
			{
				for (l in 1:(k-1))
					{
						index1 = tipNodeIndices[k]
						index2 = tipNodeIndices[l]
						x1 = cbind(tab[index1,"endLon"], tab[index1,"endLat"])
						x2 = cbind(tab[index2,"endLon"], tab[index2,"endLat"])
						distsGeo[k,l] = rdist.earth(x1, x2, miles=F, R=NULL)
						distsGeo[l,k] = distsGeo[k,l]
					}
			}	
		rP2s[i] = cor(distTree[lower.tri(distTree)],log(distsGeo[lower.tri(distsGeo)]+1), method="pearson")
	}
meanV = round(mean(rP2s), 3); hpd95 = round(hdi(rP2s)[1:2], 3)
cat("rP2 = ",meanV,", 95% HPD = [",hpd95[1],", ",hpd95[2],"]\n",sep="")

	# 9.2. Approach based on the pairwise patristic distances directly computed on trees

nberOfExtractionFiles = 100; rP2s = rep(NA, nberOfExtractionFiles)
trees = read.nexus(paste0("BEAST_DTA_analysis/With_empirical_trees/Alignment_",analysis,"_",nberOfExtractionFiles,".trees"))
for (i in 1:length(trees))
	{
		tab = read.csv(paste0("BEAST_DTA_analysis/With_empirical_trees/Alignment_",analysis,"_",nberOfExtractionFiles,"_ext/TreeExtractions_",i,".csv"), head=T)
		distTree = as.matrix(distTips(trees[[i]], method="patristic"))
		distsGeo = matrix(nrow=dim(distTree)[1], ncol=dim(distTree)[2])
		for (j in 2:dim(distsGeo)[1])
			{
				for (k in 1:(j-1))
					{
						index1 = which(tab[,"tipLabel"]==row.names(distTree)[j])
						index2 = which(tab[,"tipLabel"]==colnames(distTree)[k])
						x1 = cbind(tab[index1,"endLon"], tab[index1,"endLat"])
						x2 = cbind(tab[index2,"endLon"], tab[index2,"endLat"])
						distsGeo[j,k] = rdist.earth(x1, x2, miles=F, R=NULL)
						distsGeo[k,j] = distsGeo[j,k]
					}
			}
		rP2s[i] = cor(distTree[lower.tri(distTree)],log(distsGeo[lower.tri(distsGeo)]+1), method="pearson")
	}
meanV = round(mean(rP2s), 3); hpd95 = round(hdi(rP2s)[1:2], 3)
cat("rP2 = ",meanV,", 95% HPD = [",hpd95[1],", ",hpd95[2],"]\n",sep="")
	# rP2 = 0.111, 95% HPD = [0.100, 0.127]

