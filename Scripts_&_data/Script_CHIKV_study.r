# TO DO:
	# - R(t) estimations: analysis with the approach of Cori et al. (2013) but this time based on all case data (not only genomic)
	# - R(t) estimations: analysis with the episodic birth–death sampling (EBDS) model (based on a fixed tree or empirical trees)
	# - R(t) estimations: analysis with a coalescent-based approach following Volz & Didelot (2018) or Borchering et al. (2021) ??
# ANALYSES:
	# - running --> Virtual Machine : SA, GLM, DTA (pending); MacStudio Pro : SKG, SA, EBDS on a fixed tree, RRW, DTA, GLM
	# - to run --> EBDS on the fnial SA MCC tree, as well as EBDS on the final SA empirical trees
# FIGURES:
	# - Figure 1: sampling map with points coloured according to time + weekly epi-curve + SA + EBDS (+ R(t) coalescent-based?)
	# - Figure 2: DTA reconstruction in 4 panels; Figure 3: results of the DTA-GLM analyses (boxplot of beta*inclusion probability?)
	# - Figure S1: 4 sampling map with points coloured according to time; Figure S2: 1° BEAST analysis; and Figure 3: time-scaled tree

library(colorspace)
library(diagram)
library(EpiEstim)
library(HDInterval)
library(lubridate)
library(MetBrewer)
library(seraphim)
library(sf)

analysis = "101025"

# 1. Investigating the temporal signal associated with the alignment
# 2. Preliminary BEAST analysis to estimate the substitution rate
# 3. Some attempts to conduct a continuous phylogeographic analysis
# 4. Preparing the skygrid and discrete phylogeographic analyses
# 5. Preparing the different predictors for the discrete-GLM analysis
# 6. Evolution of Ne and R(t): skygrid and sampling aware analyses
# 7. Discrete phylogeographic analysis based on the municipalities

tab = read.table(paste0("BEAST_CTA_analysis/Alignment_",analysis,".txt"), head=T)
mostRecentSamplingDatum = max(decimal_date(ymd(tab[,"collection_date"]))) # 2025.595
admin0 = shapefile("GADM_REU_shapefile/GADM_REU_0.shp") # island borders
admins2 = shapefile("GADM_REU_shapefile/GADM_REU_2.shp") # municipalities
locations = unique(admins2@data[,"GID_2"]) # municipality IDs
elevation = mask(crop(raster("Copernicus_DEM30.tif"), admin0), admin0)
elevation[elevation[]<0] = 0; elevation[elevation[]>3000] = 3000
osm_lines = st_read("OpenStreetMap_files/OSM_25102014.pbf", layer="lines")
osm_polygons = st_read("OpenStreetMap_files/OSM_25102014.pbf", layer="multipolygons")
main_roads = subset(osm_lines, highway%in%c("trunk","primary"))
residential_areas = subset(osm_polygons, landuse=="residential")
construction_areas = subset(osm_polygons, landuse=="construction")
forest_areas = subset(osm_polygons, landuse=="forest")
elevation_cols = paste0(divergingx_hcl(14,"fall")[4:14],"BF") # 75% transparency
collection_week_cols = met.brewer(name="Hiroshige", n=60, type="continuous")[1:51]
collection_month_cols = met.brewer(name="Hiroshige", n=14, type="continuous")[1:12]

# 1. Investigating the temporal signal associated with the alignment

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

# 2. Preliminary BEAST analysis to estimate the substitution rate

tab = read.csv(paste0("BEAST_first_analysis/Alignment_",analysis,".csv"), head=T)
log = read.table(paste0("BEAST_first_analysis/Alignment_",analysis,".log"), head=T, sep="\t")
ucld_mean = log[ceiling(0.1*dim(log)[1]):dim(log)[1],"ucld.mean"]
tree = readAnnotatedNexus(paste0("BEAST_first_analysis/Alignment_",analysis,".tree"))
rootHeight = max(nodeHeights(tree)); root_time = mostRecentSamplingDatum-rootHeight
minYear = mostRecentSamplingDatum-tree$root.annotation$`height_95%_HPD`[[2]]
maxYear = mostRecentSamplingDatum; nodes_cols = c(); locations = c()
locations = unique(tab[,"location"]); locations = locations[order(locations)]
location_cols = colorRampPalette(brewer.pal(11,"Spectral"))(length(locations))
location_cols = c("#FAA521","#4676BB","#D1E5F0","#C0C0C0","#6A3D9A","#DE4327","#BF812D")

pdf(paste0("Figure_S2_",analysis,"_NEW1.pdf"), width=8, height=3.5) # dev.new(width=8, height=4)
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

pdf(paste0("Figure_S2_",analysis,"_NEW2.pdf"), width=3, height=1.5) # dev.new(width=3, height=2)
par(oma=c(0,0,0,0), mar=c(1.5,2,0.7,1), lwd=0.3, bty="o", col="gray30", col.axis="gray30", fg="gray30")
plot(density(ucld_mean), axes=F, ann=F, col=NA); polygon(stats::density(ucld_mean), col="#DE43274D", border=NA); lines(density(ucld_mean), col="#DE4327", lwd=1)
ats = seq(0.0002,0.0007,0.0001); labels = c("",expression("3x10"^"-4"),expression("4x10"^"-4"),expression("5x10"^"-4"),expression("6x10"^"-4"),"")
axis(side=1, lwd=0.5, cex.axis=0.60, mgp=c(0,0,0), lwd.tick=0.5, col.lab="gray30", col="gray30", tck=-0.03, las=1, at=ats, label=labels)
axis(side=2, lwd=0.5, cex.axis=0.60, mgp=c(0,0.1,-0.3), lwd.tick=0.5, col.lab="gray30", col="gray30", tck=-0.03, las=1, padj=0.22)
dev.off()

# 3. Some attempts to conduct a continuous phylogeographic analysis

	# 3.1. Subsampling by phylogenetic clustering and then keeping only one sample per set of coordinates

tab = read.csv(paste0("Alignment_",analysis,".csv"), head=T, sep=";")
tre = read.tree(paste0("Alignment_",analysis,".tre")); tree = tre; subTrees = subtrees(tree, wait=F)
c1 = 0; clusters1 = list() # all clusters (including the "nested" ones)
c2 = 0; clusters2 = list() # all clusters minus the "nested" clades
for (i in 2:length(subTrees)) # the 1st subtree is the entire tree
	{
		subTree = subTrees[i][[1]]; labels = subTree$tip.label
		locations = rep(NA, length(labels))
		for (j in 1:length(labels))
			{
				coordinates = tab[which(tab[,1]==labels[j]),2:3]
				locations[j] = paste0(coordinates[1],"_",coordinates[2])
			}
		if (length(unique(locations)) == 1)
			{
				c1 = c1+1; clusters1[[c1]] = labels
			}
	}
for (i in 1:length(clusters1))
	{
		nested = FALSE
		for (j in 1:length(clusters1))
			{
				if (i != j)
					{
						allSequencesIncluded = TRUE
						for (k in 1:length(clusters1[[i]]))
							{
								if (!clusters1[[i]][k]%in%clusters1[[j]]) allSequencesIncluded = FALSE
							}
						if (allSequencesIncluded == TRUE) nested = TRUE
					}
			}
		if (nested == FALSE)
			{
				c2 = c2+1; clusters2[[c2]] = clusters1[[i]]
			}	
	}
sequencesToRemove1 = c(); sequencesToRemove2 = c()
for (i in 1:length(clusters2))
	{
		sequencesToRemove1 = c(sequencesToRemove1, sample(clusters2[[i]],length(clusters2[[i]])-1,replace=F))
	}
fasta1 = scan(paste0("Alignment_",analysis,".fas"), what="", sep="\n", quiet=T); fasta2 = c(); fasta3 = c()
for (i in 1:length(fasta1))
	{
		if ((grepl(">",fasta1[i]))&&(!gsub(">","",fasta1[i])%in%sequencesToRemove1)) fasta2 = c(fasta2, fasta1[c(i,i+1)])
	}
tab = read.csv(paste0("Alignment_",analysis,".csv"), head=T, sep=";")
tab1 = read.table(paste0("Alignment_",analysis,".txt"), head=T)
tab1 = tab1[which(!tab1[,1]%in%sequencesToRemove1),]
tab2 = tab[which(!tab[,1]%in%sequencesToRemove1),]
tab = cbind(tab1, matrix(nrow=dim(tab1)[1],ncol=2))
for (i in 1:dim(tab)[1])
	{
		tab[i,3:4] = tab2[which(tab2[,1]==tab[i,1]),3:2]
	}
colnames(tab) = c("trait","collection_date","latitude","longitude")
coordinates = rep(NA, dim(tab)[1])
for (i in 1:dim(tab)[1])
	{
		coordinates[i] = paste0(tab[i,3],"_",tab[i,4])
	}
unique_coordinates = unique(coordinates) # 3946/3105
for (i in 1:length(unique_coordinates))
	{
		indices = which(coordinates==unique_coordinates[i])
		if (length(indices) > 1)
			{
				sequencesToRemove2 = c(sequencesToRemove2, sample(tab[indices,1],length(indices)-1,replace=F))
			}
	}
tab = tab[which(!tab[,1]%in%sequencesToRemove2),]
for (i in 1:length(fasta2))
	{
		if ((grepl(">",fasta2[i]))&&(!gsub(">","",fasta2[i])%in%sequencesToRemove2)) fasta3 = c(fasta3, fasta2[c(i,i+1)])
	}
write.table(tab, paste0("BEAST_CTA_analysis/Alignment_",analysis,".txt"), row.names=F, quote=F, sep="\t")
write(fasta3, paste0("BEAST_CTA_analysis/Alignment_",analysis,".fas"))

	# 3.2. Extracting the spatio-temporal information embedded in posterior trees

tab = read.table(paste0("BEAST_CTA_analysis/Alignment_",analysis,".txt"), head=T)
mostRecentSamplingDatum = max(decimal_date(ymd(tab[,"collection_date"])))
trees = readAnnotatedNexus(paste0("BEAST_CTA_analysis/With_empirical_trees/Alignment_",analysis,"_emp.trees"))
for (i in 1:length(trees))
	{
		tab = postTreeExtractions(trees[[i]], mostRecentSamplingDatum)
		write.csv(tab, paste0("BEAST_CTA_analysis/With_empirical_trees/Alignment_",analysis,"_emp_ext/TreeExtractions_",i,".csv"), row.names=F, quote=F)
	}

# 4. Preparing the skygrid and discrete phylogeographic analyses

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

# 5. Preparing the different predictors for the discrete-GLM analysis

	# 5.1. Population count at the municipality of origin and at the destination municipality

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
write.csv(population_counts, "GLM_predictors_data/Prepared_predictors/Population_counts.csv", quote=F)

	# 5.2. Climatic data at the municipality of origin and at the destination municipality

climatic_variables = matrix(nrow=dim(admins2)[1], ncol=3) # monthly averages further averaged from 2024-08 to 2025-08
colnames(climatic_variables) = c("temperature","precipitation","relative_humidity"); row.names(climatic_variables) = admins2@data[,"GID_2"]
climatic_data = read.csv("GLM_predictors_data/Original_data_files/Monthly_climate_data_2024-25.csv", head=T, sep=";")
selected_months = c("202408","202409","202410","202411","202412","202501","202502","202503","202504","202505","202506","202507","202508")
climatic_data = climatic_data[which(climatic_data[,"AAAAMM"]%in%selected_months),]
lon_lat = unique(paste0(climatic_data[,"LON"],"_",climatic_data[,"LAT"])); stations = c()
for (i in 1:length(lon_lat))
	{
		stations = rbind(stations, cbind(as.numeric(unlist(strsplit(lon_lat[i],"_"))[1]),as.numeric(unlist(strsplit(lon_lat[i],"_"))[2])))
	}
for (i in 1:dim(climatic_variables)[1])
	{
		pol_x = admins2@polygons[[i]]@Polygons[[1]]@coords[,1]; pol_y = admins2@polygons[[i]]@Polygons[[1]]@coords[,2]
		indices = which(point.in.polygon(climatic_data[,"LON"],climatic_data[,"LAT"],pol_x,pol_y) == 1)		
		if (length(indices) == 0) # no station within municipality REU.3.3_1 --> taking the closed station for associating the climatic data
			{
				minimum_distances = rep(NA, dim(stations)[1])
				for (j in 1:dim(stations)[1]) minimum_distances[j] = min(rdist.earth(cbind(pol_x,pol_y), cbind(stations[j,1],stations[j,2]), miles=F))
				index = which(minimum_distances==min(minimum_distances, na.rm=T))
				indices = which((climatic_data[,"LON"]==stations[index,1])&(climatic_data[,"LAT"]==stations[index,2]))
			}
		climatic_variables[i,"temperature"] = mean(climatic_data[indices,"TMM"], na.rm=T) # "moyenne mensuelle des températures moyennes (TM) quotidiennes (en °C et 1/10)"
		climatic_variables[i,"precipitation"] = mean(climatic_data[indices,"RR"], na.rm=T) # "cumul mensuel des hauteurs de précipitation (en mm et 1/10)"
		climatic_variables[i,"relative_humidity"] = mean(climatic_data[indices,"UMM"], na.rm=T) # "moyenne mensuelle des humidités moyennes (UM) quotidiennes (en %)"
	}
write.csv(climatic_variables[,1:2], "GLM_predictors_data/Prepared_predictors/Climatic_variables.csv", quote=F) # too many missing values for the relative humidity

	# 5.3. Pairwise great-circle geographic distance between population-weighted centroid points

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
write.csv(geographic_distances, "GLM_predictors_data/Prepared_predictors/Geographic_distances.csv", quote=F)

	# 5.4. Pairwise measure of the mobility (commuting) flux between municipalities (symetric)

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
write.csv(mobility_metric, "GLM_predictors_data/Prepared_predictors/Pairwise_mobility_metric.csv", quote=F)

# 6. Evolution of Ne and R(t): skygrid and sampling aware analyses

	# 6.1. Computation of the evolution of R(t) following the method of Cori et al. (2013)
	
tab = read.table(paste0("BEAST_DTA_analysis/Alignment_",analysis,".txt"), head=T, sep="\t") # to be replaced by all the known positive cases
total_number_of_days = interval(min(ymd(tab[,"collection_date"])),max(ymd(tab[,"collection_date"])))%/%days(1)-1
cases_day = interval(min(ymd(tab[,"collection_date"])),ymd(tab[,"collection_date"]))%/%days(1)+1
daily_cases = rep(NA, total_number_of_days)
for (i in 1:length(daily_cases)) daily_cases[i] = sum(cases_day==i)

n_MC = 1000 # number of Monte-Carlo draws
mean_range = c(10, 23); sd_range = c(4, 8) # days
t_start = seq(2, length(daily_cases)-6); t_end = seq(8, length(daily_cases)) # sliding window
all_R = matrix(NA, nrow=length(t_start), ncol=n_MC)
for(i in 1:n_MC)
	{
		mean_si_i = runif(1, mean_range[1], mean_range[2])
		sd_si_i = runif(1, sd_range[1], sd_range[2])
  		res_i = estimate_R(incid=daily_cases, method="parametric_si", config=make_config(list(
      						mean_si=mean_si_i, std_si=sd_si_i, t_start=t_start, t_end=t_end)))
		all_R[,i] = res_i$R$`Mean(R)`
	}
R_median = apply(all_R, 1, median); R_days = (t_start+t_end)/2
R_dates = decimal_date(min(ymd(tab[,"collection_date"]))+R_days)
R_lower = apply(all_R, 1, quantile, probs=0.025, na.rm=T)
R_upper = apply(all_R, 1, quantile, probs=0.975, na.rm=T)

	# 6.2. Visualisation of the evolution through time of the number of cases, R(t), and Ne

tab = read.table(paste0("BEAST_DTA_analysis/Alignment_",analysis,".txt"), head=T, sep="\t")
tab_weeks = interval(min(ymd(tab[,"collection_date"])),ymd(tab[,"collection_date"]))%/%weeks(1)+1
skg = read.csv(paste0("BEAST_DTA_analysis/Without_DTA_model/Alignment_",analysis,"_sa.csv"), head=T)[1:51,]
skg = skg[,c("time","median","lower","upper")]; colnames(skg) = c("time","median","95pHDP_lower","95pHDP_upper")
skg[,"median"] = log(skg[,"median"]+1); timeSlice = skg[1,"time"]-skg[2,"time"]; skg_weeks = 51:1
skg[,"95pHDP_lower"] = log(skg[,"95pHDP_lower"]+1); skg[,"95pHDP_upper"] = log(skg[,"95pHDP_upper"]+1)

pdf(paste0("Figure_1B_",analysis,"_NEW.pdf"), width=8/2, height=8/2) # dev.new(width=8/2, height=8/2)
par(mfrow=c(3,1), oma=c(0,0,0,0), mar=c(2.0,3.5,0.1,0.5), lwd=0.3, bty="o", col="gray30", col.axis="gray30", fg="gray30")
hist(decimal_date(ymd(tab[,"collection_date"])), breaks=51, col="gray70", border=NA, axes=F, ann=F, ylim=c(0,380), xlim=c(2024.62,2025.59))
hist(decimal_date(ymd(tab[,"collection_date"])), breaks=51, col=paste0(collection_week_cols,"BF"), add=T)
dates = c("2024-08-07","2024-09-01","2024-10-01","2024-11-01","2024-12-01","2025-01-01","2025-02-01","2025-03-01","2025-04-01","2025-05-01","2025-06-01","2025-07-01","2025-08-08")
ats1 = decimal_date(ymd(dates)); ats2 = ats1[2:(length(ats1)-1)]; labels1 = gsub("-","\\/",dates[2:(length(dates)-1)])
axis(side=1, lwd=0.5, cex.axis=0.7, mgp=c(0,0.17,0), lwd.tick=0, col.lab="gray30", col="gray30", tck=-0.04, las=1, at=ats1, label=rep("",length(ats1)))
axis(side=1, lwd=0, cex.axis=0.7, mgp=c(0,0.17,0), lwd.tick=0.5, col.lab="gray30", col="gray30", tck=-0.04, las=1, at=ats2, label=labels1)
axis(side=2, lwd=0.5, cex.axis=0.7, mgp=c(0,0.4,-0.1), lwd.tick=0.5, col.lab="gray30", col="gray30", tck=-0.04, las=1, padj=0.4, at=seq(0,300,100))
mtext("Genomic samples", side=2, col="gray30", cex=0.6, line=1.7, las=3)
plot(skg[,"time"], skg[,"median"], lwd=0.7, type="l", cex.axis=0.8, cex.lab=0.8, col="gray30", axes=F, xlab=NA, ylab=NA, ylim=c(0,7), xlim=c(2024.62,2025.59))
xx_l = c(skg[,c("time")],rev(skg[,c("time")])); yy_l = c(skg[,"95pHDP_lower"],rev(skg[,"95pHDP_upper"]))
getOption("scipen"); opt = options("scipen"=20); polygon(xx_l,yy_l,col=rgb(187/255,187/255,187/255,0.25),border=0)
for (i in 1:length(skg[,"time"]))
	{
		colour = paste0(collection_week_cols[skg_weeks[i]],"BF")
		x1 = skg[i,"time"]-(timeSlice/2); x2 = skg[i,"time"]+(timeSlice/2)
		y1 = skg[i,"95pHDP_lower"]-1; y2 = skg[i,"95pHDP_upper"]+1
		polygon(c(x1,x2,x2,x1), c(y1,y1,y2,y2), col="gray70", border=NA)
		polygon(c(x1,x2,x2,x1), c(y1,y1,y2,y2), col=colour, border=NA)
	}
getOption("scipen"); opt = options("scipen"=20); polygon(xx_l,yy_l,col=NA,border="gray30")
lines(skg[,"time"], skg[,"median"], lwd=0.3, type="l", cex.axis=0.8, cex.lab=0.8, col="gray30")
labels2 = c(1,5,10,20,50,100,200,500,2000); ats3 = c(log(labels2))
axis(side=1, lwd=0.5, cex.axis=0.7, mgp=c(0,0.17,0), lwd.tick=0, col.lab="gray30", col="gray30", tck=-0.04, las=1, at=ats1, label=rep("",length(ats1)))
axis(side=1, lwd=0, cex.axis=0.7, mgp=c(0,0.17,0), lwd.tick=0.5, col.lab="gray30", col="gray30", tck=-0.04, las=1, at=ats2, label=labels1)
axis(side=2, lwd=0.5, cex.axis=0.7, mgp=c(0,0.4,-0.1), lwd.tick=0.5, col.lab="gray30", col="gray30", tck=-0.04, las=1, padj=0.4, at=ats3, label=labels2)
mtext("Effective population size", side=2, col="gray30", cex=0.55, line=1.7, las=3)
plot(R_dates, R_median, lwd=0.7, type="l", cex.axis=0.8, cex.lab=0.8, col="gray30", axes=F, xlab=NA, ylab=NA, ylim=c(0,9.6), xlim=c(2024.62,2025.59))
xx_l = c(R_dates,rev(R_dates)); yy_l = c(R_lower,rev(R_upper)); polygon(xx_l,yy_l,col=rgb(187/255,187/255,187/255,0.25),border=0)
for (i in 1:length(R_dates))
	{
		week = interval(min(ymd(tab[,"collection_date"])),date_decimal(R_dates[i]))%/%weeks(1)+1
		colour = paste0(collection_week_cols[week],"BF")
		x1 = R_dates[i]-(timeSlice/2); x2 = R_dates[i]+(timeSlice/2)
		y1 = R_lower[i]-1; y2 = R_upper[i]+1
		polygon(c(x1,x2,x2,x1), c(y1,y1,y2,y2), col="gray70", border=NA)
		polygon(c(x1,x2,x2,x1), c(y1,y1,y2,y2), col=colour, border=NA)
	}
getOption("scipen"); opt = options("scipen"=20); polygon(xx_l,yy_l,col="NA",border="gray30")
lines(R_dates, R_median, lwd=0.3, type="l", cex.axis=0.8, cex.lab=0.8, col="gray30"); abline(h=1, lty=2, lwd=0.3, col="gray30")
axis(side=1, lwd=0.5, cex.axis=0.7, mgp=c(0,0.17,0), lwd.tick=0, col.lab="gray30", col="gray30", tck=-0.04, las=1, at=ats1, label=rep("",length(ats1)))
axis(side=1, lwd=0, cex.axis=0.7, mgp=c(0,0.17,0), lwd.tick=0.5, col.lab="gray30", col="gray30", tck=-0.04, las=1, at=ats2, label=labels1)
axis(side=2, lwd=0.5, cex.axis=0.7, mgp=c(0,0.4,-0.1), lwd.tick=0.5, col.lab="gray30", col="gray30", tck=-0.04, las=1, padj=0.4)
mtext("Effective reproduction number", side=2, col="gray30", cex=0.55, line=1.7, las=3)
dev.off()

# 7. Discrete phylogeographic analysis based on the municipalities

	# 7.1. Extracting the spatio-temporal information embedded in posterior trees

nberOfTreesToSample = 1000; nberOfExtractionFiles = nberOfTreesToSample; burnIn = XXXX
log = scan(paste0("BEAST_DTA_analysis/With_empirical_trees/Alignment_",analysis,"_emp1.log"), what="", sep="\n", quiet=T, blank.lines.skip=F)
index1 = 6+burnIn; index2 = length(log); interval = round((index2-index1)/nberOfTreesToSample)
indices = seq(index2-((nberOfTreesToSample-1)*interval),index2,interval)
write(log[c(5,indices)], paste0("BEAST_DTA_analysis/With_empirical_trees/Alignment_",analysis,"_",nberOfExtractionFiles,".log"))
trees = scan(paste0("BEAST_DTA_analysis/With_empirical_trees/Alignment_",analysis,"_emp1.trees"), what="", sep="\n", quiet=T, blank.lines.skip=F)
index1 = which(trees=="\t\t;")[length(which(trees=="\t\t;"))]; index2 = index1 + burnIn + 1
indices3 = which(grepl("tree STATE",trees)); index3 = indices3[length(indices3)]
interval = floor((index3-(index1+burnIn))/nberOfTreesToSample)
indices = seq(index3-((nberOfTreesToSample-1)*interval),index3,interval)
selected_trees = c(trees[c(1:index1,indices)],"End;")
write(selected_trees, paste0("BEAST_DTA_analysis/With_empirical_trees/Alignment_",analysis,"_",nberOfExtractionFiles,".trees"))

source("treeExtractions_DTA.r"); matrices_list = list(); matrix_means = list()
trees = readAnnotatedNexus(paste0("BEAST_DTA_analysis/With_empirical_trees/Alignment_",analysis,"_",nberOfExtractionFiles,".trees"))
# trees = readAnnotatedNexus(paste0("Alignment_",analysis,"_emp1.trees"))[[2]]
for (i in 1:nberOfExtractionFiles)
	{
		if (length(trees) == 1) dta_tab = treeExtractions_DTA(trees, mostRecentSamplingDatum)	
		if (length(trees) > 1) dta_tab = treeExtractions_DTA(trees[i[]], mostRecentSamplingDatum)	
		write.csv(dta_tab, paste0("BEAST_DTA_analysis/With_empirical_trees/Alignment_",analysis,"_1000_ext/TreeExtractions_",i,".csv"), row.names=F, quote=F)
	}

	# 7.2. Visualising the dispersal history of viral lineages among municipalities

tMRCAs = rep(NA, nberOfExtractionFiles)
for (i in 1:nberOfExtractionFiles)
	{
		dta_tab = read.csv(paste0("BEAST_DTA_analysis/With_empirical_trees/Alignment_",analysis,"_1000_ext/TreeExtractions_",i,".csv"))
		tMRCAs[i] = min(dta_tab[,"startYear"])
	}
minYear = as.numeric(HDInterval::hdi(tMRCAs)[1]); maxYear = mostRecentSamplingDatum
cutOffs = c(decimal_date(ymd(c("2025-01-01","2025-03-01","2025-05-01","2025-09-01"))))
titles1 = c("Aug. 2024 -","Jan. 2025 -","Mar. 2025 -","May 2025 -")
titles2 = c("Dec. 2024","Feb. 2025","Apr. 2025","Aug 2025")
for (h in 1:length(cutOffs))
	{
		matrices = list()
		for (i in 1:nberOfExtractionFiles)
			{
				mat = matrix(0, nrow=length(locations), ncol=length(locations))
				row.names(mat) = locations; colnames(mat) = locations
				dta_tab = read.csv(paste0("BEAST_DTA_analysis/With_empirical_trees/Alignment_",analysis,"_1000_ext/TreeExtractions_",i,".csv"), head=T)
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
saveRDS(matrices_list, paste0("BEAST_DTA_analysis/With_empirical_trees/Alignment_",analysis,"_1000.rds"))
matrices_list = readRDS(paste0("BEAST_DTA_analysis/With_empirical_trees/Alignment_",analysis,"_1000.rds"))
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

pdf(paste0("Figure_2_",analysis,".pdf"), width=8, height=7.2) # dev.new(width=8, height=7.2)
par(mfrow=c(2,2), oma=c(0,0,0,0), mar=c(0.0,0.0,0.0,0), lwd=0.3, bty="o", col="gray30", col.axis="gray30", fg="gray30", lheight=0.85)
for (h in 1:length(cutOffs))
	{
		plot(admin0, col=NA, border=NA); plot(elevation, add=T, legend=F, col=elevation_cols)
		plot(residential_areas$geometry, add=T, border=NA, col=rgb(77,77,77,120,maxColorValue=255))
		if (h == 1) plot(admins2, add=T, col=NA, border="white", lwd=0.3)
		plot(main_roads$geometry, add=T, lwd=0.7, col=rgb(222,67,39,255,maxColorValue=255))
		plot(admin0, col=NA, border="gray30", lwd=0.3, add=T)
		mtext(titles1[h], at=55.265, line=-1.6, cex=0.7, col="gray30")
		mtext(titles2[h], at=55.257, line=-2.4, cex=0.7, col="gray30")
		mat = matrix_means[[h]]; multiplier1 = 500; multiplier2 = 2; multiplier3 = 0.1
		points(centroids, cex=sqrt((multiplier1*((diag(mat)-minVals1)/(maxVals1-minVals1)))/pi), pch=16, col="#DE432780")
		for (i in 1:dim(centroids)[1])
			{
				for (j in 1:dim(centroids)[1])
					{
						if ((i!=j)&(mat[i,j]>1)&(mat[i,j]<30))
							{
								LWD = (((mat[i,j]-minVals2)/(maxVals2-minVals2))*multiplier2)+0.2; arrow = 0
								curvedarrow(centroids[i,], centroids[j,], arr.length=arrow*1.3, arr.width=arrow, lwd=LWD, lty=1,
											lcol="black", arr.col=NA, arr.pos=0.5, curve=0.15, dr=NA, endhead=F, arr.type="none")
							}
						if ((i!=j)&(mat[i,j]>1)&(mat[i,j]>=30))
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
				legend(x=55.737, y=-21.013, c("Residential","areas"), text.col="gray30", pch=16, pt.cex=1.3, col=c(rgb(77,77,77,120,maxColorValue=255),NA), box.lty=0, cex=0.65, x.intersp=0.75, y.intersp=0.8)
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

pdf(paste0("Figure_S1_",analysis,".pdf"), width=8, height=7.2) # dev.new(width=8, height=7.2)
par(mfrow=c(2,2), oma=c(0,0,0,0), mar=c(0.0,0.0,0.0,0), lwd=0.3, bty="o", col="gray30", col.axis="gray30", fg="gray30", lheight=0.85)
months = interval(min(ymd(tab[,"collection_date"])),ymd(tab[,"collection_date"]))%/%months(1)+1
dates = decimal_date(ymd(tab[,"collection_date"])); minYear = min(decimal_date(ymd("2024-08-01")))
maxYear = max(dates); endYears_colours = collection_month_cols[(((dates-minYear)/(maxYear-minYear))*100)+1]
cexs = c(1.0, 0.8, 0.6, 0.7)
for (h in 1:length(cutOffs))
	{
		plot(admin0, col=NA, border=NA); plot(elevation, add=T, legend=F, col=elevation_cols)
		plot(residential_areas$geometry, add=T, border=NA, col=rgb(77,77,77,120,maxColorValue=255))
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
				legend(x=55.737, y=-21.013, c("Residential","areas"), text.col="gray30", pch=16, pt.cex=1.3, col=c(rgb(77,77,77,120,maxColorValue=255),NA), box.lty=0, cex=0.65, x.intersp=0.75, y.intersp=0.8)
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
		plot(residential_areas$geometry, add=T, border=NA, col=rgb(77,77,77,120,maxColorValue=255))
		plot(main_roads$geometry, add=T, lwd=0.7, col=rgb(222,67,39,255,maxColorValue=255))
		plot(admin0, col=NA, border="gray30", lwd=0.3, add=T)
		sub = tab; sub = sub[order(sub[,"collection_date"]),]
		weeks = interval(min(ymd(tab[,"collection_date"])),ymd(sub[,"collection_date"]))%/%weeks(1)+1
		for (i in dim(sub)[1]:1)
			{
				points(sub[i,c("longitude","latitude")], pch=16, cex=0.6, col=collection_week_cols[weeks[i]])
				points(sub[i,c("longitude","latitude")], pch=1, cex=0.6, lwd=0.4, col="black")
			}
		labels = rev(c("0-300m","300-600m","600-900m","900-1200m","1200-1500m","1500-1800m","1800-2100m","2100-2400m","2400-2700m","2700-3050m"))
		legend(x=55.725, y=-20.867, labels, text.col="gray30", pch=15, pt.cex=1.3, col=rev(elevation_cols), box.lty=0, cex=0.65, x.intersp=0.75, y.intersp=0.80)
		legend(x=55.737, y=-21.013, c("",""), text.col=NA, pch=16, pt.cex=1.3, col=c(elevation_cols[1],NA), box.lty=0, cex=0.65, x.intersp=0.75, y.intersp=0.75)
		legend(x=55.737, y=-21.013, c("Residential","areas"), text.col="gray30", pch=16, pt.cex=1.3, col=c(rgb(77,77,77,120,maxColorValue=255),NA), box.lty=0, cex=0.65, x.intersp=0.75, y.intersp=0.8)
	}
dev.off()

system(paste0("magick -units PixelsPerInch -density 1000 Figure_2_",analysis,".pdf -background white -alpha remove -flatten Alignment_",analysis,".png"))
system(paste0("magick -units PixelsPerInch -density 1000 Figure_S1_",analysis,".pdf -background white -alpha remove -flatten Sampling_4_maps.png"))

