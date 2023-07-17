library(ade4)
library(ape)
library(blockCV)
library(diagram)
library(dismo)
library(gbm)
library(geosphere)
library(ggplot2)
library(HDInterval)
library(lubridate)
library(ncdf4)
library(ncf)
library(picante)
library(pgirmess)
library(phytools)
library(RColorBrewer)
library(raster)
library(rgdal)
library(rgeos)
library(seqinr)
library(sp)
library(vioplot)

# 1. Preparation of land cover and climatic environmental rasters
# 2. Loading and ploting all the Bombus records by continent
# 3. Defining the background area for each time period
# 4. Boosted regression trees (BRT) analyses with standard or spatial cross-validation
# 5. Computation of the prevalence-pseudoabsence-calibrated Sørensen index
# 6. Comparison of the two BRT models trained for the two time periods
# 7. Estimation of spatial sorting bias (SSB)
# 8. Comparison of the response curves for each environmental factor
# 9. Analyses of the relative influence of each environmental factor
# 10. BRT projections based on past and present environmental variables
# 11. Post hoc analyses to compare past and present ENM projections
# 12. Plotting the environmental variables and their future projections
# 13. BRT projections on historical and future scenarios
# 14. Post hoc analyses to compare present and future ENM projections

directory = "Bombus_obs_160620"; savingPlots = FALSE
timePeriods = c("1901_1974","2000_2014"); periods = list()
periods[[1]] = c(1901,1974); periods[[2]] = c(2000,2014)
species = read.csv("Bombus_species.csv", header=T)
europe1 = shapefile("Continent_shapefile/Europe.shp")
coastlines = shapefile("Continents_shapefile/Coastlines.shp")
mask = raster("Mask_raster_file.nc4")

# 1. Preparation of land cover and climatic environmental rasters

	# 1.1. Preparation of the European shapefile that will be used as a mask

polygons = list(); c = 0
for (i in 1:length(europe1@polygons[[1]]@Polygons))
	{
		if (europe1@polygons[[1]]@Polygons[[i]]@area > 1)
			{
				c = c+1; polygons[[c]] = europe1@polygons[[1]]@Polygons[[i]]
			}
	}
pols = Polygons(polygons, 1); pols_list = list(); pols_list[[1]] = pols
europe2 = SpatialPolygons(pols_list); europe3 = gSimplify(europe2, 0.1)

	# 1.2. Loading the human population, temperature, precipitation, and land cover rasters

populations = list(); temperatures = list(); precipitations = list(); land_covers = list(); buffer_land_covers = list()
for (t in 1:length(periods))
	{
		populations[[t]] = raster(paste0("Environmental_rasters/Calibration/population_histsoc_0p5deg_annual_",timePeriods[t],"_timmean.nc4"))
		temperatures[[t]] = raster(paste0("Environmental_rasters/Calibration/tas_day_GSWP3+EWEMBI_HistObs_r1i1p1_EWEMBI_",timePeriods[t],"_timmean.nc4"))
		precipitations[[t]] = raster(paste0("Environmental_rasters/Calibration/pr_day_GSWP3+EWEMBI_HistObs_r1i1p1_EWEMBI_",timePeriods[t],"_timmean.nc4"))
		land_covers[[t]] = nc_open(paste0("Environmental_rasters/Calibration/landcover_HistObs_annual_",timePeriods[t],"_timmean.nc4"))
		names(populations[[t]]) = "population"; names(temperatures[[t]]) = "temperature"; names(precipitations[[t]]) = "precipitation"
	}

	# 1.3. Preparation of distinct land cover rasters from the initial ".nc" object

variable_codes = c("croplands","pastures","urbanAreas","primaryForest","primaryNonF","secondaryForest","secondaryNonF")
variable_names = c("crops","pasture","urban land","forested primary land","non-forested primary land",
				   "potentially forested secondary land","potentially non-forested secondary land")
landCoverVariableNames = as.character(read.csv("Land_cover_vars.csv")[1:12,2])
for (t in 1:length(periods))
	{
		landCoverVariableIDs = names(land_covers[[t]]$var); cols = list()
		land_covers1 = list(); land_covers2 = list(); land_covers3 = list()
		for (i in 2:13) land_covers1[[i-1]] = brick(paste0("Environmental_rasters/Calibration/landcover_HistObs_annual_",timePeriods[t],"_timmean.nc4"),
													varname=landCoverVariableIDs[i])
		for (i in 1:length(land_covers1))
			{
				names(land_covers1[[i]]) = landCoverVariableNames[i]
				cols[[i]] = colorRampPalette(brewer.pal(9,"YlGn"))(120)[11:110]
				if (i == 1)
					{
						r_global = land_covers1[[1]]
					}	else	{
						r_global[] = r_global[]+land_covers1[[i]][]	
					}
			}
		for (i in 1:length(variable_names))
			{
				names = gsub("\\."," ",landCoverVariableNames); indices = which(landCoverVariableNames==variable_names[i])
				if (length(indices) == 0) indices = which(grepl(variable_names[i],names))
				if (variable_names[i] == "pasture") indices = c(indices, which(grepl("rangeland",names)))
				land_cover = land_covers1[[indices[1]]]; names(land_cover) = variable_codes[i]
				if (length(indices) > 1)
					{
						for (j in 2:length(indices)) land_cover[] = land_cover[]+land_covers1[[indices[j]]][]
					}
				land_covers2[[i]] = land_cover[[1]]; land_covers3[[i]] = raster::aggregate(land_cover[[1]],2)
			}
		buffer_land_covers[[t]] = land_covers3
	}
land_covers = buffer_land_covers

	# 1.4. Selecting and treating the environmental rasters for the ENM analyses

envVariables_list = list(); nullRasters = list()
for (t in 1:length(periods))
	{
		envVariables = list()
		envVariables[[1]] = temperatures[[t]]
		envVariables[[2]] = precipitations[[t]]
		envVariables[[3]] = land_covers[[t]][[4]] # primary forest areas
		envVariables[[4]] = land_covers[[t]][[5]] # primary non-forest areas
		envVariables[[5]] = land_covers[[t]][[6]] # secondary forest areas
		envVariables[[6]] = land_covers[[t]][[7]] # secondary non-forest areas
		envVariables[[7]] = land_covers[[t]][[1]] # croplands (all catergories)
		envVariables[[8]] = land_covers[[t]][[2]] # managed pasture + rangeland
		pLog = populations[[t]]; pLog[] = log10(pLog[]+1); envVariables[[9]] = pLog
		for (i in 1:length(envVariables)) envVariables[[i]][is.na(mask[])] = NA
		for (i in 1:length(envVariables)) envVariables[[i]] = crop(envVariables[[i]], europe2, snap="out")
		for (i in 1:length(envVariables)) envVariables[[i]] = mask(envVariables[[i]], europe2)
		envVariables[[1]][] = envVariables[[1]][]-273.15 # conversion to Celcius degrees
		envVariables[[2]][] = envVariables[[2]][]*60*60*24 # conversion to kg/m2/day
		for (i in c(1,2,9)) envVariables[[i]][is.na(envVariables[[3]][])] = NA
		envVariables_list[[t]] = envVariables
		nullRaster = envVariables_list[[t]][[1]]; nullRaster[!is.na(nullRaster[])] = 1
		names(nullRaster) = "nullRaster"; nullRasters[[t]] = nullRaster
	}

# 2. Loading and ploting all the Bombus records by continent

observations_list = list()
for (t in 1:length(periods))
	{
		observations = list(); c = 0
		minYear = periods[[t]][1]; maxYear = periods[[t]][2]
		for (i in 1:dim(species)[1])
			{
				tab1 = read.csv(paste0(directory,"/",species[i,1],".csv"), header=T)
				tab2 = tab1[which((tab1[,"year"]>=minYear)&(tab1[,"year"]<=maxYear)),c("longitude","latitude")]
				observations[[i]] = tab2[which(!is.na(raster::extract(nullRasters[[t]],tab2))),]
				c = c+dim(observations[[i]])[1]
			}
		observations_list[[t]] = observations; # print(c)
		if (savingPlots == TRUE)
			{
				pdf(paste0("All_the_figures_&_SI/Bombus_observations_",t,"_NEW.pdf"), width=14, height=14)
				par(mfrow=c(6,8), oma=c(0,0,0,0), mar=c(0,1.5,1.5,0), lwd=0.4, col="gray30")
				for (i in 1:dim(species)[1])
					{
						plot(nullRasters[[t]], col="gray90", ann=F, legend=F, axes=F, box=F)
						plot(coastlines, add=T, col="gray50", lwd=1.5)
						points(observations[[i]], col="gray30", pch=3, cex=0.3, lwd=0.3)
						mtext(paste0("B. ",species[i,1]), side=3, line=-2, at=0, cex=0.75, col="gray30")
					}
				dev.off()
			}
	}

# 3. Defining the background area for each time period

backgrounds = list()
for (t in 1:length(periods))
	{
		species = read.csv("Bombus_species.csv", header=T, sep=",")
		nullRaster = envVariables_list[[t]][[1]]; nullRaster[!is.na(nullRaster[])] = 1
		names(nullRaster) = "nullRaster"; allObservationsOnTheContinent = c()
		for (i in 1:dim(species)[1])
			{
				allObservationsOnTheContinent = rbind(allObservationsOnTheContinent, observations_list[[t]][[i]])
			}
		backgroundCells = unique(raster::extract(nullRaster, allObservationsOnTheContinent, cellnumbers=T))
		background = nullRaster; background[!(1:length(background[]))%in%backgroundCells] = NA; backgrounds[[t]] = background	
	}

# 4. Boosted regression trees (BRT) analyses with standard or spatial cross-validation

samplingPtsMinDist = function(observations, minDist=500, nberOfPoints=5)
	{
		indices = rep(NA, nberOfPoints)
		selection_list = list(1:nrow(observations)) 
  		indices[1] = sample(1:dim(observations)[1], 1)
		dists = list(spDistsN1(as.matrix(observations), as.matrix(observations[indices[1],]), longlat=T))
		for (i in 2:nberOfPoints)
			{
    			selection = which(dists[[(i-1)]] > minDist)
    			if (length(selection) == 0)
    				{
    					stop("Restarts the function with a smaller minimum distance")
					}
    			selection_list[[i]] = selection
    			test = table(unlist(selection_list))
    			indices_minDist = as.numeric(names(which(test==i)))
    			indices[i] = sample(indices_minDist, 1)   
				dists[[i]] = spDistsN1(as.matrix(observations), as.matrix(observations[indices[i],]), longlat=T)
			}
		return(indices)
	}
foldSelection = function(observations, selectedPoints)
	{
		fold_selection = sapply(1:nrow(observations), function(i) which.min(spDistsN1(as.matrix(selectedPoints), as.matrix(observations[i,]), longlat=T)))
		return(fold_selection)
	}

newAnalyses = FALSE; spatialCrossValidation1 = TRUE; spatialCrossValidation2 = TRUE; all_data = list(); datas = list()
occurrence_data_summary = matrix(nrow=dim(species)[1], ncol=4); row.names(occurrence_data_summary) = species[,1]
colnames(occurrence_data_summary) = c("n_1901_1974","n_1901_1974_filtered","n_2000_2014","n_2000_2014_filtered")
if (newAnalyses == TRUE) { for (t in 1:length(periods)) { for (i in 1:dim(species)[1]) {
		minYear = periods[[t]][1]; maxYear = periods[[t]][2]
		rasters_stack = stack(envVariables_list[[t]])
		observations = observations_list[[t]][[i]]
		probas = values(backgrounds[[t]])[!is.na(values(backgrounds[[t]]))]; n = 1000
		if (n > sum(!is.na(values(backgrounds[[t]])))) n = sum(!is.na(values(backgrounds[[t]])))
		pseudo_absences = xyFromCell(backgrounds[[t]], sample(which(!is.na(values(backgrounds[[t]]))), n, prob=probas, replace=F))
		presences = cbind(observations, rep(1,dim(observations)[1]), raster::extract(rasters_stack, observations))
		absences = cbind(pseudo_absences, rep(0,dim(pseudo_absences)[1]), raster::extract(rasters_stack, pseudo_absences))
		colnames(absences)[1] = "longitude"; colnames(absences)[2] = "latitude"; colnames(absences)[3] = "response"
		colnames(presences)[3] = "response"; data = rbind(presences,absences); data_to_discard = c()
		for (j in 1:length(rasters_stack@layers))
			{
				data_to_discard = c(data_to_discard, which(is.na(raster::extract(rasters_stack[[j]],data[,1:2]))))
			}
		data_to_discard = data_to_discard[order(data_to_discard)]
		data = data[which(!c(1:dim(data)[1])%in%data_to_discard),]
		occurrence_data_summary[i,((t-1)*2)+1] = sum(data[,"response"])
		cellIDs = cellFromXY(rasters_stack[[1]], data[,1:2]); buffer = c()
		for (j in 1:length(unique(cellIDs)))
			{	# Keeping only one presence or pseudo-absence point per raster cell (priority = presence points):
				if (sum(cellIDs==unique(cellIDs)[j]) > 1)
					{
						tmp = data[which(cellIDs==unique(cellIDs)[j]),]
						if (sum(tmp[,"response"]==1) == 0)
							{
								buffer = rbind(buffer, tmp[sample(1:dim(tmp)[1],1),])
							}
						if (sum(tmp[,"response"]==1) == 1)
							{
								buffer = rbind(buffer, tmp[which(tmp[,"response"]==1),])
							}
						if (sum(tmp[,"response"]==1) >= 2)
							{
								indices = which(tmp[,"response"]==1)
								buffer = rbind(buffer, tmp[sample(indices,1),])
							}
					}	else	{
						buffer = rbind(buffer, data[which(cellIDs==unique(cellIDs)[j]),])
					}
			}
		data = buffer; datas[[i]] = data; all_data[[t]] = datas
		occurrence_data_summary[i,((t-1)*2)+2] = sum(data[,"response"]) # }}
		if (!file.exists(paste0("BRT_projection_files/BRT_models/B_",species[i,1],"_",t,"_CCV_SCV_AUCs.csv"))) {
		plottingCorrelogram = FALSE
		if (plottingCorrelogram == TRUE)
			{
				correlogram = ncf::correlog(data[,"longitude"], data[,"latitude"], data[,"response"], na.rm=T, increment=10, resamp=0, latlon=T)
				dev.new(width=4.5, height=3); par(mar=c(2.2,2.2,1.5,1.5))
				plot(correlogram$mean.of.class[-1], correlogram$correlation[-1], ann=F, axes=F, lwd=0.2, cex=0.5, col=NA, ylim=c(-0.4,1.0), xlim=c(0,8500))
				abline(h=0, lwd=0.5, col="red", lty=2)
				points(correlogram$mean.of.class[-1], correlogram$correlation[-1], lwd=0.2, cex=0.35, col="gray30")
				lines(correlogram$mean.of.class[-1], correlogram$correlation[-1], lwd=0.2, col="gray30")
				axis(side=1, pos=-0.4, lwd.tick=0.2, cex.axis=0.6, lwd=0.2, tck=-0.015, col.axis="gray30", mgp=c(0,-0.05,0), at=seq(0,9000,1000))
				axis(side=2, pos=0, lwd.tick=0.2, cex.axis=0.6, lwd=0.2, tck=-0.015, col.axis="gray30", mgp=c(0,0.18,0), at=seq(-0.4,1,0.2))
				title(xlab="distance (km2)", cex.lab=0.7, mgp=c(0.3,0,0), col.lab="gray30")
				title(ylab="correlation", cex.lab=0.7, mgp=c(0.4,0,0), col.lab="gray30")
			}
		theRanges = c(500,500)*1000 # distance in meters
		nberOfReplicates = 10 # one replicate = one folds partition
		gbm.x = names(rasters_stack)
		gbm.y = colnames(data)[3]
		offset = NULL
		tree.complexity = 5 # "tc" = number of nodes in the trees
		learning.rate = 0.005 # "lr" = contribution of each tree to the growing model
		bag.fraction = 0.80 # proportion of data used to train a given tree
		site.weights = rep(1, dim(data)[1])
		var.monotone = rep(0, length(gbm.x))
		n.folds = 5
		prev.stratify = TRUE
		family = "bernoulli"
		n.trees = 10 # initial number of trees
		step.size = 5 # interval at which the predictive deviance is computed and logged
					  # (at each interval, the folds are successively used as test data set
					  # nd the remaining folds as training data sets to compute the deviance)
		max.trees = 10000 # maximum number of trees that will be considered
		tolerance.method = "auto"
		tolerance = 0.001
		plot.main = TRUE
		plot.folds = FALSE
		verbose = TRUE
		silent = FALSE
		keep.fold.models = FALSE
		keep.fold.vector = FALSE
		keep.fold.fit = FALSE
		showingFoldsPlot = FALSE
		brt_model_ccvs = list() # classic cross-validations (CCVs)
		brt_model_scv1 = list() # spatial cross-validations 1 (SCV1)
		brt_model_scv2 = list() # spatial cross-validations 2 (SCV2)
		if (spatialCrossValidation2 == TRUE)	
			{
				AUCs = matrix(nrow=nberOfReplicates, ncol=3); colnames(AUCs) = c("CCV_AUC","SCV1_AUC","SCV2_AUC")
			}	else	{
				if (spatialCrossValidation2 == TRUE)
					{
						AUCs = matrix(nrow=nberOfReplicates, ncol=2); colnames(AUCs) = c("CCV_AUC","SCV_AUC")
					}	else	{
						AUCs = matrix(nrow=nberOfReplicates, ncol=1); colnames(AUCs) = c("AUC")
					}
			}
		for (j in 1:nberOfReplicates)
			{
				# BRT with classic (standard) cross-validation (CCV):
				pdf(file=paste0("BRT_projection_files/BRT_models/B_",species[i,1],"_",t,"_CCV_replicate_",j,".pdf"))
				n.trees = 10; learning.rate = 0.005; step.size = 5; fold.vector = NULL; worked = FALSE
				while (worked == FALSE)
					{
						trycatch = tryCatch(
							{
								brt_model_ccvs[[j]] = gbm.step(data, gbm.x, gbm.y, offset, fold.vector, tree.complexity, learning.rate, bag.fraction, site.weights,
									var.monotone, n.folds, prev.stratify, family, n.trees, step.size, max.trees, tolerance.method, tolerance, plot.main, plot.folds,
									verbose, silent, keep.fold.models, keep.fold.vector, keep.fold.fit); # summary(brt_model_scv) # gbm.plot(brt_model_scv, plot.layout=c(4,4))
							},	error = function(cond) {
							},	warning = function(cond) {
							},	finally = {
							})
						if (length(brt_model_ccvs) == j) worked = TRUE
					}
				dev.off()
				AUCs[j,1] = brt_model_ccvs[[j]]$cv.statistics$discrimination.mean # Mean test AUC (from the AUCs computed on each fold tested as test data in the CCV)
				object = brt_model_ccvs[[j]]; df = as.data.frame(rasters_stack)
				not_NA = which(!is.na(rowMeans(df))); newdata = df[not_NA,]
				n.trees = brt_model_ccvs[[j]]$gbm.call$best.trees; type = "response"; single.tree = FALSE
				prediction = predict.gbm(object, newdata, n.trees, type, single.tree)
				rast = rasters_stack[[1]]; rast[!is.na(rast[])] = prediction
				if (spatialCrossValidation1 == TRUE)
					{
						# BRT with spatial (geographic) cross-validation (SCV) based on the folds generation of Dhingra, Artois et al. (2016, eLife):
						folds_with_similar_sizes = FALSE; c = 0
						while (folds_with_similar_sizes == FALSE) # while loop to select a partition where the x folds gather at least
							{									  #  proportion = (1/(x+1)) of the total number of presence points
								data_presence = data[which(data[,3]==1),]; c = c+1; # print(c)
								fivePoints = samplingPtsMinDist(data_presence[,1:2], minDist=200, nberOfPoints=n.folds)
								fold.vector = foldSelection(data[,1:2], selectedPoints=data_presence[fivePoints,1:2])
								fold.vector_presences = fold.vector[which(data[,3]==1)]
								counts = hist(fold.vector_presences, plot=F)$counts
								props = counts[which(counts > 0)]/sum(counts); print(round(props,2))
								if (min(props) > (1/(n.folds*2))) folds_with_similar_sizes = TRUE
							}
						if (showingFoldsPlot == TRUE)
							{
								par(mar=c(0,0,0,0), oma=c(0.0,3.6,0.0,0.0), mgp=c(0,0.4,0), lwd=0.2, bty="o")
								cols = c("olivedrab3","tan3","steelblue3","orange1","tomato2","mediumseagreen")[fold.vector]
								plot(backgrounds[[t]], col="gray90", useRaster=T, colNA=NA, box=F, axes=F, legend=F)
								pchs = c(16,3)[data[,3]+1]; cexs = c(0.25,0.5)[data[,3]+1]
								points(data[,1:2], col=cols, pch=pchs, cex=cexs, lwd=0.7)
							}
						if (spatialCrossValidation2 == TRUE) pdf(file=paste0("BRT_projection_files/BRT_models/B_",species[i,1],"_",t,"_SCV1_replicate_",j,".pdf"))
						if (spatialCrossValidation2 == FALSE) pdf(file=paste0("BRT_projection_files/BRT_models/B_",species[i,1],"_",t,"_SCV_replicate_",j,".pdf"))
						n.trees = 10; learning.rate = 0.001; step.size = 2; worked = FALSE
						while (worked == FALSE)
							{
								trycatch = tryCatch(
									{
										brt_model_scv1[[j]] = gbm.step(data, gbm.x, gbm.y, offset, fold.vector, tree.complexity, learning.rate, bag.fraction, site.weights,
											var.monotone, n.folds, prev.stratify, family, n.trees, step.size, max.trees, tolerance.method, tolerance, plot.main, plot.folds,
											verbose, silent, keep.fold.models, keep.fold.vector, keep.fold.fit)
									},	error = function(cond) {
									},	warning = function(cond) {
									},	finally = {
									})
								if (length(brt_model_scv1) == j) worked = TRUE
							}
						dev.off()
						if (spatialCrossValidation2 == TRUE) AUCs[j,"SCV1_AUC"] = brt_model_scv1[[j]]$cv.statistics$discrimination.mean
						if (spatialCrossValidation2 == FALSE) AUCs[j,"SCV_AUC"] = brt_model_scv1[[j]]$cv.statistics$discrimination.mean
							# Mean test AUC (from the AUCs computed on each fold tested as test data in the SCV)		
						object = brt_model_scv1[[j]]; df = as.data.frame(rasters_stack)
						not_NA = which(!is.na(rowMeans(df))); newdata = df[not_NA,]
						n.trees = brt_model_scv1[[j]]$gbm.call$best.trees; type = "response"; single.tree = FALSE
						prediction = predict.gbm(object, newdata, n.trees, type, single.tree)
						rast = rasters_stack[[1]]; rast[!is.na(rast[])] = prediction
					}
				if (spatialCrossValidation2 == TRUE)
					{
						# BRT with spatial (geographic) cross-validation (SCV) based on the blocks generation of Valavi et al. (2019, MEE):
						spdf = SpatialPointsDataFrame(data[c("longitude","latitude")], data[,3:dim(data)[2]], proj4string=crs(nullRasters[[t]]))
						worked = FALSE
						while (worked == FALSE)
							{
								trycatch = tryCatch(
									{
										myblocks = NULL
										myblocks = spatialBlock(spdf, species="response", rasterLayer=nullRasters[[t]], k=n.folds, theRange=theRanges[1], selection="random")
									},	error = function(cond) {
									},	finally = {
									})
								if (!is.null(myblocks)) worked = TRUE
							}
						fold.vector = myblocks$foldID; n.trees = 100
						pdf(file=paste0("BRT_projection_files/BRT_models/B_",species[i,1],"_",t,"_SCV2_replicate_",j,".pdf"))
						brt_model_scv2[[j]] = gbm.step(data, gbm.x, gbm.y, offset, fold.vector, tree.complexity, learning.rate, bag.fraction, site.weights,
							var.monotone, n.folds, prev.stratify, family, n.trees, step.size, max.trees, tolerance.method, tolerance, plot.main, plot.folds,
							verbose, silent, keep.fold.models, keep.fold.vector, keep.fold.fit); # summary(brt_model_scv) # gbm.plot(brt_model_scv, plot.layout=c(4,4))
						dev.off()
						AUCs[j,"SCV2_AUC"] = brt_model_scv2[[j]]$cv.statistics$discrimination.mean
							# Mean test AUC (from the AUCs computed on each fold tested as test data in the SCV)
						object = brt_model_scv2[[j]]; df = as.data.frame(rasters_stack)
						not_NA = which(!is.na(rowMeans(df))); newdata = df[not_NA,]
						n.trees = brt_model_scv2[[j]]$gbm.call$best.trees; type = "response"; single.tree = FALSE
						prediction = predict.gbm(object, newdata, n.trees, type, single.tree)
						rast = rasters_stack[[1]]; rast[!is.na(rast[])] = prediction
					}
			}
		saveRDS(brt_model_ccvs, paste0("BRT_projection_files/BRT_models/B_",species[i,1],"_",t,"_models_CCV.rds"))
		if (spatialCrossValidation1 == TRUE)	
			{
				if (spatialCrossValidation2 == TRUE)	
					{
						saveRDS(brt_model_scv1, paste0("BRT_projection_files/BRT_models/B_",species[i,1],"_",t,"_models_SCV1.rds"))
					}	else		{
						saveRDS(brt_model_scv1, paste0("BRT_projection_files/BRT_models/B_",species[i,1],"_",t,"_models_SCV.rds"))
					}
			}
		if (spatialCrossValidation2 == TRUE) saveRDS(brt_model_scv2, paste0("BRT_projection_files/BRT_models/B_",species[i,1],"_",t,"_models_SCV2.rds"))
		if (spatialCrossValidation1 == TRUE)	
			{
				write.csv(AUCs, paste0("BRT_projection_files/BRT_models/B_",species[i,1],"_",t,"_CCV_SCV_AUCs.csv"), row.names=F, quote=F)
			}	else	{
				write.csv(AUCs, paste0("BRT_projection_files/BRT_models/B_",species[i,1],"_",t,"_CCV_AUC_vals.csv"), row.names=F, quote=F)
			}
	}}}}
if (!file.exists(paste0("Occurrence_data.csv")))
	{
		write.csv(occurrence_data_summary, "Occurrence_data.csv", quote=F)
	}
AUC_values = matrix(nrow=dim(species)[1], ncol=6); row.names(AUC_values) = species[,"species"]
colnames(AUC_values) = c("CCV_t1","SCV1_t1","SCV2_t1","CCV_t2","SCV1_t2","SCV2_t2")
for (t in 1:length(periods))
	{
		for (i in 1:dim(species)[1])
			{
				tab = read.csv(paste0("BRT_projection_files/BRT_models/B_",species[i,"species"],"_",t,"_CCV_SCV_AUCs.csv"), head=T)
				for (j in 1:dim(tab)[2])
					{
						AUC_values[i,((t-1)*3)+j] = paste0(round(mean(tab[,j]),3)," (",round(sd(tab[,j]),3),")")
					}
			}
	}
if (!file.exists(paste0("All_AUC_values.csv")))
	{
		write.csv(AUC_values, "All_AUC_values.csv", quote=F)
	}
AUC_values = read.csv("All_AUC_values.csv", head=T)

# 5. Computation of the prevalence-pseudoabsence-calibrated Sørensen index

	# Sources:
		# - computation performed according to the formulas of Leroi et al. (2018, J. Biogeography)
		# - optimisation of the threshold with a 0.01 step increment according to Li & Guo (2013, Ecography)

if (!file.exists(paste0("All_SIppc_values.csv")))
	{
		tab = matrix(nrow=dim(species)[1], ncol=4); row.names(tab) = species[,1]; tabs_list1 = list()
		colnames(tab) = c("ppac_sorensenIndex_t1","optimisedThreshold_t1","ppac_sorensenIndex_t2","optimisedThreshold_t2")
		for (t in 1:length(periods))
			{
				rasters_stack = stack(envVariables_list[[t]]); tabs_list2 = list()
				background_cells = sum(!is.na(envVariables_list[[t]][[1]][]))
				for (i in 1:dim(species)[1])
					{
						brt_model_scv2 = readRDS(paste0("BRT_projection_files/BRT_models/B_",species[i,1],"_",t,"_models_SCV2.rds"))
						sorensen_ppcs = rep(NA, length(brt_model_scv2)); thresholds = rep(NA, length(brt_model_scv2)); tabs_list3 = list()
						for (j in 1:length(brt_model_scv2))
							{
								tmp = matrix(nrow=101, ncol=2); tmp[,1] = seq(0,1,0.01)
								df = brt_model_scv2[[j]]$gbm.call$dataframe
								responses = df$response; data = df[,4:dim(df)[2]]
								n.trees = brt_model_scv2[[j]]$gbm.call$best.trees; type = "response"; single.tree = FALSE
								prediction = predict.gbm(brt_model_scv2[[j]], data, n.trees, type, single.tree)		
								N = background_cells
								P = sum(responses==1); A = sum(responses==0)
								prev = P/(P+A) # proportion of recorded sites where the species is present
								x = (P/A)*((1-prev)/prev)
								sorensen_ppc = 0
								for (threshold in seq(0,1,0.01))
									{
										TP = length(which((responses==1)&(prediction>=threshold))) # true positives
										FN = length(which((responses==1)&(prediction<threshold))) # false negatives
										# UPR = FN/(TP+FN); FFP_total = (1-UPR)*(prev*N-P)
										# FP_sample = length(which((responses==0)&(prediction>=threshold)))
										# FP_prim = (FP_sample/A)*(N-P)-FFP_total
										# FP_pa = A*(FP_prim/((1-prev)*N))
										FP_pa = length(which((responses==0)&(prediction>=threshold))) # false positives
										sorensen_ppc_tmp = (2*TP)/((2*TP)+(x*FP_pa)+(FN))
										tmp[which(tmp[,1]==threshold),2] = sorensen_ppc_tmp
										if (sorensen_ppc < sorensen_ppc_tmp)
											{
												sorensen_ppc = sorensen_ppc_tmp
												optimised_threshold = threshold
											}
									}
								tabs_list3[[j]] = tmp
								sorensen_ppcs[j] = sorensen_ppc
								thresholds[j] = optimised_threshold
							}
						tabs_list2[[i]] = tabs_list3
						tab[i,((t-1)*2)+1] = paste0(round(median(sorensen_ppcs),2)," [",round(min(sorensen_ppcs),2),"-",round(max(sorensen_ppcs),2),"]")
						tab[i,((t-1)*2)+2] = paste0(round(median(thresholds),2)," [",round(min(thresholds),2),"-",round(max(thresholds),2),"]")
					}
				tabs_list1[[t]] = tabs_list2
			}
		write.csv(tab, "All_SIppc_values.csv", quote=F)		
		pdf(paste0("All_the_figures_&_SI/SI_ppc_species_curves_NEW.pdf"), width=10, height=12)
		par(mfrow=c(8,6), oma=c(0,0,0,0), mar=c(2.5,2.5,0.5,0.5), lwd=0.4, col="gray30")
		for (i in 1:dim(species)[1])
			{
				plot(tabs_list1[[1]][[i]][[1]], col=NA, ann=F, axes=F, xlim=c(0,1), ylim=c(0,1))
				for (j in 1:length(tabs_list1[[1]][[i]])) lines(tabs_list1[[1]][[i]][[j]], lwd=0.3, col="gray80", lty=1)
				for (j in 1:length(tabs_list1[[2]][[i]])) lines(tabs_list1[[2]][[i]][[j]], lwd=0.3, col="gray30", lty=1)
				axis(side=1, lwd.tick=0.2, cex.axis=0.7, lwd=0, tck=-0.030, col.axis="gray30", mgp=c(0,0.07,0))
				axis(side=2, lwd.tick=0.2, cex.axis=0.7, lwd=0, tck=-0.030, col.axis="gray30", mgp=c(0,0.30,0))
				if (i %in% c(1,7,13,19,25,31,37,43)) title(ylab=expression("SI"["ppc"]), cex.lab=0.9, mgp=c(1.3,0,0), col.lab="gray30")
				if (i %in% c(41,42,43,44,45,46)) title(xlab="threshold", cex.lab=0.9, mgp=c(1.1,0,0), col.lab="gray30")
				box(lwd=0.2, col="gray30"); mtext(paste0("B. ",species[i,1]), side=3, line=-1.3, at=0.98, cex=0.55, col="gray30", adj=1)
			}
		dev.off()
	}

# 6. Comparison of the two BRT models trained for the two time periods

if (!file.exists(paste0("BRT_comparison.csv")))
	{
		tab = matrix(nrow=dim(species)[1], ncol=4); row.names(tab) = species[,1]
		colnames(tab) = c("ppac_sorensenIndex_t1_model_projected_on_t2","optimisedThreshold_t1_model_projected_on_t2",
						  "ppac_sorensenIndex_t2_model_projected_on_t1","optimisedThreshold_t2_model_projected_on_t1")
		for (t in 1:length(periods))
			{
				rasters_stack = stack(envVariables_list[[t]])
				background_cells = sum(!is.na(envVariables_list[[t]][[1]][]))
				for (i in 1:dim(species)[1])
					{
						brt_model_scv2_1 = readRDS(paste0("BRT_projection_files/BRT_models/B_",species[i,1],"_",t,"_models_SCV2.rds"))
						if (t == 1) brt_model_scv2_2 = readRDS(paste0("BRT_projection_files/BRT_models/B_",species[i,1],"_2_models_SCV2.rds"))
						if (t == 2) brt_model_scv2_2 = readRDS(paste0("BRT_projection_files/BRT_models/B_",species[i,1],"_1_models_SCV2.rds"))
						sorensen_ppcs = rep(NA, length(brt_model_scv2_1)); thresholds = rep(NA, length(brt_model_scv2_1))
						for (j in 1:length(brt_model_scv2_1))
							{
								df = brt_model_scv2_2[[1]]$gbm.call$dataframe
								responses = df$response; data = df[,4:dim(df)[2]]
								n.trees = brt_model_scv2_1[[j]]$gbm.call$best.trees; type = "response"; single.tree = FALSE
								prediction = predict.gbm(brt_model_scv2_1[[j]], data, n.trees, type, single.tree)		
								N = background_cells
								P = sum(responses==1); A = sum(responses==0)
								prev = P/(P+A)
								x = (P/A)*((1-prev)/prev)
								sorensen_ppc = 0
								for (threshold in seq(0,1,0.01))
									{
										TP = length(which((responses==1)&(prediction>=threshold))) # true positives
										FN = length(which((responses==1)&(prediction<threshold))) # false negatives
										FP_pa = length(which((responses==0)&(prediction>=threshold))) # false positives
										sorensen_ppc_tmp = (2*TP)/((2*TP)+(x*FP_pa)+(FN))
										if (sorensen_ppc < sorensen_ppc_tmp)
											{
												sorensen_ppc = sorensen_ppc_tmp
												optimised_threshold = threshold
											}
									}
								sorensen_ppcs[j] = sorensen_ppc; thresholds[j] = optimised_threshold
							}
						tab[i,((t-1)*2)+1] = paste0(round(median(sorensen_ppcs),2)," [",round(min(sorensen_ppcs),2),"-",round(max(sorensen_ppcs),2),"]")
						tab[i,((t-1)*2)+2] = paste0(round(median(thresholds),2)," [",round(min(thresholds),2),"-",round(max(thresholds),2),"]")
					}
			}
		write.csv(tab, "BRT_comparison.csv", quote=F)
		tab = read.csv("All_SIppc_values.csv", head=T)
		v1as = rep(NA, dim(tab)[1]); v2as = rep(NA, dim(tab)[1])
		for (i in 1:length(v1as))
			{
				v1as[i] = as.numeric(unlist(strsplit(tab[i,"ppac_sorensenIndex_t1"]," \\["))[1])
				v2as[i] = as.numeric(unlist(strsplit(tab[i,"ppac_sorensenIndex_t2"]," \\["))[1])
			}
		tab = read.csv("BRT_comparison.csv", head=T)
		v1bs = rep(NA, dim(tab)[1]); v2bs = rep(NA, dim(tab)[1])
		for (i in 1:length(v1bs))
			{
				v1bs[i] = as.numeric(unlist(strsplit(tab[i,"ppac_sorensenIndex_t1_model_projected_on_t2"]," \\["))[1])
				v2bs[i] = as.numeric(unlist(strsplit(tab[i,"ppac_sorensenIndex_t2_model_projected_on_t1"]," \\["))[1])
			}
		diffs1 = v1as-v1bs; cat(round(median(diffs1),2),", 95% CI = [",round(quantile(diffs1,0.025),2),", ",round(quantile(diffs1,0.975),2),"]",sep="") # 0.19, 95% CI = [-0.03, 0.44]
		diffs2 = v2as-v2bs; cat(round(median(diffs2),2),", 95% CI = [",round(quantile(diffs2,0.025),2),", ",round(quantile(diffs2,0.975),2),"]",sep="") # 0.09, 95% CI = [-0.09, 0.31]	
	}

# 7. Estimation of spatial sorting bias (SSB)

	# SSB = Dp/Da (Hijsmans 2012, Ecology), where:
		# Dp = mean distance between testing presence sites and nearest training-presence sites
		# Da = mean distance between testing absence sites and nearest training-presence sites
		# --> SSB = 1 suggests there is no spatial sorting bias
		# --> SSB = 0 suggests extreme spatial sorting bias

newAnalyses = FALSE
if (newAnalyses == TRUE) { for (t in 1:length(periods))
	{
		n.folds = 5; theRanges = c(500,500)*1000
		SSBs = matrix(ncol=n.folds, nrow=dim(species)[1])
		row.names(SSBs) = species[,"species"]
		colnames(SSBs) = c("fold1","fold2","fold3","fold4","fold5")
		for (i in 1:dim(species)[1])
			{
				data = all_data[[t]][[i]]
				fold.vector = rep(NA, dim(data)[1])
				for (j in 1:dim(data)[1]) fold.vector[j] = sample(1:n.folds,1)
				for (j in 1:n.folds)
					{
						p = data[which((data[,"response"]==1)&(fold.vector!=j)),1:2]
						a = data[which((data[,"response"]==0)&(fold.vector!=j)),1:2]
						reference = data[which((data[,"response"]==1)&(fold.vector==j)),1:2]
						SSB = ssb(p, a, reference); SSBs[i,j] = SSB[1,"p"]/SSB[1,"a"]
					}
			}
		write.csv(round(SSBs,2), paste0("BRT_projection_files/BRT_models/SSBs_CCV_for_time_period_",t,".csv"), quote=F)
		SSBs = matrix(ncol=n.folds, nrow=dim(species)[1])
		row.names(SSBs) = species[,"species"]
		for (i in 1:dim(species)[1])
			{
				colnames(SSBs) = c("fold1","fold2","fold3","fold4","fold5")
				folds_with_similar_sizes = FALSE; c = 0
				while (folds_with_similar_sizes == FALSE)
					{
						data_presence = data[which(data[,3]==1),]; c = c+1
						fivePoints = samplingPtsMinDist(data_presence[,1:2], minDist=200, nberOfPoints=n.folds)
						fold.vector = foldSelection(data[,1:2], selectedPoints=data_presence[fivePoints,1:2])
						fold.vector_presences = fold.vector[which(data[,3]==1)]
						counts = hist(fold.vector_presences, plot=F)$counts
						props = counts[which(counts > 0)]/sum(counts); print(round(props,2))
						if (min(props) > (1/(n.folds*2))) folds_with_similar_sizes = TRUE
					}
				for (j in 1:n.folds)
					{
						p = data[which((data[,"response"]==1)&(fold.vector!=j)),1:2]
						a = data[which((data[,"response"]==0)&(fold.vector!=j)),1:2]
						reference = data[which((data[,"response"]==1)&(fold.vector==j)),1:2]
						SSB = ssb(p, a, reference); SSBs[i,j] = SSB[1,"p"]/SSB[1,"a"]
					}
			}
		write.csv(round(SSBs,2), paste0("BRT_projection_files/BRT_models/SSBs_SCV1_for_time_period_",t,".csv"), quote=F)
		SSBs = matrix(ncol=n.folds, nrow=dim(species)[1])
		row.names(SSBs) = species[,"species"]
		colnames(SSBs) = c("fold1","fold2","fold3","fold4","fold5")
		for (i in 1:dim(species)[1])
			{
				folds_with_similar_sizes = FALSE; c = 0
				spdf = SpatialPointsDataFrame(data[c("longitude","latitude")], data[,3:dim(data)[2]], proj4string=crs(nullRaster))
				worked = FALSE
				while (worked == FALSE)
					{
						trycatch = tryCatch(
							{
								myblocks = NULL
								myblocks = spatialBlock(spdf, species="response", rasterLayer=nullRasters[[t]], k=n.folds, theRange=theRanges[1], selection="random")
							},	error = function(cond) {
							},	finally = {
							})
						if (!is.null(myblocks)) worked = TRUE
					}
				for (j in 1:n.folds)
					{
						fold.vector = myblocks$foldID
						p = data[which((data[,"response"]==1)&(fold.vector!=j)),1:2]
						a = data[which((data[,"response"]==0)&(fold.vector!=j)),1:2]
						reference = data[which((data[,"response"]==1)&(fold.vector==j)),1:2]
						if (dim(reference)[1]>0)
							{
								SSB = ssb(p, a, reference); SSBs[i,j] = SSB[1,"p"]/SSB[1,"a"]
							}
					}
			}
		write.csv(round(SSBs,2), paste0("BRT_projection_files/BRT_models/SSBs_SCV2_for_time_period_",t,".csv"), quote=F)
	}}

# 8. Comparison of the response curves for each environmental factor

envVariableNames1 = rep(NA, length(envVariables))
for (i in 1:length(envVariables)) envVariableNames1[i] = names(envVariables[[i]])
envVariableNames2 = c("Temperature","Precipitation","Primary forest areas","Primary non-forest areas",
					  "Secondary forest areas","Secondary non-forest areas","Croplands","Pastures","Human population")
newAnalyses = FALSE
if (newAnalyses == TRUE)
	{
		selectedModel = "SCV2"
		envVariableValues_list1 = list()
		for (t in 1:length(periods))
			{
				envVariableValues_list2 = list()
				for (i in 1:dim(species)[1])
					{
						data = all_data[[t]][[i]]; data = data[which(data[,"response"]==1),]
						envVariableValues = matrix(nrow=3, ncol=length(envVariables))
						row.names(envVariableValues) = c("median","minV","maxV")
						colnames(envVariableValues) = envVariableNames1
						for (j in 1:length(envVariables))
							{
								minV = min(data[,envVariableNames1[j]], na.rm=T)
								maxV = max(data[,envVariableNames1[j]], na.rm=T)
								medianV = median(data[,envVariableNames1[j]], na.rm=T)
								envVariableValues[,j] = cbind(medianV, minV, maxV)
							}
						envVariableValues_list2[[i]] = envVariableValues
					}
				envVariableValues_list1[[t]] = envVariableValues_list2
			}
		for (tt in 1:3)
			{
				pdf(paste0("All_the_figures_&_SI/All_response_curves_",tt,"_NEW.pdf"), width=7.5, height=4)
				par(mfrow=c(3,3), oma=c(1,1,1,1), mar=c(2,1.3,0.5,0.5), lwd=0.2, col="gray30"); t = tt
				if (tt == 3) t = tt-1
				for (i in 1:length(envVariables))
					{
						predictions = list(); dfs = list()
						for (j in 1:dim(species)[1])
							{
								valuesInterval = 0.1; valuesInterval = (envVariableValues_list1[[t]][[j]]["maxV",i]-envVariableValues_list1[[t]][[j]]["minV",i])/100
								df = data.frame(matrix(nrow=length(seq(envVariableValues_list1[[t]][[j]]["minV",i],envVariableValues_list1[[t]][[j]]["maxV",i],valuesInterval)),
													   ncol=length(envVariables))); colnames(df) = envVariableNames1
								for (k in 1:length(envVariables))
									{
										valuesInterval = 0.1; valuesInterval = (envVariableValues_list1[[t]][[j]]["maxV",k]-envVariableValues_list1[[t]][[j]]["minV",k])/100
										if (i == k) df[,envVariableNames1[k]] = seq(envVariableValues_list1[[t]][[j]]["minV",k],envVariableValues_list1[[t]][[j]]["maxV",k],valuesInterval)
										if (i != k) df[,envVariableNames1[k]] = rep(envVariableValues_list1[[t]][[j]]["median",k],dim(df)[1])
									}
								dfs[[j]] = df
								brt_model = readRDS(paste0("BRT_projection_files/BRT_models/B_",species[j,"species"],"_",t,"_models_",selectedModel,".rds"))
								AUC_values = read.csv(paste0("BRT_projection_files/BRT_models/B_",species[j,"species"],"_",t,"_CCV_SCV_AUCs.csv"))[,paste0(selectedModel,"_AUC")]
								index = which(AUC_values==max(AUC_values))[1]
								n.trees = brt_model[[index]]$gbm.call$best.trees; type = "response"; single.tree = FALSE
								prediction = predict.gbm(brt_model[[index]], newdata=df, n.trees, type, single.tree)
								if (j == 1)
									{
										minX = min(df[,envVariableNames1[i]]); maxX = max(df[,envVariableNames1[i]])
										minY = min(prediction); maxY = max(prediction)
									}	else	{
										if (minX > min(df[,envVariableNames1[i]])) minX = min(df[,envVariableNames1[i]])
										if (maxX < max(df[,envVariableNames1[i]])) maxX = max(df[,envVariableNames1[i]])
										if (minY > min(prediction)) minY = min(prediction)
										if (maxY < max(prediction)) maxY = max(prediction)
									}
								predictions[[j]] = prediction
							}
						if ((tt <= 2)|((tt == 3)&(!file.exists("ESD_2070SSP3t0.csv"))))
							{
								for (j in 1:dim(species)[1])
									{
										if (j == 1)
											{
												plot(dfs[[j]][,envVariableNames1[i]],predictions[[j]],col=cols[t],ann=F,axes=F,lwd=0.2,type="l",xlim=c(minX,maxX),ylim=c(minY,maxY))
											}	else	{
												lines(dfs[[j]][,envVariableNames1[i]],predictions[[j]],col=cols[t],lwd=0.2)
											}
									}
							}	else		{
								ESD_2070_SSP3_t0 = read.csv("ESD_2070SSP3t0.csv", head=T)
								for (j in 1:dim(species)[1])
									{
										colour = ESD_2070_SSP3_t0[j,"colour"]
										if (j == 1)
											{
												plot(dfs[[j]][,envVariableNames1[i]],predictions[[j]],col=colour,ann=F,axes=F,lwd=0.2,type="l",xlim=c(minX,maxX),ylim=c(minY,maxY))
											}	else	{
												lines(dfs[[j]][,envVariableNames1[i]],predictions[[j]],col=colour,lwd=0.2)
											}
									}
							}
						axis(side=1, lwd.tick=0.2, cex.axis=0.6, lwd=0, tck=-0.030, col.axis="gray30", mgp=c(0,0.07,0))
						axis(side=2, lwd.tick=0.2, cex.axis=0.6, lwd=0, tck=-0.030, col.axis="gray30", mgp=c(0,0.30,0))
						title(ylab="predicted values", cex.lab=0.8, mgp=c(1.3,0,0), col.lab="gray30")
						title(xlab=envVariableNames2[i], cex.lab=0.8, mgp=c(0.9,0,0), col.lab="gray30")
						box(lwd=0.2, col="gray30")
					}
				dev.off()
			}
	}

# 9. Analyses of the relative influence of each environmental factor

	# 9.1. Computation of the relative influence of each environmental factor

selectedModel = "SCV2"
for (t in 1:length(periods))
	{
		if (t == 1) fileName = "RI_BRT_1901-74.csv"
		if (t == 2) fileName = "RI_BRT_2000-14.csv"
		if (!file.exists(paste0(fileName)))
			{
				relativeInfluences = matrix(0, nrow=dim(species)[1], ncol=length(envVariables))
				row.names(relativeInfluences) = species[,"species"]
				for (i in 1:dim(species)[1])
					{
						brt_models = readRDS(paste0("BRT_projection_files/BRT_models/B_",species[i,"species"],"_",t,"_models_",selectedModel,".rds"))
						for (j in 1:length(brt_models))
							{
								if ((i == 1)&(j == 1)) envVariableNames1 = rep(NA, length(envVariables))
								for (k in 1:length(envVariables))
									{
										if ((i == 1)&(j == 1)) envVariableNames1[k] = names(envVariables[[k]])
										relativeInfluences[i,k] = relativeInfluences[i,k] + summary(brt_models[[j]])[names(envVariables[[k]]),"rel.inf"]
									}
							}
						if (i == 1) colnames(relativeInfluences) = envVariableNames1
						relativeInfluences[i,] = relativeInfluences[i,]/length(brt_models)
					}
				write.table(round(relativeInfluences,2), fileName, quote=F, sep=",")
			}
	}
if (!file.exists("RIs_comparison.csv"))
	{
		fileNames = c("RI_BRT_1901-74.csv","RI_BRT_2000-14.csv")
		meanRelativeInfluences = matrix(0, nrow=length(envVariables), ncol=3)
		colnames(meanRelativeInfluences) = c("RI_1901-74", "RI_2000-14", "mean_abs_dif")
		for (t in 1:2)
			{
				relativeInfluences = read.csv(fileNames[t], header=T)
				if (t == 1)
					{
						row.names(meanRelativeInfluences) = colnames(relativeInfluences)
					}
				for (i in 1:dim(relativeInfluences)[2])
					{
						RIs = relativeInfluences[,i]
						median = round(median(RIs),1)
						HPD = round(HDInterval::hdi(RIs)[1:2],1)
						CR = round(quantile(RIs,c(0.025,0.975)),1)
						meanRelativeInfluences[i,t] = paste0(median," [",CR[1],"-",CR[2],"]")
					}
			}
		for (i in 1:dim(meanRelativeInfluences)[1])
			{
				diffs = rep(NA, dim(relativeInfluences)[1])
				for (j in 1:dim(relativeInfluences)[1])
					{
						relativeInfluence1 = read.csv(fileNames[1], header=T)[j,i]
						relativeInfluence2 = read.csv(fileNames[2], header=T)[j,i]
						diffs[j] = abs(relativeInfluence1-relativeInfluence2)
					}
				median = round(median(diffs),1)
				HPD = round(HDInterval::hdi(diffs)[1:2],1)
				CR = round(quantile(RIs,c(0.025,0.975)),1)
				meanRelativeInfluences[i,3] = paste0(median," [",CR[1],"-",CR[2],"]")
			}
		write.csv(meanRelativeInfluences, "RIs_comparison.csv", quote=F)
	}

	# 9.2. PCA and investigation of the relatioship with IUCN status
	
relativeInfluences = read.csv("RI_BRT_2000-14.csv", header=T)
IUCN_status = c("least_concern", "near_threatened", "vulnerable"); colours = rep(NA,dim(relativeInfluences)[1])
IUCN_colours = c("chartreuse3", rgb(250,165,33,255,maxColorValue=255), rgb(222,67,39,255,maxColorValue=255))
pca = dudi.pca(relativeInfluences, scannf=F, nf=9); lis = pca$li[,1:2]; cos = pca$co
for (i in 1:dim(lis)[1])
	{
		index = which(species[,"species"]==row.names(lis)[i])
		status = species[index,"IUCN_status"]
		colours[i] = IUCN_colours[which(IUCN_status==status)]
	}
if (savingPlots == TRUE)
	{
		pdf("All_the_figures_&_SI/BRT_relative_influences_1.pdf", width=4.5, height=4); par(mar=c(3,3,0,1.5))
		plot(lis, col=colours, cex=0.7, pch=16, ann=F, axes=F, xlim=c(-4.0,3.5), ylim=c(-4.0,3.5))
		# s.corcircle(2*cos,xax=1,yax=2,box=F,sub="",csub=1,clabel=1.5,possub="topleft",grid=T,cgrid=1,full=T,add.plot=T)
		text(lis[,1], lis[,2], labels=gsub("B_","",row.names(lis)), cex=0.6, col="gray50", pos=4, offset=0.25)
		points(lis, col=colours, cex=0.9, pch=16); points(lis, col="gray30", cex=0.9, pch=1, lwd=0.3)
		axis(side=1, lwd.tick=0.2, cex.axis=0.6, lwd=0.2, tck=-0.017, col.axis="gray30", mgp=c(0,0.05,0), at=seq(-8,5,1))
		axis(side=2, lwd.tick=0.2, cex.axis=0.6, lwd=0.2, tck=-0.017, col.axis="gray30", mgp=c(0,0.30,0), at=seq(-5,9,1))
		title(xlab="PCA axis 1", cex.lab=0.7, mgp=c(1.0,0,0), col.lab="gray30")
		title(ylab="PCA axis 2", cex.lab=0.7, mgp=c(1.3,0,0), col.lab="gray30")
		dev.off()
	}

# 10. BRT projections based on past and present environmental variables

predictions1a = list() # past projections based on past data, and present projections based on present data
predictions1b = list() # past projections based on present data, and present projections based on past data
selectedModel = "SCV2" # selectedModel = "CCV"
for (i in 1:dim(species)[1])
	{
		predictions2 = list()
		for (t in 1:length(periods))
			{
				rasters_stack = stack(envVariables_list[[t]]); replicates = list()
				brt_model = readRDS(paste0("BRT_projection_files/BRT_models/B_",species[i,"species"],"_",t,"_models_",selectedModel,".rds"))
				for (j in 1:length(brt_model))
					{
						df = as.data.frame(rasters_stack); not_NA = which(!is.na(rowMeans(df))); newdata = df[not_NA,]
						n.trees = brt_model[[j]]$gbm.call$best.trees; type = "response"; single.tree = FALSE
						prediction = predict.gbm(brt_model[[j]], newdata, n.trees, type, single.tree)
						rast = rasters_stack[[1]]; rast[!is.na(rast[])] = prediction; replicates[[j]] = rast
					}
				rasts = stack(replicates); predictions2[[t]] = mean(rasts)
			}
		predictions2[[3]] = predictions2[[2]]-predictions2[[1]]
		for (t in 1:3)
			{
				if (t < 3)
					{
						speciesName = as.character(species[i,"species"])
						writeRaster(predictions2[[t]], paste0("BRT_projection_files/BRT_",selectedModel,"/B_",speciesName,"_",t,"_projection.asc"), overwrite=T)
					}
			}
		predictions1a[[i]] = predictions2
		predictions2 = list()
		for (t in 1:length(periods))
			{
				if (t == 1) tt = 2
				if (t == 2) tt = 1
				rasters_stack = stack(envVariables_list[[tt]]); replicates = list()
				brt_model = readRDS(paste0("BRT_projection_files/BRT_models/B_",species[i,"species"],"_",t,"_models_",selectedModel,".rds"))
				for (j in 1:length(brt_model))
					{
						df = as.data.frame(rasters_stack); not_NA = which(!is.na(rowMeans(df))); newdata = df[not_NA,]
						n.trees = brt_model[[j]]$gbm.call$best.trees; type = "response"; single.tree = FALSE
						prediction = predict.gbm(brt_model[[j]], newdata, n.trees, type, single.tree)
						rast = rasters_stack[[1]]; rast[!is.na(rast[])] = prediction; replicates[[j]] = rast
					}
				rasts = stack(replicates); predictions2[[t]] = mean(rasts)
			}
		predictions2[[3]] = predictions2[[2]]-predictions2[[1]]
		predictions1b[[i]] = predictions2
	}
if (savingPlots == TRUE)
	{
		useCutOff = TRUE; cutOff = 0.10 # not used anymore
		for (t in 1:length(periods))
			{
				cols = c(rep("#E5E5E5",10),rev(colorRampPalette(brewer.pal(11,"RdYlBu"))(120))[21:110])
				legend = raster(as.matrix(c(0,1))); # plot(legend, col=cols, ann=F, legend=T, axes=F, box=F)
				pdf(paste0("All_the_figures_&_SI/Projections_for_",selectedModel,"_",t,"_NEW.pdf"), width=14, height=14)
				par(mfrow=c(6,8), oma=c(0,0,0,0), mar=c(0,1.5,1.5,0), lwd=0.4, col="gray30")
				for (i in 1:dim(species)[1])
					{
						prediction = predictions1a[[i]][[t]]
						plot(prediction, col=cols[1:(max(prediction[],na.rm=T)*100)], ann=F, legend=F, axes=F, box=F)
						plot(coastlines, add=T, col="gray50", lwd=1.5)
						mtext(paste0("B. ",species[i,1]), side=3, line=-2, at=0, cex=0.75, col="gray30")
					}
				dev.off()
			}
	}

# 11. Post hoc analyses to compare past and present ENM projections

ESIs = list(); SRIs = list()
counterRasters1a = list(); bufferRasters1a = list()
counterRasters2a = list(); bufferRasters2a = list()
counterRasters1b = list(); bufferRasters1b = list()
counterRasters2b = list(); bufferRasters2b = list()
cutOff = 0.10 # a fixed cut-off value previously used
SIppcs = read.csv("All_SIppc_values.csv", head=T)
for (i in 1:dim(species)[1])
	{
		txt = SIppcs[i,"optimisedThreshold_t1"]
		cutOff = as.numeric(unlist(strsplit(txt," \\["))[1])
		bufferRasters1a[[i]] = predictions1a[[i]][[1]]; c1 = predictions1a[[i]][[1]]
		c1[c1[]<cutOff] = 0; c1[c1[]>cutOff] = 1; counterRasters1a[[i]] = c1
		txt = SIppcs[i,"optimisedThreshold_t2"]
		cutOff = as.numeric(unlist(strsplit(txt," \\["))[1])
		bufferRasters2a[[i]] = predictions1a[[i]][[2]]; c2 = predictions1a[[i]][[2]]
		c2[c2[]<cutOff] = 0; c2[c2[]>cutOff] = 1; counterRasters2a[[i]] = c2
		txt = SIppcs[i,"optimisedThreshold_t1"]
		cutOff = as.numeric(unlist(strsplit(txt," \\["))[1])
		bufferRasters1b[[i]] = predictions1b[[i]][[1]]; c1 = predictions1b[[i]][[1]]
		c1[c1[]<cutOff] = 0; c1[c1[]>cutOff] = 1; counterRasters1b[[i]] = c1
		txt = SIppcs[i,"optimisedThreshold_t2"]
		cutOff = as.numeric(unlist(strsplit(txt," \\["))[1])
		bufferRasters2b[[i]] = predictions1b[[i]][[2]]; c2 = predictions1b[[i]][[2]]
		c2[c2[]<cutOff] = 0; c2[c2[]>cutOff] = 1; counterRasters2b[[i]] = c2
	}
minNberOfRasterCells = length(bufferRasters1a[[1]][]); index = 1
for (i in 2:length(bufferRasters1a))
	{
		if (minNberOfRasterCells > length(bufferRasters1a[[i]][]))
			{
				minNberOfRasterCells > length(bufferRasters1a[[i]][]); index = i
			}
	}
for (i in 1:length(bufferRasters1a))
	{
		counterRasters1a[[i]] = crop(counterRasters1a[[i]], extent(counterRasters1a[[index]]))
		bufferRasters1a[[i]] = crop(bufferRasters1a[[i]], extent(bufferRasters1a[[index]]))
	}
for (i in 1:length(bufferRasters2a))
	{
		counterRasters2a[[i]] = crop(counterRasters2a[[i]], extent(counterRasters1a[[index]]))
		bufferRasters2a[[i]] = crop(bufferRasters2a[[i]], extent(bufferRasters1a[[index]]))
	}
for (i in 1:length(bufferRasters1b))
	{
		counterRasters1b[[i]] = crop(counterRasters1b[[i]], extent(counterRasters1b[[index]]))
		bufferRasters1b[[i]] = crop(bufferRasters1b[[i]], extent(bufferRasters1b[[index]]))
	}
for (i in 1:length(bufferRasters2b))
	{
		counterRasters2b[[i]] = crop(counterRasters2b[[i]], extent(counterRasters1b[[index]]))
		bufferRasters2b[[i]] = crop(bufferRasters2b[[i]], extent(bufferRasters1b[[index]]))
	}
bufferRasters1a = stack(bufferRasters1a); bufferRasters2a = stack(bufferRasters2a)
counterRasters1a = stack(counterRasters1a); counterRasters2a = stack(counterRasters2a)
bufferRaster1a = mean(bufferRasters1a); bufferRaster2a = mean(bufferRasters2a)
counterRaster1a = sum(counterRasters1a); counterRaster2a = sum(counterRasters2a)
bufferRasters1b = stack(bufferRasters1b); bufferRasters2b = stack(bufferRasters2b)
counterRasters1b = stack(counterRasters1b); counterRasters2b = stack(counterRasters2b)
bufferRaster1b = mean(bufferRasters1b); bufferRaster2b = mean(bufferRasters2b)
counterRaster1b = sum(counterRasters1b); counterRaster2b = sum(counterRasters2b)
ESIs[[1]] = bufferRaster1a; ESIs[[2]] = bufferRaster2a # ESI: ecological suitability index (mean ecological suitability among species)
SRIs[[1]] = counterRaster1a; SRIs[[2]] = counterRaster2a # SRI: species richness index
if (savingPlots == TRUE)
	{
		pdf(paste0("All_the_figures_&_SI/ESI_and_SRI_differences.pdf"), width=10, height=2.0)
		par(mfrow=c(1,6), oma=c(0,0,0,0), mar=c(0.1,0.75,0.1,2.5), lwd=0.4, col="gray30")
		bufferRaster3a = bufferRaster1a; bufferRaster3a[] = bufferRaster2a[]-bufferRaster1a[]
		counterRaster3a = counterRaster1a; counterRaster3a[] = counterRaster2a[]-counterRaster1a[]
		cols = rev(hcl.colors(100,palette="terrain2")[1:100])
		colsDiff = colorRampPalette(brewer.pal(11,"RdBu"))(101)
		maxV1 = max(c(max(bufferRaster1a[],na.rm=T),max(bufferRaster2a[],na.rm=T)))
		minV2 = min(bufferRaster3a[],na.rm=T); maxV2 = max(bufferRaster3a[],na.rm=T)
		if (abs(minV2) < abs(maxV2)) minV2 = -maxV2
		if (abs(maxV2) < abs(minV2)) maxV2 = -minV2
		legend1 = raster(as.matrix(c(0,maxV1))); legend2 = raster(as.matrix(c(minV2,maxV2)))
		plot(bufferRaster1a, col=cols[1:((max(bufferRaster1a[],na.rm=T)/maxV1)*100)], ann=F, legend=F, axes=F, box=F)
		plot(coastlines, add=T, col="gray50", lwd=1.0)
		mtext("1901 - 1974", side=3, line=-1.5, cex=0.75, col="gray30", at=1.5)
		mtext("ESI", side=3, line=-3.1, cex=1.0, col="gray30", at=-5)
		plot(legend1, col=cols, legend.only=T, add=T, legend.width=0.5, legend.shrink=0.3, 
			 smallplot=c(0.83,0.86,0.03,0.96), adj=3, axis.args=list(cex.axis=0.9, lwd=0, lwd.tick=0.4, col.tick="gray30", tck=-0.8, col="gray30", 
			 col.axis="gray30", line=0, mgp=c(0,0.5,0)), alpha=1, side=3, horizontal=F)
		plot(bufferRaster2a, col=cols[1:((max(bufferRaster2a[],na.rm=T)/maxV1)*100)], ann=F, legend=F, axes=F, box=F)
		plot(coastlines, add=T, col="gray50", lwd=1.0)
		mtext("2000 - 2014", side=3, line=-1.5, cex=0.75, col="gray30", at=1.5)
		mtext("ESI", side=3, line=-3.1, cex=1.0, col="gray30", at=-5)
		plot(legend1, col=cols, legend.only=T, add=T, legend.width=0.5, legend.shrink=0.3, 
			 smallplot=c(0.83,0.86,0.03,0.96), adj=3, axis.args=list(cex.axis=0.9, lwd=0, lwd.tick=0.4, col.tick="gray30", tck=-0.8, col="gray30", 
			 col.axis="gray30", line=0, mgp=c(0,0.5,0)), alpha=1, side=3, horizontal=F)
		index1 = round(((min(bufferRaster3a[],na.rm=T)-minV2)/(maxV2-minV2))*100)
		index2 = round(((max(bufferRaster3a[],na.rm=T)-minV2)/(maxV2-minV2))*100)
		plot(bufferRaster3a, col=colsDiff[index1:index2], ann=F, legend=F, axes=F, box=F)
		plot(coastlines, add=T, col="gray50", lwd=1.0)
		mtext("Difference", side=3, line=-1.5, cex=0.75, col="gray30", at=0)
		mtext("ESI", side=3, line=-3.1, cex=1.0, col="gray30", at=-5)
		plot(legend2, col=colsDiff, legend.only=T, add=T, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.83,0.86,0.03,0.96), adj=3,
			 axis.args=list(cex.axis=0.9, lwd=0, lwd.tick=0.4, col.tick="gray30", tck=-0.8, col="gray30", col.axis="gray30", line=0, 
			 mgp=c(0,0.5,0)), alpha=1, side=3, horizontal=F)		
		minV1 = min(c(min(counterRaster1a[],na.rm=T),min(counterRaster2a[],na.rm=T)))
		maxV1 = max(c(max(counterRaster1a[],na.rm=T),max(counterRaster2a[],na.rm=T)))
		minV2 = min(counterRaster3a[],na.rm=T); maxV2 = max(counterRaster3a[],na.rm=T)
		if (abs(minV2) < abs(maxV2)) minV2 = -maxV2
		if (abs(maxV2) < abs(minV2)) maxV2 = -minV2
		legend1 = raster(as.matrix(c(0,maxV1))); legend2 = raster(as.matrix(c(minV2,maxV2)))
		plot(counterRaster1a, col=cols[1:((max(counterRaster1a[],na.rm=T)/maxV1)*100)], ann=F, legend=F, axes=F, box=F)
		plot(coastlines, add=T, col="gray50", lwd=1.0)
		mtext("1901 - 1974", side=3, line=-1.5, cex=0.75, col="gray30", at=1.5)
		mtext("SRI", side=3, line=-3.1, cex=1.0, col="gray30", at=-5)
		plot(legend1, col=cols, legend.only=T, add=T, legend.width=0.5, legend.shrink=0.3, 
			 smallplot=c(0.83,0.86,0.03,0.96), adj=3, axis.args=list(cex.axis=0.9, lwd=0, lwd.tick=0.4, col.tick="gray30", tck=-0.8, col="gray30", 
			 col.axis="gray30", line=0, mgp=c(0,0.5,0)), alpha=1, side=3, horizontal=F)
		plot(counterRaster2a, col=cols[1:((max(counterRaster2a[],na.rm=T)/maxV1)*100)], ann=F, legend=F, axes=F, box=F)
		plot(coastlines, add=T, col="gray50", lwd=1.0)
		mtext("2000 - 2014", side=3, line=-1.5, cex=0.75, col="gray30", at=1.5)
		mtext("SRI", side=3, line=-3.1, cex=1.0, col="gray30", at=-5)
		plot(legend1, col=cols, legend.only=T, add=T, legend.width=0.5, legend.shrink=0.3, 
			 smallplot=c(0.83,0.86,0.03,0.96), adj=3, axis.args=list(cex.axis=0.9, lwd=0, lwd.tick=0.4, col.tick="gray30", tck=-0.8, col="gray30", 
			 col.axis="gray30", line=0, mgp=c(0,0.5,0)), alpha=1, side=3, horizontal=F)
		index1 = round(((min(counterRaster3a[],na.rm=T)-minV2)/(maxV2-minV2))*100)
		index2 = round(((max(counterRaster3a[],na.rm=T)-minV2)/(maxV2-minV2))*100)
		plot(counterRaster3a, col=colsDiff[index1:index2], ann=F, legend=F, axes=F, box=F)
		plot(coastlines, add=T, col="gray50", lwd=1.0)
		mtext("Difference", side=3, line=-1.5, cex=0.75, col="gray30", at=0)
		mtext("SRI", side=3, line=-3.1, cex=1.0, col="gray30", at=-5)
		plot(legend2, col=colsDiff, legend.only=T, add=T, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.83,0.86,0.03,0.96), adj=3,
			 axis.args=list(cex.axis=0.9, lwd=0, lwd.tick=0.4, col.tick="gray30", tck=-0.8, col="gray30", col.axis="gray30", line=0, 
			 mgp=c(0,0.5,0)), alpha=1, side=3, horizontal=F)		
		dev.off()
	}

# 12. Plotting the environmental variables and their future projections

years = c(2030, 2050, 2070)
year_intervals = c("2021_2040","2041_2060","2061_2080") 
year_names = c("2021-2040","2041-2060","2061-2080") 
scenarios = c("RCP_26","RCP_60","RCP_85")
scenario_names = c("SSP1","SSP3","SSP5")
models = c("GFDL-ESM2M","HadGEM2-ES","IPSL-CM5A-LR","MIROC5")
variables = c("population","pr","tas","landcover")
for (y in 1:length(years))
	{
		buffer1 = list()
		for (s in 1:length(scenarios))
			{
				buffer2 = list()
				for (m in 1:length(models))
					{
						files = list.files(paste0("Environmental_rasters/",scenarios[s],"/",models[m],"/"))
						if (s != 1) files = files[grepl(year_intervals[y],files)]
						files = files[!grepl("5min",files)]
						index_population = which(grepl("population",files)); index_temperature = which(grepl("tas_day",files))
						index_precipitation = which(grepl("pr_day",files)); index_land_cover = which(grepl("landcover",files))
						population = raster(paste0("Environmental_rasters/",scenarios[s],"/",models[m],"/",files[index_population]))
						temperature = raster(paste0("Environmental_rasters/",scenarios[s],"/",models[m],"/",files[index_temperature]))
						precipitation = raster(paste0("Environmental_rasters/",scenarios[s],"/",models[m],"/",files[index_precipitation]))
						land_cover = nc_open(paste0("Environmental_rasters/",scenarios[s],"/",models[m],"/",files[index_land_cover]))
						names(population) = "population"; names(temperature) = "temperature"; names(precipitation) = "precipitation"
						landCoverVariableIDs = names(land_cover$var); cols = list()
						landCoverVariableNames = as.character(read.csv("Land_cover_vars.csv")[1:12,2])
						land_covers1 = list(); land_covers2 = list(); land_covers3 = list()
						for (i in 2:13)
							{
								land_covers1[[i-1]] = brick(paste0("Environmental_rasters/",scenarios[s],"/",models[m],"/",
															files[index_land_cover]), varname=landCoverVariableIDs[i])
							}
						for (i in 1:length(land_covers1)) names(land_covers1[[i]]) = landCoverVariableNames[i]
						variable_codes = c("croplands","pastures","urbanAreas","primaryForest","primaryNonF","secondaryForest","secondaryNonF")
						variable_names = c("crops","pasture","urban land","forested primary land","non-forested primary land",
										   "potentially forested secondary land","potentially non-forested secondary land")
						for (i in 1:length(variable_names))
							{
								names = gsub("\\."," ",landCoverVariableNames); indices = which(landCoverVariableNames==variable_names[i])
								if (length(indices) == 0) indices = which(grepl(variable_names[i],names))
								if (variable_names[i] == "pasture") indices = c(indices, which(grepl("rangeland",names)))
								land_cover = land_covers1[[indices[1]]]; names(land_cover) = variable_codes[i]
								if (length(indices) > 1)
									{
										for (j in 2:length(indices)) land_cover[] = land_cover[]+land_covers1[[indices[j]]][]
									}
								land_covers2[[i]] = land_cover[[1]]; land_covers3[[i]] = raster::aggregate(land_cover[[1]],2)
							}
						envVariables = list()
						envVariables[[1]] = temperature; envVariables[[2]] = precipitation
						envVariables[[3]] = land_covers3[[4]] # primary forest areas
						envVariables[[4]] = land_covers3[[5]] # primary non-forest areas
						envVariables[[5]] = land_covers3[[6]] # secondary forest areas
						envVariables[[6]] = land_covers3[[7]] # secondary non-forest areas
						envVariables[[7]] = land_covers3[[1]] # croplands (all catergories)
						envVariables[[8]] = land_covers3[[2]] # managed pasture + rangeland
						pLog = population; pLog[] = log10(pLog[]+1); envVariables[[9]] = pLog
						for (i in 1:length(envVariables)) envVariables[[i]][is.na(mask[])] = NA
						for (i in 1:length(envVariables)) envVariables[[i]] = crop(envVariables[[i]], europe2, snap="out")
						for (i in 1:length(envVariables)) envVariables[[i]] = mask(envVariables[[i]], europe2)
						envVariables[[1]][] = envVariables[[1]][]-273.15 # conversion to Celcius degrees
						envVariables[[2]][] = envVariables[[2]][]*60*60*24 # conversion to kg/m2/day
						for (i in c(1,2,9)) envVariables[[i]][is.na(envVariables[[3]][])] = NA
						buffer2[[m]] = envVariables
					}
				buffer1[[s]] = buffer2
			}
		envVariables_list[[2+y]] = buffer1
	}
rasters_toPlot = list(); rasters_types = list()
rasters_titles_1 = list(); rasters_titles_2 = list()
for (i in 1:length(envVariables))
	{
		buffer = list(); types = list(); titles_1 = list(); titles_2 = list()
		buffer[[1]] = envVariables_list[[1]][[i]]; types[[1]] = "values"; titles_1[[1]] = "1901-1970"; titles_2[[1]] = "(past)"
		buffer[[2]] = envVariables_list[[2]][[i]]; types[[2]] = "values"; titles_1[[2]] = "2000-2014"; titles_2[[2]] = "(t0)"
		buffer[[9]] = NULL; types[[9]] = "empty"
		buffer[[10]] = buffer[[2]]-buffer[[1]]; types[[10]] = "diffs"; titles_1[[10]] = "2000-2014"; titles_2[[10]] = "(t0 - past)"
		buffer[[17]] = NULL; types[[17]] = "empty"
		buffer[[18]] = NULL; types[[18]] = "empty"
		for (s in 1:length(scenarios))
			{
				for (y in 1:length(years))
					{
						rasts = list()
						for (m in 1:length(envVariables_list[[2+y]][[s]]))
							{
								rasts[[m]] = envVariables_list[[2+y]][[s]][[m]][[i]]
							}
						buffer[[((s-1)*8)+(2+y)]] = mean(stack(rasts))
						types[[((s-1)*8)+(2+y)]] = "values"
						titles_1[[((s-1)*8)+(2+y)]] = paste0(year_names[y])
						titles_2[[((s-1)*8)+(2+y)]] = paste0("(",scenario_names[s],")")
					}
			}
		for (s in 1:length(scenarios))
			{
				for (y in 1:length(years))
					{
						buffer[[((s-1)*8)+(5+y)]] = buffer[[((s-1)*8)+(2+y)]]-buffer[[2]]
						types[[((s-1)*8)+(5+y)]] = "diffs"
						titles_1[[((s-1)*8)+(5+y)]] = paste0(year_names[y])
						titles_2[[((s-1)*8)+(5+y)]] = paste0("(",scenario_names[s]," - t0)")
					}
			}
		rasters_toPlot[[i]] = buffer; rasters_types[[i]] = types
		rasters_titles_1[[i]] = titles_1; rasters_titles_2[[i]] = titles_2
	}
if (savingPlots == TRUE)
	{
		envVariableNames = list(); cols = list(); colsDiff = colorRampPalette(brewer.pal(11,"BrBG"))(100) # 1° choice
		envVariableNames = list(); cols = list(); colsDiff = colorRampPalette(brewer.pal(11,"RdBu"))(100) # 2° choice
		envVariableNames[[1]] = "Temperature near surface (°C)"; cols[[1]] = colorRampPalette(brewer.pal(9,"YlOrRd"))(150)[1:100]
		envVariableNames[[2]] = "Precipitation (kg/m2/day)"; cols[[2]] = colorRampPalette(brewer.pal(9,"YlGnBu"))(100)
		envVariableNames[[3]] = "Forested primary land"; cols[[3]] = colorRampPalette(c("gray97","chartreuse4"),bias=1)(100)
		envVariableNames[[4]] = "Non-forested primary land"; cols[[4]] = colorRampPalette(c("gray97","orange3"),bias=1)(100)
		envVariableNames[[5]] = "Forested secondary land"; cols[[5]] = colorRampPalette(c("gray97","olivedrab3"),bias=1)(100)
		envVariableNames[[6]] = "Non-forested secondary land"; cols[[6]] = colorRampPalette(c("gray97","darkseagreen3"),bias=1)(100)
		envVariableNames[[7]] = "Croplands (all categories)"; cols[[7]] = colorRampPalette(c("gray97","navajowhite4"),bias=1)(100)
		envVariableNames[[8]] = "Pastures and rangeland"; cols[[8]] = colorRampPalette(c("gray97","burlywood3"),bias=1)(100)
		envVariableNames[[9]] = "Human population (log10)"; cols[[9]] = colorRampPalette(brewer.pal(9,"BuPu"))(150)[1:100]
		for (i in 1:length(rasters_toPlot))
			{
				colsDiffI = colsDiff
				if (i == 1)
					{
						colsDiffI = rev(colsDiff)
					}
				minV1 = 9999; maxV1 = -9999; minV2 = 9999; maxV2 = -9999
				for (j in 1:length(rasters_toPlot[[i]]))
					{
						if (rasters_types[[i]][[j]] == "values")
							{
								if (minV1 > min(rasters_toPlot[[i]][[j]][],na.rm=T)) minV1 = min(rasters_toPlot[[i]][[j]][],na.rm=T)
								if (maxV1 < max(rasters_toPlot[[i]][[j]][],na.rm=T)) maxV1 = max(rasters_toPlot[[i]][[j]][],na.rm=T)
							}
						if (rasters_types[[i]][[j]] == "diffs")
							{
								if (minV2 > min(rasters_toPlot[[i]][[j]][],na.rm=T)) minV2 = min(rasters_toPlot[[i]][[j]][],na.rm=T)
								if (maxV2 < max(rasters_toPlot[[i]][[j]][],na.rm=T)) maxV2 = max(rasters_toPlot[[i]][[j]][],na.rm=T)
							}
					}
				if (abs(minV2) < abs(maxV2)) minV2 = -maxV2
				if (abs(maxV2) < abs(minV2)) maxV2 = -minV2
				pdf(paste0("All_the_figures_&_SI/Environmental_variable_",i,".pdf"), width=14, height=6.3)
				par(mfrow=c(3,8), oma=c(0,0,0,0.6), mar=c(0.1,0.75,0.1,2.5), lwd=0.6, col="gray30")
				for (j in 1:length(rasters_toPlot[[i]]))
					{
						if (rasters_types[[i]][[j]] == "values")
							{
								index1 = round(((min(rasters_toPlot[[i]][[j]][],na.rm=T)-minV1)/(maxV1-minV1))*100)
								index2 = round(((max(rasters_toPlot[[i]][[j]][],na.rm=T)-minV1)/(maxV1-minV1))*100)
								plot(rasters_toPlot[[i]][[j]], col=cols[[i]][index1:index2], ann=F, legend=F, axes=F, box=F)
								plot(coastlines, add=T, col="gray50", lwd=1.5)
								mtext(rasters_titles_1[[i]][[j]], side=3, line=-1.5, at=0, cex=0.8, col="gray30")
								mtext(rasters_titles_2[[i]][[j]], side=3, line=-2.8, at=0, cex=0.8, col="gray30")
								rastLegend = raster(t(as.matrix(c(minV1,maxV1))))
								plot(rastLegend, col=cols[[i]], legend.only=T, add=T, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.83,0.86,0.03,0.965), adj=3,
									 axis.args=list(cex.axis=1.0,lwd=0,lwd.tick=0.6,col.tick="gray30",tck=-0.8,col="gray30",col.axis="gray30",line=0,
									 mgp=c(0,0.6,0)),alpha=1, side=3, horizontal=F)						
							}
						if (rasters_types[[i]][[j]] == "diffs")
							{
								index1 = round(((min(rasters_toPlot[[i]][[j]][],na.rm=T)-minV2)/(maxV2-minV2))*100)
								index2 = round(((max(rasters_toPlot[[i]][[j]][],na.rm=T)-minV2)/(maxV2-minV2))*100)
								plot(rasters_toPlot[[i]][[j]], col=colsDiffI[index1:index2], ann=F, legend=F, axes=F, box=F)
								plot(coastlines, add=T, col="gray50", lwd=1.5)
								mtext(rasters_titles_1[[i]][[j]], side=3, line=-1.5, at=0, cex=0.8, col="gray30")
								mtext(rasters_titles_2[[i]][[j]], side=3, line=-2.6, at=0, cex=0.7, col="gray30")
								rastLegend = raster(t(as.matrix(c(minV2,maxV2))))
								plot(rastLegend, col=colsDiffI, legend.only=T, add=T, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.83,0.86,0.03,0.965), adj=3,
									 axis.args=list(cex.axis=1.0,lwd=0,lwd.tick=0.6,col.tick="gray30",tck=-0.8,col="gray30",col.axis="gray30",line=0,
									 mgp=c(0,0.6,0)),alpha=1, side=3, horizontal=F)						
							}
						if (rasters_types[[i]][[j]] == "empty") plot.new()
					}
				dev.off()		
			}
	 }

# 13. BRT projections on historical and future scenarios

selectedModel = "SCV2"; predictions1 = list(); predictions2 = list(); std_deviations = list()
rasters_toPlot = list(); rasters_types = list(); rasters_titles_1 = list(); rasters_titles_2 = list()
cols = c(rep("#F2F4F4",10),rev(colorRampPalette(brewer.pal(11,"RdYlBu"))(120))[21:110]) # 1° choice
cols = c(rep("#F2F4F4",10),rev(hcl.colors(100,palette="terrain2")[1:90])) # 2° choice
colsDiff = colorRampPalette(brewer.pal(11,"BrBG"))(100) # 1° choice
colsDiff = colorRampPalette(brewer.pal(11,"RdBu"))(100) # 2° choice
for (i in 1:dim(species)[1])
	{
		buffer1a = list(); buffer1b = list(); buffer1c = list()
		for (t in 1:2)
			{
				rasters_stack = stack(envVariables_list[[t]]); rast_models = list() # provious error: envVariables_list[[1]]
				brt_model = readRDS(paste0("BRT_projection_files/BRT_models/B_",species[i,"species"],"_",t,"_models_",selectedModel,".rds"))
				for (j in 1:length(brt_model))
					{
						df = as.data.frame(rasters_stack); not_NA = which(!is.na(rowMeans(df))); newdata = df[not_NA,]
						n.trees = brt_model[[j]]$gbm.call$best.trees; type = "response"; single.tree = FALSE
						prediction = predict.gbm(brt_model[[j]], newdata, n.trees, type, single.tree)
						rast = rasters_stack[[1]]; rast[!is.na(rast[])] = prediction; rast_models[[j]] = rast
					}
				buffer1a[[t]] = rast_models; rasts = stack(rast_models); buffer1a[[t]] = mean(rasts)
				buffer1b[[t]] = rast_models; rasts = stack(rast_models); buffer1b[[t]] = mean(rasts)
				buffer1c[[t]] = calc(stack(rasts), fun=sd)
			}
		brt_model = readRDS(paste0("BRT_projection_files/BRT_models/B_",species[i,"species"],"_2_models_",selectedModel,".rds"))
		for (y in 1:length(years))
			{
				buffer2a = list(); buffer2b = list(); buffer2c = list()
				for (s in 1:length(scenarios))
					{
						buffer3 = list(); all_predictions = list()
						for (m in 1:length(models))
							{
								rasters_stack = stack(envVariables_list[[2+y]][[s]][[m]]); replicates = list()
								for (j in 1:length(brt_model))
									{
										df = as.data.frame(rasters_stack); not_NA = which(!is.na(rowMeans(df))); newdata = df[not_NA,]
										n.trees = brt_model[[j]]$gbm.call$best.trees; type = "response"; single.tree = FALSE
										prediction = predict.gbm(brt_model[[j]], newdata, n.trees, type, single.tree)
										rast = rasters_stack[[1]]; rast[!is.na(rast[])] = prediction; replicates[[j]] = rast
										all_predictions[[length(all_predictions)+1]] = rast
									}
								buffer3[[m]] = mean(stack(replicates))
							}
						buffer2a[[s]] = buffer3
						buffer2b[[s]] = mean(stack(buffer3))
						buffer2c[[s]] = calc(stack(all_predictions), fun=sd)
					}
				buffer1a[[2+y]] = buffer2a
				buffer1b[[2+y]] = buffer2b
			}	
		predictions1[[i]] = buffer1a; predictions2[[i]] = buffer1b
		buffer = list(); types = list(); titles_1 = list(); titles_2 = list()
		buffer[[1]] = predictions2[[i]][[1]]; types[[1]] = "values"; titles_1[[1]] = "1901-1970"; titles_2[[1]] = "(past)"
		buffer[[2]] = predictions2[[i]][[2]]; types[[2]] = "values"; titles_1[[2]] = "2000-2014"; titles_2[[2]] = "(t0)"
		buffer[[9]] = NULL; types[[9]] = "observations1"; titles_1[[9]] = "1901-1970"; titles_2[[9]] = "(past)"
		buffer[[10]] = NULL; types[[10]] = "observations2"; titles_1[[10]] = "2000-2014"; titles_2[[10]] = "(t0)"
		buffer[[17]] = NULL; types[[17]] = "empty"
		buffer[[18]] = buffer[[2]]-buffer[[1]]; types[[18]] = "diffs"; titles_1[[18]] = "2000-2014"; titles_2[[18]] = "(t0 - past)"
		for (s in 1:length(scenarios))
			{
				for (y in 1:length(years))
					{
						buffer[[((s-1)*8)+(2+y)]] = predictions2[[i]][[2+y]][[s]]
						types[[((s-1)*8)+(2+y)]] = "values"
						titles_1[[((s-1)*8)+(2+y)]] = paste0(year_names[y])
						titles_2[[((s-1)*8)+(2+y)]] = paste0("(",scenario_names[s],")")
					}
			}
		for (s in 1:length(scenarios))
			{
				for (y in 1:length(years))
					{
						buffer[[((s-1)*8)+(5+y)]] = buffer[[((s-1)*8)+(2+y)]]-buffer[[2]]
						types[[((s-1)*8)+(5+y)]] = "diffs"
						titles_1[[((s-1)*8)+(5+y)]] = paste0(year_names[y])
						titles_2[[((s-1)*8)+(5+y)]] = paste0("(",scenario_names[s]," - t0)")
					}
			}
		rasters_toPlot[[i]] = buffer; rasters_types[[i]] = types
		rasters_titles_1[[i]] = titles_1; rasters_titles_2[[i]] = titles_2
		colsDiffI = colsDiff
		minV1 = 9999; maxV1 = -9999; minV2 = 9999; maxV2 = -9999
		for (j in 1:length(rasters_toPlot[[i]]))
			{
				if (rasters_types[[i]][[j]] == "values")
					{
						if (minV1 > min(rasters_toPlot[[i]][[j]][],na.rm=T)) minV1 = min(rasters_toPlot[[i]][[j]][],na.rm=T)
						if (maxV1 < max(rasters_toPlot[[i]][[j]][],na.rm=T)) maxV1 = max(rasters_toPlot[[i]][[j]][],na.rm=T)
					}
				if (rasters_types[[i]][[j]] == "diffs")
					{
						if (minV2 > min(rasters_toPlot[[i]][[j]][],na.rm=T)) minV2 = min(rasters_toPlot[[i]][[j]][],na.rm=T)
						if (maxV2 < max(rasters_toPlot[[i]][[j]][],na.rm=T)) maxV2 = max(rasters_toPlot[[i]][[j]][],na.rm=T)
					}
			}
		minV1 = 0; maxV1 = 1
		if (abs(minV2) < abs(maxV2)) minV2 = -maxV2
		if (abs(maxV2) < abs(minV2)) maxV2 = -minV2
		if (savingPlots == TRUE)
			{
				pdf(paste0("All_the_figures_&_SI/All_the_species_projections/B_",species[i,"species"],"_projections.pdf"), width=14, height=6.3)
				par(mfrow=c(3,8), oma=c(0,0,0,0.6), mar=c(0.1,0.75,0.1,2.5), lwd=0.6, col="gray30")
				for (j in 1:length(rasters_toPlot[[i]]))
					{
						if (rasters_types[[i]][[j]] == "observations1")
							{
								plot(nullRasters[[1]], col="#F2F4F4", ann=F, legend=F, axes=F, box=F)
								plot(coastlines, add=T, col="gray50", lwd=1.5)
								points(observations_list[[1]][[i]], col="gray30", pch=3, cex=0.6, lwd=0.3)
								mtext(rasters_titles_1[[i]][[j]], side=3, line=-1.5, at=0, cex=0.8, col="gray30")
								mtext(rasters_titles_2[[i]][[j]], side=3, line=-2.8, at=0, cex=0.8, col="gray30")
							}
						if (rasters_types[[i]][[j]] == "observations2")
							{
								plot(nullRasters[[2]], col="#F2F4F4", ann=F, legend=F, axes=F, box=F)
								plot(coastlines, add=T, col="gray50", lwd=1.5)
								points(observations_list[[2]][[i]], col="gray30", pch=3, cex=0.6, lwd=0.3)
								mtext(rasters_titles_1[[i]][[j]], side=3, line=-1.5, at=0, cex=0.8, col="gray30")
								mtext(rasters_titles_2[[i]][[j]], side=3, line=-2.8, at=0, cex=0.8, col="gray30")
							}
						if (rasters_types[[i]][[j]] == "values")
							{
								index1 = round(((min(rasters_toPlot[[i]][[j]][],na.rm=T)-minV1)/(maxV1-minV1))*100)
								index2 = round(((max(rasters_toPlot[[i]][[j]][],na.rm=T)-minV1)/(maxV1-minV1))*100)
								plot(rasters_toPlot[[i]][[j]], col=cols[index1:index2], ann=F, legend=F, axes=F, box=F)
								plot(coastlines, add=T, col="gray50", lwd=1.5)
								mtext(rasters_titles_1[[i]][[j]], side=3, line=-1.5, at=0, cex=0.8, col="gray30")
								mtext(rasters_titles_2[[i]][[j]], side=3, line=-2.8, at=0, cex=0.8, col="gray30")
								rastLegend = raster(t(as.matrix(c(minV1,maxV1))))
								plot(rastLegend, col=cols, legend.only=T, add=T, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.83,0.86,0.03,0.965), adj=3,
									 axis.args=list(cex.axis=1.0,lwd=0,lwd.tick=0.6,col.tick="gray30",tck=-0.8,col="gray30",col.axis="gray30",line=0,
									 mgp=c(0,0.6,0)),alpha=1, side=3, horizontal=F)						
							}
						if (rasters_types[[i]][[j]] == "diffs")
							{
								index1 = round(((min(rasters_toPlot[[i]][[j]][],na.rm=T)-minV2)/(maxV2-minV2))*100)
								index2 = round(((max(rasters_toPlot[[i]][[j]][],na.rm=T)-minV2)/(maxV2-minV2))*100)
								plot(rasters_toPlot[[i]][[j]], col=colsDiffI[index1:index2], ann=F, legend=F, axes=F, box=F)
								plot(coastlines, add=T, col="gray50", lwd=1.5)
								mtext(rasters_titles_1[[i]][[j]], side=3, line=-1.5, at=0, cex=0.8, col="gray30")
								mtext(rasters_titles_2[[i]][[j]], side=3, line=-2.6, at=0, cex=0.7, col="gray30")
								rastLegend = raster(t(as.matrix(c(minV2,maxV2))))
								plot(rastLegend, col=colsDiffI, legend.only=T, add=T, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.83,0.86,0.03,0.965), adj=3,
									 axis.args=list(cex.axis=1.0,lwd=0,lwd.tick=0.6,col.tick="gray30",tck=-0.8,col="gray30",col.axis="gray30",line=0,
									 mgp=c(0,0.6,0)),alpha=1, side=3, horizontal=F)						
							}
						if (rasters_types[[i]][[j]] == "empty") plot.new()
					}
				dev.off()
			}
	}
if (savingPlots == TRUE)
	{
		c = 0
		for (i in 1:dim(species)[1])
			{
				c = c+1
cat(paste0("\\renewcommand{\\thefigure}{S",9+c,"}"),"\n",sep="")
cat("\\begin{figure}[H]","\n",sep="")
cat("\\centering","\n",sep="")
cat("\\makebox[\\textwidth][c]{\\includegraphics[width=1.00\\textwidth]{All_the_species_projections/B_",as.character(species[i,"species"]),"_projections.pdf}}","\n",sep="")
cat(paste0("\\caption{\\small{\\textbf{Ecological niche modelling for \\emph{Bombus ",as.character(species[i,"species"]),"}.}"),"\n",sep="")
cat("We report the projection for the current period as well as differences between future and current projections.}}","\n",sep="")
cat("\\label{FigureS",9+c,"}","\n",sep="")
cat("\\end{figure}","\n",sep="")
			}
	}

# 14. Post hoc analyses to compare present and future ENM projections

rasters_ESI_SRI_list = list(); suffixes = list(); suffixes[[1]] = "ESI"; suffixes[[2]] = "SRI"
ESD_2070_SSP3_t0_list = list() # for 2070 (SSP3) - t0, "ESD" = "ecological suitability difference"
EVD_2070_SSP3_t0_list = list() # for 2070 (SSP3) - t0, "EVD" = "environmental value differences"
for (m in 1:length(models))
	{
		rasters_ESI_SRI = list(); rasters_ESI = list(); rasters_SRI = list()
		counterRasters1 = list(); counterRasters2 = list()
		bufferRasters1 = list(); bufferRasters2 = list()
		ESD_2070_SSP3_t0 = list(); EVD_2070_SSP3_t0 = list()
		for (i in 1:dim(species)[1])
			{
				txt = SIppcs[i,"optimisedThreshold_t1"]; cutOff = as.numeric(unlist(strsplit(txt," \\["))[1])
				bufferRasters1[[i]] = predictions1[[i]][[1]]; c1 = predictions1[[i]][[1]]; c1[c1[]<cutOff] = 0; c1[c1[]>cutOff] = 1; counterRasters1[[i]] = c1
				txt = SIppcs[i,"optimisedThreshold_t2"]; cutOff = as.numeric(unlist(strsplit(txt," \\["))[1])
				bufferRasters2[[i]] = predictions1[[i]][[2]]; c2 = predictions1[[i]][[2]]; c2[c2[]<cutOff] = 0; c2[c2[]>cutOff] = 1; counterRasters2[[i]] = c2
				ESD_2070_SSP3_t0[[i]] = predictions1[[i]][[2]]
			}
		rasters_ESI[[1]] = mean(stack(bufferRasters1)); rasters_ESI[[2]] = mean(stack(bufferRasters2))
		rasters_SRI[[1]] = sum(stack(counterRasters1)); rasters_SRI[[2]] = sum(stack(counterRasters2))
		for (y in 1:length(years))
			{
				buffer1 = list(); buffer2 = list()
				for (s in 1:length(scenarios))
					{
						counterRasters = list(); bufferRasters = list()
						for (i in 1:dim(species)[1])
							{
								EVD_2070_SSP3_t0_sp = list()
								txt = SIppcs[i,"optimisedThreshold_t2"]
								cutOff = as.numeric(unlist(strsplit(txt," \\["))[1])
								bufferRasters[[i]] = predictions1[[i]][[y+2]][[s]][[m]]
								c = predictions1[[i]][[y+2]][[s]][[m]]
								c[c[]<cutOff] = 0; c[c[]>cutOff] = 1; counterRasters[[i]] = c
								if ((years[y] == 2070)&(scenarios[s] == "RCP_60"))
									{
										r1 = ESD_2070_SSP3_t0[[i]]; r2 = bufferRasters[[i]]
										r3 = r1; r3[] = r2[]-r1[]; r3[(r1[]<cutOff)&(r2[]<cutOff)] = NA
										ESD_2070_SSP3_t0[[i]] = r3
										for (j in 1:length(envVariables))
											{
												r4 = envVariables_list[[1]][[j]]
												r5 = envVariables_list[[2+y]][[s]][[m]][[j]]
												r6 = r4; r6[] = r5[]-r4[]; r6[(r1[]<cutOff)&(r2[]<cutOff)] = NA
												EVD_2070_SSP3_t0_sp[[j]] = r6
											}
										EVD_2070_SSP3_t0[[i]] = EVD_2070_SSP3_t0_sp
									}
							}
						buffer1[[s]] = mean(stack(bufferRasters)); buffer2[[s]] = sum(stack(counterRasters))
					}
				rasters_ESI[[2+y]] = buffer1; rasters_SRI[[2+y]] = buffer2
			}
		rasters_ESI_SRI[[1]] = rasters_ESI; rasters_ESI_SRI[[2]] = rasters_SRI
		rasters_ESI_SRI_list[[m]] = rasters_ESI_SRI
		ESD_2070_SSP3_t0_list[[m]] = ESD_2070_SSP3_t0
		EVD_2070_SSP3_t0_list[[m]] = EVD_2070_SSP3_t0
	}
ESD_2070_SSP3_t0 = matrix(nrow=dim(species)[1], ncol=2)
EVD_2070_SSP3_t0 = matrix(nrow=dim(species)[1], ncol=length(envVariables))
colnames(ESD_2070_SSP3_t0) = c("ESD_2070_SSP3_t0", "colour")
colnames(EVD_2070_SSP3_t0) = paste0("EVD_2070_SSP3_t0_",envVariableNames1)
row.names(ESD_2070_SSP3_t0) = species[,"species"]
row.names(EVD_2070_SSP3_t0) = species[,"species"]
for (i in 1:dim(species)[1])
	{
		buffer = list()
		for (m in 1:length(models)) buffer[[m]] = ESD_2070_SSP3_t0_list[[m]][[i]]
		ESD_2070_SSP3_t0[i,1] = mean(mean(stack(buffer), na.rm=T)[], na.rm=T)
		for (j in 1:length(envVariables))
			{
				buffer = list()
				for (m in 1:length(models)) buffer[[m]] = EVD_2070_SSP3_t0_list[[m]][[i]][[j]]
				EVD_2070_SSP3_t0[i,j] = mean(mean(stack(buffer), na.rm=T)[], na.rm=T)
			}
	}
rasters_ESI_SRI = list()
for (i in 1:length(rasters_ESI_SRI_list[[1]]))
	{
		bufferRasters1 = list()
		for (j in 1:length(rasters_ESI_SRI_list[[1]][[i]]))
			{
				if (!is.list(rasters_ESI_SRI_list[[1]][[i]][[j]]))
					{
						bufferRasters2 = list()
						for (m in 1:length(models)) bufferRasters2[[m]] = rasters_ESI_SRI_list[[m]][[i]][[j]]
						bufferRasters1[[j]] = mean(stack(bufferRasters2))
					}	else	{
						bufferRasters2 = list()
						for (k in 1:length(rasters_ESI_SRI_list[[1]][[i]][[j]]))
							{
								bufferRasters3 = list()
								for (m in 1:length(models)) bufferRasters3[[m]] = rasters_ESI_SRI_list[[m]][[i]][[j]][[k]]
								bufferRasters2[[k]] = mean(stack(bufferRasters3))
							}
						bufferRasters1[[j]] = bufferRasters2
					}
			}
		rasters_ESI_SRI[[i]] = bufferRasters1
	}
difference1 = rasters_ESI_SRI[[1]][[2]]-rasters_ESI_SRI[[1]][[1]] # ESI difference between "t0" and the "past"
difference2 = (difference1/rasters_ESI_SRI[[1]][[1]])*100 # ESI difference between "t0" and the "past"
cat("\tMean ESI difference between t0 and the past = ",round(mean(difference2[],na.rm=T),4),sep="") # -4.5 %
cat("\tMaximum ESI difference between t0 and the past = ",round(min(difference2[],na.rm=T),4),sep="") # -33.1 %
cols = rev(colorRampPalette(brewer.pal(11,"RdYlBu"))(131))[21:131] # 1° choice
cols = rev(hcl.colors(100,palette="terrain2")[1:100]) # 2° choice
colsDiff = colorRampPalette(brewer.pal(11,"BrBG"))(100) # 1° choice
colsDiff = colorRampPalette(brewer.pal(11,"RdBu"))(100) # 2° choice
for (i in 1:length(rasters_ESI_SRI))
	{
		if (i == 1) cexAxisLegend = 0.8
		if (i == 2) cexAxisLegend = 1.0
		buffer = list(); types = list(); titles_1 = list(); titles_2 = list()
		buffer[[1]] = rasters_ESI_SRI[[i]][[1]]; types[[1]] = "values"; titles_1[[1]] = "1901-1970"; titles_2[[1]] = "(past)"
		buffer[[2]] = rasters_ESI_SRI[[i]][[2]]; types[[2]] = "values"; titles_1[[2]] = "2000-2014"; titles_2[[2]] = "(t0)"
		buffer[[9]] = NULL; types[[9]] = "empty"
		buffer[[10]] = buffer[[2]]-buffer[[1]]; types[[10]] = "diffs"; titles_1[[10]] = "2000-2014"; titles_2[[10]] = "(t0 - past)"
		buffer[[17]] = NULL; types[[17]] = "empty"
		buffer[[18]] = NULL; types[[18]] = "empty"
		for (s in 1:length(scenarios))
			{
				for (y in 1:length(years))
					{
						buffer[[((s-1)*8)+(2+y)]] = rasters_ESI_SRI[[i]][[2+y]][[s]]
						types[[((s-1)*8)+(2+y)]] = "values"
						titles_1[[((s-1)*8)+(2+y)]] = paste0(year_names[y])
						titles_2[[((s-1)*8)+(2+y)]] = paste0("(",scenario_names[s],")")
					}
			}
		for (s in 1:length(scenarios))
			{
				for (y in 1:length(years))
					{
						buffer[[((s-1)*8)+(5+y)]] = buffer[[((s-1)*8)+(2+y)]]-buffer[[2]]
						types[[((s-1)*8)+(5+y)]] = "diffs"
						titles_1[[((s-1)*8)+(5+y)]] = paste0(year_names[y])
						titles_2[[((s-1)*8)+(5+y)]] = paste0("(",scenario_names[s]," - t0)")
					}
			}
		rasters_toPlot[[i]] = buffer; rasters_types[[i]] = types
		rasters_titles_1[[i]] = titles_1; rasters_titles_2[[i]] = titles_2
		colsDiffI = colsDiff
		minV1 = 9999; maxV1 = -9999; minV2 = 9999; maxV2 = -9999
		for (j in 1:length(rasters_toPlot[[i]]))
			{
				if (rasters_types[[i]][[j]] == "values")
					{
						if (minV1 > min(rasters_toPlot[[i]][[j]][],na.rm=T)) minV1 = min(rasters_toPlot[[i]][[j]][],na.rm=T)
						if (maxV1 < max(rasters_toPlot[[i]][[j]][],na.rm=T)) maxV1 = max(rasters_toPlot[[i]][[j]][],na.rm=T)
					}
				if (rasters_types[[i]][[j]] == "diffs")
					{
						if (minV2 > min(rasters_toPlot[[i]][[j]][],na.rm=T)) minV2 = min(rasters_toPlot[[i]][[j]][],na.rm=T)
						if (maxV2 < max(rasters_toPlot[[i]][[j]][],na.rm=T)) maxV2 = max(rasters_toPlot[[i]][[j]][],na.rm=T)
					}
			}
		if (abs(minV2) < abs(maxV2)) minV2 = -maxV2
		if (abs(maxV2) < abs(minV2)) maxV2 = -minV2
		pdf(paste0("All_the_figures_&_SI/Bombus_",suffixes[[i]],"_diffrences_NEW.pdf"), width=14, height=6.3)
		par(mfrow=c(3,8), oma=c(0,0,0,0.6), mar=c(0.1,0.75,0.1,2.5), lwd=0.6, col="gray30")
		for (j in 1:length(rasters_toPlot[[i]]))
			{
				if (rasters_types[[i]][[j]] == "values")
					{
						index1 = round(((min(rasters_toPlot[[i]][[j]][],na.rm=T)-minV1)/(maxV1-minV1))*100)
						index2 = round(((max(rasters_toPlot[[i]][[j]][],na.rm=T)-minV1)/(maxV1-minV1))*100)
						plot(rasters_toPlot[[i]][[j]], col=cols[index1:index2], ann=F, legend=F, axes=F, box=F)
						plot(coastlines, add=T, col="gray50", lwd=1.5)
						mtext(rasters_titles_1[[i]][[j]], side=3, line=-1.5, at=0, cex=0.8, col="gray30")
						mtext(rasters_titles_2[[i]][[j]], side=3, line=-2.8, at=0, cex=0.8, col="gray30")
						rastLegend = raster(t(as.matrix(c(minV1,maxV1))))
						plot(rastLegend, col=cols, legend.only=T, add=T, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.83,0.86,0.03,0.965), adj=3,
							 axis.args=list(cex.axis=cexAxisLegend,lwd=0,lwd.tick=0.6,col.tick="gray30",tck=-0.8,col="gray30",col.axis="gray30",line=0,
							 mgp=c(0,0.6,0)),alpha=1, side=3, horizontal=F)						
					}
				if (rasters_types[[i]][[j]] == "diffs")
					{
						index1 = round(((min(rasters_toPlot[[i]][[j]][],na.rm=T)-minV2)/(maxV2-minV2))*100)
						index2 = round(((max(rasters_toPlot[[i]][[j]][],na.rm=T)-minV2)/(maxV2-minV2))*100)
						plot(rasters_toPlot[[i]][[j]], col=colsDiffI[index1:index2], ann=F, legend=F, axes=F, box=F)
						plot(coastlines, add=T, col="gray50", lwd=1.5)
						mtext(rasters_titles_1[[i]][[j]], side=3, line=-1.5, at=0, cex=0.8, col="gray30")
						mtext(rasters_titles_2[[i]][[j]], side=3, line=-2.6, at=0, cex=0.7, col="gray30")
						rastLegend = raster(t(as.matrix(c(minV2,maxV2))))
						plot(rastLegend, col=colsDiffI, legend.only=T, add=T, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.83,0.86,0.03,0.965), adj=3,
							 axis.args=list(cex.axis=cexAxisLegend,lwd=0,lwd.tick=0.6,col.tick="gray30",tck=-0.8,col="gray30",col.axis="gray30",line=0,
							 mgp=c(0,0.6,0)),alpha=1, side=3, horizontal=F)						
					}
				if (rasters_types[[i]][[j]] == "empty") plot.new()
			}
		dev.off()
	}
if (savingPlots == TRUE)
	{
		relativeInfluences = read.csv("RI_BRT_2000-14.csv", header=T); colours = c()
		pca = dudi.pca(relativeInfluences, scannf=F, nf=9); lis = pca$li[,1:2]; cos = pca$co
		colourScale1 = colorRampPalette(brewer.pal(9,"Greens"))(53); maxV = 0.071
		colourScale2 = rev(colorRampPalette(brewer.pal(9,"Reds"))(364)); minV = -0.51
		legendRast = raster(as.matrix(cbind(minV,maxV)))
		colourScale = c(colourScale2,colourScale1)
		for (i in 1:dim(lis)[1])
			{
				index1 = which(species[,"species"]==row.names(lis)[i])
				v = as.numeric(ESD_2070_SSP3_t0[index1,"ESD_2070_SSP3_t0"])
				index2 = (((v-minV)/(maxV-minV))*length(colourScale))+1
				ESD_2070_SSP3_t0[index1,"colour"] = colourScale[index2]
			}
		pdf("All_the_figures_&_SI/BRT_relative_influences_2_NEW.pdf", width=4.5, height=4); par(mar=c(3,3,0,1.5))
		plot(lis, col=ESD_2070_SSP3_t0[,"colour"], cex=0.7, pch=16, ann=F, axes=F, xlim=c(-4.0,3.5), ylim=c(-4.0,3.5))
		s.corcircle(4*cos,xax=1,yax=2,box=F,sub="",csub=1,clabel=0.5,possub="topleft",grid=T,cgrid=1,full=T,add.plot=T)
		text(lis[,1], lis[,2], labels=gsub("B_","",row.names(lis)), cex=0.6, col="gray50", pos=4, offset=0.25)
		points(lis, col=ESD_2070_SSP3_t0[,"colour"], cex=0.9, pch=16); points(lis, col="gray30", cex=0.9, pch=1, lwd=0.3)
		axis(side=1, lwd.tick=0.2, cex.axis=0.6, lwd=0.2, tck=-0.017, col.axis="gray30", mgp=c(0,0.05,0), at=seq(-8,5,1))
		axis(side=2, lwd.tick=0.2, cex.axis=0.6, lwd=0.2, tck=-0.017, col.axis="gray30", mgp=c(0,0.30,0), at=seq(-5,9,1))
		title(xlab="PCA axis 1", cex.lab=0.7, mgp=c(1.0,0,0), col.lab="gray30")
		title(ylab="PCA axis 2", cex.lab=0.7, mgp=c(1.3,0,0), col.lab="gray30")
		dev.off()
		pdf("All_the_figures_&_SI/BRT_relative_influences_3_NEW.pdf", width=1, height=3); par(mar=c(1,1,1,1), col="gray30", lwd=0.2)
		plot(legendRast, legend.only=T, add=T, col=colourScale, legend.width=0.5, alpha=1, legend.shrink=0.3, smallplot=c(0.32,0.4,0.4,0.8),
			 legend.args=list(text="", cex=0.6, line=0.5, col.lab="gray30", col.axis="gray30", col="gray30"),
			 axis.args=list(cex.axis=0.7, lwd=0, lwd.tick=0.2, tck=-0.5, line=0, mgp=c(0,0.4,0), col.lab="gray30", col.axis="gray30", col="gray30"))
		dev.off()
		pca = dudi.pca(EVD_2070_SSP3_t0, scannf=F, nf=9); lis = pca$li[,1:2]; cos = pca$co
		colourScale1 = colorRampPalette(brewer.pal(9,"Greens"))(53); maxV = 0.048
		colourScale2 = rev(colorRampPalette(brewer.pal(9,"Reds"))(364)); minV = -0.406
		legendRast = raster(as.matrix(cbind(minV,maxV)))
		colourScale = c(colourScale2,colourScale1); colours = c()
		for (i in 1:dim(lis)[1])
			{
				index = which(species[,"species"]==row.names(lis)[i])
				status = species[index,"IUCN_status"]
				colours = c(colours, IUCN_colours[which(IUCN_status==status)])
			}
		for (i in 1:dim(lis)[1])
			{
				index1 = which(species[,"species"]==row.names(lis)[i])
				v = as.numeric(ESD_2070_SSP3_t0[index1,"ESD_2070_SSP3_t0"])
				index2 = (((v-minV)/(maxV-minV))*length(colourScale))+1
				ESD_2070_SSP3_t0[index1,"colour"] = colourScale[index2]
			}
		pdf("All_the_figures_&_SI/BRT_EVD_2070_SSP3_t0_NEW.pdf", width=9, height=3.5); par(mfrow=c(1,2), mar=c(2.0,2.0,0,0.5))
		plot(lis, col=colours, cex=0.7, pch=16, ann=F, axes=F, xlim=c(-6.3,4.0), ylim=c(-2.15,3.05))
		text(lis[,1], lis[,2], labels=gsub("B_","",row.names(lis)), cex=0.6, col="gray50", pos=4, offset=0.25)
		points(lis, col=colours, cex=0.9, pch=16); points(lis, col="gray30", cex=0.9, pch=1, lwd=0.3)
		axis(side=1, lwd.tick=0.2, cex.axis=0.6, lwd=0.2, tck=-0.017, col.axis="gray30", mgp=c(0,0.05,0), at=seq(-8,5,1))
		axis(side=2, lwd.tick=0.2, cex.axis=0.6, lwd=0.2, tck=-0.017, col.axis="gray30", mgp=c(0,0.30,0), at=seq(-5,9,1))
		title(xlab="PCA axis 1", cex.lab=0.7, mgp=c(1.0,0,0), col.lab="gray30")
		title(ylab="PCA axis 2", cex.lab=0.7, mgp=c(1.3,0,0), col.lab="gray30")
		plot(lis, col=ESD_2070_SSP3_t0[,"colour"], cex=0.7, pch=16, ann=F, axes=F, xlim=c(-6.3,4.0), ylim=c(-2.15,3.05))
		s.corcircle(2.5*cos,xax=1,yax=2,box=F,sub="",csub=1,clabel=0.5,possub="topleft",grid=T,cgrid=1,full=T,add.plot=T)
		text(lis[,1], lis[,2], labels=gsub("B_","",row.names(lis)), cex=0.6, col="gray50", pos=4, offset=0.25)
		points(lis, col=ESD_2070_SSP3_t0[,"colour"], cex=0.9, pch=16); points(lis, col="gray30", cex=0.9, pch=1, lwd=0.3)
		axis(side=1, lwd.tick=0.2, cex.axis=0.6, lwd=0.2, tck=-0.017, col.axis="gray30", mgp=c(0,0.05,0), at=seq(-8,5,1))
		axis(side=2, lwd.tick=0.2, cex.axis=0.6, lwd=0.2, tck=-0.017, col.axis="gray30", mgp=c(0,0.30,0), at=seq(-5,9,1))
		title(xlab="PCA axis 1", cex.lab=0.7, mgp=c(1.0,0,0), col.lab="gray30")
		title(ylab="PCA axis 2", cex.lab=0.7, mgp=c(1.3,0,0), col.lab="gray30")
		dev.off()
	}

relativeInfluences = read.csv("RI_BRT_2000-14.csv", header=T)
df = data.frame(as.matrix(cbind(ESD_2070_SSP3_t0[,"ESD_2070_SSP3_t0"], relativeInfluences))); colnames(df)[1] = "ESD_2070_SSP3_t0"
for (i in 1:dim(df)[2]) df[,i] = as.numeric(df[,i])
for (i in 2:dim(df)[2])
	{
		cat(colnames(df)[i],": ",round(cor(df[,1],df[,i],method="spearman"),2),"\n",sep="")
	}
i = 1 # to plot successive scatterplots:
i = i+1; plot(df[,i], df[,"ESD_2070_SSP3_t0"], main=colnames(df)[i])

df = data.frame(as.matrix(cbind(ESD_2070_SSP3_t0[,"ESD_2070_SSP3_t0"], EVD_2070_SSP3_t0)))
colnames(df) = c("ESD_2070_SSP3_t0","temperature","precipitation","primaryForest","primaryNonF","secondaryForest","secondaryNonF","croplands","pastures","population")
for (i in 1:dim(df)[2]) df[,i] = as.numeric(df[,i])
for (i in 2:dim(df)[2])
	{
		cat(colnames(df)[i],": ",round(cor(df[,1],df[,i],method="spearman"),2),"\n",sep="")
	}
i = 1 # to plot successive scatterplots:
i = i+1; plot(df[,i], df[,"ESD_2070_SSP3_t0"], main=colnames(df)[i])

	# ANOVA tests on IUCN data:

y = as.numeric(ESD_2070_SSP3_t0[,"ESD_2070_SSP3_t0"])
x = species[,"IUCN_status"]
boxplot(y ~ x); anova = aov(y ~ x)
shapiro.test(residuals(anova)) # p-value < 0.001
bartlett.test(y ~ x) # p-value = 0.142
kruskal.test(y ~ x) # p-value = 0.762

	# Computing the ESR values:

ESR_comparison = matrix(nrow=dim(species)[1], ncol=length(scenarios))
ESR_meanValues = matrix(nrow=dim(species)[1], ncol=length(scenarios))
row.names(ESR_comparison) = paste0("B. ",species[,1])
colnames(ESR_comparison) = scenario_names
row.names(ESR_meanValues) = paste0("B. ",species[,1])
colnames(ESR_meanValues) = scenario_names
for (i in 1:dim(species)[1]) # ESR = "ecological suitability ratio"
	{
		txt = SIppcs[i,"optimisedThreshold_t2"]
		cutOff = as.numeric(unlist(strsplit(txt," \\["))[1])
		for (j in 1:length(scenarios))
			{
				r1 = predictions1[[i]][[2]]
				r1[r1[]<cutOff] = NA; s1 = sum(!is.na(r1[]))
				if (s1 != 0)
					{
						ESRs = rep(NA, length(models))
						for (m in 1:length(models))
							{
								r2 = predictions1[[i]][[5]][[j]][[m]] # 2070
								r2[r2[]<cutOff] = NA; s2 = sum(!is.na(r2[]))
								ESRs[m] = s2/s1
							}
						meanV = round(mean(ESRs),2)
						minV = round(min(ESRs),2)
						maxV = round(max(ESRs),2)
						ESR_comparison[i,j] = paste0(meanV," [",minV,"-",maxV,"]")
						ESR_meanValues[i,j] = meanV
					}	else	{
						ESR_comparison[i,j] = mean(ESRs)
					}
			}
	}
write.csv(ESR_comparison, "ESR_comparison.csv", quote=F)

ESR_meanValues_LC = ESR_meanValues[which(species["IUCN_status"]=="least_concern"),]
sum(ESR_meanValues_LC[,"SSP1"]<0.70)/dim(ESR_meanValues_LC)[1] # 32%
sum(ESR_meanValues_LC[,"SSP3"]<0.70)/dim(ESR_meanValues_LC)[1] # 57%
sum(ESR_meanValues_LC[,"SSP5"]<0.70)/dim(ESR_meanValues_LC)[1] # 76%
sum(ESR_meanValues[,"SSP1"]<0.70)/dim(ESR_meanValues)[1] # 35%
sum(ESR_meanValues[,"SSP3"]<0.70)/dim(ESR_meanValues)[1] # 57%
sum(ESR_meanValues[,"SSP5"]<0.70)/dim(ESR_meanValues)[1] # 76%

