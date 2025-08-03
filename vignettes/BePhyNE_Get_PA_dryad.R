
{
  library(rgbif)
  library(ggplot2)
  library(maptools)
  library(CoordinateCleaner)
  library(raster)
  library(dismo)
  library(dplyr)
  library(sp)
  library(countrycode)
  library(taxize)
  library(Imap)
  library(tidyjson)
  library(R.utils)
  library(rlist)
  #library(prob)
  library(rnaturalearthdata)
  #library(bazar)
  library(rgdal)
  library(rgeos)
  library(stringi)
  library(rlang)
  library(sf)
  #library(splus2R)
  library(treeplyr)
  library(devtools)
  library(dismo)
  library(geiger)
  #library(datelife)
  library(ape)
  library(castor)
  library(spThin)

  library(RColorBrewer)
  library(maptools)
  library(maps)
  library(mapdata)
  library(rangeBuilder)
  #devtools::install_git("https://github.com/SantanderMetGroup/mopa.git")
  #library(mopa)
}

library(BePhyNE)



pleth_pruned<-read.tree("ENA_98_sp.tre")



##unfortunately when pulling from ti[s we still have some taxa that are a challenge to get data for specifically, E. troglodytes and E.tynerensis lineages
## to resolve this for now I will edit and provide and
df_sp_list <- as.data.frame(pleth_pruned$tip.label)




df_sp_list <- as.data.frame(pleth_pruned$tip.label)


names(df_sp_list)<-"scientificName"

class(df_sp_list)






#sp_dat <- readRDS("Pletho_gbif_raw_list.RDS")




#### Match names to gbif taxon keys (more convenient) #######

gbif_taxon_keys<-list()

#changed this to a loop so we know which taxa we dont get a taxon key for, these will show up at "integer(o)"

for (i in 1:nrow(df_sp_list)){

  gbif_taxon_keys[[i]] <- df_sp_list[i,] %>%
    #pull(scientificName) %>%
    taxize::get_gbifid_(method = "backbone") %>% # match names to the GBIF backbone to get taxonkeys
    bind_rows() %>% # combine all data.frames into one
    filter(matchtype =="EXACT" & status =="ACCEPTED") %>% # get only accepted and matched names
    filter(family =="Plethodontidae") %>% # remove anything that might have matched to non-animalia
    pull(usagekey) # get gbif taxonkeys

  if(length(gbif_taxon_keys[[i]])==0){
    gbif_taxon_keys[[i]]=NA
  }

}
# error with filter (matchtype...) ?

gbif_id_df<-data.frame(unlist(gbif_taxon_keys),df_sp_list)

length(gbif_taxon_keys)
length(unique(gbif_taxon_keys))
nrow(df_sp_list)

#### Get occurrences######
#### ***There is a limit of 500 species for occ_search***

lapply(1:length(gbif_taxon_keys), function(i) occ_count(gbif_taxon_keys[[i]], georeferenced = T))

max(unlist)

values_pres<-occ_search(gbif_taxon_keys, hasCoordinate=TRUE, limit = 10000, start = 0)



{

#values_pres <- readRDS("GBIF_clim_pres.RDS")


names(values_pres)<-gsub(" ","_",names(values_pres)  )

values_pres<-values_pres[names(values_pres)%in%pleth_pruned$tip.label]

lapply(values_pres, nrow)

#calculate scale from occurences

values_pres_concat<-do.call(rbind, values_pres)

values_pres_concat[,4]

scaled_occurrences<-scale(do.call(rbind, values_pres))


#saveRDS(scaled_occurrences, file = "scaled_GBIF_clim_pres.RDS")



#View(values_pres)
#trim to only the bioclim variables we want (1 and 12)





#values_pres[[sp]][,c(1,2,3,4,15)]

#values_pres <- lapply(1:length(values_pres), function(sp) na.omit(values_pres[[sp]][,c(1,2,3,4,15)]))

#filter all occurrences for other species as absences....gotta filter from here though otherwise you will get biased garbage

values_abs<- lapply(1:length(values_pres), function(sp) cbind(do.call(rbind,values_pres[-sp]), seq(0, by=0, length=nrow(do.call(rbind,values_pres[-sp])) ) ) )




#names(values_pres)<-names(sp_occs)

#lass(expert_filtered_list[[1]])

#class(values_pres[[sp]][,2:3])



sp_occs_exp_filtered <- lapply(1:length(values_pres), function(sp)   try(SpatialPoints(coords = values_pres[[sp]][,2:3],
                                                                                       proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 ")) )
)


#one has super few points (like 1) so the dimensions are differemt subset and try again for this one (this is also why we have the SpatialPoint fn above in the try() function)


sp_occs_exp_filtered[[40]]<-SpatialPoints(coords = rbind(values_pres[[43]][2:3],values_pres[[43]][2:3]),
                                          proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 "))

sp_occs_abs <- lapply(1:length(values_abs), function(sp)   SpatialPoints(coords = values_abs[[sp]][,2:3],
                                                                         proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 "))
)




names(sp_occs_exp_filtered)<-names(values_pres)
names(sp_occs_abs)<-names(sp_occs_exp_filtered)

#######setup spatial object for plotting and parsing##########

library(RColorBrewer)
library(maptools)
library(maps)
library(mapdata)
library(rangeBuilder)


max.lat<- 66
min.lat<- 27

max.lon<- (-67)
min.lon<- (-105)



xlim<-c(min.lon, max.lon)
ylim<-c(min.lat, max.lat)

pleth_shp <- shapefile('data_0.shp')
pleth_shp$BINOMIAL<-gsub(" ","_",pleth_shp$BINOMIAL  )


ENA_pleth_shp<- pleth_shp[pleth_shp$BINOMIAL%in%names(values_pres),]

nrow(ENA_pleth_shp)

ENA_pleth_shp_trans<-spTransform(ENA_pleth_shp, "+init=epsg:3857")


buffers<-unlist(lapply(1:nrow(ENA_pleth_shp_trans), function(x) sqrt(area(ENA_pleth_shp_trans[x,])/100)))

buffers/1000

ENA_pleth_shp_trans_buff<-gBuffer(ENA_pleth_shp_trans, width=buffers, byid = T)

ENA_pleth_shp_buff<-spTransform(ENA_pleth_shp_trans_buff, CRS(proj4string(ENA_pleth_shp)))





{

  #based on this plot the absences totally overwhelm presences for every species..yikes!!!

  pdf(file= paste("PA_IUCN_DS_maps_71.pdf"))

  par(mfrow = c(2, 2))

  cols <- rainbow(115)
  for( i in 1:length(sp_occs_exp_filtered)){
    map('world', fill = F, col = 1, xlim =xlim, ylim = ylim);
    plot(sp_occs_abs[[i]],add=T,col=2, pch=9, cex=.3);
    plot(sp_occs_exp_filtered[[i]],add=T,col=1, pch=9, cex=.3);
    if(length(grep(names(sp_occs_exp_filtered)[[i]], ENA_pleth_shp$BINOMIAL) )==0){
      title(paste(i," ",names(sp_occs_exp_filtered)[[i]], " NA IUCN" ))
    } else{
      plot(ENA_pleth_shp_buff[(ENA_pleth_shp_buff$BINOMIAL==names(sp_occs_exp_filtered)[[i]]),], add=T, col=transparentColor(namedColor = "red", alpha=.1) )
      plot(ENA_pleth_shp[(ENA_pleth_shp$BINOMIAL==names(sp_occs_exp_filtered)[[i]]),], add=T, col=transparentColor(namedColor = "blue", alpha=.1) )

      title(paste(i," ", names(sp_occs_exp_filtered)[[i]]))

    }

  }
  dev.off()
}

values_abs_out_range<-list()
values_abs_in_range<-list()
values_abs_all<-list()
#sp=8

#pres_test<-values_pres

#sample down presences so we dont have an obscene number


for (sp in 1: length(values_pres)){

  if(is.null(nrow(values_pres[[sp]]))==F){
    if(nrow(values_pres[[sp]])>1000){

      values_pres[[sp]]<-values_pres[[sp]][sample(1:nrow(values_pres[[sp]]), size=1000 ),]

    }
  }


}

#sample(names(values_pres)[sp], replace=T, size = 10)

#values_pres<-lapply(1:length(values_pres), function(sp) cbind(values_pres[[sp]], species=gsub(" ", "_", names(values_pres)[[sp]])) )

#colnames(values_pres[[1]][,c(1,6)])<-c("y","species")


#View(values_pres)


#do.call(rbind, values_pres)

#View(values_pres)

#spThin::thin(loc.data = values_pres[[1]], lat.col = "decimallatitude", long.col = "decimallongitude", spec.col = "species" )

#how about we downsample occurrences for ultra common species with a ton of overlap and account for that spatial autocorrelation

#values_pres[[1]][[which(names(values_pres[1])=="species")]]



#na.omit(values_abs[[sp]])
#set maximum absences out of environmental range
#filter absences so that each species has equal amount of absences and presences within environmental range, and then set max number of absences outside of range

#add the max out as the sample size for values_abs_out_range if want to set number, im just going to chose double the # of presences.. but this is still a lot.. I dunno
#max_out_abs=500
#or
#2*nrow(values_pres[[sp]])

out_range_size<-list()

for (sp in 1:length(values_pres)) {

  print(2*nrow(values_pres[[sp]]))

  if(is.null(nrow(values_pres[[sp]]))){

    btwn_logical_filter<- (between(x = values_abs[[sp]][,4], left = min(values_pres[[sp]][4])-10, right = max(values_pres[[sp]][4])+10 ) &
                             between(x = values_abs[[sp]][,5], left = min(values_pres[[sp]][5])-10, right = max(values_pres[[sp]][5])+10 ))

    #pulling all absences within range and
    values_abs_in_range[[sp]]<- values_abs[[sp]][btwn_logical_filter, ]
    values_abs_out_range[[sp]]<-values_abs[[sp]][!btwn_logical_filter,]
    values_abs_all[[sp]]<-rbind(values_abs_in_range[[sp]][sample(nrow(values_abs_in_range[[sp]]),size=1), ],
                                values_abs_out_range[[sp]][sample(nrow(values_abs_out_range[[sp]]),size=100),])

  }else{
    #filter to be between climatic range of occurrences
    btwn_logical_filter<- (between(x = values_abs[[sp]][,4], left = min(values_pres[[sp]][,4])-10, right = max(values_pres[[sp]][,4])+10 ) &
                             between(x = values_abs[[sp]][,5], left = min(values_pres[[sp]][,5])-10, right = max(values_pres[[sp]][,5])+10 ))

    values_abs_in_range[[sp]]<-values_abs[[sp]][btwn_logical_filter,]
    values_abs_out_range[[sp]]<-values_abs[[sp]][!btwn_logical_filter,]


    if(nrow(values_abs_out_range[[sp]])<2000 ) {
      out_range_size[[sp]]=nrow(values_abs_out_range[[sp]])

    }else{
      out_range_size[[sp]]=2*nrow(values_pres[[sp]])
    }


    values_abs_all[[sp]]<-rbind(values_abs_in_range[[sp]][sample(nrow(values_abs_in_range[[sp]]),size=nrow(values_pres[[sp]])), ],
                                values_abs_out_range[[sp]][sample(nrow(values_abs_out_range[[sp]]),size=out_range_size[[sp]]),])

  }
}                  #     &
#      between(x = values_abs[[sp]][,5], left = min(values_pres[[sp]][,5])-10, right = max(values_pres[[sp]][,5])+10 ))

unlist(lapply(values_abs_in_range, nrow))
unlist(lapply(values_abs_out_range, nrow))
unlist(lapply(values_abs_all, nrow))+unlist(lapply(values_pres, nrow))



unlist(lapply(values_pres, nrow))


unlist(lapply(data_new_final, function(x) length(x$y==0) ))

names(values_pres)[[40]]

sp_occs_abs <- lapply(1:length(values_abs_all), function(sp)   SpatialPoints(coords = values_abs_all[[sp]][,2:3],
                                                                             proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 "))
)

sp_occs_pres <- lapply(1:length(values_pres), function(sp)   try(SpatialPoints(coords = values_pres[[sp]][,2:3],
                                                                           proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 "))
))

sp_occs_pres[[40]] <- SpatialPoints(coords = rbind(values_pres[[40]][2:3],values_pres[[40]][2:3]),
                                          proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 "))


{

  #okay lets try now with absences separated by within and outside of climatic range of occurrences, with those within the climatic range environemntally filtered

  pdf(file= paste("filtered_PA_IUCN_maps_71.pdf"))

  par(mfrow = c(2, 2))

  cols <- rainbow(115)
  lapply(1:length(sp_occs_exp_filtered), function(i){
    map('world', fill = F, col = 1, xlim =xlim, ylim = ylim);
    plot(sp_occs_abs[[i]],add=T,col=transparentColor(namedColor = "red", alpha=.5), pch=9, cex=.3);
    plot(sp_occs_exp_filtered[[i]],add=T,col=transparentColor(namedColor = "grey", alpha=.7), pch=9, cex=.3);

    if(length(grep(names(sp_occs_exp_filtered)[[i]], ENA_pleth_shp$BINOMIAL) )==0){
      title(paste(i," ",names(sp_occs_exp_filtered)[[i]], " NA IUCN" ))
    } else{
      plot(ENA_pleth_shp_buff[(ENA_pleth_shp_buff$BINOMIAL==names(sp_occs_exp_filtered)[[i]]),], add=T, col=transparentColor(namedColor = "green", alpha=.1) )
      plot(ENA_pleth_shp[(ENA_pleth_shp$BINOMIAL==names(sp_occs_exp_filtered)[[i]]),], add=T, col=transparentColor(namedColor = "blue", alpha=.1) )

      title(paste(i," ", names(sp_occs_exp_filtered)[[i]]))

    }

  })
  dev.off()
}


values_pres_all<-lapply(1:length(values_abs_all), function(x)  try(values_pres[[x]][,c(2:5,1)] ) )
values_pres_all[[40]]<-values_pres[[40]][c(2:5,1)]

unlist(lapply(values_pres_all, function(x) length(x[,2]) ))



values_abs_all<-lapply(1:length(values_abs_all), function(x)  values_abs_all[[x]][,2:6] )

unlist(lapply(values_pres_abs_all, function(x) length(x[,5]) ))

values_pres_abs_all[[1]][,5]


values_pres_abs_all<-lapply(1:length(values_abs_out_range), function(x){print(x);
  try(rbind(values_pres[[x]][,c(2:5,1)],
            values_abs_all[[x]]) )})

names(values_pres_abs_out)<-names(values_pres)

values_pres_abs_all[[40]]<-rbind(values_pres[[40]][c(2:5,1)],values_abs_all[[40]])

data_new_final_unscaled<-list()
for (i in 1:length(values_pres_abs_out)){
  data_new_final_unscaled[[i]]=list("Lat"=values_pres_abs_all[[i]][,1],
                                    "Lon"=values_pres_abs_all[[i]][,2],
                                    species=gsub(" ", "_", names(values_pres_abs_out)[[i]]),
                                    "y"=values_pres_abs_all[[i]][,5],
                                    "X1"=values_pres_abs_all[[i]][,4],
                                    "X2"=values_pres_abs_all[[i]][,3])
}

#scaled data final by scaled presences (I was bum an didnt scale them originally by presences and thought i could scale it by pres and absences...but then you would be scaling them by the absences too which is strange..didnt want to be a bum and did it by presences, ideally you determine scaling yb full worlclim data but i was having trouble extracting the full values from the raster)
data_new_final<-data_new_final_unscaled


for(sp in 1:length(data_final)){

  data_new_final[[sp]]$X1<-(data_new_final_unscaled[[sp]]$X1-scale_atr$center[[1]])/scale_atr$scale[[1]]
  data_new_final[[sp]]$X2<-(data_new_final_unscaled[[sp]]$X2-scale_atr$center[[2]])/scale_atr$scale[[2]]

}



{

  pdf(file= paste("data_new_final_Espace.pdf"))

  par(mfrow=c(3,3))
  lapply(1:length(data_new_final), function(sp){
    plot(c(-4,4),c(4,-7), xlab="precip",ylab="temp",main=data_new_final[[sp]]$species, col=4);
    points(data_new_final[[sp]]$X1[data_new_final[[sp]]$y==0],data_new_final[[sp]]$X2[data_new_final[[sp]]$y==0],col=2);
    points(data_new_final[[sp]]$X1[data_new_final[[sp]]$y==1],data_new_final[[sp]]$X2[data_new_final[[sp]]$y==1]);
  })
  dev.off()
}
#reorder list of lists to match tree

{
  pdf(file= paste("data_new_final_Espace_hist.pdf"))

  par(mfrow=c(3,4))
  lapply(1:length(data_new_final), function(sp){
    hist(data_new_final[[sp]]$X1[data_new_final[[sp]]$y==0],col=2, xlab="precip",ylab="temp",main=data_new_final[[sp]]$species);
    hist(data_new_final[[sp]]$X1[data_new_final[[sp]]$y==1], xlab="precip",ylab="temp",main=data_new_final[[sp]]$species);
  })
  dev.off()
}



}

ENA_Pleth_Tree<-pleth_pruned

cave_sals<-c(extract.clade(pleth_pruned, node=115)$tip.label, "Eurycea_wallacei",  "Eurycea_subfluvicola", "Gyrinophilus_subterraneus", "Urspelerpes_brucei", "Gyrinophilus_porphyriticus", "Hemidactylium_scutatum" , "Plethodon_albagula"  )

ENA_Pleth_Tree=drop.tip(ENA_Pleth_Tree, cave_sals)



phylo$edge.length<-phylo$edge.length/(max(branching.times(phylo)))



ENA_Pletho_PA<-do.call(rbind, lapply(1:length(data_new_final), function(sp) do.call(cbind,data_new_final[[sp]]) ))

#saveRDS(ENA_Pleth_PA,"ENA_Pleth_PA" )

ENA_Pleth_PA<-as.data.frame(readRDS("ENA_Pletho_PA.RDS"),row.names = F)
as.numeric(ENA_Pleth_PA)



ENA_Pleth_PA<-ENA_Pleth_PA[!ENA_Pleth_PA$species%in%cave_sals,]


#save(ENA_Pleth_PA, file="ENA_Pleth_PA.RData")

#save(ENA_Pleth_Tree, file="ENA_Pleth_Tree.RData")


##########load old data


#pull occurrence data from GBIF with worldclim data concatenated and scaled
pres_data_scaled<-readRDS("scaled_GBIF_clim_pres.RDS")

pres_data_scaled[,colnames(pres_data_scaled)=="bio1"]

attr(pres_data_scaled, "scaled:center")[colnames(pres_data_scaled)=="bio1"]



scale_atr <- list(scale=list(attr(pres_data_scaled, "scaled:scale")[colnames(pres_data_scaled)=="bio12"],
                             attr(pres_data_scaled, "scaled:scale")[colnames(pres_data_scaled)=="bio1"]) ,
                  center= list(attr(pres_data_scaled, "scaled:center")[colnames(pres_data_scaled)=="bio12"],
                               attr(pres_data_scaled, "scaled:center")[colnames(pres_data_scaled)=="bio1"]))

#pull full data with filtered target group absences
data_final_unscaled<-readRDS("GBIF_clim_PA_DS.RDS")




# data_final_unscaled<-readRDS("GBIF_clim_PA_abs_out_DS.RDS")



####



#scaled data final by scaled presences (I was bum an didnt scale them originally by presences and thought i could scale it by pres and absences...but then you would be scaling them by the absences too which is strange..didnt want to be a bum and did it by presences, ideally you determine scaling yb full worlclim data but i was having trouble extracting the full values from the raster)

for(sp in 1:length(data_final)){

  data_final[[sp]]$X1<-(data_final_unscaled[[sp]]$X1-scale_atr$center[[1]])/scale_atr$scale[[1]]
  data_final[[sp]]$X2<-(data_final_unscaled[[sp]]$X2-scale_atr$center[[2]])/scale_atr$scale[[2]]

}




{

  pdf(file= paste("Data_final_Espace.pdf"))

  par(mfrow=c(3,3))
  lapply(1:length(data_final), function(sp){
    plot(c(-4,4),c(4,-7), xlab="precip",ylab="temp",main=data_final[[sp]]$species, col=4);
    points(data_final[[sp]]$X1[data_final[[sp]]$y==0],data_final[[sp]]$X2[data_final[[sp]]$y==0],col=2);
    points(data_final[[sp]]$X1[data_final[[sp]]$y==1],data_final[[sp]]$X2[data_final[[sp]]$y==1]);
  })
  dev.off()
}
#reorder list of lists to match tree

{
  pdf(file= paste("Data_final_Espace_hist.pdf"))

  par(mfrow=c(3,4))
  lapply(1:length(data_final), function(sp){
    hist(data_final[[sp]]$X1[data_final[[sp]]$y==0],col=2, xlab="precip",ylab="temp",main=data_final[[sp]]$species);
    hist(data_final[[sp]]$X1[data_final[[sp]]$y==1], xlab="precip",ylab="temp",main=data_final[[sp]]$species);
  })
  dev.off()
}



