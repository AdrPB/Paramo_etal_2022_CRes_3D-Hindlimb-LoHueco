##########
## Three-dimensional analysis of the titanosaurian limb skeleton: 
## implications for systematic analysis
##
## Paramo, A.*, Mocho, P., Ortega, F.
##  
## This study was carried under R version 1.6.1
## The code for image reproduction have not been included, only the
## lines for the analyses.
## 
##########

## 0. Set the workspace --------------------------------------------------

# The seed for this study is set to 500
set.seed(500)

# These are the initial libraries for conducting the main anayses
# but they may require additional libraries that are already installed
# and loaded with the next packages.

library(geomorph)
library(Morpho)
library(rgl)
library(MASS)

# Setting the working directory with a folder "data", "img", "plot",
# "meshes" and "results" help to streamline all the analyses.

wdir <- getwd()
datadir <- paste(wdir, "/data", sep= "")
imgdir <- file.path(wdir,"img")
plotdir <- file.path(wdir,"plot")
meshdir <- file.path(wdir,"meshes")
resultdir <- file.path(wdir,"results")

## 0.1. Custom functions -------------------------------------------------
#
# The current study also uses several custom functions, mostly wrappers of
# combined functions already available or streamlined analyses.
# They are a separated in thir own section so it can be minimized for
# R-studio or similar gui users, as some of them are long functions.
#
# These routines are also available from A. Páramo github repository
#
##

# This function computes an intralandmark distance matrix from
# two procrustes coordinates arrays.

ldkdist <- function (x,y)
{
  m <- matrix(ncol=1, nrow= dim(x)[1])
  for (i in 1:dim(x)[1])
  {
    m[i,] <- sqrt(((x[i,1]-y[i,1])^2) + 
                    ((x[i,2]-y[i,2])^2) + 
                    ((x[i,3]-y[i,3])^2))
  }
  
  d <- data.frame(ldk_distance= m, row.names = rownames(x))
  return(d)
}

# This is a function to obtain the configurations at the extreme
# values of each LD from a LDA. It is a modified version of
# Adams & Sherratt "plotTangentSpace", In: Adams et al. (2019)
# "geomorph" v3.1.2 package
#
# It also include a transformation from the within group variance
# matrix following Claude (2008) proposed code.

plotTangent.LDA <- function(LDA.object, LDA.data, orpdata, fac, atlas.lm, mesh)
{
  atlas <- atlas.lm
  A <- orpdata
  LDA.data <- LDA.data
  k <- dim(A)[2]
  p <- dim(A)[1]
  n <- dim(A)[3]
  
  ref <- mshape(A)
  mshap <- aperm(ref,c(2,1))
  dim(mshap) <- p*k
  refmesh <- rotmesh.onto(mesh= mesh, refmat= atlas, tarmat= ref, scale= TRUE)
  refmesh <- tps3d(x= refmesh$mesh, refmat= refmesh$yrot, tarmat= ref)
  shapes <- shape.names <- NULL
  
  mod <- lm(two.d.array(A)~fac)
  dfw <- p - length(levels(fac))
  SSw <- var(mod$residuals) * (p-1)
  VCVw <- SSw/dfw
  
  LD <- LDA.object$scaling
  LDs <- VCVw %*% LD
  
  ldadata <- data.frame(predict(object= LDA.object, newdata= LDA.data)$x)
  
  for (i in 1:ncol(ldadata)) {
    ldaxis.min <- min(ldadata[,i])
    ldaxis.max <- max(ldadata[,i])
    lda.min <- lda.max <- matrix(NA,nrow= dim(ldadata)[2], ncol= dim(mshap))
    lda.min[i,] <- as.matrix(ldaxis.min * LDs[,i]) + as.vector(mshap)
    lda.max[i,] <- as.matrix(ldaxis.max * LDs[,i]) + as.vector(mshap)
    shapes <- rbind(shapes, lda.min[i,], lda.max[i,])
    
    shape.names <- c(shape.names, paste("LD", i, "min", sep= ""), 
                     paste("LD", i, "max", sep=""))
  }
  
  
  shapes <- arrayspecs(shapes, p, k)
  dimnames(shapes)[c(1:2)] <- dimnames(A)[c(1:2)]
  dimnames(shapes)[[3]] <- unlist(shape.names)
  
  
  tps.mesh <- list(vector(length= dim(shapes)[3]))
  
  for (i in 1:dim(shapes)[3])
  {
    RotAtlas <- rotmesh.onto(mesh= refmesh, refmat= ref, tarmat= shapes[,,i], scale= TRUE)
    tps.mesh[[i]] <- tps3d(x= RotAtlas$mesh, refmat= RotAtlas$yrot, tarmat= shapes[,,i])
  }
  
  names(tps.mesh) <- unlist(dimnames(shapes)[3])
  
  out <- list(lda.shapes <- shapes, lda.mesh <- tps.mesh, ref <- ref)
  names(out) <- c("lda.shapes", "lda.mesh", "ref")
  
  return(out)
}

# This function is used for mass reconstruction of a complete
# specimen mesh from a set of estimated landmarks.
#
# For cases when the specimen mesh is incomplete and the 3D
# landmarks have been estimated, this function warps the template
# mesh to each new configuration producing a complete set of
# specimen meshes that can be used for sliding of semilandmarks,
# surface analyses, etc.
# The resultant meshes are also batch exported to a folder.
#

Restore.Virt <- function(Spc_list, Ldks, Atlas_lm, Atlas_mesh, folder)
{
  restore.spcs <- Ldks[,,unlist(Spc_list)]
  mesh_length <- length(unlist(Spc_list))
  mesh_list <- ls()
  for (i in 1:mesh_length)
  {
    tarmat <- restore.spcs[,,unlist(Spc_list)[i]]
    rot_Atlas <- rotmesh.onto(mesh= Atlas_mesh, refmat= Atlas_lm, tarmat= tarmat)
    new.mesh <- tps3d(rot_Atlas$mesh, rot_Atlas$yrot, tarmat)
    vcgObjWrite(new.mesh,  
                filename= paste(folder,"/",unlist(Spc_list)[i],".obj", sep= ""), 
                writeNormals = TRUE)
  }
}

# The virtual restoration may produce errors in the
# meshes by any sample bias, specially in small samples.
# For comparison between the meshes, here we present a 
# wrapping function for assessing mesh differences between
# the original specimen 3D mesh representation and the virtual
# restoration of th specimen mesh.
#
# Both mshes are compared via vcgMetro algorithm included
# in the package Morpho. Following XXX (XX)
# The algorithm needs three folders. One folder were the initial
# meshes are placed. Other folder where all the meshes for
# comparison are stored. Lastly a folder where all the
# metro resultant meshes will be stored. These resulting meshes
# are the original specimens plus a texture equal to the differences
# in mm from the comparison mesh to the original mesh.
#

metro.polys <- function(Landmarks, folder.metro, folder.mesh, folder.mesh2) 
{
  metro.dist <- vector("list", dim(Landmarks)[3])
  metro.df <- data.frame(matrix(NA, nrow= dim(Landmarks)[3], ncol= 7), row.names= unlist(dimnames(Landmarks)[3][[1]]))
  
  for (i in 1:dim(Landmarks)[3])
  {
    
    mesh1 <- file2mesh(paste(folder.mesh, "/", unlist(dimnames(Landmarks)[3][[1]][i]),".obj",sep=""))
    mesh2 <- file2mesh(paste(folder.mesh2, "/", unlist(dimnames(Landmarks)[3][[1]][i]),".obj",sep=""))
    
    test.metro <- vcgMetro(mesh1, 
                           mesh2,
                           edgeSamp= FALSE,
                           faceSamp= FALSE)
    
    metro.df[i,] <- test.metro$ForwardSampling
    
    metro.dist[[i]] <- meshDist(mesh1, distvec= test.metro$distances1, 
                                file= paste(folder.metro,"/", unlist(dimnames(Landmarks)[3][[1]][i]),sep=""), save= TRUE, plot= FALSE)
    
    
  }
  
  listA <- vector("list")
  listA <- list(Metro.distances = metro.df, Mesh.dist = metro.dist)
  return(listA)
}

# A small function to generate the gradients from the
# interlandmark distance matrices.
#
# From: Daniel Hoop (Feb.23, 2016) as response to "Gradient of n colors ranging
# from color 1 and color 2". https://stackoverflow.com/q/13353213
#

color.gradient <- function(x, colors=c("red","yellow","green"), colsteps=100) 
{
  return( colorRampPalette(colors) (colsteps) [ findInterval(x, seq(min(x),max(x), length.out=colsteps)) ] )
}

# Mann-Whitney U Test for each PCA/LDA shape variable

MannU <- function(coords,groups,ro= 999)
{
  if (!ro)
  {
    ro <- 999
  }
  else
  {
    ro <- ro # The number of decimals after round()
  }
  
  coords <- data.frame(coords)
  coords$factor <- groups
  lvl <- as.character(levels(groups))
  
  if (length(lvl) > 2)
  {
    comb_names <- list()
    combs <- combn(lvl,2)
    
    for (j in 1:dim(combs)[2])
    {
      comb_names[j] <- as.character(paste(unlist(combs[1,j]),"_",unlist(combs[2,j]),sep=""))
    }
    
    results <- setNames(data.frame(matrix(ncol=(dim(coords)[2]-1), nrow=dim(combs)[2]), row.names= unlist(comb_names)),
                        nm= c(paste("PC",seq(from= 1, to= (dim(coords)[2]-1)),sep="")))
    for (l in 1:dim(results)[1])
    {
      for (i in 1:dim(results)[2])
      {
        results[l,i] <-  wilcox.test(coords[coords$factor %in% combs[,l],i]~coords$factor[coords$factor %in% combs[,l]])$p.value
      }
      results[l,] <- round(results[l,],ro)
    }
    return(results)
    
  }
  
  else
  {
    
    results <- setNames(data.frame(matrix(ncol=(dim(coords)[2]-1), nrow=1), row.names= paste(lvl[1],"_",lvl[2],sep="")),
                        nm= c(paste("PC",seq(from= 1, to= (dim(coords)[2]-1)),sep="")))
    for (i in 1:dim(results)[2])
    {
      results[,i] <-  wilcox.test(coords[,i]~coords$factor)$p.value
    }
    results <- round(results,ro)
    
    return(results)
  }
}

 
## 1. Load the data ------------------------------------------------------
#
# In this section all the data with the specimen labels
# and the landmark configurations will be loaded.
#
# The Atlas for each specimen will be also defined, as
# well as the vectors for each curve and the indices
# of the surface semilandmarks.
#
##

## 1.1. Specimen dataframes ----------------------------------------------

Femur.db <- data.frame(read.csv(paste(datadir, "/", "Femur_specimens", ".csv", sep=""),
                                sep =";", header= T, row.names= "Femur"))
Femur.db$Morphotype <- as.factor(Femur.db$Morphotype)

Tibia.db <- data.frame(read.csv(paste(datadir, "/", "Tibia_specimens", ".csv", sep=""),
                                sep =";", header= T, row.names= "Tibia"))
Tibia.db$Morphotype <- as.factor(Tibia.db$Morphotype)

Fibula.db <- data.frame(read.csv(paste(datadir, "/", "Fibula_specimens", ".csv", sep=""),
                                 sep =";", header= T, row.names= "Fibula"))
Fibula.db$Morphotype <- as.factor(Fibula.db$Morphotype)

Humerus.db <- data.frame(read.csv(paste(datadir, "/", "Humerus_specimens", ".csv", sep=""),
                                  sep =";", header= T, row.names= "Humerus"))
Humerus.db$Morphotype <- as.factor(Humerus.db$Morphotype)

Ulna.db <- data.frame(read.csv(paste(datadir, "/", "Ulna_specimens", ".csv", sep=""),
                               sep =";", header= T, row.names= "Ulna"))
Ulna.db$Morphotype <- as.factor(Ulna.db$Morphotype)

Radius.db <- data.frame(read.csv(paste(datadir, "/", "Radius_specimens", ".csv", sep=""),
                                 sep =";", header= T, row.names = "Radius"))
Radius.db$Morphotype <- as.factor(Radius.db$Morphotype)


## 1.2. Setting up the landmark definitions for each element type --------

# Femur landmark definitions 
       
Femur_model <- file2mesh(paste(meshdir,"/","Atlas_femur", ".obj", sep= ""))
Femur_Atlas.lm <- read.pts(paste(datadir, "/", "Atlas_femur", ".pts", sep= ""),na= 9999)
colnames(Femur_Atlas.lm) <- c(" x","y","z")

Femur_Atlas_surface <- file2mesh(paste(meshdir,"/","Atlas_femur_smldk mesh",".obj", sep= ""))
Fm_atlas_surf <- vert2points(Femur_Atlas_surface)

Fm.fix <- c(1:24)
Fm.supp <- c(25,26,27) # supp = supporting landmarks 

# Some elements have supporting landmarks.
# These landmarks are defined only for
# improving estimation and surface restoration.
# However, they are deemed type-III and
# therefore not useful for the analyses.
#
# They are removed after landamark estimation
# and mesh restoration.   

Fm.fix.supp <- c(1:27)
Fm.cur1 <- c(28:77)
Fm.cur2 <- c(78:97)
Fm.cur3 <- c(98:137)
Fm.cur4 <- c(138:197)
Fm.curves <- list(Fm.cur1,Fm.cur2,Fm.cur3,Fm.cur4)
Fm.cur <- c(Fm.cur1,Fm.cur2,Fm.cur3,Fm.cur4)
Fm.surf <- c((length(c(Fm.fix.supp,unlist(Fm.curves)))+1):length(c(Fm.fix.supp,Fm.cur,Fm_atlas_surf[,1])))

Femur_Atlas.complete.lm <- rbind(Femur_Atlas.lm, Fm_atlas_surf)
colnames(Femur_Atlas.complete.lm) <- c(" x","y","z")
rownames(Femur_Atlas.complete.lm) <- Fm.row.cmp <- c(paste("s",seq(1,length(Fm.fix.supp)),sep=""),
                                                     paste("c1_",seq(1,length(Fm.cur1)),sep=""),
                                                     paste("c2_",seq(1,length(Fm.cur2)),sep=""),
                                                     paste("c3_",seq(1,length(Fm.cur3)),sep=""),
                                                     paste("c4_",seq(1,length(Fm.cur4)),sep=""),
                                                     paste("sur",seq(1,(dim(Femur_Atlas.complete.lm)[1]-length(c(Fm.fix.supp,Fm.cur)))),sep=""))

# Tibia landmark definitions

Tibia_model <- file2mesh(paste(meshdir,"/","Atlas_Tibia", ".ply", sep= ""))
Tibia_Atlas.lm <- read.pts(paste(datadir, "/", "Atlas_Tibia", ".pts", sep= ""),na= 9999)
colnames(Tibia_Atlas.lm) <- c(" x","y","z")

Tibia_Atlas_surface <- file2mesh(paste(meshdir,"/","Atlas_Tibia_smldk mesh",".obj", sep= ""))
Tb_atlas_surf <- vert2points(Tibia_Atlas_surface)

Tb.fix <- c(1:14)
Tb.fix.supp <- c(1:16)
Tb.supp <- c(15,16)
Tb.cur1 <- c(17:57)
Tb.cur2 <- c(58:76)
Tb.cur3 <- c(77:116)
Tb.cur4 <- c(117:166)
Tb.curves <- list(Tb.cur1,Tb.cur2,Tb.cur3,Tb.cur4)
Tb.cur <- c(Tb.cur1,Tb.cur2,Tb.cur3,Tb.cur4)
Tb.surf <- c((length(c(Tb.fix.supp,unlist(Tb.curves)))+1):length(c(Tb.fix.supp,Tb.cur,Tb_atlas_surf[,1])))
 
Tibia_Atlas.complete.lm <- rbind(Tibia_Atlas.lm, Tb_atlas_surf)
colnames(Tibia_Atlas.complete.lm) <- c(" x","y","z")
rownames(Tibia_Atlas.complete.lm) <- Tb.row.cmp <- c(paste("s",seq(1,length(Tb.fix.supp)),sep=""),
                                                     paste("c1_",seq(1,length(Tb.cur1)),sep=""),
                                                     paste("c2_",seq(1,length(Tb.cur2)),sep=""),
                                                     paste("c3_",seq(1,length(Tb.cur3)),sep=""),
                                                     paste("c4_",seq(1,length(Tb.cur4)),sep=""),
                                                     paste("sur",seq(1,(dim(Tibia_Atlas.complete.lm)[1]-length(c(Tb.fix.supp,Tb.cur)))),sep=""))

# Fibula landmark definitions

Fibula_model <- file2mesh(paste(meshdir,"/","Atlas_Fibula", ".ply", sep= ""))
Fibula_Atlas.lm <- read.pts(paste(datadir, "/", "Atlas_Fibula", ".pts", sep= ""),na= 9999)
colnames(Fibula_Atlas.lm) <- c(" x","y","z")

Fibula_Atlas_surface <- file2mesh(paste(meshdir,"/","Atlas_Fibula_smldk mesh",".ply", sep= ""))
Fb_atlas_surf <- vert2points(Fibula_Atlas_surface)

Fb.fix <- c(1:10)
Fb.fix.supp <- c(1:12)
Fb.supp <- c(11,12)
Fb.cur1 <- c(13:52)
Fb.cur2 <- c(53:92)
Fb.cur3 <- c(93:136)
Fb.cur4 <- c(137:146)
Fb.cur5 <- c(147:186)
Fb.curves <- list(Fb.cur1,Fb.cur2,Fb.cur3,Fb.cur4,Fb.cur5)
Fb.cur <- c(Fb.cur1,Fb.cur2,Fb.cur3,Fb.cur4,Fb.cur5)
Fb.surf <- c((length(c(Fb.fix.supp,unlist(Fb.curves)))+1):length(c(Fb.fix.supp,Fb.cur,Fb_atlas_surf[,1])))

Fibula_Atlas.complete.lm <- rbind(Fibula_Atlas.lm, Fb_atlas_surf)
colnames(Fibula_Atlas.complete.lm) <- c(" x","y","z")
rownames(Fibula_Atlas.complete.lm) <- Fb.row.cmp <- c(paste("s",seq(1,length(Fb.fix.supp)),sep=""),
                                                      paste("c1_",seq(1,length(Fb.cur1)),sep=""),
                                                      paste("c2_",seq(1,length(Fb.cur2)),sep=""),
                                                      paste("c3_",seq(1,length(Fb.cur3)),sep=""),
                                                      paste("c4_",seq(1,length(Fb.cur4)),sep=""),
                                                      paste("c4_",seq(1,length(Fb.cur5)),sep=""),
                                                      paste("sur",seq(1,(dim(Fibula_Atlas.complete.lm)[1]-length(c(Fb.fix.supp,Fb.cur)))),sep=""))

# Humerus landmark definitions

Humerus_model <- file2mesh(paste(meshdir,"/","Atlas_Humerus", ".ply", sep= ""))
Humerus_Atlas.lm <- read.pts(paste(datadir, "/", "Atlas_Humerus", ".pts", sep= ""),na= 9999)
colnames(Humerus_Atlas.lm) <- c(" x","y","z")

Humerus_Atlas_surface <- file2mesh(paste(meshdir,"/","Atlas_Humerus_smldk mesh",".ply", sep= ""))
Hm_atlas_surf <- vert2points(Humerus_Atlas_surface)

Hm.fix <- c(1:18)
Hm.cur1 <- c(19:48)
Hm.cur2 <- c(49:68)
Hm.cur3 <- c(69:88)
Hm.cur4 <- c(89:128)
Hm.cur5 <- c(129:198)
Hm.curves <- list(Hm.cur1,Hm.cur2,Hm.cur3,Hm.cur4,Hm.cur5)
Hm.cur <- c(Hm.cur1,Hm.cur2,Hm.cur3,Hm.cur4,Hm.cur5)
Hm.surf <- c((length(c(Hm.fix,unlist(Hm.curves)))+1):length(c(Hm.fix,Hm.cur,Hm_atlas_surf[,1])))

Humerus_Atlas.complete.lm <- rbind(Humerus_Atlas.lm, Hm_atlas_surf)
colnames(Humerus_Atlas.complete.lm) <- c(" x","y","z")
rownames(Humerus_Atlas.complete.lm) <- Hm.row.cmp <- c(paste("s",seq(1,length(Hm.fix)),sep=""),
                                                       paste("c1_",seq(1,length(Hm.cur1)),sep=""),
                                                       paste("c2_",seq(1,length(Hm.cur2)),sep=""),
                                                       paste("c3_",seq(1,length(Hm.cur3)),sep=""),
                                                       paste("c4_",seq(1,length(Hm.cur4)),sep=""),
                                                       paste("c5_",seq(1,length(Hm.cur5)),sep=""),
                                                       paste("sur",seq(1,(dim(Humerus_Atlas.complete.lm)[1]-length(c(Hm.fix,Hm.cur)))),sep=""))

# Ulna landmark definitions

Ulna_model <- file2mesh(paste(meshdir,"/","Atlas_Ulna", ".ply", sep= ""))
Ulna_Atlas.lm <- read.pts(paste(datadir, "/", "Atlas_Ulna", ".pts", sep= ""),na= 9999)
colnames(Ulna_Atlas.lm) <- c(" x","y","z")

Ulna_Atlas_surface <- file2mesh(paste(meshdir,"/","Atlas_Ulna_smldk mesh",".obj", sep= ""))
Ul_atlas_surf <- vert2points(Ulna_Atlas_surface)
 
Ul.fix <- c(1:6,8:13)
Ul.fix.supp <- c(1:13)
Ul.cur1 <- c(14:73)
Ul.cur2 <- c(74:103)
Ul.cur3 <- c(104:133)
Ul.curves <- list(Ul.cur1,Ul.cur2,Ul.cur3)
Ul.cur <- c(Ul.cur1,Ul.cur2,Ul.cur3)
Ul.surf <- c((length(c(Ul.fix.supp,unlist(Ul.curves)))+1):length(c(Ul.fix.supp,Ul.cur,Ul_atlas_surf[,1])))
Ul.supp <- (7)

Ulna_Atlas.complete.lm <- rbind(Ulna_Atlas.lm, Ul_atlas_surf)
colnames(Ulna_Atlas.complete.lm) <- c(" x","y","z")
rownames(Ulna_Atlas.complete.lm) <- Ul.row.cmp <- c(paste("s",seq(1,length(Ul.fix.supp)),sep=""),
                                                    paste("c1_",seq(1,length(Ul.cur1)),sep=""),
                                                    paste("c2_",seq(1,length(Ul.cur2)),sep=""),
                                                    paste("c3_",seq(1,length(Ul.cur3)),sep=""),
                                                    paste("sur",seq(1,(dim(Ulna_Atlas.complete.lm)[1]-length(c(Ul.fix.supp,Ul.cur)))),sep=""))

# Radius landmark definitions

Radius_model <- file2mesh(paste(meshdir,"/","Atlas_Radius", ".ply", sep= ""))
Radius_Atlas.lm <- read.pts(paste(datadir, "/", "Atlas_Radius", ".pts", sep= ""),na= 9999)
colnames(Radius_Atlas.lm) <- c(" x","y","z")

Radius_Atlas_surface <- file2mesh(paste(meshdir,"/","Atlas_Radius_smldk mesh",".obj", sep= ""))
Rd_atlas_surf <- vert2points(Radius_Atlas_surface)

Rd.fix <- c(1:7)
Rd.cur2 <- c(8:35)
Rd.cur3 <- c(36:63)
Rd.cur1 <- c(64:109)
Rd.cur4 <- c(110:150)
Rd.curves <- list(Rd.cur1,Rd.cur2,Rd.cur3,Rd.cur4)
Rd.cur <- c(Rd.cur1,Rd.cur2,Rd.cur3,Rd.cur4)
Rd.surf <- c((length(c(Rd.fix,unlist(Rd.curves)))+1):length(c(Rd.fix,Rd.cur,Rd_atlas_surf[,1])))

Radius_Atlas.complete.lm <- rbind(Radius_Atlas.lm, Rd_atlas_surf)
colnames(Radius_Atlas.complete.lm) <- c(" x","y","z")
rownames(Radius_Atlas.complete.lm) <- Rd.row.cmp <- c(paste("s",seq(1,length(Rd.fix)),sep=""),
                                                      paste("c1_",seq(1,length(Rd.cur1)),sep=""),
                                                      paste("c2_",seq(1,length(Rd.cur2)),sep=""),
                                                      paste("c3_",seq(1,length(Rd.cur3)),sep=""),
                                                      paste("c4_",seq(1,length(Rd.cur4)),sep=""),
                                                      paste("sur",seq(1,(dim(Radius_Atlas.complete.lm)[1]-length(c(Rd.fix,Rd.cur)))),sep=""))

## 1.3. Import the landmark configurations -------------------------------

Humerus.ldk <- readland.tps(paste(datadir,"/Humerus_landmarks.tps"), specID = "ID")
Ulna.ldk <- readland.tps(paste(datadir,"/Ulna_landmarks.tps"), specID = "ID")
Radius.ldk <- readland.tps(paste(datadir,"/Radius_landmarks.tps"), specID = "ID")

Femur.ldk <- readland.tps(paste(datadir,"/Femur_landmarks.tps"), specID = "ID")
Tibia.ldk <- readland.tps(paste(datadir,"/Tibia_landmarks.tps"), specID = "ID")
Fibula.ldk <- readland.tps(paste(datadir,"/Fibula_landmarks.tps"), specID = "ID")

## 2. Estimate Missing Landmarks -----------------------------------------
#
# Here we estimate the missing landmarks in each configuration.
#
# This is possible with the use of Thin Plate Spline algorithm
# implemented via geomorph package.
#
##

Humerus.ildk <- estimate.missing(Humerus.ldk, method= "TPS")
Ulna.ildk <- estimate.missing(Ulna.ldk, method= "TPS")
Radius.ildk <- estimate.missing(Radius.ldk, method= "TPS")

Femur.ildk <- estimate.missing(Femur.ldk, method= "TPS")
Tibia.ildk <- estimate.missing(Tibia.ldk, method= "TPS")
Fibula.ildk <- estimate.missing(Fibula.ldk, method= "TPS")

## 2.1. Virtual statitical restoration of the specimen meshes ------------
#
# This procedure will warp the template mesh into the 
# estimated landmark configuration.
# The resulting complete 3D mesh restored for each
# specimen will be saved into a new folder.
#
##

Restore.Virt(Spc_list= rownames(Humerus.ldk), Ldks= Humerus.ildk, 
             Atlas_lm= Humerus_Atlas.lm, Atlas_mesh= Humerus_model, 
             folder= paste(meshdir,"/humerus_warp", sep=""))

Restore.Virt(Spc_list= rownames(Ulna.ldk), Ldks= Ulna.ildk, 
             Atlas_lm= Ulna_Atlas.lm[-Ul.supp,], Atlas_mesh= Ulna_model, 
             folder= paste(meshdir,"/ulna_warp", sep=""))

Restore.Virt(Spc_list= rownames(Radius.ldk), Ldks= Radius.ildk, 
             Atlas_lm= Radius_Atlas.lm, Atlas_mesh= Radius_model, 
             folder= paste(meshdir,"/radius_warp", sep=""))

Restore.Virt(Spc_list= rownames(Femur.ldk), Ldks= Femur.ildk, 
             Atlas_lm= Femur_Atlas.lm[-Fm.supp,], Atlas_mesh= Femur_model, 
             folder= paste(meshdir,"/femur_warp", sep=""))

Restore.Virt(Spc_list= rownames(Tibia.ldk), Ldks= Tibia.ildk, 
             Atlas_lm= Tibia_Atlas.lm[-Tb.supp,], Atlas_mesh= Tibia_model, 
             folder= paste(meshdir,"/tibia_warp", sep=""))

Restore.Virt(Spc_list= rownames(Fibula.ldk), Ldks= Tibia.ildk, 
             Atlas_lm= Fibula_Atlas.lm[-Fb.supp,], Atlas_mesh= Fibula_model, 
             folder= paste(meshdir,"/fibula_warp", sep=""))

## 2.2. Remove the supporting landmarks -------------------------------------
#
# The supporting landmarks are removed as they are
# no longer needed.
# The actual Atlas is defined for the surface semilandmark
# projection.
#
##


Fm.cur1.s <- c(25:74)
Fm.cur2.s <- c(75:94)
Fm.cur3.s <- c(95:134)
Fm.cur4.s <- c(135:194)

Fm.curves.s <- list(Fm.cur1.s,Fm.cur2.s,Fm.cur3.s,Fm.cur4.s)
Fm.cur.s <- c(Fm.cur1.s,Fm.cur2.s,Fm.cur3.s,Fm.cur4.s)
Fm.surf.s <- c((min(Fm.surf)-length(Fm.supp)):(max(Fm.surf)-length(Fm.supp)))

Tb.cur1.s <- c(15:55)
Tb.cur2.s <- c(56:74)
Tb.cur3.s <- c(75:114)
Tb.cur4.s <- c(115:164)

Tb.curves.s <- list(Tb.cur1.s,Tb.cur2.s,Tb.cur3.s,Tb.cur4.s)
Tb.cur.s <- c(Tb.cur1.s,Tb.cur2.s,Tb.cur3.s,Tb.cur4.s)
Tb.surf.s <- c((min(Tb.surf)-length(Tb.supp)):(max(Tb.surf)-length(Tb.supp)))

Fb.cur1.s <- c(11:50)
Fb.cur2.s <- c(51:90)
Fb.cur3.s <- c(91:134)
Fb.cur4.s <- c(135:144)
Fb.cur5.s <- c(145:184)

Fb.curves.s <- list(Fb.cur1.s,Fb.cur2.s,Fb.cur3.s,Fb.cur4.s,Fb.cur5.s)
Fb.cur.s <- c(Fb.cur1.s,Fb.cur2.s,Fb.cur3.s,Fb.cur4.s,Fb.cur5.s)
Fb.surf.s <- c((min(Fb.surf)-length(Fb.supp)):(max(Fb.surf)-length(Fb.supp)))

Ul.fix.s <- c(1:12)
Ul.cur1.s <- c(13:72)
Ul.cur2.s <- c(73:102)
Ul.cur3.s <- c(103:132)

Ul.curves.s <- list(Ul.cur1.s,Ul.cur2.s,Ul.cur3.s)
Ul.cur.s <- c(Ul.cur1.s,Ul.cur2.s,Ul.cur3.s)
Ul.surf.s <- c((min(Ul.surf)-length(Ul.supp)):(max(Ul.surf)-length(Ul.supp)))

## 2.3. Estimated landmarks and atlas creation ---------------------------
#
# Alternatvely, the set of estimated landmarks is also 
# provided and can be loaded into the working space.
#
##

Humerus.ildk <- readland.tps(paste(datadir,"/Humerus_estimated_landmarks.tps"), specID = "ID")
Ulna.ildk <- readland.tps(paste(datadir,"/Ulna_estimated_landmarks.tps"), specID = "ID")
Radius.ildk <- readland.tps(paste(datadir,"/Radius_estimated_landmarks.tps"), specID = "ID")

Femur.ildk <- readland.tps(paste(datadir,"/Femur_estimated_landmarks.tps"), specID = "ID")
Tibia.ildk <- readland.tps(paste(datadir,"/Tibia_estimated_landmarks.tps"), specID = "ID")
Fibula.ildk <- readland.tps(paste(datadir,"/Fibula_estimated_landmarks.tps"), specID = "ID")

## Atlas creation

Humerus_Atlas <- createAtlas(Humerus_model, 
                             landmarks= Humerus_Atlas.lm, 
                             patch= Humerus_Atlas.complete.lm[Hm.surf,],
                             corrCurves= Hm.curves,
                             keep.fix= Hm.fix)

Ulna_Atlas <- createAtlas(Ulna_model, 
                          landmarks= Ulna_Atlas.lm[-Ul.supp,], 
                          patch= Ulna_Atlas.complete.lm[Ul.surf,],
                          corrCurves= Ul.curves.s,
                          keep.fix= Ul.fix.s)

Radius_Atlas <- createAtlas(Radius_model, 
                            landmarks= Radius_Atlas.lm, 
                            patch= Radius_Atlas.complete.lm[Rd.surf,],
                            corrCurves= Rd.curves,
                            keep.fix= Rd.fix) 

Femur_Atlas <- createAtlas(Femur_model, 
                           landmarks= Femur_Atlas.lm[-Fm.supp,], 
                           patch= Femur_Atlas.complete.lm[Fm.surf,],
                           corrCurves= Fm.curves.s,
                           keep.fix= Fm.fix)

Tibia_Atlas <- createAtlas(Tibia_model, 
                           landmarks= Tibia_Atlas.lm[-Tb.supp,], 
                           patch= Tibia_Atlas.complete.lm[Tb.surf,],
                           keep.fix= c(Tb.fix, Tb.cur.s))

Fibula_Atlas <- createAtlas(Fibula_model, 
                            landmarks= Fibula_Atlas.complete.lm[c(Fb.fix, Fb.cur),], 
                            patch= Fibula_Atlas.complete.lm[Fb.surf.s,],
                            corrCurves= Fb.curves.s,
                            keep.fix= Fb.fix)
                            
## 2.3. High density surface semilandmark projection ---------------------

Humerus.ildk.complete <- placePatch(atlas= Humerus_Atlas, dat.array= Humerus.ildk, keep.fix= Hm.fix,
                                    path= paste(meshdir,"/","humerus_warp/",sep=""), fileext= ".obj")
Ulna.ildk.complete <- placePatch(atlas= Humerus_Atlas, dat.array= Ulna.ildk, keep.fix= Hm.fix,
                                 path= paste(meshdir,"/","humerus_warp/",sep=""), fileext= ".obj")
Radius.ildk.complete <- placePatch(atlas= Humerus_Atlas, dat.array= Radius.ildk, keep.fix= Hm.fix,
                                   path= paste(meshdir,"/","humerus_warp/",sep=""), fileext= ".obj")

Femur.ildk.complete <- placePatch(atlas= Humerus_Atlas, dat.array= Femur.ildk, keep.fix= Hm.fix,
                                  path= paste(meshdir,"/","humerus_warp/",sep=""), fileext= ".obj")
Tibia.ildk.complete <- placePatch(atlas= Humerus_Atlas, dat.array= Tibia.ildk, keep.fix= Hm.fix,
                                  path= paste(meshdir,"/","humerus_warp/",sep=""), fileext= ".obj")                                   
Fibula.ildk.complete <- placePatch(atlas= Humerus_Atlas, dat.array= Fibula.ildk, keep.fix= Hm.fix,
                                   path= paste(meshdir,"/","humerus_warp/",sep=""), fileext= ".obj")
                                   
## 3. Generalized Procrustes Analysis ------------------------------------                                                                      
Hm.gpa <- procSym(Humerus.ildk, SMvector= Hm.fix, outlines= Hm.curves, 
                  deselect= TRUE, recursive= TRUE, iterations= 10)
Rd.gpa <- procSym(Radius.ildk, SMvector= Rd.fix, outlines= Rd.curves, 
                  deselect= TRUE, recursive= TRUE, iterations= 10)
Ul.gpa <- procSym(Ulna.ildk, SMvector= Ul.fix.s, outlines= Ul.curves.s, 
                  deselect= TRUE, recursive= TRUE, iterations= 10)

# Also GPA aligment of the surface semilandmark datasets    
          
Hm.gpa.complete <- procSym(Humerus.ildk.complete, 
                           SMvector= Hm.fix, outlines= list(unlist(Hm.curves), Hm.surf),
                           deselect= TRUE, recursive= TRUE, iterations= 10)
Ul.gpa.complete <- procSym(Ulna.ildk.complete, 
                           SMvector= Ul.fix.s, outlines= list(unlist(Ul.curves.s), Ul.surf.s),
                           deselect= TRUE, recursive= TRUE, iterations= 10)
Rd.gpa.complete <- procSym(Radius.ildk.complete, 
                           SMvector= Rd.fix, outlines= list(unlist(Rd.curves), Rd.surf),
                           deselect= TRUE, recursive= TRUE, iterations= 10)
                           
Fm.gpa <- procSym(Femur.ildk, SMvector= Fm.fix, outlines= Fm.curves.s, 
                  deselect= TRUE, recursive= TRUE, iterations= 10)
Tb.gpa <- procSym(Tibia.ildk, outlines= Tb.curves.s, 
                  recursive= TRUE, iterations= 10)
Fb.gpa <- procSym(Fibula.ildk, SMvector= Fb.fix, outlines= Fb.curves.s, 
                  deselect= TRUE, recursive= TRUE, iterations= 10)                               

Fm.gpa.complete <- procSym(Femur.ildk.complete, 
                           SMvector= Fm.fix, outlines= list(unlist(Fm.curves.s), Fm.surf.s),
                           deselect= TRUE, recursive= TRUE, iterations= 10)
Tb.gpa.complete <- procSym(Tibia.ildk.complete, 
                           SMvector= Tb.fix, outlines= list(unlist(Tb.curves.s), Tb.surf.s),
                           deselect= TRUE, recursive= TRUE, iterations= 10)
Fb.gpa.complete <- procSym(Fibula.ildk.complete, 
                           SMvector= Fb.fix, outlines= list(unlist(Fb.curves.s), Fb.surf.s),
                           deselect= TRUE, recursive= TRUE, iterations= 10)   
                           
# PCA and LDA needs that the coordinates are in a p*k x m format,
# and the LDA especifically needs the "Clade" factor.
# We will create a dataframe with the transposed coordinates
# following two.d.array() function from package geomorpho.
# Also in each data frame we will include the "Clade" and "OTU"
# classification so they can be called from the same object. 

Hm.2darr.df <- data.frame(OTU= Humerus.db$Genus, 
                          Clade= Humerus.db$Clade,
                          two.d.array(Hm.gpa$orpdata),
                          row.names=dimnames(Hm.gpa$orpdata)[[3]])
Ul.2darr.df <- data.frame(OTU= Ulna.db$Genus, 
                          Clade= Ulna.db$Clade,
                          two.d.array(Ul.gpa$orpdata),
                          row.names=dimnames(Ul.gpa$orpdata)[[3]])
Rd.2darr.df <- data.frame(OTU= Radius.db$Genus, 
                          Clade= Radius.db$Clade,
                          two.d.array(Rd.gpa$orpdata),
                          row.names=dimnames(Rd.gpa$orpdata)[[3]])

Fm.2darr.df <- data.frame(OTU= Femur.db$Genus, 
                          Clade= Femur.db$Clade,
                          two.d.array(Fm.gpa$orpdata),
                          row.names=dimnames(Fm.gpa$orpdata)[[3]])
Tb.2darr.df <- data.frame(OTU= Tibia.db$Genus, 
                          Clade= Tibia.db$Clade,
                          two.d.array(Tb.gpa$orpdata),
                          row.names=dimnames(Tb.gpa$orpdata)[[3]])
Fb.2darr.df <- data.frame(OTU= Fibula.db$Genus, 
                          Clade= Fibula.db$Clade,
                          two.d.array(Fb.gpa$orpdata),
                          row.names=dimnames(Fb.gpa$orpdata)[[3]])

Hm.complete.2darr.df <- data.frame(OTU= Humerus.db$Genus, 
                                   Clade= Humerus.db$Clade,
                                   two.d.array(Hm.gpa.complete$orpdata),
                                   row.names=dimnames(Hm.gpa.complete$orpdata)[[3]])
Ul.complete.2darr.df <- data.frame(OTU= Ulna.db$Genus, 
                                   Clade= Ulna.db$Clade,
                                   two.d.array(Ul.gpa.complete$orpdata),
                                   row.names=dimnames(Ul.gpa.complete$orpdata)[[3]])
Rd.complete.2darr.df <- data.frame(OTU= Radius.db$Genus, 
                                   Clade= Radius.db$Clade,
                                   two.d.array(Rd.gpa.complete$orpdata),
                                   row.names=dimnames(Rd.gpa.complete$orpdata)[[3]])
                          
Fm.complete.2darr.df <- data.frame(OTU= Femur.db$Genus, 
                                   Clade= Femur.db$Clade,
                                   two.d.array(Fm.gpa.complete$orpdata),
                                   row.names=dimnames(Fm.gpa.complete$orpdata)[[3]]) 
Tb.complete.2darr.df <- data.frame(OTU= Tibia.db$Genus, 
                                   Clade= Tibia.db$Clade,
                                   two.d.array(Tb.gpa.complete$orpdata),
                                   row.names=dimnames(Tb.gpa.complete$orpdata)[[3]])                                                      
Fb.complete.2darr.df <- data.frame(OTU= Fibula.db$Genus, 
                                   Clade= Fibula.db$Clade,
                                   two.d.array(Fb.gpa.complete$orpdata),
                                   row.names=dimnames(Fb.gpa.complete$orpdata)[[3]])                           
                           
## 4. PCA and LDA of the procrustess aligned coordinates -----------------
#
# Now that all the configurations are aligned, rotated and scaled
# we can set-up the different analyses and visualizations.
#
##

Hm.PCA.ldk <- prcomp(Hm.2darr.df[,c(-1:-2)])
Hm.LDA.ldk <- lda(Clade ~., data= Hm.2darr.df[,-1], CV=FALSE) 

#
# In the LDA we need to obtain the specimen coordinates with pred()
#
Hm.LDA.pred <- data.frame(predict(object= Hm.LDA.ldk, newdata= Hm.2darr.df[-1])$x)

#
# We can visualize the results of the PCA and LDA.
# The PCA results as meshes are obtained with plotTangentSpace()
# from the geomorph package (Adams et al. 2019).
# The LDA results as meshes wre obtained with our custom function
# but not used in the current study.
# However, they help to visualize the actual plots.
#
Hm.LDA.projections <- plotTangent.LDA(LDA.object= Hm.LDA.ldk, LDA.data= Hm.2darr.df[,-1], fac= Hm.2darr.df$Clade,
                                      orpdata= Hm.gpa$orpdata, atlas= Humerus_Atlas.lm, 
                                      mesh= Humerus_model)
Hm.PCA.projections <- plotTangentSpace(Hm.gpa$orpdata)

Ul.PCA.ldk <- prcomp(Ul.2darr.df[,c(-1:-2)])
Ul.LDA.ldk <- lda(Clade ~., data= Ul.2darr.df[,-1], CV=FALSE) 
Ul.LDA.pred <- data.frame(predict(object= Ul.LDA.ldk, newdata= Ul.2darr.df[-1])$x)

Ul.LDA.projections <- plotTangent.LDA(LDA.object= Ul.LDA.ldk, LDA.data= Ul.2darr.df[,-1], fac= Ul.2darr.df$Clade,
                                      orpdata= Ul.gpa$orpdata, atlas= Ulna_Atlas.lm[-Ul.supp,], 
                                      mesh= Ulna_model)
Ul.PCA.projections <- plotTangentSpace(Ul.gpa$orpdata)

Rd.PCA.ldk <- prcomp(Rd.2darr.df[,c(-1:-2)])
Rd.LDA.ldk <- lda(Clade ~., data= Rd.2darr.df[,-1], CV=FALSE) 
Rd.LDA.pred <- data.frame(predict(object= Rd.LDA.ldk, newdata= Rd.2darr.df[-1])$x)

Rd.LDA.projections <- plotTangent.LDA(LDA.object= Rd.LDA.ldk, LDA.data= Rd.2darr.df[,-1], fac= Rd.2darr.df$Clade,
                                      orpdata= Rd.gpa$orpdata, atlas= Radius_Atlas.lm, 
                                      mesh= Radius_model)
Rd.PCA.projections <- plotTangentSpace(Rd.gpa$orpdata)

Fm.PCA.ldk <- prcomp(Fm.2darr.df[,c(-1:-2)])
Fm.LDA.ldk <- lda(Clade ~., data= Fm.2darr.df[,-1], CV=FALSE) 
Fm.LDA.pred <- data.frame(predict(object= Fm.LDA.ldk, newdata= Fm.2darr.df[-1])$x)

Fm.LDA.projections <- plotTangent.LDA(LDA.object= Fm.LDA.ldk, LDA.data= Fm.2darr.df[,-1], fac= Fm.2darr.df$Clade,
                                      orpdata= Fm.gpa$orpdata, atlas= Femur_Atlas.lm[-Fm.supp,], 
                                      mesh= Femur_model)
Fm.PCA.projections <- plotTangentSpace(Fm.gpa$orpdata)

Tb.PCA.ldk <- prcomp(Tb.2darr.df[,c(-1:-2)])
Tb.LDA.ldk <- lda(Clade ~., data= Tb.2darr.df[,-1], CV=FALSE) 
Tb.LDA.pred <- data.frame(predict(object= Tb.LDA.ldk, newdata= Tb.2darr.df[-1])$x)

Tb.LDA.projections <- plotTangent.LDA(LDA.object= Tb.LDA.ldk, LDA.data= Tb.2darr.df[,-1], fac= Tb.2darr.df$Clade,
                                      orpdata= Tb.gpa$orpdata, atlas= Tibia_Atlas.lm[-Tb.supp,], 
                                      mesh= Tibia_model)
Tb.PCA.projections <- plotTangentSpace(Tb.gpa$orpdata)

Fb.PCA.ldk <- prcomp(Fb.2darr.df[,c(-1:-2)])
Fb.LDA.ldk <- lda(Clade ~., data= Fb.2darr.df[,-1], CV=FALSE) 
Fb.LDA.pred <- data.frame(predict(object= Fb.LDA.ldk, newdata= Fb.2darr.df[-1])$x)

Fb.LDA.projections <- plotTangent.LDA(LDA.object= Fb.LDA.ldk, LDA.data= Fb.2darr.df[,-1], fac= Fb.2darr.df$Clade,
                                      orpdata= Fb.gpa$orpdata, atlas= Fibula_Atlas.lm[-Fb.supp,], 
                                      mesh= Fibula_model)
Fb.PCA.projections <- plotTangentSpace(Fb.gpa$orpdata)

## 4.1. Kruskal Wallis and Mann Whitney U's tests ------------------------
#
# Here we will conduct Kruskal Wallis test over the
# PCA and LDA shape variables in order to assess which
# axes plot significant differences between the clades
# and the different operative taxonomic units analysed.
#
# Also, a Mann Whitney U pairwise test allow us to
# assess which clades and OTUs are different in each axe.
# For this analysis it is necessary the MannU custom
# function.
#

# 4.1.1 Group differences pairwise in the PCA -----------------------------------------------------

Hm.ManU.PCA.OTU <- MannU(Hm.PCA.ldk$x, Humerus.db$OTU, ro=3)
Hm.ManU.PCA.Clade <- MannU(Hm.PCA.ldk$x, Humerus.db$Clade, ro=3)

Hm.KWallis.PCA.Clade <- setNames(data.frame(matrix(ncol=3,nrow=dim(Hm.PCA.ldk$x)[2]),
                                  row.names= paste("PC",seq(from= 1, to= dim(Hm.PCA.ldk$x)[2]),sep="")),
                             nm= c("Chi-sq","p-value","p-adjusted"))
for (i in 1:dim(Hm.PCA.ldk$x)[2])
{
  Hm.KWallis.PCA.Clade[i,] <- kruskal.test(Hm.PCA.ldk$x[,i], Humerus.db$Clade)[c(1,3)]
  Hm.KWallis.PCA.Clade[i,3] <- p.adjust(Hm.KWallis.PCA.Clade[i,2], method= "bonferroni", length(levels(Humerus.db$Clade))) # a post-hoc Bonferroni test allow to adjust the p-value
  Hm.KWallis.PCA.Clade[i,] <- round(Hm.KWallis.PCA.Clade[i,],3)
}

Hm.KWallis.PCA.OTU <- setNames(data.frame(matrix(ncol=3,nrow=dim(Hm.PCA.ldk$x)[2]),
                                      row.names= paste("PC",seq(from= 1, to= dim(Hm.PCA.ldk$x)[2]),sep="")),
                             nm= c("Chi-sq","p-value"))
for (i in 1:dim(Hm.PCA.ldk$x)[2])
{
  Hm.KWallis.PCA.OTU[i,] <- kruskal.test(Hm.PCA.ldk$x[,i], Humerus.db$Morphotype)[c(1,3)]
  Hm.KWallis.PCA.OTU[i,3] <- p.adjust(Hm.KWallis.PCA.OTU[i,2], method= "bonferroni", length(levels(Humerus.db$Morphotype)))
  Hm.KWallis.PCA.OTU[i,] <- round(Hm.KWallis.PCA.OTU[i,],3)
}


Hm.report.Clade <- data.frame(Hm.KWallis.PCA.Clade,t(Hm.ManU.PCA.Clade))
write.table(Hm.report.Clade, file= paste(resultdir,"/","Humerus PCA report by Clade.csv", sep=""),sep=";")

Hm.report.OTU <- data.frame(Hm.KWallis.PCA.OTU,t(Hm.ManU.PCA.OTU))
write.table(Hm.report.OTU, file= paste(resultdir,"/","Humerus PCA report by OTU.csv", sep=""),sep=";")

Ul.ManU.PCA.OTU <- MannU(Ul.PCA.ldk$x, Ulna.db[dimnames(Ul.PCA.ldk$x)[[1]],"Morphotype"], ro=3)
Ul.ManU.PCA.Clade <- MannU(Ul.PCA.ldk$x, Ulna.db[dimnames(Ul.PCA.ldk$x)[[1]],"Clade"], ro=3)

Ul.KWallis.PCA.Clade <- setNames(data.frame(matrix(ncol=3,nrow=dim(Ul.PCA.ldk$x)[2]),
                                        row.names= paste("PC",seq(from= 1, to= dim(Ul.PCA.ldk$x)[2]),sep="")),
                             nm= c("Chi-sq","p-value"))
for (i in 1:dim(Ul.PCA.ldk$x)[2])
{
  Ul.KWallis.PCA.Clade[i,] <- kruskal.test(Ul.PCA.ldk$x[,i], Ulna.db[dimnames(Ul.PCA.ldk$x)[[1]],"Clade"])[c(1,3)]
  Ul.KWallis.PCA.Clade[i,3] <- p.adjust(Ul.KWallis.PCA.Clade[i,2], method= "bonferroni", length(levels(Ulna.db$Clade)))
  Ul.KWallis.PCA.Clade[i,] <- round(Ul.KWallis.PCA.Clade[i,],3)
}

Ul.KWallis.PCA.OTU <- setNames(data.frame(matrix(ncol=3,nrow=dim(Ul.PCA.ldk$x)[2]),
                                      row.names= paste("PC",seq(from= 1, to= dim(Ul.PCA.ldk$x)[2]),sep="")),
                           nm= c("Chi-sq","p-value"))
for (i in 1:dim(Ul.PCA.ldk$x)[2])
{
  Ul.KWallis.PCA.OTU[i,] <- kruskal.test(Ul.PCA.ldk$x[,i], Ulna.db[dimnames(Ul.PCA.ldk$x)[[1]],"Morphotype"])[c(1,3)]
  Ul.KWallis.PCA.OTU[i,3] <- p.adjust(Ul.KWallis.PCA.OTU[i,2], method= "bonferroni", length(levels(Ulna.db$Morphotype)))
  Ul.KWallis.PCA.OTU[i,] <- round(Ul.KWallis.PCA.OTU[i,],3)
}


Ul.report.Clade <- data.frame(Ul.KWallis.PCA.Clade,t(Ul.ManU.PCA.Clade))
write.table(Ul.report.Clade, file= paste(resultdir,"/","Ulna PCA report by Clade.csv", sep=""),sep=";")

Ul.report.OTU <- data.frame(Ul.KWallis.PCA.OTU,t(Ul.ManU.PCA.OTU))
write.table(Ul.report.OTU, file= paste(resultdir,"/","Ulna PCA report by OTU.csv", sep=""),sep=";")


Rd.ManU.PCA.OTU <- MannU(Rd.PCA.ldk$x, Radius.db[dimnames(Rd.PCA.ldk$x)[[1]],"Morphotype"], ro=3)
Rd.ManU.PCA.Clade <- MannU(Rd.PCA.ldk$x, Radius.db[dimnames(Rd.PCA.ldk$x)[[1]],"Clade"], ro=3)

Rd.KWallis.PCA.Clade <- setNames(data.frame(matrix(ncol=3,nrow=dim(Rd.PCA.ldk$x)[2]),
                                        row.names= paste("PC",seq(from= 1, to= dim(Rd.PCA.ldk$x)[2]),sep="")),
                             nm= c("Chi-sq","p-value"))
for (i in 1:dim(Rd.PCA.ldk$x)[2])
{
  Rd.KWallis.PCA.Clade[i,] <- kruskal.test(Rd.PCA.ldk$x[,i], Radius.db[dimnames(Rd.PCA.ldk$x)[[1]],"Clade"])[c(1,3)]
  Rd.KWallis.PCA.Clade[i,3] <- p.adjust(Rd.KWallis.PCA.Clade[i,2], method= "bonferroni", length(levels(Radius.db$Clade)))
  Rd.KWallis.PCA.Clade[i,] <- round(Rd.KWallis.PCA.Clade[i,],3)
}

Rd.KWallis.PCA.OTU <- setNames(data.frame(matrix(ncol=3,nrow=dim(Rd.PCA.ldk$x)[2]),
                                      row.names= paste("PC",seq(from= 1, to= dim(Rd.PCA.ldk$x)[2]),sep="")),
                           nm= c("Chi-sq","p-value"))
for (i in 1:dim(Rd.PCA.ldk$x)[2])
{
  Rd.KWallis.PCA.OTU[i,] <- kruskal.test(Rd.PCA.ldk$x[,i], Radius.db[dimnames(Rd.PCA.ldk$x)[[1]],"Morphotype"])[c(1,3)]
  Rd.KWallis.PCA.OTU[i,3] <- p.adjust(Rd.KWallis.PCA.OTU[i,2], method= "bonferroni", length(levels(Radius.db$Morphotype)))
  Rd.KWallis.PCA.OTU[i,] <- round(Rd.KWallis.PCA.OTU[i,],3)
}


Rd.report.Clade <- data.frame(Rd.KWallis.PCA.Clade,t(Rd.ManU.PCA.Clade))
write.table(Rd.report.Clade, file= paste(resultdir,"/","Radius PCA report by Clade.csv", sep=""),sep=";")

Rd.report.OTU <- data.frame(Rd.KWallis.PCA.OTU,t(Rd.ManU.PCA.OTU))
write.table(Rd.report.OTU, file= paste(resultdir,"/","Radius PCA report by OTU.csv", sep=""),sep=";")

Fm.ManU.PCA.OTU <- MannU(Fm.PCA.ldk$x, Femur.db[dimnames(Fm.PCA.ldk$x)[[1]],"Morphotype"], ro=3)
Fm.ManU.PCA.Clade <- MannU(Fm.PCA.ldk$x, Femur.db[dimnames(Fm.PCA.ldk$x)[[1]],"Clade"], ro=3)

Fm.KWallis.PCA.Clade <- setNames(data.frame(matrix(ncol=3,nrow=dim(Fm.PCA.ldk$x)[2]),
                                        row.names= paste("PC",seq(from= 1, to= dim(Fm.PCA.ldk$x)[2]),sep="")),
                             nm= c("Chi-sq","p-value"))
for (i in 1:dim(Fm.PCA.ldk$x)[2])
{
  Fm.KWallis.PCA.Clade[i,] <- kruskal.test(Fm.PCA.ldk$x[,i], Femur.db[dimnames(Fm.PCA.ldk$x)[[1]],"Clade"])[c(1,3)]
  Fm.KWallis.PCA.Clade[i,3] <- p.adjust(Fm.KWallis.PCA.Clade[i,2], method= "bonferroni", length(levels(Femur.db$Clade)))
  Fm.KWallis.PCA.Clade[i,] <- round(Fm.KWallis.PCA.Clade[i,],3)
}

Fm.KWallis.PCA.OTU <- setNames(data.frame(matrix(ncol=3,nrow=dim(Fm.PCA.ldk$x)[2]),
                                      row.names= paste("PC",seq(from= 1, to= dim(Fm.PCA.ldk$x)[2]),sep="")),
                           nm= c("Chi-sq","p-value"))
for (i in 1:dim(Fm.PCA.ldk$x)[2])
{
  Fm.KWallis.PCA.OTU[i,] <- kruskal.test(Fm.PCA.ldk$x[,i], Femur.db[dimnames(Fm.PCA.ldk$x)[[1]],"Morphotype"])[c(1,3)]
  Fm.KWallis.PCA.OTU[i,3] <- p.adjust(Fm.KWallis.PCA.OTU[i,2], method= "bonferroni", length(levels(Femur.db$Morphotype)))
  Fm.KWallis.PCA.OTU[i,] <- round(Fm.KWallis.PCA.OTU[i,],3)
}


Fm.report.Clade <- data.frame(Fm.KWallis.PCA.Clade,t(Fm.ManU.PCA.Clade))
write.table(Fm.report.Clade, file= paste(resultdir,"/","Femur PCA report by Clade.csv", sep=""),sep=";")

Fm.report.OTU <- data.frame(Fm.KWallis.PCA.OTU,t(Fm.ManU.PCA.OTU))
write.table(Fm.report.OTU, file= paste(resultdir,"/","Femur PCA report by OTU.csv", sep=""),sep=";")


Tb.ManU.PCA.OTU <- MannU(Tb.PCA.ldk$x, Tibia.db[dimnames(Tb.PCA.ldk$x)[[1]],"Morphotype"], ro=3)
Tb.ManU.PCA.Clade <- MannU(Tb.PCA.ldk$x, Tibia.db[dimnames(Tb.PCA.ldk$x)[[1]],"Clade"], ro=3)

Tb.KWallis.PCA.Clade <- setNames(data.frame(matrix(ncol=3,nrow=dim(Tb.PCA.ldk$x)[2]),
                                        row.names= paste("PC",seq(from= 1, to= dim(Tb.PCA.ldk$x)[2]),sep="")),
                             nm= c("Chi-sq","p-value"))
for (i in 1:dim(Tb.PCA.ldk$x)[2])
{
  Tb.KWallis.PCA.Clade[i,] <- kruskal.test(Tb.PCA.ldk$x[,i], Tibia.db[dimnames(Tb.PCA.ldk$x)[[1]],"Clade"])[c(1,3)]
  Tb.KWallis.PCA.Clade[i,3] <- p.adjust(Tb.KWallis.PCA.Clade[i,2], method= "bonferroni", length(levels(Tibia.db$Morphotype)))
  Tb.KWallis.PCA.Clade[i,] <- round(Tb.KWallis.PCA.Clade[i,],3)
}

Tb.KWallis.PCA.OTU <- setNames(data.frame(matrix(ncol=3,nrow=dim(Tb.PCA.ldk$x)[2]),
                                      row.names= paste("PC",seq(from= 1, to= dim(Tb.PCA.ldk$x)[2]),sep="")),
                           nm= c("Chi-sq","p-value"))
for (i in 1:dim(Tb.PCA.ldk$x)[2])
{
  Tb.KWallis.PCA.OTU[i,] <- kruskal.test(Tb.PCA.ldk$x[,i], Tibia.db[dimnames(Tb.PCA.ldk$x)[[1]],"Morphotype"])[c(1,3)]
  Tb.KWallis.PCA.OTU[i,3] <- p.adjust(Tb.KWallis.PCA.OTU[i,2], method= "bonferroni", length(levels(Tibia.db$Morphotype)))
  Tb.KWallis.PCA.OTU[i,] <- round(Tb.KWallis.PCA.OTU[i,],3)
}


Tb.report.Clade <- data.frame(Tb.KWallis.PCA.Clade,t(Tb.ManU.PCA.Clade))
write.table(Tb.report.Clade, file= paste(resultdir,"/","Tibia PCA report by Clade.csv", sep=""),sep=";")

Tb.report.OTU <- data.frame(Tb.KWallis.PCA.OTU,t(Tb.ManU.PCA.OTU))
write.table(Tb.report.OTU, file= paste(resultdir,"/","Tibia PCA report by OTU.csv", sep=""),sep=";")


Fb.ManU.PCA.OTU <- MannU(Fb.PCA.ldk$x, Fibula.db[dimnames(Fb.PCA.ldk$x)[[1]],"Morphotype"], ro=3)
Fb.ManU.PCA.Clade <- MannU(Fb.PCA.ldk$x, Fibula.db[dimnames(Fb.PCA.ldk$x)[[1]],"Clade"], ro=3)

Fb.KWallis.PCA.Clade <- setNames(data.frame(matrix(ncol=3,nrow=dim(Fb.PCA.ldk$x)[2]),
                                        row.names= paste("PC",seq(from= 1, to= dim(Fb.PCA.ldk$x)[2]),sep="")),
                             nm= c("Chi-sq","p-value"))
for (i in 1:dim(Fb.PCA.ldk$x)[2])
{
  Fb.KWallis.PCA.Clade[i,] <- kruskal.test(Fb.PCA.ldk$x[,i], Fibula.db[dimnames(Fb.PCA.ldk$x)[[1]],"Clade"])[c(1,3)]
  Fb.KWallis.PCA.Clade[i,3] <- p.adjust(Fb.KWallis.PCA.Clade[i,2], method= "bonferroni", length(levels(Fibula.db$Morphotype)))
  Fb.KWallis.PCA.Clade[i,] <- round(Fb.KWallis.PCA.Clade[i,],3)
}

Fb.KWallis.PCA.OTU <- setNames(data.frame(matrix(ncol=3,nrow=dim(Fb.PCA.ldk$x)[2]),
                                      row.names= paste("PC",seq(from= 1, to= dim(Fb.PCA.ldk$x)[2]),sep="")),
                           nm= c("Chi-sq","p-value"))
for (i in 1:dim(Fb.PCA.ldk$x)[2])
{
  Fb.KWallis.PCA.OTU[i,] <- kruskal.test(Fb.PCA.ldk$x[,i], Fibula.db[dimnames(Fb.PCA.ldk$x)[[1]],"Morphotype"])[c(1,3)]
  Fb.KWallis.PCA.OTU[i,3] <- p.adjust(Fb.KWallis.PCA.OTU[i,2], method= "bonferroni", length(levels(Fibula.db$Morphotype)))
  Fb.KWallis.PCA.OTU[i,] <- round(Fb.KWallis.PCA.OTU[i,],3)
}


Fb.report.Clade <- data.frame(Fb.KWallis.PCA.Clade,t(Fb.ManU.PCA.Clade))
write.table(Fb.report.Clade, file= paste(resultdir,"/","Fibula PCA report by Clade.csv", sep=""),sep=";")

Fb.report.OTU <- data.frame(Fb.KWallis.PCA.OTU,t(Fb.ManU.PCA.OTU))
write.table(Fb.report.OTU, file= paste(resultdir,"/","Fibula PCA report by OTU.csv", sep=""),sep=";")

# 4.1.2. Group differences pairwise in the LDA -----------------------------------------------------

Hm.ManU.LDA.OTU <- MannU(Hm.LDA.pred, Humerus.db$OTU, ro=3)
Hm.ManU.LDA.Clade <- MannU(Hm.LDA.pred, Humerus.db$Clade, ro=3)

Hm.KWallis.LDA.Clade <- setNames(data.frame(matrix(ncol=3,nrow=dim(Hm.LDA.pred)[2]),
                                        row.names= paste("PC",seq(from= 1, to= dim(Hm.LDA.pred)[2]),sep="")),
                             nm= c("Chi-sq","p-value","p-adjusted"))
for (i in 1:dim(Hm.LDA.pred)[2])
{
  Hm.KWallis.LDA.Clade[i,] <- kruskal.test(Hm.LDA.pred[,i], Humerus.db$Clade)[c(1,3)]
  Hm.KWallis.LDA.Clade[i,3] <- p.adjust(Hm.KWallis.LDA.Clade[i,2], method= "bonferroni", length(levels(Humerus.db$Clade)))
  Hm.KWallis.LDA.Clade[i,] <- round(Hm.KWallis.LDA.Clade[i,],3)
}

Hm.KWallis.LDA.OTU <- setNames(data.frame(matrix(ncol=3,nrow=dim(Hm.LDA.pred)[2]),
                                      row.names= paste("PC",seq(from= 1, to= dim(Hm.LDA.pred)[2]),sep="")),
                           nm= c("Chi-sq","p-value"))
for (i in 1:dim(Hm.LDA.pred)[2])
{
  Hm.KWallis.LDA.OTU[i,] <- kruskal.test(Hm.LDA.pred[,i], Humerus.db$Morphotype)[c(1,3)]
  Hm.KWallis.LDA.OTU[i,3] <- p.adjust(Hm.KWallis.LDA.OTU[i,2], method= "bonferroni", length(levels(Humerus.db$Morphotype)))
  Hm.KWallis.LDA.OTU[i,] <- round(Hm.KWallis.LDA.OTU[i,],3)
}


Hm.report.Clade <- data.frame(Hm.KWallis.LDA.Clade,t(Hm.ManU.LDA.Clade))
write.table(Hm.report.Clade, file= paste(resultdir,"/","Humerus LDA report by Clade.csv", sep=""),sep=";")

Hm.report.OTU <- data.frame(Hm.KWallis.LDA.OTU,t(Hm.ManU.LDA.OTU))
write.table(Hm.report.OTU, file= paste(resultdir,"/","Humerus LDA report by OTU.csv", sep=""),sep=";")

Ul.ManU.LDA.OTU <- MannU(Ul.LDA.pred, Ulna.db[dimnames(Ul.LDA.pred)[[1]],"Morphotype"], ro=3)
Ul.ManU.LDA.Clade <- MannU(Ul.LDA.pred, Ulna.db[dimnames(Ul.LDA.pred)[[1]],"Clade"], ro=3)

Ul.KWallis.LDA.Clade <- setNames(data.frame(matrix(ncol=3,nrow=dim(Ul.LDA.pred)[2]),
                                        row.names= paste("PC",seq(from= 1, to= dim(Ul.LDA.pred)[2]),sep="")),
                             nm= c("Chi-sq","p-value"))
for (i in 1:dim(Ul.LDA.pred)[2])
{
  Ul.KWallis.LDA.Clade[i,] <- kruskal.test(Ul.LDA.pred[,i], Ulna.db[dimnames(Ul.LDA.pred)[[1]],"Clade"])[c(1,3)]
  Ul.KWallis.LDA.Clade[i,3] <- p.adjust(Ul.KWallis.LDA.Clade[i,2], method= "bonferroni", length(levels(Ulna.db$Clade)))
  Ul.KWallis.LDA.Clade[i,] <- round(Ul.KWallis.LDA.Clade[i,],3)
}

Ul.KWallis.LDA.OTU <- setNames(data.frame(matrix(ncol=3,nrow=dim(Ul.LDA.pred)[2]),
                                      row.names= paste("PC",seq(from= 1, to= dim(Ul.LDA.pred)[2]),sep="")),
                           nm= c("Chi-sq","p-value"))
for (i in 1:dim(Ul.LDA.pred)[2])
{
  Ul.KWallis.LDA.OTU[i,] <- kruskal.test(Ul.LDA.pred[,i], Ulna.db[dimnames(Ul.LDA.pred)[[1]],"Morphotype"])[c(1,3)]
  Ul.KWallis.LDA.OTU[i,3] <- p.adjust(Ul.KWallis.LDA.OTU[i,2], method= "bonferroni", length(levels(Ulna.db$Morphotype)))
  Ul.KWallis.LDA.OTU[i,] <- round(Ul.KWallis.LDA.OTU[i,],3)
}


Ul.report.Clade <- data.frame(Ul.KWallis.LDA.Clade,t(Ul.ManU.LDA.Clade))
write.table(Ul.report.Clade, file= paste(resultdir,"/","Ulna LDA report by Clade.csv", sep=""),sep=";")

Ul.report.OTU <- data.frame(Ul.KWallis.LDA.OTU,t(Ul.ManU.LDA.OTU))
write.table(Ul.report.OTU, file= paste(resultdir,"/","Ulna LDA report by OTU.csv", sep=""),sep=";")


Rd.ManU.LDA.OTU <- MannU(Rd.LDA.pred, Radius.db[dimnames(Rd.LDA.pred)[[1]],"Morphotype"], ro=3)
Rd.ManU.LDA.Clade <- MannU(Rd.LDA.pred, Radius.db[dimnames(Rd.LDA.pred)[[1]],"Clade"], ro=3)

Rd.KWallis.LDA.Clade <- setNames(data.frame(matrix(ncol=3,nrow=dim(Rd.LDA.pred)[2]),
                                        row.names= paste("PC",seq(from= 1, to= dim(Rd.LDA.pred)[2]),sep="")),
                             nm= c("Chi-sq","p-value"))
for (i in 1:dim(Rd.LDA.pred)[2])
{
  Rd.KWallis.LDA.Clade[i,] <- kruskal.test(Rd.LDA.pred[,i], Radius.db[dimnames(Rd.LDA.pred)[[1]],"Clade"])[c(1,3)]
  Rd.KWallis.LDA.Clade[i,3] <- p.adjust(Rd.KWallis.LDA.Clade[i,2], method= "bonferroni", length(levels(Radius.db$Clade)))
  Rd.KWallis.LDA.Clade[i,] <- round(Rd.KWallis.LDA.Clade[i,],3)
}

Rd.KWallis.LDA.OTU <- setNames(data.frame(matrix(ncol=3,nrow=dim(Rd.LDA.pred)[2]),
                                      row.names= paste("PC",seq(from= 1, to= dim(Rd.LDA.pred)[2]),sep="")),
                           nm= c("Chi-sq","p-value"))
for (i in 1:dim(Rd.LDA.pred)[2])
{
  Rd.KWallis.LDA.OTU[i,] <- kruskal.test(Rd.LDA.pred[,i], Radius.db[dimnames(Rd.LDA.pred)[[1]],"Morphotype"])[c(1,3)]
  Rd.KWallis.LDA.OTU[i,3] <- p.adjust(Rd.KWallis.LDA.OTU[i,2], method= "bonferroni", length(levels(Radius.db$Morphotype)))
  Rd.KWallis.LDA.OTU[i,] <- round(Rd.KWallis.LDA.OTU[i,],3)
}


Rd.report.Clade <- data.frame(Rd.KWallis.LDA.Clade,t(Rd.ManU.LDA.Clade))
write.table(Rd.report.Clade, file= paste(resultdir,"/","Radius LDA report by Clade.csv", sep=""),sep=";")

Rd.report.OTU <- data.frame(Rd.KWallis.LDA.OTU,t(Rd.ManU.LDA.OTU))
write.table(Rd.report.OTU, file= paste(resultdir,"/","Radius LDA report by OTU.csv", sep=""),sep=";")

Fm.ManU.LDA.OTU <- MannU(Fm.LDA.pred, Femur.db[dimnames(Fm.LDA.pred)[[1]],"Morphotype"], ro=3)
Fm.ManU.LDA.Clade <- MannU(Fm.LDA.pred, Femur.db[dimnames(Fm.LDA.pred)[[1]],"Clade"], ro=3)

Fm.KWallis.LDA.Clade <- setNames(data.frame(matrix(ncol=3,nrow=dim(Fm.LDA.pred)[2]),
                                        row.names= paste("PC",seq(from= 1, to= dim(Fm.LDA.pred)[2]),sep="")),
                             nm= c("Chi-sq","p-value"))
for (i in 1:dim(Fm.LDA.pred)[2])
{
  Fm.KWallis.LDA.Clade[i,] <- kruskal.test(Fm.LDA.pred[,i], Femur.db[dimnames(Fm.LDA.pred)[[1]],"Clade"])[c(1,3)]
  Fm.KWallis.LDA.Clade[i,3] <- p.adjust(Fm.KWallis.LDA.Clade[i,2], method= "bonferroni", length(levels(Femur.db$Clade)))
  Fm.KWallis.LDA.Clade[i,] <- round(Fm.KWallis.LDA.Clade[i,],3)
}

Fm.KWallis.LDA.OTU <- setNames(data.frame(matrix(ncol=3,nrow=dim(Fm.LDA.pred)[2]),
                                      row.names= paste("PC",seq(from= 1, to= dim(Fm.LDA.pred)[2]),sep="")),
                           nm= c("Chi-sq","p-value"))
for (i in 1:dim(Fm.LDA.pred)[2])
{
  Fm.KWallis.LDA.OTU[i,] <- kruskal.test(Fm.LDA.pred[,i], Femur.db[dimnames(Fm.LDA.pred)[[1]],"Morphotype"])[c(1,3)]
  Fm.KWallis.LDA.OTU[i,3] <- p.adjust(Fm.KWallis.LDA.OTU[i,2], method= "bonferroni", length(levels(Femur.db$Morphotype)))
  Fm.KWallis.LDA.OTU[i,] <- round(Fm.KWallis.LDA.OTU[i,],3)
}


Fm.report.Clade <- data.frame(Fm.KWallis.LDA.Clade,t(Fm.ManU.LDA.Clade))
write.table(Fm.report.Clade, file= paste(resultdir,"/","Femur LDA report by Clade.csv", sep=""),sep=";")

Fm.report.OTU <- data.frame(Fm.KWallis.LDA.OTU,t(Fm.ManU.LDA.OTU))
write.table(Fm.report.OTU, file= paste(resultdir,"/","Femur LDA report by OTU.csv", sep=""),sep=";")


Tb.ManU.LDA.OTU <- MannU(Tb.LDA.pred, Tibia.db[dimnames(Tb.LDA.pred)[[1]],"Morphotype"], ro=3)
Tb.ManU.LDA.Clade <- MannU(Tb.LDA.pred, Tibia.db[dimnames(Tb.LDA.pred)[[1]],"Clade"], ro=3)

Tb.KWallis.LDA.Clade <- setNames(data.frame(matrix(ncol=3,nrow=dim(Tb.LDA.pred)[2]),
                                        row.names= paste("PC",seq(from= 1, to= dim(Tb.LDA.pred)[2]),sep="")),
                             nm= c("Chi-sq","p-value"))
for (i in 1:dim(Tb.LDA.pred)[2])
{
  Tb.KWallis.LDA.Clade[i,] <- kruskal.test(Tb.LDA.pred[,i], Tibia.db[dimnames(Tb.LDA.pred)[[1]],"Clade"])[c(1,3)]
  Tb.KWallis.LDA.Clade[i,3] <- p.adjust(Tb.KWallis.LDA.Clade[i,2], method= "bonferroni", length(levels(Tibia.db$Morphotype)))
  Tb.KWallis.LDA.Clade[i,] <- round(Tb.KWallis.LDA.Clade[i,],3)
}

Tb.KWallis.LDA.OTU <- setNames(data.frame(matrix(ncol=3,nrow=dim(Tb.LDA.pred)[2]),
                                      row.names= paste("PC",seq(from= 1, to= dim(Tb.LDA.pred)[2]),sep="")),
                           nm= c("Chi-sq","p-value"))
for (i in 1:dim(Tb.LDA.pred)[2])
{
  Tb.KWallis.LDA.OTU[i,] <- kruskal.test(Tb.LDA.pred[,i], Tibia.db[dimnames(Tb.LDA.pred)[[1]],"Morphotype"])[c(1,3)]
  Tb.KWallis.LDA.OTU[i,3] <- p.adjust(Tb.KWallis.LDA.OTU[i,2], method= "bonferroni", length(levels(Tibia.db$Morphotype)))
  Tb.KWallis.LDA.OTU[i,] <- round(Tb.KWallis.LDA.OTU[i,],3)
}


Tb.report.Clade <- data.frame(Tb.KWallis.LDA.Clade,t(Tb.ManU.LDA.Clade))
write.table(Tb.report.Clade, file= paste(resultdir,"/","Tibia LDA report by Clade.csv", sep=""),sep=";")

Tb.report.OTU <- data.frame(Tb.KWallis.LDA.OTU,t(Tb.ManU.LDA.OTU))
write.table(Tb.report.OTU, file= paste(resultdir,"/","Tibia LDA report by OTU.csv", sep=""),sep=";")


Fb.ManU.LDA.OTU <- MannU(Fb.LDA.pred, Fibula.db[dimnames(Fb.LDA.pred)[[1]],"Morphotype"], ro=3)
Fb.ManU.LDA.Clade <- MannU(Fb.LDA.pred, Fibula.db[dimnames(Fb.LDA.pred)[[1]],"Clade"], ro=3)

Fb.KWallis.LDA.Clade <- setNames(data.frame(matrix(ncol=3,nrow=dim(Fb.LDA.pred)[2]),
                                        row.names= paste("PC",seq(from= 1, to= dim(Fb.LDA.pred)[2]),sep="")),
                             nm= c("Chi-sq","p-value"))
for (i in 1:dim(Fb.LDA.pred)[2])
{
  Fb.KWallis.LDA.Clade[i,] <- kruskal.test(Fb.LDA.pred[,i], Fibula.db[dimnames(Fb.LDA.pred)[[1]],"Clade"])[c(1,3)]
  Fb.KWallis.LDA.Clade[i,3] <- p.adjust(Fb.KWallis.LDA.Clade[i,2], method= "bonferroni", length(levels(Fibula.db$Morphotype)))
  Fb.KWallis.LDA.Clade[i,] <- round(Fb.KWallis.LDA.Clade[i,],3)
}

Fb.KWallis.LDA.OTU <- setNames(data.frame(matrix(ncol=3,nrow=dim(Fb.LDA.pred)[2]),
                                      row.names= paste("PC",seq(from= 1, to= dim(Fb.LDA.pred)[2]),sep="")),
                           nm= c("Chi-sq","p-value"))
for (i in 1:dim(Fb.LDA.pred)[2])
{
  Fb.KWallis.LDA.OTU[i,] <- kruskal.test(Fb.LDA.pred[,i], Fibula.db[dimnames(Fb.LDA.pred)[[1]],"Morphotype"])[c(1,3)]
  Fb.KWallis.LDA.OTU[i,3] <- p.adjust(Fb.KWallis.LDA.OTU[i,2], method= "bonferroni", length(levels(Fibula.db$Morphotype)))
  Fb.KWallis.LDA.OTU[i,] <- round(Fb.KWallis.LDA.OTU[i,],3)
}


Fb.report.Clade <- data.frame(Fb.KWallis.LDA.Clade,t(Fb.ManU.LDA.Clade))
write.table(Fb.report.Clade, file= paste(resultdir,"/","Fibula LDA report by Clade.csv", sep=""),sep=";")

Fb.report.OTU <- data.frame(Fb.KWallis.LDA.OTU,t(Fb.ManU.LDA.OTU))
write.table(Fb.report.OTU, file= paste(resultdir,"/","Fibula LDA report by OTU.csv", sep=""),sep=";")

## 5. High density surface semilandmark LDA ------------------------------
#
# As well as the LDA of the landmark and curve semilandmarks
# we proceed to do the LDA of the complete surfaces.
# We will use the predicted results with the custom function
# for obtaining the configurations along each LD.
# 
# Then we will calculate the procrustes distance between
# these extreme configurations.
#
##

Hm.complete.LDA.ldk <- lda(Clade ~., data= Hm.complete.2darr.df[,-1], CV=FALSE) 
Hm.complete.LDA.pred <- data.frame(predict(object= Hm.complete.LDA.ldk, 
                                           newdata= Hm.complete.2darr.df[-1])$x)
Hm.complete.LDA.projections <- plotTangent.LDA(LDA.object= Hm.complete.LDA.ldk, LDA.data= Hm.complete.2darr.df[,-1], fac= Hm.2darr.df$Clade,
                                      orpdata= Hm.gpa.complete$orpdata, atlas= Humerus_Atlas.complete.lm, 
                                      mesh= Humerus_model)

Hm.LD1surf.ldkdist <- ldkdist(Hm.complete.LDA.projections$lda.shapes[,,"LD1min"],Hm.complete.LDA.projections$lda.shapes[,,"LD1max"])
Hm.LD1surf.ldkdist <- unlist(Hm.LD1surf.ldkdist)
Hm.LD2surf.ldkdist <- ldkdist(Hm.complete.LDA.projections$lda.shapes[,,"LD2min"],Hm.complete.LDA.projections$lda.shapes[,,"LD2max"])
Hm.LD2surf.ldkdist <- unlist(Hm.LD2surf.ldkdist)

#
# The color scale for the surface semilandmarks is easy to create
# from the intralandmark distances.
# We will use the color.gradient() custom function.
# The palette selected was "viridis" from ViridisLite package
#
# Simon Garnier (2018). viridisLite: Default Color Maps from 'matplotlib' (Lite Version). R package version 0.3.0.
# https://CRAN.R-project.org/package=viridisLite
#                                                                 
color.Hm.LD1surf.ldkdist <- color.gradient(t(Hm.LD1surf.ldkdist[Hm.surf]), colors= viridisLite::viridis(3))
#
# This is added as example but omited for the Humerus LD2 and other
# results.
#

Ul.complete.LDA.ldk <- lda(Clade ~., data= Ul.complete.2darr.df[,-1], CV=FALSE) 
Ul.complete.LDA.pred <- data.frame(predict(object= Ul.complete.LDA.ldk, 
                                           newdata= Ul.complete.2darr.df[-1])$x)
Ul.complete.LDA.projections <- plotTangent.LDA(LDA.object= Ul.complete.LDA.ldk, LDA.data= Ul.complete.2darr.df[,-1], fac= Ul.2darr.df$Clade,
                                               orpdata= Ul.gpa.complete$orpdata, atlas= Ulna_Atlas.complete.lm, 
                                               mesh= Ulna_model)
                                               
Ul.LD1surf.ldkdist <- ldkdist(Ul.complete.LDA.projections$lda.shapes[,,"LD1min"],Ul.complete.LDA.projections$lda.shapes[,,"LD1max"])
Ul.LD1surf.ldkdist <- unlist(Ul.LD1surf.ldkdist)
Ul.LD2surf.ldkdist <- ldkdist(Ul.complete.LDA.projections$lda.shapes[,,"LD2min"],Ul.complete.LDA.projections$lda.shapes[,,"LD2max"])
Ul.LD2surf.ldkdist <- unlist(Ul.LD2surf.ldkdist)

Rd.complete.LDA.ldk <- lda(Clade ~., data= Rd.complete.2darr.df[,-1], CV=FALSE) 
Rd.complete.LDA.pred <- data.frame(predict(object= Rd.complete.LDA.ldk, 
                                           newdata= Rd.complete.2darr.df[-1])$x)
Rd.complete.LDA.projections <- plotTangent.LDA(LDA.object= Rd.complete.LDA.ldk, LDA.data= Rd.complete.2darr.df[,-1], fac= Rd.2darr.df$Clade,
                                               orpdata= Rd.gpa.complete$orpdata, atlas= Radius_Atlas.complete.lm, 
                                               mesh= Radius_model)
                                               
Rd.LD1surf.ldkdist <- ldkdist(Rd.complete.LDA.projections$lda.shapes[,,"LD1min"],Rd.complete.LDA.projections$lda.shapes[,,"LD1max"])
Rd.LD1surf.ldkdist <- unlist(Rd.LD1surf.ldkdist)
Rd.LD2surf.ldkdist <- ldkdist(Rd.complete.LDA.projections$lda.shapes[,,"LD2min"],Rd.complete.LDA.projections$lda.shapes[,,"LD2max"])
Rd.LD2surf.ldkdist <- unlist(Rd.LD2surf.ldkdist)

Fm.complete.LDA.ldk <- lda(Clade ~., data= Fm.complete.2darr.df[,-1], CV=FALSE) 
Fm.complete.LDA.pred <- data.frame(predict(object= Fm.complete.LDA.ldk, 
                                           newdata= Fm.complete.2darr.df[-1])$x)
Fm.complete.LDA.projections <- plotTangent.LDA(LDA.object= Fm.complete.LDA.ldk, LDA.data= Fm.complete.2darr.df[,-1], fac= Fm.2darr.df$Clade,
                                               orpdata= Fm.gpa.complete$orpdata, atlas= Femur_Atlas.complete.lm[-Fm.supp,], 
                                               mesh= Femur_model)
                                               
Fm.LD1surf.ldkdist <- ldkdist(Fm.complete.LDA.projections$lda.shapes[,,"LD1min"],Fm.complete.LDA.projections$lda.shapes[,,"LD1max"])
Fm.LD1surf.ldkdist <- unlist(Fm.LD1surf.ldkdist)
Fm.LD2surf.ldkdist <- ldkdist(Fm.complete.LDA.projections$lda.shapes[,,"LD2min"],Fm.complete.LDA.projections$lda.shapes[,,"LD2max"])
Fm.LD2surf.ldkdist <- unlist(Fm.LD2surf.ldkdist)


Tb.complete.LDA.ldk <- lda(Clade ~., data= Tb.complete.2darr.df[,-1], CV=FALSE) 
Tb.complete.LDA.pred <- data.frame(predict(object= Tb.complete.LDA.ldk, 
                                           newdata= Tb.complete.2darr.df[-1])$x)
Tb.complete.LDA.projections <- plotTangent.LDA(LDA.object= Tb.complete.LDA.ldk, LDA.data= Tb.complete.2darr.df[,-1], fac= Tb.2darr.df$Clade,
                                               orpdata= Tb.gpa.complete$orpdata, atlas= Tibia_Atlas.complete.lm, 
                                               mesh= Tibia_model)  
                                                
Tb.LD1surf.ldkdist <- ldkdist(Tb.complete.LDA.projections$lda.shapes[,,"LD1min"],Tb.complete.LDA.projections$lda.shapes[,,"LD1max"])
Tb.LD1surf.ldkdist <- unlist(Tb.LD1surf.ldkdist) 
Tb.LD2surf.ldkdist <- ldkdist(Tb.complete.LDA.projections$lda.shapes[,,"LD2min"],Tb.complete.LDA.projections$lda.shapes[,,"LD2max"])
Tb.LD2surf.ldkdist <- unlist(Tb.LD2surf.ldkdist)  


Fb.complete.2darr.df <- data.frame(OTU= Fibula.db$Genus, 
                                   Clade= Fibula.db$Clade,
                                   two.d.array(Fb.gpa.complete$orpdata),
                                   row.names=dimnames(Fb.gpa.complete$orpdata)[[3]])                                                        
Fb.complete.LDA.ldk <- lda(Clade ~., data= Fb.complete.2darr.df[,-1], CV=FALSE) 
Fb.complete.LDA.pred <- data.frame(predict(object= Fb.complete.LDA.ldk, 
                                           newdata= Fb.complete.2darr.df[-1])$x)
Fb.complete.LDA.projections <- plotTangent.LDA(LDA.object= Fb.complete.LDA.ldk, LDA.data= Fb.complete.2darr.df[,-1], fac= Fb.2darr.df$Clade,
                                               orpdata= Fb.gpa.complete$orpdata, atlas= Fibula_Atlas.complete.lm, 
                                               mesh= Fibula_model) 
                                               
Fb.LD1surf.ldkdist <- ldkdist(Fb.complete.LDA.projections$lda.shapes[,,"LD1min"],Fb.complete.LDA.projections$lda.shapes[,,"LD1max"])
Fb.LD1surf.ldkdist <- unlist(Fb.LD1surf.ldkdist)
Fb.LD2surf.ldkdist <- ldkdist(Fb.complete.LDA.projections$lda.shapes[,,"LD2min"],Fb.complete.LDA.projections$lda.shapes[,,"LD2max"])
Fb.LD2surf.ldkdist <- unlist(Fb.LD2surf.ldkdist)

## 6. Assess estimation errors -------------------------------------------
#
# Here we present as example how we calculated the mesh error
# and the landmark errors as distances for each specimens.
# We used the original 3D specimen meshes as reference against
# each of the TPS virtually restored 3D specimen mesh.
#
# To collect the results we created a folder named "*element_type*_metro"
# were the algorithm export a original mesh with a texture which
# colors the areas according to the distances between the original
# and the compared meshes.
#

Humerus.metro <- metro.polys(Humerus.ldk, 
                             folder.metro= paste(meshdir, "/humerus_metro",sep= ""),
                             folder.mesh= paste(meshdir, "/humerus",sep= ""),
                             folder.mesh2= paste(meshdir, "/humerus_warp",sep= ""))

#
# The error of the imputation method over the landmarks were
# assessed with our custom function and 50 bootstrap samples.
#
                             
Humerus.PDSE <- LdkPDSE(Humerus.ldk, boots= 50, method= "TPS")
                                  