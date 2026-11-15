# Atmospheric gas exploratory analysis, imputation, PCA and clustering
# Input: data/MaunaLoa_miss.csv (authorised copy only; see data/README.md).
# This is exploratory multivariate analysis, not time-series forecasting.
required <- c("readr","dplyr","naniar","mice","FactoMineR","factoextra","cluster")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly=TRUE)]
if(length(missing)) stop("Install packages: ", paste(missing, collapse=", "))
suppressPackageStartupMessages({library(readr);library(dplyr);library(naniar);library(mice);library(FactoMineR);library(factoextra);library(cluster)})
input <- file.path("data","MaunaLoa_miss.csv")
if(!file.exists(input)) stop("Missing ",input,"; review data/README.md.")
out <- "outputs"; dir.create(out,showWarnings=FALSE,recursive=TRUE)
raw <- read_csv(input,show_col_types=FALSE)
if("Date" %in% names(raw)) raw <- raw |> mutate(Date=as.Date(Date))
x <- raw |> select(where(is.numeric))
if(ncol(x)<2) stop("At least two numeric gas variables are required.")
# IQR flags are diagnostic. Replacing flagged observations with NA reproduces
# the coursework sensitivity choice; it is not asserted to be scientifically optimal.
flags <- lapply(x,function(v){q<-quantile(v,c(.25,.75),na.rm=TRUE); d<-IQR(v,na.rm=TRUE); !is.na(v)&(v<q[1]-1.5*d|v>q[2]+1.5*d)})
write.csv(data.frame(variable=names(flags),flagged_n=vapply(flags,sum,integer(1))),file.path(out,"iqr_outlier_flags.csv"),row.names=FALSE)
for(nm in names(flags)) x[[nm]][flags[[nm]]] <- NA_real_
write.csv(data.frame(variable=names(x),missing_n=colSums(is.na(x)),missing_pct=round(colMeans(is.na(x))*100,2)),file.path(out,"missingness_summary.csv"),row.names=FALSE)
# A non-significant MCAR test does not prove MCAR.
write.csv(as.data.frame(naniar::mcar_test(x)),file.path(out,"mcar_test.csv"),row.names=FALSE)
set.seed(123); imp <- mice(x,m=10,method="pmm",seed=123,printFlag=FALSE)
# Keep first completed dataset for coursework comparability; this is not pooled inference.
completed <- complete(imp,1)
write.csv(completed,file.path(out,"imputed_dataset_1.csv"),row.names=FALSE)
pca <- FactoMineR::PCA(completed,scale.unit=TRUE,graph=FALSE)
eig <- as.data.frame(pca$eig); eig$component <- seq_len(nrow(eig))
write.csv(eig,file.path(out,"pca_eigenvalues.csv"),row.names=FALSE)
png(file.path(out,"pca_scree_plot.png"),1200,800,res=140); print(factoextra::fviz_eig(pca,addlabels=TRUE)); dev.off()
ncomp <- which(cumsum(pca$eig[,2])>=85)[1]; if(is.na(ncomp)) ncomp <- ncol(pca$ind$coord)
coords <- as.data.frame(pca$ind$coord[,seq_len(ncomp),drop=FALSE])
# Coursework selected k=3 visually; silhouette is reported as a diagnostic.
set.seed(123); k <- 3L; km <- kmeans(coords,centers=k,nstart=25)
sil <- cluster::silhouette(km$cluster,dist(coords))
write.csv(data.frame(method="kmeans",k=k,between_within_ratio_pct=100*km$betweenss/km$totss,mean_silhouette=mean(sil[,"sil_width"])),file.path(out,"kmeans_summary.csv"),row.names=FALSE)
write.csv(data.frame(observation=seq_len(nrow(coords)),cluster=km$cluster),file.path(out,"kmeans_assignments.csv"),row.names=FALSE)
hc <- hclust(dist(coords),method="ward.D2"); groups <- cutree(hc,k=k)
write.csv(data.frame(observation=seq_len(nrow(coords)),cluster=groups),file.path(out,"hierarchical_assignments.csv"),row.names=FALSE)
png(file.path(out,"hierarchical_dendrogram.png"),1200,800,res=140); plot(hc,labels=FALSE,main=paste("Ward.D2 clustering, k =",k)); rect.hclust(hc,k=k,border=2:4); dev.off()
message("Outputs saved to ",out,". PCA/clustering use imputation 1 only; temporal dependence is not modelled.")
