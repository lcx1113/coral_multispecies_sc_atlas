find_pca_elbow = function(sdev){

	x<-1:length(sdev)
	eb<-data.frame(x,sdev)

	# piecewise linear model
	sse<-NULL
	for (k in 1:length(sdev)) {
	D=ifelse(x<=k,0,1)
	x2<-(x-k)*D
	sse[k]=sum(lm(sdev~x+x2)$residuals^2)
	}
	dim_plm<-which.min(sse)
	names(dim_plm)<-'Piecewise linear model'

	# first derivative
	df1<-diff(sdev,differences=1)
	den<-density(df1)
	rlev<-rle(diff(den$y)>0)$lengths
	if(length(rlev)<=2) {dim_fid<-max(which(df1<mean(df1)))+1
	} else {
	cutoff<-sum(rlev[-((length(rlev)-1):length(rlev))])
	dim_fid<-max(which(df1<den$x[cutoff]))+1 }
	names(dim_fid)<-'First derivative'

	# second derivative
	df2<-diff(sdev,differences=2)
	df2p<-df2[df2>0]
	den<-density(df2p)
	rlev<-rle(diff(den$y)>0)$lengths
	if(length(rlev)<=2) {dim_sed<-which.max(df2)+1
	} else {
	cutoff<-sum(rlev[1:2])
	dim_sed<-max(which(df2>den$x[cutoff]))+1 }
	names(dim_sed)<-'Second derivative'

	# preceding residual
	fit<-NULL;res<-NULL
	for (i in 1:(length(sdev)-2)) {
	fit[[i]]<-lm(sdev~x,data=eb[(i+1):length(sdev),])
	res[i]<-sdev[i]-predict(fit[[i]],newdata=data.frame(x=i))
	}
	den<-density(res)
	rlev<-rle(diff(den$y)>0)$lengths
	if(length(rlev)<=2) {dim_pr<-max(which(res>mean(res)))+1
	} else {
	cutoff<-sum(rlev[1:2])
	dim_pr<-max(which(res>den$x[cutoff]))+1 }
	names(dim_pr)<-'Preceding residual'

	# perpendicular line
	A<-c(1,sdev[1]);B<-c(length(sdev),sdev[length(sdev)]);Dist<-NULL
	for (i in 2:(length(sdev)-1)) {
	C<-c(i,sdev[i]);D<-cbind(rbind(A,B,C),rep(1,3))
	S<-1/2*abs(det(D));Dist[i]<-2*S/dist(rbind(A,B))
	}
	dim_perl<-which.max(Dist)
	names(dim_perl)<-'Perpendicular line'

	# k-means clustering
	set.seed(2022)
	dim_clu<-min(kmeans(sdev,2)$size)+1
	names(dim_clu)<-'K-means clustering'

	dim_all<-c(dim_plm,dim_fid,dim_sed,dim_pr,dim_perl,dim_clu)
	return(dim_all)
}
