# Inter-chromosomal normalization, taken from Hi-C norm https://academic.oup.com/bioinformatics/article/28/23/3131/192582
normalize_mat_trans<-function(u, v_a, v_b){
#change matrix into vector
 u_vec<-c(as.matrix(u))
#get cov matrix
 len_m<-as.matrix(log(v_a[,eff_length]%o%v_b[,eff_length]))
 gcc_m<-as.matrix(log(v_a[,gc]%o%v_b[,gc]))
 map_m<-as.matrix(log(v_a[,cov]%o%v_b[,cov]))
#centralize cov matrix of enz, gcc
 len_m<-(len_m-mean(c(len_m)))/sd(c(len_m))
 gcc_m<-(gcc_m-mean(c(gcc_m)))/sd(c(gcc_m))
#change matrix into vector
 len_vec<-c(len_m)
 gcc_vec<-c(gcc_m)
 map_vec<-c(map_m)
#fit Poisson regression: u~len+gcc+offset(map)
 fit<-glm(u_vec~len_vec+gcc_vec+offset(map_vec),family="poisson")
#summary(fit)
 coeff <- round(fit$coeff,4)
 res <- round(u/exp(coeff[1]+coeff[2]*len_m+coeff[3]*gcc_m+map_m), 4)
 data.table(id1=v_a$id, res)->res
 melt(res, id.var="id1", value.name="nlinks_norm", variable.name="id2")->res
 res[nlinks_norm > 0]
}