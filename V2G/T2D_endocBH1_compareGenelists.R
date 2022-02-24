
library(LDlinkR)
library(tidyverse)

t2d_mah = read.csv("/mnt/isilon/sfgi/pahlm/projects/T2D_v2g/leadSNPs/leadSNPs_T2D_McCarthy.csv")
t2d_mah  = t2d_mah %>% mutate(sentinel = rsid, study = "Mahajan2018") %>% 
select(sentinel, chr, P, study, Ancestry)

t2d_vuj = read.csv("/mnt/isilon/sfgi/pahlm/projects/T2D_v2g/leadSNPs/leadSNPs_T2D_Voight.csv")
t2d_vuj  = t2d_vuj %>% mutate(sentinel = rsid, study = "Vujkovic2020") %>% 
select(sentinel, chr, P, study, Ancestry)

sent_list = rbind(t2d_vuj, t2d_mah)

sent_list = sent_list[grepl("rs",sent_list$sentinel),]

sent_list.chr = split(sent_list, sent_list$chr) 

out = lapply(sent_list.chr, function(variant_list){
	if(nrow(variant_list)>1){
	x= LDmatrix(variant_list$sentinel, pop="EUR", r2d = "r2",  token = "7d6ca3508279" )
	Sys.sleep(5)
	return(x)
	}else{
	y = data.frame(RS_number = variant_list$sentinel, temp = 1)
	row.names(y) = variant_list$sentinel 
	y
	}
})

out = lapply(out, function(x){
	x %>% gather(key = RS_number_2, value = R2, names(x)[-1])
})

out = do.call("rbind", out)

out = out[is.na(out$R2)==F,]

out = out[(out$RS_number != out$RS_number_2) & out$R2> 0.6,]

