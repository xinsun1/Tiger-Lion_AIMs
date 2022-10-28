#### PCAngsd tiger/lion AIMs ####
setwd('~/Documents/Projects/Tger:Lion_AIMs/1.pca/')

library(tidyverse)

C = as.matrix(read_table('./gl_tv_maf05_mis50.fwm_d100k_20k.cov', col_names = FALSE))
e = eigen(C)

spp = read_table('list.id', col_names = c("ID"))



p = data.frame(id=spp$ID,
               PC1=e$vectors[,1],
               PC2=e$vectors[,2],
               stringsAsFactors=F)


p %>%
    mutate(sp=substr(id,1,4)) %>%
    
    mutate(spp=case_when(sp=="Pleo" ~ "Lion",
                         sp=="Ptig" ~ "Tiger",
                         sp=="Pspa" ~ "Cave lion",
                         )) %>%
    
    
    ggplot(aes(x=PC1, y=PC2)) +
    geom_point(aes(color= spp, shape=spp), size =5, alpha = 0.8) +
   
    
    theme(panel.background = element_rect(fill = 'white', colour = "black", linetype = "solid"),
          # plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "right",
          aspect.ratio=1) +   # set main plot ratio
    xlab('PC1') + ylab('PC2')

ggsave('pca_gl_raw_8x8.png',width = 8, height = 8, device = "png", dpi = 500)  
ggsave('pca_gl_20K_8x8.png',width = 8, height = 8, device = "png", dpi = 500)  



#### snp weight filter ####
library(RcppCNPy)

snpw = npyLoad("./pca.gl_tv_maf05_mis50.snp_weight.weights.npy")
snpw = as.data.frame(snpw)

snpw %>%
    ggplot() +
    geom_density(aes(x=V1))

# take 99.9%
# use 0.00651
quantile(abs(snpw$V1), 0.6)
quantile(abs(snpw$V2), 0.99)
 
snp_index = snpw %>%
    rownames_to_column(var="index") %>%
    filter(abs(V1) > 0.001 | abs(V2) > 0.001) %>%
    select(index)
    
    # ggplot() +
    # geom_density(aes(x=V2))


snp_raw = read_table("./gl_tv_maf05_mis50.bed",
                     col_names = FALSE)
snp_f = read_table("./gl_tv_maf05_mis50.rep.trf.numt_f100.bed",
                   col_names = FALSE)
    
snp_f = snp_f %>%
    unite("id", c(X1,X2), remove = F)

# 1. get the top 0.01 SNPs with PCA weights
snp_f_w = snp_raw[snp_index$index,] %>%
    unite("id", c(X1,X2), remove = F)

# 2. masked SNP
snp_f_w_m = snp_f_w %>%
    filter(id %in% snp_f$id)

# 3. check distance
snp_f_w_m %>% 
    group_by(X1) %>%
    count() 


snpw %>%
    rownames_to_column(var="index") %>%
    filter(abs(V1) > 0.001 | abs(V2) > 0.001) %>%
    mutate(X1=snp_raw[index,]$X1) %>%
    group_by(X1) %>%
    count()


index_keep = as.matrix(rep(0,nrow(snp_f_w_m)))
index_keep[1,1] =1
t=snp_f_w_m[1,]
for(i in 2:nrow(snp_f_w_m)){
    # current line
    li = snp_f_w_m[i,]
    # distance to previous line
    add=0
    if(t$X1 == li$X1){
        if(li$X2 - t$X2 >= 100000){
            add=1
        }
    }else{
        add=1
    }
    
    if(add==1){
        index_keep[i,1] = 1
        t=li
    }
}

snp_fwm_d10k = snp_f_w_m[index_keep[,1]==1,]

write_tsv(snp_fwm_d10k, file='snp.fwm_d100k')
