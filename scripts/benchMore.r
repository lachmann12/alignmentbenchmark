
lines = readLines("supportfiles/avi.txt")


library("plyr")
library("stringr")

res = load("supportfiles/gene_mapping.rda")



getSalmon <- function(sra, mapping){
    file = paste0("quant/salmon/",sra,"/quant.sf")
    input = read.table(file, sep="\t", stringsAsFactors=F, header=T)[,c(1,5)]
    input[,1] = gsub("\\.[0-9]$","", input[,1])
    
    gene_match = match(input[,1], mapping[,1])
    
    transcript_counts = data.frame(mapping[gene_match,2], round(as.numeric(input[,2])))
    colnames(transcript_counts) = c("gene", "value")
    
    dd = ddply(transcript_counts,.(gene),summarize,sum=sum(value),number=length(gene))
    ge = dd[,2]
    names(ge) = dd[,1]
    return(ge)
}


getSalmonk15 <- function(sra, mapping){
    file = paste0("quant/salmon/",sra,"_k15/quant.sf")
    input = read.table(file, sep="\t", stringsAsFactors=F, header=T)[,c(1,5)]
    input[,1] = gsub("\\.[0-9]$","", input[,1])
    
    gene_match = match(input[,1], mapping[,1])
    
    transcript_counts = data.frame(mapping[gene_match,2], round(as.numeric(input[,2])))
    colnames(transcript_counts) = c("gene", "value")
    
    dd = ddply(transcript_counts,.(gene),summarize,sum=sum(value),number=length(gene))
    ge = dd[,2]
    names(ge) = dd[,1]
    return(ge)
}


getKallistok15 <- function(sra, mapping){
    file = paste0("quant/kallisto/",sra,"_k15/abundance.tsv")
    input = read.table(file, sep="\t", stringsAsFactors=F, header=T)[,c(1,4)]
    input[,1] = gsub("\\.[0-9]$","", input[,1])
    
    gene_match = match(input[,1], mapping[,1])
    
    transcript_counts = data.frame(mapping[gene_match,2], round(as.numeric(input[,2])))
    colnames(transcript_counts) = c("gene", "value")
    
    dd = ddply(transcript_counts,.(gene),summarize,sum=sum(value),number=length(gene))
    ge = dd[,2]
    names(ge) = dd[,1]
    return(ge)
}

exp = list()

for(i in lines){
    exp[[length(exp)+1]] = getKallistok15(i, mouse_map)
}

exp3 = do.call(cbind,exp)


exps15 = list()
for(i in lines){
    exps15[[length(exps15)+1]] = getSalmonk15(i, mouse_map)
}
exp15 = do.call(cbind,exps15)


exps31 = list()
for(i in lines){
    exps31[[length(exps31)+1]] = getSalmon(i, mouse_map)
}
exp31 = do.call(cbind,exps31)


inter = intersect(rownames(exp3), rownames(exp4))

diag(cor(exp3[inter,], exp4[inter,], method="spearman"))

diag(cor(exp31, exp4, method="spearman"))

diag(cor(exp2[inter,], exp4[inter,], method="spearman"))


diag(cor(exp2, exp3))

diag(cor(exp2, exp3, method="spearman"))

exp2 = do.call(cbind,exp)
fn = unlist(fr[,1])
rownames(exp) = fn
colnames(exp) = lines

colnames(exp2) = c("GSM3439006", "GSM3439007", "GSM3439008", 
"GSM3439009", "GSM3439010", "GSM3439011", "GSM3439012", "GSM3439013", 
"GSM3439014", "GSM3439015", "GSM3439016", "GSM3439017", "GSM3439018", 
"GSM3439019")

write.table(exp2, file="GSM3439006.tsv", quote=F, sep="\t")



kallisto = fread("supportfiles/kallisto_expression.tsv")
genes = unlist(kallisto[,1])
kallisto = as.matrix(kallisto[,-1])
rownames(kallisto) = genes




signatures = read.table("supportfiles/signatures.tsv", stringsAsFactors=F, sep="\t", fill=T, quote="", header=T)

whuman = which(signatures[,14] == "homo sapiens")

controlList = list()
perturbationList = list()

for(co in whuman){
    controlList[[length(controlList)+1]] = sort(unlist(c(strsplit(signatures[co, 6], ","), strsplit(signatures[co,10], ","))))
    perturbationList[[length(perturbationList)+1]] = sort(unlist(c(strsplit(signatures[co,7], ","), strsplit(signatures[co,9], ","))))
}
ww = which((lapply(perturbationList,length) == 3) & (lapply(controlList,length)  == 3))
controlList = controlList[ww]
perturbationList = perturbationList[ww]


imatch <- function(expr){
    innerMatchControl = list()
    for(i in 1:length(controlList)){
        inter = intersect(controlList[[i]], colnames(expr))
        inter2 = intersect(perturbationList[[i]], colnames(expr))
        
        if(length(inter) > 1 && length(inter2) > 1){
            cc = cor(expr[,inter])
            
            if(length(cc) > 1){
                ut = mean(unlist(cc[upper.tri(cc)]))
                innerMatchControl[[length(innerMatchControl)+1]] = ut
            }
        
            cc = cor(expr[,inter2])
            
            if(length(cc) > 1){
                ut = mean(unlist(cc[upper.tri(cc)]))
                #innerMatchControl[[length(innerMatchControl)+1]] = ut
            }
        }
    }
    return(unlist(innerMatchControl))
}


exmatch <- function(expr){
    innerMatchControl = list()
    for(i in 1:length(controlList)){
        inter = intersect(controlList[[i]], colnames(expr))
        inter2 = intersect(perturbationList[[i]], colnames(expr))
        
        if(length(inter) > 1 && length(inter2) > 1){
            cc = mean(unlist(cor(expr[,inter], expr[,inter2])))
            cc1 = mean(cor(expr[,inter]))
            innerMatchControl[[length(innerMatchControl)+1]] = cc1/cc
        }
        
    }
    return(unlist(innerMatchControl))
}




kallisto = fread("supportfiles/kallisto_expression.tsv")
genes = unlist(kallisto[,1])
kallisto = as.matrix(kallisto[,-1])
rownames(kallisto) = genes
ikallisto = imatch(kallisto)
exkallisto = exmatch(kallisto)

star = fread("supportfiles/star_expression.tsv")
genes = unlist(star[,1])
star = as.matrix(star[,-1])
rownames(star) = genes
istar = imatch(star)
exstar = exmatch(star)

hisat2 = fread("supportfiles/hisat2_expression.tsv")
genes = unlist(hisat2[,1])
hisat2 = as.matrix(hisat2[,-1])
rownames(hisat2) = genes
ihisat2 = imatch(hisat2)
exhisat2 = exmatch(hisat2)

bwa = fread("supportfiles/bwa_expression.tsv")
genes = unlist(bwa[,1])
bwa = as.matrix(bwa[,-1])
rownames(bwa) = genes
ibwa = imatch(bwa)
exbwa = exmatch(bwa)

salmon = fread("supportfiles/salmon_expression.tsv")
genes = unlist(salmon[,1])
salmon = as.matrix(salmon[,-1])
rownames(salmon) = genes
isalmon = imatch(salmon)
exsalmon = exmatch(salmon)

exs = list(exkallisto, exsalmon, exbwa, exhisat2, exstar)

par(mar=c(1,1,1,1))
par(mfrow=c(5,5))
for(i in 1:length(exs)){
    for(j in 1:length(exs)){
        tt = t.test(exs[[i]], exs[[j]])$p.value
        plot(exs[[i]], exs[[j]], pch=20, main=tt)
        abline(0,1, col="red", lwd=2)
    }
}

tt = t.test(exkallisto, exsalmon)$p.value
plot(exkallisto, exsalmon, pch=20, main=tt)
abline(0,1, col="red", lwd=2)

tt = t.test(exkallisto, exstar)$p.value
plot(exkallisto, exstar, pch=20, main=tt)
abline(0,1, col="red", lwd=2)

tt = t.test(exkallisto, exbwa)$p.value
plot(exkallisto, exbwa, pch=20, main=tt)
abline(0,1, col="red", lwd=2)

tt = t.test(exkallisto, exhisat2)$p.value
plot(exkallisto, exhisat2, pch=20, main=tt)
abline(0,1, col="red", lwd=2)


t.test(exkallisto, exhisat2)


plot(isalmon, istar, pch=20)
abline(0,1, col="red", lwd=2)

plot(isalmon, ibwa, pch=20)
abline(0,1, col="red", lwd=2)

plot(isalmon, ihisat2, pch=20)
abline(0,1, col="red", lwd=2)


plot(istar, ihisat2, pch=20)
abline(0,1, col="red", lwd=2)

plot(istar, ibwa, pch=20)
abline(0,1, col="red", lwd=2)

plot(istar, ikallisto, pch=20)
abline(0,1, col="red", lwd=2)

plot(istar, isalmon, pch=20)
abline(0,1, col="red", lwd=2)




