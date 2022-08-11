

args = commandArgs(trailingOnly=TRUE)
ld_score_file_name = args[1]
output_base_name = args[2]

lds_seg = read.table(ld_score_file_name,header=T,colClasses=c("character",rep("numeric",8)))
quartiles=summary(lds_seg$ldscore_SNP)

lb1 = which(lds_seg$ldscore_SNP <= quartiles[2])
lb2 = which(lds_seg$ldscore_SNP > quartiles[2] & lds_seg$ldscore_SNP <= quartiles[3])
lb3 = which(lds_seg$ldscore_SNP > quartiles[3] & lds_seg$ldscore_SNP <= quartiles[5])
lb4 = which(lds_seg$ldscore_SNP > quartiles[5])

lb1_snp = lds_seg$SNP[lb1]
lb2_snp = lds_seg$SNP[lb2]
lb3_snp = lds_seg$SNP[lb3]
lb4_snp = lds_seg$SNP[lb4]

write.table(lb1_snp, paste(output_base_name, "_snp_group_1.txt", sep = ""), row.names=F, quote=F, col.names=F)
write.table(lb1_snp, paste(output_base_name, "_snp_group_2.txt", sep = ""), row.names=F, quote=F, col.names=F)
write.table(lb1_snp, paste(output_base_name, "_snp_group_3.txt", sep = ""), row.names=F, quote=F, col.names=F)
write.table(lb1_snp, paste(output_base_name, "_snp_group_4.txt", sep = ""), row.names=F, quote=F, col.names=F)
