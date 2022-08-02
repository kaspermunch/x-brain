
https://yanglab.westlake.edu.cn/software/gcta/#GREMLinWGSorimputeddata

## Step 1: segment based LD score

    gcta64 --bfile test --ld-score-region 200 --out test

## Step 2: stratify the SNPs by LD scores of individual SNPs in R

```R
lds_seg = read.table("test.score.ld",header=T,colClasses=c("character",rep("numeric",8)))
quartiles=summary(lds_seg$ldscore_SNP)

lb1 = which(lds_seg$ldscore_SNP <= quartiles[2])
lb2 = which(lds_seg$ldscore_SNP > quartiles[2] & lds_seg$ldscore_SNP <= quartiles[3])
lb3 = which(lds_seg$ldscore_SNP > quartiles[3] & lds_seg$ldscore_SNP <= quartiles[5])
lb4 = which(lds_seg$ldscore_SNP > quartiles[5])

lb1_snp = lds_seg$SNP[lb1]
lb2_snp = lds_seg$SNP[lb2]
lb3_snp = lds_seg$SNP[lb3]
lb4_snp = lds_seg$SNP[lb4]

write.table(lb1_snp, "snp_group1.txt", row.names=F, quote=F, col.names=F)
write.table(lb2_snp, "snp_group2.txt", row.names=F, quote=F, col.names=F)
write.table(lb3_snp, "snp_group3.txt", row.names=F, quote=F, col.names=F)
write.table(lb4_snp, "snp_group4.txt", row.names=F, quote=F, col.names=F)
```

## Step 3: making GRMs using SNPs stratified into different groups

    gcta64 --bfile test --extract snp_group1.txt --make-grm --out test_group1
    gcta64 --bfile test --extract snp_group2.txt --make-grm --out test_group2
    gcta64 --bfile test --extract snp_group3.txt --make-grm --out test_group3
    gcta64 --bfile test --extract snp_group4.txt --make-grm --out test_group4


or 

    gcta64 --bfile test --extract snp_group1.txt --chr 1 --make-grm --out test_chr1 --thread-num 10
    gcta64 --bfile test --extract snp_group1.txt --chr 2 --make-grm --out test_chr2 --thread-num 10
    ...

    gcta64 --bfile test --extract snp_group1.txt --chr 1 --make-grm --out test_chr1 --thread-num 10
    gcta64 --bfile test --extract snp_group1.txt --chr 2 --make-grm --out test_chr2 --thread-num 10
    ...



## Step 4: REML analysis with multiple GRMs

    gcta64 --reml --mgrm multi_GRMs.txt --pheno phen.txt --out test


Format of multi_GRMs.tx:
test_group1
test_group2
...    