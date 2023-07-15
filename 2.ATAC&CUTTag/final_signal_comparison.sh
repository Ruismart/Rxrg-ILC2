#### some results processed and generated in shell
## heatmap of Rxrg-bound peak regions        
# a. save Rxrg-CUTTag macs3-q0.25 diffbind-DESeq2-p0.05 differential peak regions (first three columns) as Rxrg_new.q0.25_p0.05.bed
# b. generate peak midpoint.bed
# c. deeptools::computeMatrix reference-point for selected samples
# d. deoptools::plotHeatmap

cat ../Rxrg_new.q0.25_p0.05.bed |\
awk 'OFS="\t"{printf "%s\t%d\t%d\n" ,$1,$2+($3-$2)/2,$2+1+($3-$2)/2}' \
> Rxrg_new.q0.25_p0.05.midpoint.bed

# only 2nd_replicate used
computeMatrix reference-point -S \
${local_1104}/ILC2_cuttag.Rxrg_WT_1_30k.1104.nodup.rpkm.bw \
${local_1104}/ILC2_cuttag.Rxrg_KO_1_30k.1104.nodup.rpkm.bw \
${local_bw}/ILC2_ATAC.WT_2_10k.1124.nodup.rpkm.bw \
${local_bw}/ILC2_ATAC.KO_2_10k.1124.nodup.rpkm.bw \
${local_1104}/ILC2_cuttag.H3K27ac_WT_10k.1104.nodup.rpkm.bw \
${local_1104}/ILC2_cuttag.H3K27ac_KO_10k.1104.nodup.rpkm.bw \
${local_1104}/ILC2_cuttag.H3K27me3_WT_10k.1104.nodup.rpkm.bw \
${local_1104}/ILC2_cuttag.H3K27me3_KO_10k.1104.nodup.rpkm.bw \
-R ../Heatmap/Rxrg_new.q0.25_p0.05.midpoint.bed \
--skipZeros \
-o ./Rxrg_2nd_rep.q0.25_p0.05.mat.gz \
-p 6 \
-a 3000 -b 3000 \
--referencePoint center \
--samplesLabel \
Rxrg_WT_30k \
Rxrg_KO_30k \
ATAC_WT_10k \
ATAC_KO_10k \
ac_WT_10k \
ac_KO_10k \
me3_WT_10k \
me3_KO_10k

#
plotHeatmap \
-m ./Rxrg_2nd_rep.q0.25_p0.05.mat.gz \
-out ./final_cuttag_Rxrg.diff_peak.plotHeatmap.2nd_rep.noline.pdf \
--sortUsing sum \
--startLabel "Peak Start" \
--endLabel "Peak End" \
--xAxisLabel "" \
--regionsLabel "Peaks" \
--colorList 'white, #EE2124' \
--heatmapWidth 3 \
--heatmapHeight 12 \
--dpi 200 \
--whatToShow 'heatmap and colorbar' 

plotProfile \
-m ./Rxrg_2nd_rep.q0.25_p0.05.mat.gz \
-out ./final_cuttag_Rxrg.diff_peak.plotHeatmap.2nd_rep.noheat.pdf \
--colors red red red red red red red red \
--regionsLabel "" \
--legendLocation "upper-right" \
--plotWidth 4 \
--plotHeight 5.4

## whole genome correlation
# Rxrg
multiBigwigSummary bins \
-p 4 \
-bs 5000 \
-b \
${local_1014}/ILC2_cuttag.Rxrg_WT_30k.1014.nodup.rpkm.bw \
${local_1104}/ILC2_cuttag.Rxrg_WT_1_30k.1104.nodup.rpkm.bw \
${local_1014}/ILC2_cuttag.Rxrg_KO_30k.1014.nodup.rpkm.bw \
${local_1104}/ILC2_cuttag.Rxrg_KO_1_30k.1104.nodup.rpkm.bw \
-l \
Rxrg_WT0_30k_1014 \
Rxrg_WT1_30k_1104 \
Rxrg_KO0_30k_1014 \
Rxrg_KO1_30k_1104 \
-o ./nodup.bin5kb.rpkm.final_cuttag_Rxrg.npz 

plotCorrelation \
-in nodup.bin5kb.rpkm.final_cuttag_Rxrg.npz \
-o nodup.bin5kb.rpkm.final_cuttag_Rxrg.pdf \
--corMethod pearson \
-p scatterplot \
--removeOutliers \
--outFileCorMatrix pearson.final_cuttag_Rxrg.txt \
--plotHeight 2.15 \
--plotWidth 2.25

# ATAC
multiBigwigSummary bins \
-p 4 \
-bs 5000 \
-b \
${local_bw}/ILC2_ATAC.WT_1_10k.1124.nodup.rpkm.bw \
${local_bw}/ILC2_ATAC.WT_2_10k.1124.nodup.rpkm.bw \
${local_bw}/ILC2_ATAC.KO_1_10k.1124.nodup.rpkm.bw \
${local_bw}/ILC2_ATAC.KO_2_10k.1124.nodup.rpkm.bw \
-l \
ATAC_WT1_10k \
ATAC_WT2_10k \
ATAC_KO1_10k \
ATAC_KO2_10k \
-o ./nodup.bin5kb.rpkm.final_ATAC_Rxrg.npz  

plotCorrelation \
-in nodup.bin5kb.rpkm.final_ATAC_Rxrg.npz \
-o nodup.bin5kb.rpkm.final_ATAC_Rxrg.pdf \
--corMethod pearson \
-p scatterplot \
--removeOutliers \
--outFileCorMatrix pearson.final_ATAC_Rxrg.txt \
--plotHeight 2.15 \
--plotWidth 2.25

# ac
multiBigwigSummary bins \
-p 4 \
-bs 5000 \
-b \
${local_1014}/ILC2_cuttag.H3K27ac_WT_10k.1014.nodup.rpkm.bw \
${local_1104}/ILC2_cuttag.H3K27ac_WT_10k.1104.nodup.rpkm.bw \
${local_1014}/ILC2_cuttag.H3K27ac_KO_10k.1014.nodup.rpkm.bw \
${local_1104}/ILC2_cuttag.H3K27ac_KO_10k.1104.nodup.rpkm.bw \
-l \
ac_WT_10k_1014 \
ac_WT_10k_1104 \
ac_KO_10k_1014 \
ac_KO_10k_1104 \
-o ./nodup.bin5kb.rpkm.final_cuttag_ac.npz  

plotCorrelation \
-in nodup.bin5kb.rpkm.final_cuttag_ac.npz \
-o nodup.bin5kb.rpkm.final_cuttag_ac.pdf \
--corMethod pearson \
-p scatterplot \
--removeOutliers \
--outFileCorMatrix pearson.final_cuttag_ac.txt \
--plotHeight 2.15 \
--plotWidth 2.25

# me3
multiBigwigSummary bins \
-p 4 \
-bs 5000 \
-b \
${local_1104}/ILC2_cuttag.H3K27me3_WT_10k.1104.nodup.rpkm.bw \
${local_1104}/ILC2_cuttag.H3K27me3_KO_10k.1104.nodup.rpkm.bw \
-l \
me3_WT_10k_1104 \
me3_KO_10k_1104 \
-o ./nodup.bin5kb.rpkm.final_cuttag_me3.npz 

plotCorrelation \
-in nodup.bin5kb.rpkm.final_cuttag_me3.npz \
-o nodup.bin5kb.rpkm.final_cuttag_me3.pdf \
--corMethod pearson \
-p scatterplot \
--removeOutliers \
--outFileCorMatrix pearson.final_cuttag_me3.txt \
--plotHeight 2.15 \
--plotWidth 2.25

## motif enrichment
# homer
findMotifsGenome.pl \
../../../Rxrg_new.q0.25_p0.05.bed \
mm10 \
./motif_Rxrg_p0.05 \
-size 200 \
-mask \
-S 100 \
-p 6

# meme
fa_file='/Shared/genomics/mouse/GRCm38_vM25/GRCm38.p6.genome.fa'
bedtools getfasta -fi ${fa_file} -bed \
../../../Rxrg_new.q0.25_p0.05.bed \
-fo Rxrg_new.q0.25_p0.05.fa

# submit to https://meme-suite.org/meme/doc/meme-chip.html
#   Maximum width: 30
#   database select:  
#       Human and Mouse(HOCOMOCO v11 FULL)
#       HOCOMOCO Mouse(v11 CORE)
#       JASPAR CORE and UniPROBE Mouse
# summarize the output results
#   could see 'AGGTCA' sequence always in top motif enrichment

## IGV check signal regions
# vector-image output using FIlE-Save SVG Image

## GREAT analysis
# peak bed file was uploaded for online GREAT analysis(http://great.stanford.edu) 
#     with default parameters
#         (Basal plus extension for finding peak-associated genes: Proximal 5kb upstream, 1kb downstream, plus Distal up to 1000kb). 
# Then 259 putative Rxrg-bound peaks got total 387 associated genes, 
# peak-gene/gene-peak tables and GO enrichment of those associated genes were downloaded from the GREAT website.

## GSEA of +1000 random genesets
#    output rownames of filtered expression matrix as: gene_expressed.txt
#    output GREAT-associated 387 genes as: GREAT387.txt
#  also try only adding the geneset into hallmark or go.bp

# remap to human symbol
remap="/Shared/projects/RNA_normal/GSEA_db/Mouse_Gene_Symbol_Remapping_Human_Orthologs_MSigDB.v7.3.chip"

:>GREAT387.remap.txt
for GG in $(cat GREAT387.txt)
do
cat ${remap} |gawk '$1=="'${GG}'"' |cut -f 1,2 >> GREAT387.remap.txt
done

# add into hallmark as new .gmt
cp h.all.v7.3.symbols.gmt h.all.v7.3.symbols.addGREAT387.gmt
echo -e \
"RxrgCUTTagPeaks_GREAT387_Genes local_added " \
$(cat GREAT387.remap.txt |cut -f 2 |tr "\n" " ") |\
tr " " "\t" \
>> h.all.v7.3.symbols.addGREAT387.gmt

# add into go.bp
cp c5.go.bp.v7.3.symbols.gmt c5.go.bp.v7.3.symbols.addGREAT387.gmt
echo -e \
"RxrgCUTTagPeaks_GREAT387_Genes local_added " \
$(cat GREAT387.remap.txt |cut -f 2 |tr "\n" " ") |\
tr " " "\t" \
>> c5.go.bp.v7.3.symbols.addGREAT387.gmt

# add into custom random genesets
:>gene_expressed.remap.txt
for GG in $(cat gene_expressed.txt)
do
cat ${remap} |gawk '$1=="'${GG}'"' |cut -f 1,2 >> gene_expressed.remap.txt
done

:> GREAT387.random1000genesets.gmt
echo -e \
"RxrgCUTTagPeaks_GREAT387_Genes local_added " \
$(cat ../GREAT387.remap.txt |cut -f 2 |tr "\n" " ") |\
tr " " "\t" \
>> GREAT387.random1000genesets.gmt

# only 270 of 387 expressed, so just select 270 as one geneset
for NN in {1..1000}
do
echo -e \
"Remap270_random"${NN} " local_added " \
$(cat gene_expressed.remap.txt |cut -f 2 |shuf |\
head -n 270 |tr "\n" " ") |tr " " "\t" \
>> GREAT387.random1000genesets.gmt
done

# use this .gmt to run GSEA enrichment




