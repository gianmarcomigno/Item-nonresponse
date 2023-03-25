########################################################
#            CELL TYPE SPECIFIC ANALYSES               #
########################################################
# Reference: https://github.com/bulik/ldsc/wiki/Cell-type-specific-analyses

first by functional annotations, second only by brain
————————————————————————————————————————————————
————————
gwas_factscores_PNA_g.tsv.gz

#Choose which set of gene sets to analyze. Options include Multi_tissue_gene_expr, Multi_tissue_chromatin, GTEx_brain, Cahoy, ImmGen, or Corces_ATAC

cts_name=GTEx_brain 
cts_name=Multi_tissue_gene_expr

#Download the LD scores
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/LDSC_SEG_ldscores/${cts_name}_1000Gv3_ldscores.tgz
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_baseline_ldscores.tgz
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/weights_hm3_no_hla.tgz
tar -xvzf ${cts_name}_1000Gv3_ldscores.tgz
tar -xvzf 1000G_Phase3_baseline_ldscores.tgz
tar -xvzf 1000G_Phase3_weights_hm3_no_MHC.tar
tar -xvzf 1000G_Phase3_cell_type_groups.tgz

### working!
###### RUN LDSCORE-REGRESSION (H2) - STRATIFIED BY FUNCTIONAL ANNOTATIONS - GTEX BRAIN ######
1)
 python ldsc.py \
 --h2-cts /Users/gianmarcomignogna/ldsc/gwas_IDK.tsv.gz \
 --ref-ld-chr /Users/gianmarcomignogna/ldsc/1000G_EUR_Phase3_baseline/baseline. \
 --ref-ld-chr-cts /Users/gianmarcomignogna/ldsc/GTEx_brain.ldcts \
 --w-ld-chr /Users/gianmarcomignogna/ldsc/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
 --print-coefficients \
 --out IDK_new1

2)
 python ldsc.py \
 --h2-cts /Users/gianmarcomignogna/ldsc/gwas_IDK.tsv.gz \
 --ref-ld-chr /Users/gianmarcomignogna/ldsc/1000G_EUR_Phase3_baseline/baseline. \
 --ref-ld-chr-cts /Users/gianmarcomignogna/ldsc/Multi_tissue_gene_expr.ldcts \
 --w-ld-chr /Users/gianmarcomignogna/ldsc/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
 --print-coefficients \
 --out IDK_new2


3)
### cell type groups analysis (no .ldcts file here)
for i in `seq 1 10`; do 
python ldsc.py --h2 gwas_IDK.tsv.gz \
      --w-ld-chr /Users/gianmarcomignogna/ldsc/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
      --ref-ld-chr /Users/gianmarcomignogna/ldsc/1000G_EUR_Phase3_baseline/baseline.,/Users/gianmarcomignogna/ldsc/1000G_Phase3_cell_type_groups/cell_type_group.${i}. \
      --overlap-annot \
      --frqfile-chr /Users/gianmarcomignogna/ldsc/1000G_Phase3_frq/1000G.EUR.QC. \
      --print-coefficients \
      --out /Users/gianmarcomignogna/ldsc/results_IDK_new3/IDK_new3_tissue${i}_baseline; 
done

cat <(for i in `seq 1 10`; do awk -v ii=$i 'NR==FNR{if($1==ii){cell=$2}}NR!=FNR && FNR==1 && ii==1{$1=$1; print $0}NR!=FNR && $1=="L2_1"{$1=cell; print $0}' /Users/gianmarcomignogna/ldsc/1000G_Phase3_cell_type_groups/names /Users/gianmarcomignogna/ldsc/results_IDK_new3/IDK_new3_tissue${i}_baseline.results; done) > IDK_new3_final_output_tissues.results

 

