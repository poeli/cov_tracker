ec19_varviz tsv2pkl \
    --tsv    test/metadata.tsv \
    --prefix test/metadata

ec19_varviz gisaid_stats \
    --meta-pkl test/metadata.pkl \
    --mut-pkl  test/metadata.mutation.pkl \
    --output   test/metadata_global.html

ec19_varviz gisaid_stats \
    --meta-pkl test/metadata.pkl \
    --mut-pkl  test/metadata.mutation.pkl \
    --geo-type country \
    --country  USA \
    --output   test/metadata_USA.html

ec19_varviz gisaid_stats \
    --meta-pkl test/metadata.pkl \
    --mut-pkl  test/metadata.mutation.pkl \
    --geo-type state \
    --country  USA \
    --state    California \
    --output   test/metadata_CA.html

ec19_varviz gisaid_stats \
    --meta-pkl test/metadata.pkl \
    --mut-pkl  test/metadata.mutation.pkl \
    --geo-type state \
    --country  USA \
    --state    "New Mexico" \
    --output   test/metadata_NM.html

ec19_varviz project \
    --meta-pkl test/metadata.pkl \
    --mut-pkl  test/metadata.mutation.pkl \
    --sample   2251_127 \
    --snps     test/NC_045512.2_consensus.SNPs_report.txt \
    --pango    test/NC_045512.2_consensus_lineage.txt \
    --metadata test/metadata_gisaid_ncbi.txt \
    --geo-type country \
    --country  USA \
    --output   test/ec19_project_2251_127.html

ec19_varviz report \
    --snps     test/lanl_project_list.SNP.tsv \
    --gaps     test/lanl_project_list.gaps.tsv \
    --alnstats test/lanl_project_list.alnstats.tsv \
    --pango    test/lanl_project_list.lineage_report.csv \
    --metadata test/lanl_project_list.metadata.tsv \
    --output   test/lanl_project_list_ec19.html
