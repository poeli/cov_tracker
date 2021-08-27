# SARS-CoV-2 variant tracking

## Installation
```
pip install -r requirements.txt
python setup.py install
```

## Test
```
# convert GISAID metadata.tsv file to .pkl files (metadata.pkl and metadata.mutation.pkl)
ec19_varviz tsv2pkl \
    --tsv    /path/to/gisaid/metadata.tsv \
    --prefix test/metadata

# generate global stats
ec19_varviz gisaid_stats \
    --meta-pkl test/metadata.pkl \
    --mut-pkl  test/metadata.mutation.pkl \
    --output   test/metadata_global.html

# generate US stats and display lineage distributions for each state
ec19_varviz gisaid_stats \
    --meta-pkl test/metadata.pkl \
    --mut-pkl  test/metadata.mutation.pkl \
    --geo-type country \
    --country  USA \
    --output   test/metadata_USA.html

# generate California stats and display lineage distributions for each county
ec19_varviz gisaid_stats \
    --meta-pkl test/metadata.pkl \
    --mut-pkl  test/metadata.mutation.pkl \
    --geo-type state \
    --country  USA \
    --state    California \
    --output   test/metadata_CA.html

# generate lineage tracking info for an EC19 project (USA only)
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

# generate report for EC19 projects
ec19_varviz report \
    --snps     test/lanl_project_list.SNP.tsv \
    --gaps     test/lanl_project_list.gaps.tsv \
    --alnstats test/lanl_project_list.alnstats.tsv \
    --pango    test/lanl_project_list.lineage_report.csv \
    --metadata test/lanl_project_list.metadata.tsv \
    --output   test/lanl_project_list_ec19.html

```
