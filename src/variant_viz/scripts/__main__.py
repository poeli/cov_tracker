"""
Commands for variant_viz
"""
import click
import sys
import logging
from variant_viz import __version__
from variant_viz.scripts.utility import EC19_data
from variant_viz.scripts.utility import GISAID_stats

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(module)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M',
)

@click.group(help=f"""EDGE Covid19 variant viz v{__version__}.""")
def vizcli():
    pass

@vizcli.command('tsv2pkl')
@click.option('-t', '--tsv',
              help='input GISAID metadata.tsv',
              required=True,
              type=click.File(mode='r'))
@click.option('-x', '--prefix',
              help='output filename for python .pkl file',
              required=True,
              type=str)
@click.option('-c', '--country',
              help='specify a country if available',
              required=False,
              type=str)
@click.option('-d', '--state',
              help='specify a state if available',
              required=False,
              type=str)
@click.option('-cp', '--complete',
              help='complete genomes only (>29Kbp))',
              default=True,
              type=bool)
@click.option('-hc', '--high-coverage',
              help='high coverage genomes only (<1% Ns)',
              default=True,
              type=bool)


def tsv2pkl(tsv, prefix, country, state, complete, high_coverage):
    """
    Convert GISAID metadata.tsv file to .pkl files
    """
    print(f"Running tsv2pkl command...")
    gisaid = GISAID_stats(gisaid_tsv=tsv, country=country, state=state, complete_only=complete, high_coverage_only=high_coverage)
    gisaid.tsv2pkl(prefix)


@vizcli.command('gisaid_stats')
@click.option('-mk', '--meta-pkl',
              help='GISAID metadata pre-parsed .pkl file',
              required=True,
              type=click.File(mode='r'))
@click.option('-tk', '--mut-pkl',
              help='GISAID mutation pre-parsed .pkl file',
              required=True,
              type=click.File(mode='r'))
@click.option('-gt', '--geo-type',
              help='Scale of geo visualization',
              default='global',
              type=click.Choice(['global', 'country', 'state'], case_sensitive=False))
@click.option('-c', '--country',
              help='specify a country',
              required=False,
              type=str)
@click.option('-d', '--state',
              help='specify a state (US only)',
              required=False,
              type=str)
@click.option('-o', '--output',
              help='output filename for stats html',
              required=True,
              type=str)

def gisaid_stats(meta_pkl, mut_pkl, geo_type, country, state, output):
    """
    Generate a stats html file from GISAID data
    """
    print(f"Running gisaid_stats command...")
    if geo_type=='state' and state==None:
        print(f"ERROR: --geo-type is set to 'state' but --state is not specified.")
        sys.exit(1)
    gisaid = GISAID_stats(gisaid_pkl=meta_pkl, gisaid_mutation_pkl=mut_pkl, country=country, state=state)
    gisaid.generate_stats_html(geo_type, output)


@vizcli.command('project')
@click.option('-mk', '--meta-pkl',
              help='GISAID metadata pre-parsed .pkl file',
              required=True,
              type=click.File(mode='r'))
@click.option('-tk', '--mut-pkl',
              help='GISAID mutation pre-parsed .pkl file',
              required=True,
              type=click.File(mode='r'))
@click.option('-n', '--sample',
              help='EC19 sample name of the project',
              required=True,
              type=str)
@click.option('-s', '--snps',
              help='EC19 SNPs tsv file (NC_045512.2_consensus.SNPs_report.txt)',
              required=True,
              type=click.File(mode='r'))
@click.option('-p', '--pango',
              help='EC19 Pangolin lineage tsv file (NC_045512.2_consensus_lineage.txt)',
              required=True,
              type=click.File(mode='r'))
@click.option('-m', '--metadata',
              help='EC19 metadata file (metadata_gisaid_ncbi.txt)',
              required=True,
              type=click.File(mode='r'))
@click.option('-o', '--output',
              help='output filename for stats html',
              required=True,
              type=str)
@click.option('-gt', '--geo-type',
              help='Scale of geo visualization',
              default='global',
              type=click.Choice(['global', 'country', 'state'], case_sensitive=False))
@click.option('-c', '--country',
              help='specify a country',
              required=False,
              type=str)
@click.option('-d', '--state',
              help='specify an US state in fullname (like: "New Mexico")',
              required=False,
              type=str)

def project(meta_pkl, mut_pkl, sample, snps, pango, metadata, output, geo_type, country, state):
    """
    Generate a stats html for a particular EC19 project
    """
    if geo_type=='state' and state==None:
        print(f"ERROR: --geo-type is set to 'state' but --state is not specified.")
        sys.exit(1)

    print(f"Running project command...")
    virus_name=None
    collection_date=None
    location=None

    for line in metadata:
        line.strip()
        f = line.split('=')
        if line.startswith('virus_name'):
            virus_name = f[1]
        if line.startswith('collection_date'):
            collection_date = f[1]
        if line.startswith('location'):
            location = f[1]

    gisaid = GISAID_stats(gisaid_pkl=meta_pkl, gisaid_mutation_pkl=mut_pkl, country=country, state=state)
    gisaid.generate_ec19_sample_html(sample, snps, pango, output, geo_type, virus_name, collection_date, location)

@vizcli.command('report')
@click.option('-s', '--snps',
              help='EC19 SNPs tsv file',
              required=True,
              type=click.File(mode='r'))
@click.option('-g', '--gaps',
              help='EC19 gaps tsv file',
              required=True,
              type=click.File(mode='r'))
@click.option('-a', '--alnstats',
              help='EC19 gaps alnstats file',
              required=True,
              type=click.File(mode='r'))
@click.option('-p', '--pango',
              help='Pangolin lineage tsv file',
              required=True,
              type=click.File(mode='r'))
@click.option('-m', '--metadata',
              help='Metadata',
              required=True,
              type=click.File(mode='r'))
@click.option('-o', '--output',
              help='output filename for EC19 report html',
              required=True,
              type=str)

def report(snps, gaps, alnstats, pango, metadata, output):
    """
    Generate a visualization report for EDGE-Covid19 workflow
    """
    print(f"Running report command...")
    report = EC19_data(snps, gaps, alnstats, pango, metadata, output)
    report.generate_EC_repot()

if __name__ == '__main__':
    vizcli()
