"""
Commands for variant_viz
"""
import click
from variant_viz import __version__
from variant_viz.scripts.utility import EC19_Report
from variant_viz.scripts.utility import GISAID_stats

@click.group(help=f"""EDGE Covid19 variant viz v{__version__}.""")
def vizcli():
    pass

@vizcli.command('report')
@click.option('--snps',
              help='EC19 SNPs tsv file',
              required=True,
              type=click.File(mode='r'))
@click.option('--gaps',
              help='EC19 gaps tsv file',
              required=True,
              type=click.File(mode='r'))
@click.option('--alnstats',
              help='EC19 gaps alnstats file',
              required=True,
              type=click.File(mode='r'))
@click.option('--pango',
              help='Pangolin lineage tsv file',
              required=True,
              type=click.File(mode='r'))
@click.option('--metadata',
              help='Metadata',
              required=True,
              type=click.File(mode='r'))
@click.option('--output',
              help='output filename for EC19 report html',
              required=True,
              type=str)

def report(snps, gaps, alnstats, pango, metadata, output):
    """
    Generate a visualization report for EDGE-Covid19 workflow
    """
    report = EC19_Report(snps, gaps, alnstats, pango, metadata, output)
    report.generate_EC_repot()

@vizcli.command('gisaid_stats')
@click.option('--meta-pkl',
              help='GISAID metadata pre-parsed .pkl file',
              required=True,
              type=click.File(mode='r'))
@click.option('--mut-pkl',
              help='GISAID mutation pre-parsed .pkl file',
              required=True,
              type=click.File(mode='r'))
@click.option('--geo-type',
              help='Scale of geo visualization',
              default='global',
              type=click.Choice(['global', 'country', 'state'], 
              case_sensitive=False))
@click.option('--country',
              help='specify a country',
              required=False,
              type=str)
@click.option('--state',
              help='specify a state (US only)',
              required=False,
              type=str)
@click.option('--output',
              help='output filename for stats html',
              required=True,
              type=str)

def gisaid_stats(meta_pkl, mut_pkl, geo_type, country, state, output):
    """
    Generate a stats html file from GISAID data
    """
    report = GISAID_stats(gisaid_pkl=meta_pkl, gisaid_mutation_pkl=mut_pkl, country=country, state=state)
    report.generate_stats_html(geo_type, output)


@vizcli.command('tsv2pkl')
@click.option('--tsv',
              help='input GISAID metadata.tsv',
              required=True,
              type=click.File(mode='r'))
@click.option('--prefix',
              help='output filename for python .pkl file',
              required=True,
              type=str)

def tsv2pkl(tsv, prefix):
    """
    Convert GISAID metadata.tsv file to .pkl files
    """
    stats = GISAID_stats(gisaid_tsv=tsv)
    stats.tsv2pkl(prefix)


if __name__ == '__main__':
    vizcli()
