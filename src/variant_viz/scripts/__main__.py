"""
Commands for variant_viz
"""
import click
from variant_viz import __version__
from variant_viz.scripts.utility import EC19_data
from variant_viz.scripts.utility import GISAID_stats

@click.group(help=f"""EDGE Covid19 variant viz v{__version__}.""")
def vizcli():
    pass

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
    report = EC19_data(snps, gaps, alnstats, pango, metadata, output)
    report.generate_EC_repot()

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
    gisaid = GISAID_stats(gisaid_pkl=meta_pkl, gisaid_mutation_pkl=mut_pkl, country=country, state=state)
    gisaid.generate_stats_html(geo_type, output)


@vizcli.command('tsv2pkl')
@click.option('-t', '--tsv',
              help='input GISAID metadata.tsv',
              required=True,
              type=click.File(mode='r'))
@click.option('-x', '--prefix',
              help='output filename for python .pkl file',
              required=True,
              type=str)

def tsv2pkl(tsv, prefix):
    """
    Convert GISAID metadata.tsv file to .pkl files
    """
    gisaid = GISAID_stats(gisaid_tsv=tsv)
    gisaid.tsv2pkl(prefix)


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
              type=click.File(mode='r'))
@click.option('-s', '--snps',
              help='EC19 SNPs tsv file (NC_045512.2_consensus.SNPs_report.txt)',
              required=True,
              type=click.File(mode='r'))
@click.option('-p', '--pango',
              help='EC19 Pangolin lineage tsv file (NC_045512.2_consensus_lineage.txt)',
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

def project(meta_pkl, mut_pkl, sample, snps, pango, output, geo_type, country, state):
    """
    Generate a stats html for a particular EC19 project
    """
    gisaid = GISAID_stats(gisaid_pkl=meta_pkl, gisaid_mutation_pkl=mut_pkl, country=country, state=state)
    gisaid.generate_ec19_sample_html(meta_pkl, mut_pkl, sample, snps, pango, output, geo_type, country, state)


if __name__ == '__main__':
    vizcli()
