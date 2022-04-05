"""
Commands for variant_viz
"""
import click
import sys
import logging
import datetime
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
              help='complete genomes only (>29Kbp) - is_complete==True (default: True)',
              default=True,
              type=bool)
@click.option('-hc', '--high_coverage_only',
              help='high_coverage_only - is_high_coverage==True (default: False)',
              default=False,
              type=bool)
@click.option('-n', '--n-content',
              help='remove genomes if percentage of Ns excess this parameter (default: 0.01)',
              default=0.01,
              type=float)
@click.option('-ds', '--date-start',
              help='select data from this date in YYYY-MM-DD format',
              type=str,
              required=False,
              default=None)
@click.option('-de', '--date-end',
              help='select data ending this date in YYYY-MM-DD format',
              type=str,
              required=False,
              default=None)
@click.option('--debug',
              help='logging at debug level',
              is_flag=True)


def tsv2pkl(tsv, prefix, country, state, complete, high_coverage_only, n_content, date_start, date_end, debug):
    """
    Convert GISAID metadata.tsv file to .pkl files
    """
    print(f"Running CovTracker v{__version__} tsv2pkl command...")
    if debug:
        logging.basicConfig(level=logging.DEBUG)
        logging.debug('Logging level set to DEBUG.')

    gisaid = GISAID_stats(gisaid_tsv=tsv, 
                          country=country, 
                          state=state, 
                          complete_only=complete, 
                          high_coverage_only=high_coverage_only, 
                          n_content=n_content, 
                          date_start=date_start, 
                          date_end=date_end, 
                          merge_meta_to_mut=False)
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
@click.option('-ds', '--date-start',
              help='select data from this date in YYYY-MM-DD format',
              type=str,
              required=False,
              default=None)
@click.option('-de', '--date-end',
              help='select data ending this date in YYYY-MM-DD format',
              type=str,
              required=False,
              default=None)
@click.option('--debug',
              help='logging at debug level',
              is_flag=True)


def gisaid_stats(meta_pkl, mut_pkl, geo_type, country, state, output, date_start, date_end, debug):
    """
    Generate a stats html file from GISAID data
    """
    print(f"Running CovTracker v{__version__} gisaid_stats command...")
    if debug:
        logging.basicConfig(level=logging.DEBUG)
        logging.debug('Logging level set to DEBUG.')
    if geo_type=='state' and state==None:
        print(f"ERROR: --geo-type is set to 'state' but --state is not specified.")
        sys.exit(1)
    gisaid = GISAID_stats(gisaid_pkl=meta_pkl, gisaid_mutation_pkl=mut_pkl, country=country, state=state, date_start=date_start, date_end=date_end)
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
@click.option('-ds', '--date-start',
              help='select data from this date in YYYY-MM-DD format',
              type=str,
              required=False,
              default=None)
@click.option('-de', '--date-end',
              help='select data ending this date in YYYY-MM-DD format',
              type=str,
              required=False,
              default=None)
@click.option('--debug',
              help='logging at debug level',
              is_flag=True)


def project(meta_pkl, mut_pkl, sample, snps, pango, metadata, output, geo_type, country, state, date_start, date_end, debug):
    """
    Generate a stats html for a particular EC19 project
    """
    if debug:
        logging.getLogger().setLevel(logging.DEBUG)
        logging.debug('Logging level set to DEBUG.')
    if geo_type=='state' and state==None:
        print(f"ERROR: --geo-type is set to 'state' but --state is not specified.")
        sys.exit(1)

    print(f"Running CovTracker v{__version__} project command...")
    virus_name=None
    collection_date=None
    location=None

    for line in metadata:
        line.strip()
        f = line.split('=')
        if line.startswith('virus_name'):
            virus_name = f[1].strip()
        if line.startswith('collection_date'):
            collection_date = f[1].strip()
        if line.startswith('location'):
            location = f[1].strip()

    gisaid = GISAID_stats(gisaid_pkl=meta_pkl, gisaid_mutation_pkl=mut_pkl, country=country, state=state, date_start=date_start, date_end=date_end)
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
@click.option('--debug',
              help='logging at debug level',
              is_flag=True)


def report(snps, gaps, alnstats, pango, metadata, output, debug):
    """
    Generate a visualization report for EDGE-Covid19 workflow
    """
    print(f"Running CovTracker v{__version__} report command...")
    if debug:
        logging.basicConfig(level=logging.DEBUG)
        logging.debug('Logging level set to DEBUG.')

    report = EC19_data(snps, gaps, alnstats, pango, metadata, output)
    report.generate_EC_repot()


@vizcli.command('lanl_summary')
@click.option('-mk', '--meta-pkl',
              help='GISAID metadata pre-parsed .pkl file',
              required=True,
              type=click.File(mode='r'))
@click.option('-o', '--output',
              help='output filename for stats html',
              required=True,
              type=str)
@click.option('--debug',
              help='logging at debug level',
              is_flag=True)


def lanl_summary(meta_pkl, output, debug):
    """
    Generate a lanl_summary html file from GISAID data
    """
    print(f"Running CovTracker v{__version__} lanl_summary command...")
    if debug:
        logging.basicConfig(level=logging.DEBUG)
        logging.debug('Logging level set to DEBUG.')
    gisaid = GISAID_stats(gisaid_pkl=meta_pkl, merge_meta_to_mut=False)
    gisaid.generate_lanl_summary_html(output)


@vizcli.command('lanl_ec19_summary')
# @click.option('-s', '--snps',
#               help='EC19 SNPs tsv file',
#               required=True,
#               type=click.File(mode='r'))
# @click.option('-g', '--gaps',
#               help='EC19 gaps tsv file',
#               required=True,
#               type=click.File(mode='r'))
# @click.option('-a', '--alnstats',
#               help='EC19 gaps alnstats file',
#               required=True,
#               type=click.File(mode='r'))
@click.option('-p', '--pango',
              help='Pangolin lineage tsv file',
              required=True,
              type=click.File(mode='r'))
@click.option('-m', '--metadata',
              help='Metadata',
              required=True,
              type=click.File(mode='r'))
@click.option('-u', '--url',
              help='EDGE-COVID19 report url',
              required=True,
              type=str)
@click.option('-o', '--output',
              help='output filename for EC19 report html',
              required=True,
              type=str)
@click.option('--debug',
              help='logging at debug level',
              is_flag=True)


def lanl_ec19_summary(#snps, gaps, alnstats, 
                      pango, metadata, url, output, debug):
    """
    Generate a visualization report for EDGE-Covid19 workflow
    """
    print(f"Running CovTracker v{__version__} report command...")
    if debug:
        logging.basicConfig(level=logging.DEBUG)
        logging.debug('Logging level set to DEBUG.')

    report = EC19_data(#snps, gaps, alnstats, 
                       pango=pango, metadata=metadata, output=output)
    report.generate_lanl_ec19_summary(url)




if __name__ == '__main__':
    vizcli()
