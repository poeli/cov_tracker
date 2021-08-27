import pandas as pd
import numpy as np
import os
import logging
import variant_viz.surv_viz.us_state_abbrev as us

class CovidMetadata(object):
    """
    This is the class for loading GISAID metadata
    
    :param filename_meta_tsv:
    :type  filename_meta_tsv: str
    :param filename_cnty_tsv:
    :type  filename_cnty_tsv: str
    :param filename_v_info_tsv:
    :type  filename_v_info_tsv: str
    :param filename_v_muta_tsv:
    :type  filename_v_muta_tsv: str
    :param filename_anno_tsv:
    :type  filename_anno_tsv: str
    :param filename_meta_pkl:
    :type  filename_meta_pkl: str
    :param filename_mutation_pkl:
    :type  filename_mutation_pkl: str
    :param save_pkl:
    :type  save_pkl: boolean
    """
    
    def __init__(self, filename_cnty_tsv=None,
                       filename_v_info_tsv=None,
                       filename_v_muta_tsv=None,
                       filename_anno_tsv=None,
                       filename_meta_pkl=None,
                       filename_mutation_pkl=None,
                       filename_meta_tsv=None,
                       country=None,
                       division=None):
        
        self.df_meta_orig     = pd.DataFrame()
        self.df_meta          = pd.DataFrame()
        self.df_mutation_orig = pd.DataFrame()
        self.df_mutation      = pd.DataFrame()
        self.df_country       = pd.DataFrame()
        self.df_v_info        = pd.DataFrame()
        self.df_v_muta        = pd.DataFrame()
        self.df_v_anno        = pd.DataFrame()

        # Loading variant information and country coordinate files
        if filename_cnty_tsv:
            self.df_country  = pd.read_csv(filename_cnty_tsv, sep="\t")
        if filename_v_info_tsv:
            self.df_v_info   = pd.read_csv(filename_v_info_tsv, sep="\t")
        if filename_v_muta_tsv:
            self.df_v_muta   = pd.read_csv(filename_v_muta_tsv, sep="\t")
        if filename_anno_tsv:
            self.df_v_anno   = pd.read_csv(filename_anno_tsv, sep="\t")
        
        # Use pre-processed pickle files first if they are provided, otherwise --
        # read the raw metadata tsv file, process and save them to pickle files
        if filename_meta_pkl and filename_mutation_pkl:
            logging.info(f'Loading {filename_meta_pkl}...')
            self.df_meta_orig     = pd.read_pickle(filename_meta_pkl)
            logging.info(f'Loading {filename_mutation_pkl}...')
            self.df_mutation_orig = pd.read_pickle(filename_mutation_pkl)
        elif filename_meta_tsv:
            logging.info(f'Parsing {filename_meta_tsv} file...')
            (self.df_meta_orig, self.df_mutation_orig) = self._prepare_metadata(filename_meta_tsv)
        
        if country==None and division==None:
            self.df_meta     = self.df_meta_orig
            self.df_mutation = self.df_mutation_orig
        else:
            logging.info(f'Specifying country to {country}...')
            self.set_country(country)
            logging.info(f'Specifying division to {division}...')
            self.set_division(division)

        
    def _prepare_metadata(self, filename_meta):
        """
        Loading and prepare metadata tsv file
        
        :param filename_meta:
        :type  filename_meta: str
        """
        import pandas as pd
        import numpy as np

        # loading metadata.tsv file
        df_meta = pd.read_csv(filename_meta, sep="\t", low_memory=False)
        
        # accommodate new GISAID metadata
        df_meta = df_meta.rename(columns={
            'Virus name': 'name',
            'Type': 'type',
            'Accession ID': 'acc',
            'Collection date': 'date',
            'Additional location information': 'Location_add',
            'Sequence length': 'length',
            'Host': 'host',
            'Patient age': 'age',
            'Gender': 'gender',
            'Clade': 'clade',
            'Pango lineage': 'lineage',
            'Pangolin version': 'pangolin_ver',
            'Variant': 'var', #e.g. VUI202012/01 GRY (B.1.1.7) GH/501Y.v2 (B.1.351)
            'AA Substitutions': 'mutation',
            'Submission date': 'submission_date',
            'Is reference?': 'is_ref',
            'Is complete?': 'is_complete',
            'Is high coverage?': 'is_high_cov',
            'Is low coverage?': 'is_low_cov',
            'N-Content': 'n_content',
            'GC-Content': 'gc_content'
        })
        
        # Remove records with no full date
        df_meta = df_meta.drop(df_meta[df_meta.date.str.contains('-XX')].index)
        # df_meta = df_meta.drop(df_meta[df_meta.date.str.len()!=10].index) # remove records with 'year' only
        df_meta = df_meta[df_meta.host=='Human']
        # Remove mislabeled collection date
        df_meta = df_meta.drop(df_meta[df_meta.date<'2019-12'].index)
        
        # for GISAID
        # df_meta['date'] = df_meta['date'].astype('datetime64[ns]')
        df_meta['date'] = pd.to_datetime(df_meta['date'], errors = 'coerce')
        df_meta['week'] = df_meta['date'].dt.strftime('%Y-%Uw')
        df_meta[['region', 'country', 'division', 'location']] = df_meta['Location'].str.split(' / ', expand=True, n=3)
        df_meta['name'] = df_meta['name'].str.replace('hCoV-19/', '')
        df_meta = df_meta.drop(columns=['type', 'Location', 'Location_add'])

        # deleete "-00w" rows
        df_meta = df_meta.drop(df_meta[df_meta.week.str.contains('-00w')].index)

        # for Genbank
        #df_meta['acc'] = df_meta['genbank_accession']
        #df_meta = df_meta[~df_meta.genbank_accession.str.startswith('LC')]

        # cleaning location
        df_meta['location'] = df_meta['location'].str.replace(' County', '')

        # variant info extraction; VOC Gamma GR/501Y.V3 (P.1) first detected in B..
        df_meta[['v_type', 'lineage_greek']] = df_meta['var'].str.extract(r'^(VO\w) (\w+)', expand=True)

        # use WHO label to replace linage name if possible
        df_meta['pango_lineage'] = df_meta['lineage']
        idx = df_meta['lineage_greek'].notnull()
        df_meta.loc[idx, 'lineage'] = df_meta.loc[idx, 'lineage_greek']        
        
        # fine touching data
        df_meta.loc[df_meta.lineage=='B.1.427/429', 'lineage'] = 'Epsilon'
        
        # preparing aa_sub in gisaid
        df_mutation = df_meta[['acc', 'mutation']].copy()
        
        # cleaning parentheses, unwanted strings in aa_sub
        df_mutation['mutation'] = df_mutation['mutation'].str.replace('\(|\)', '')
        df_mutation['mutation'] = df_mutation['mutation'].str.replace(r'NS[\w\d]+,?', '') # NSP\d+ and NS\d+ are removed for now until we can map them to the protein
        df_mutation['mutation'] = df_mutation['mutation'].str.replace(r',[^,]+_ins[\w\d]+', '') # remove insertion
        df_mutation['mutation'] = df_mutation['mutation'].str.replace(r'Spike', 'S')
        df_mutation['mutation'] = df_mutation['mutation'].str.replace('del', '*')
        df_mutation['mutation'] = df_mutation['mutation'].str.replace('_', ':')

        # split variants, explode into rows, and extract genes/positions
        df_mutation = df_mutation.assign(variant=df_mutation['mutation'].str.split(',')).explode('variant')
        df_mutation = df_mutation.drop(columns=['mutation'])
        df_mutation = df_mutation.rename(columns={'variant': 'mutation'})
        df_mutation[['gene','pos']] = df_mutation['mutation'].str.extract(r'(\w+):\w(\d+)')

        # dropping reference genomes
        df_mutation = df_mutation.dropna()
        df_mutation['pos'] = df_mutation['pos'].astype(int)

        cols = ['acc', 'lineage', 'date', 'week', 'country', 'division']
        df_mutation = df_mutation.merge(df_meta[cols], on='acc', how='left')

        return df_meta, df_mutation

    def save_pkl(self, outfile_prefix):
        """
        Saving df_meta and df_mutation to .pkl files respectively

        :param outfile_prefix: output filename prefix
        :type  outfile_prefix: str
        """
        logging.info(f'Saving metadata to {outfile_prefix}.pkl...')
        self.df_meta.to_pickle(f'{outfile_prefix}.pkl')
        logging.info(f'Saving metadata to {outfile_prefix}.mutation.pkl...')
        self.df_mutation.to_pickle(f'{outfile_prefix}.mutation.pkl')
        
    def set_country(self, country):
        """
        Specific country
        
        :param country: country
        :type  country: str
        """
        if country:
            self.df_meta     = self.df_meta_orig[self.df_meta_orig['country']==country].copy()
            self.df_mutation = self.df_mutation_orig[self.df_mutation_orig['country']==country].copy()
    
    
    def set_division(self, division):
        """
        Specific division/state
        
        :param division: state
        :type  division: str
        """
        if division:
            self.df_meta     = self.df_meta_orig[self.df_meta_orig['division']==division].copy()
            self.df_mutation = self.df_mutation_orig[self.df_mutation_orig['division']==division].copy()