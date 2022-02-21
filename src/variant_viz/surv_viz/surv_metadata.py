import pandas as pd
import logging

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
                       division=None, 
                       complete_only=True,
                       n_content=0.01):
        
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
            self.df_v_muta[['gene','pos']] = self.df_v_muta['mutation'].str.extract(r'(\w+):\w(\d+)')
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
            (self.df_meta_orig, self.df_mutation_orig) = self._prepare_metadata(filename_meta_tsv, complete_only, n_content)
        
        if country==None and division==None:
            self.df_meta     = self.df_meta_orig
            self.df_mutation = self.df_mutation_orig
        else:
            logging.info(f'Specifying country to {country}...')
            self.set_country(country)
            logging.info(f'Specifying division to {division}...')
            self.set_division(division)

        ## add additional info to df_mutation
        # if len(self.df_mutation):
        #     logging.info(f'Adding metadata to mutations...')
        #     cols = ['acc', 'lineage', 'date', 'week', 'country', 'division']
        #     self.df_mutation = self.df_mutation.merge(self.df_meta[cols], on='acc', how='left')
        #     self.df_mutation[['gene','pos']] = self.df_mutation['mutation'].str.extract(r'(\w+):\w(\d+)')
        #     self.df_mutation['pos'] = self.df_mutation['pos'].astype(int)
        
    def _prepare_metadata(self, filename_meta, complete_only=True, n_content=0.01):
        """
        Loading and prepare metadata tsv file
        
        :param filename_meta: GISAID metadata.tsv
        :type  filename_meta: str
        :param complete_only: filter in the complete_only filed
        :type  complete_only: bool
        :param high_coverage_only: filter in the high_coverage_only filed
        :type  high_coverage_only: bool
        """
        import pandas as pd
        
        # loading metadata.tsv file
        df_meta = pd.DataFrame()
        with pd.read_csv(filename_meta, chunksize=1000000, sep="\t") as reader:
            for chunk in reader:
                df_meta = df_meta.append(chunk)
        
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
        # remove NOT high_coverage genomes
        # if high_coverage_only:
        #     df_meta = df_meta[df_meta.is_high_cov==True]
        # remove genomes with percentage of Ns exess n_content
        if n_content:
            df_meta = df_meta[df_meta.n_content<=n_content]
        # remove NOT complete genomes
        if complete_only:
            df_meta = df_meta[df_meta.is_complete==True]
        
        # for GISAID
        # df_meta['date'] = df_meta['date'].astype('datetime64[ns]')
        df_meta['date'] = pd.to_datetime(df_meta['date'], errors = 'coerce')
        df_meta['submission_date'] = pd.to_datetime(df_meta['submission_date'], errors = 'coerce')
        df_meta['week'] = df_meta['date'].dt.strftime('%Y-%Uw')
        df_meta[['region', 'country', 'division', 'location']] = df_meta['Location'].str.split(' / ', expand=True, n=3)
        df_meta['name'] = df_meta['name'].str.replace('hCoV-19/', '')
        df_meta = df_meta.drop(columns=['type', 'Location', 'Location_add', 'gc_content'])

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
        
        ## fine touching data
        #df_meta.loc[df_meta.lineage=='B.1.427/429', 'lineage'] = 'Epsilon'
        
        logging.info(f'Extracting mutations...')

        # preparing aa_sub in gisaid
        df_mutation = df_meta[['acc', 'mutation']].copy()
        
        # cleaning parentheses, unwanted strings in aa_sub
        df_mutation['mutation'] = df_mutation['mutation'].str.replace(r'\(|\)', '', regex=True)
        df_mutation['mutation'] = df_mutation['mutation'].str.replace(r'[^,]+_ins[^,]+', '', regex=True) # remove insertion
        df_mutation['mutation'] = df_mutation['mutation'].str.replace(r'Spike', 'S')
        df_mutation['mutation'] = df_mutation['mutation'].str.replace('del', '*')
        df_mutation['mutation'] = df_mutation['mutation'].str.replace('_', ':')

        # split variants, explode into rows, and extract genes/positions
        df_mutation = df_mutation.assign(variant=df_mutation['mutation'].str.split(',')).explode('variant')
        df_mutation = df_mutation.drop(columns=['mutation'])
        df_mutation = df_mutation.rename(columns={'variant': 'mutation'})

        # dropping reference genomes
        df_mutation = df_mutation[df_mutation.mutation!=""]

        logging.info(f'Converting genes and positions...')

        # converting GISAID mutations to EC19's format
        gene_prod = {
            'NSP1':  {'gene': 'ORF1a',  'offset':    0},
            'NSP2':  {'gene': 'ORF1a',  'offset':  180},
            'NSP3':  {'gene': 'ORF1a',  'offset':  818},
            'NSP4':  {'gene': 'ORF1a',  'offset': 2763},
            'NSP5':  {'gene': 'ORF1a',  'offset': 3263},
            'NSP6':  {'gene': 'ORF1a',  'offset': 3569},
            'NSP7':  {'gene': 'ORF1a',  'offset': 3859},
            'NSP8':  {'gene': 'ORF1a',  'offset': 3942},
            'NSP9':  {'gene': 'ORF1a',  'offset': 4140},
            'NSP10': {'gene': 'ORF1a',  'offset': 4253},
            'NSP11': {'gene': 'ORF1a',  'offset': 4392},
            'NSP12': {'gene': 'ORF1ab', 'offset':   -9},
            'NSP13': {'gene': 'ORF1ab', 'offset':  923},
            'NSP14': {'gene': 'ORF1ab', 'offset': 1524},
            'NSP15': {'gene': 'ORF1ab', 'offset': 2051},
            'NSP16': {'gene': 'ORF1ab', 'offset': 2397},
            'NS3':   {'gene': 'ORF3a',  'offset':    0},
            'NS6':   {'gene': 'ORF6',   'offset':    0},
            'NS7a':  {'gene': 'ORF7a',  'offset':    0},
            'NS7b':  {'gene': 'ORF7b',  'offset':    0},
            'NS8':   {'gene': 'ORF8',   'offset':    0},
            'NS10':  {'gene': 'ORF10',  'offset':    0},
        }

        df_mutation[['aa_Ref', 'aa_Sub']] = df_mutation['mutation'].str.extract(r':(\D+)\d+(\D+)')
        df_mutation[['gene', 'pos']] = df_mutation['mutation'].str.extract(r'(\w+):\w(\d+)')
        df_mutation = df_mutation.dropna(subset=['gene', 'pos'])
        df_mutation['pos'] = df_mutation['pos'].astype(int)

        for prod in gene_prod:
            idx = df_mutation.gene==prod
            df_mutation.loc[idx, 'gene'] = gene_prod[prod]['gene']
            df_mutation.loc[idx, 'pos'] = df_mutation.loc[idx, 'pos'] + gene_prod[prod]['offset']

        # remove rows with invalid positions
        idx = df_mutation.pos<=0
        df_mutation = df_mutation.loc[~idx,:]

        # update mutations to new gene names and positions
        df_mutation['mutation'] = df_mutation['gene']+":"+df_mutation['aa_Ref']+df_mutation['pos'].astype(str)+df_mutation['aa_Sub']

        logging.info(f'Cleaning dataframe...')

        # dropping columns
        df_mutation = df_mutation.drop(columns=['aa_Ref', 'aa_Sub', 'gene', 'pos'])

        return df_meta, df_mutation

    def save_pkl(self, outfile_prefix):
        """
        Saving original (before filtering country and ste) df_meta and df_mutation to .pkl files respectively

        :param outfile_prefix: output filename prefix
        :type  outfile_prefix: str
        """
        logging.info(f'Saving metadata to {outfile_prefix}.pkl...')
        self.df_meta_orig.to_pickle(f'{outfile_prefix}.pkl')
        logging.info(f'Saving metadata to {outfile_prefix}.mutation.pkl...')
        self.df_mutation_orig.to_pickle(f'{outfile_prefix}.mutation.pkl')
        
    def set_country(self, country):
        """
        Specific country and return 
        
        :param country: country
        :type  country: str
        """
        if country:
            self.df_meta     = self.df_meta_orig[self.df_meta_orig['country']==country].copy()
            self.df_mutation = self.df_mutation_orig[self.df_mutation_orig.acc.isin(self.df_meta.acc)].copy()
    
    def set_division(self, division):
        """
        Specific division/state
        
        :param division: state
        :type  division: str
        """
        if division:
            self.df_meta     = self.df_meta_orig[self.df_meta_orig['division']==division].copy()
            self.df_mutation = self.df_mutation_orig[self.df_mutation_orig.acc.isin(self.df_meta.acc)].copy()
