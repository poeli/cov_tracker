from math import pi
from bokeh.models import ColumnDataSource
import numpy as np
import pandas as pd
import logging

class PlotDataSource(object):
    """
    This is the class for generating bokeh datasource
    
    color_mapper(num, pt='Category20')
    voc_voi_line_color(lineage_list)
    voc_hatch_pattern(lineage_list)
    get_size_bin(min_s=5, max_s=50)
    prepare_week_lineage_bar_ds()
    prepare_geo_country()
    prepare_geo_state()
    prepare_geo_county()
    prepare_spike_protein_hotspot()
    prepare_protein_hotspot()
    prepare_week_spike_mutation()
    """
    
    def __init__(self, data, target_lineage=None):

        self.data = data

        # init target_lineage, voc_list, voi_list
        self.target_lineage = target_lineage
        self.voc_list = self.data.df_v_info.loc[self.data.df_v_info['type']=='VOC', 'lineage'].to_list()
        self.voi_list = self.data.df_v_info.loc[self.data.df_v_info['type']=='VOI', 'lineage'].to_list()
        
        # find top 10 lineages
        top10_list = []
        if len(self.data.df_meta):
            top10_list = self.data.df_meta.groupby('lineage').size().reset_index().sort_values(0, ascending=False).lineage[:10].to_list()
        logging.debug(f"top10_list: {top10_list}")

        # display top 10 and VOC/VOI lineages, rest of them go to "others"
        self.dis_lineage = [x for x in top10_list if x not in self.voc_list+self.voi_list]
        self.dis_lineage = ["Others"] + self.voc_list+ self.voi_list + self.dis_lineage

        # adding target lineage to list of displaying lineages
        if target_lineage and not target_lineage in self.dis_lineage:
            self.dis_lineage = self.dis_lineage + [target_lineage]

        # set lineage colors
        self.lineage_colors = self.color_mapper(len(self.dis_lineage))

        # set lineage patterns
        self.lineage_patterns = self.variant_hatch_pattern(self.dis_lineage)
        
        # assign colors to subtitutions
        df_variant = self.data.df_v_info.merge(self.data.df_v_muta, on='lineage', how='left')    
        df_variant['color'] = '#AAAAAA'
        df_variant.loc[(df_variant['type']=='VOC') & (df_variant['not_all']=='N'), 'color'] = '#D03A49'
        df_variant.loc[(df_variant['type']=='VOI') & (df_variant['not_all']=='N'), 'color'] = '#FDE724'
        df_variant.loc[(df_variant['type']=='EDGE') & (df_variant['not_all']=='N'), 'color'] = '#F87748'
        self.data.df_variant = df_variant
        
        # axis labels
        self.week_ticks = None
        
        # dataScource
        self.ds_variant = ColumnDataSource(df_variant)
        self.ds_week_variant = None
        self.ds_week_variant_freq = None
        self.ds_week_lineage = None
        self.ds_week_lineage_freq = None
        self.ds_geo_country = None
        self.ds_geo_states = None
        self.ds_geo_county = None
        self.ds_protein_hotspot = None
        self.ds_trend_lineage = None
        self.ds_trend_position = None
        self.trend_lineage_10 = None
        self.trend_position_10 = None
        

    def color_mapper(self, num, pt='Category20'):
        from bokeh.palettes import all_palettes
        try:
            return ['#AAAAAA']+list(all_palettes[pt][num-1])
        except:
            return ['#AAAAAA']+list(all_palettes['Turbo256'][num-1])

        
    def variant_color(self, lineage_list):
        colors = []
        for lin in lineage_list:
            if lin in self.voc_list:
                colors.append('#D44A58')
            elif lin in self.voi_list:
                colors.append('#FEDC85')
            else:
                colors.append('#FFFFFF')

        return colors

    
    def variant_hatch_pattern(self, lineage_list):
        pattern = []
        for lin in lineage_list:
            if lin in self.voc_list:
                pattern.append('x')
            elif lin in self.voi_list:
                pattern.append('/')
            else:
                pattern.append(' ')

        return pattern
    
    
    def get_size_bin(self, ser_val, min_s=5, max_s=50):
        intv = int(np.ceil(ser_val.max()/(max_s-min_s)))
        bins = list(range(0, ser_val.max(), intv))
        bins[-1] = ser_val.max()
        labels = list(range(min_s, max_s))[-len(bins)+1:] #labels must be one fewer than the number of bin edges

        if len(ser_val) > 1:
            s = pd.cut(x=ser_val, 
                       bins=bins,
                       labels=labels)
        else:
            s = max_s

        return s
    
    
    def prepare_week_lineage_bar_ds(self):
        logging.info(f'Preparing datasource for week lineage bar...')
        df = self.data.df_meta.copy()
        col = ~df.lineage.isin(self.dis_lineage)
        df.loc[col, 'lineage'] = "Others"

        ds_week_variant = ColumnDataSource(pd.crosstab(df.week, df.lineage))
        ds_week_variant_freq = ColumnDataSource(pd.crosstab(df.week, df.lineage, normalize='index'))

        ds_week_lineage = ColumnDataSource(pd.crosstab(df.week, df.pango_lineage))
        ds_week_lineage_freq = ColumnDataSource(pd.crosstab(df.week, df.pango_lineage, normalize='index'))

        self.week_ticks = df.week.unique()
        self.week_ticks.sort()
        self.ds_week_variant = ds_week_variant
        self.ds_week_variant_freq = ds_week_variant_freq

        self.ds_week_lineage = ds_week_lineage
        self.ds_week_lineage_freq = ds_week_lineage_freq    
    
    def prepare_geo_country(self):
        logging.info(f'Preparing datasource for geographical plot in country-scale...')
        df = self.data.df_meta
        df_group_country = df.groupby('country').count()[['acc']].reset_index()
        df_group_country = df_group_country.merge(self.data.df_country, on='country')
        df_group_country['size'] = self.get_size_bin(df_group_country.acc, min_s=10, max_s=200)
        ds_geo_country = ColumnDataSource(df_group_country)
    
        self.ds_geo_country = ds_geo_country
    
    
    def prepare_geo_state(self, periods=90):
        logging.info(f'Preparing datasource geographical plot in state-scale...')
        
        import bokeh.sampledata
        bokeh.sampledata.download()
        # import pyproj for geo coordinates
        from pyproj import Transformer, CRS
        from bokeh.transform import cumsum
        from bokeh.sampledata.us_states import data as states
        from datetime import datetime

        df_meta = self.data.df_meta

        # prepare data
        df_meta['division'] = df_meta['division'].str.replace(' /', '')
        df_periods = df_meta[df_meta['date'].isin(pd.date_range(pd.Timedelta(-periods, unit='d')+datetime.today().date(), periods=periods))].copy()

        lineage_list = df_periods.groupby('pango_lineage').size().reset_index().sort_values(0, ascending=False).pango_lineage.to_list()
        top10_list = lineage_list
        if len(lineage_list) > 10:
            top10_list = lineage_list[:10]
        if self.target_lineage and not self.target_lineage in top10_list:
            top10_list.append(self.target_lineage)

        logging.debug(f'toplist: {top10_list}')
        
        df_periods['lineage_type'] = df_periods['pango_lineage']
        df_periods.loc[~df_periods['pango_lineage'].isin(top10_list), 'lineage_type'] = 'Others'
        
        # sequence count for each state (division)
        df_periods_cnt = df_periods.groupby(['division']).count()[['name']].rename(columns={'name': 'count'})
        df_periods_states = df_periods.groupby(['division', 'lineage_type']).count()[['name']].rename(columns={'name': 'count'})
        df_periods_states['prop'] = df_periods_states.div(df_periods_cnt)
        
        # generate states datasource
        state_lons = []
        state_lats = []
        state_names = []
        state_count = []
        state_xs = []
        state_ys = []

        for code in states:
            ## skip HI, AK
            #if code=='AK':
            #    continue

            transformer = Transformer.from_crs("EPSG:4326", "EPSG:3857", always_xy=True)
            xs, ys = transformer.transform(states[code]["lons"], states[code]["lats"])
            name = states[code]['name']
            state_lons.append(xs)
            state_lats.append(ys)
            state_names.append(name)
            state_xs.append((max(xs)+min(xs))/2)
            state_ys.append((max(ys)+min(ys))/2)

            if name in df_periods_cnt.index:
                state_count.append(df_periods_cnt.loc[name, 'count'])
            else:
                state_count.append(0)

        source = ColumnDataSource(data=dict(
            lons=state_lons,
            lats=state_lats,
            x=state_xs,
            y=state_ys,
            name=state_names,
            count=state_count,
        ))

        # create circle to indicate the proportion of VOC/VOI
        df_state_xy = pd.DataFrame({'division': state_names, 'x': state_xs, 'y': state_ys})

        df_periods_states = df_periods_states.reset_index()
        df_periods_states = df_periods_states.merge(df_state_xy, on='division')
        df_periods_states          = df_periods_states.rename(columns={'division': 'name'})
        df_periods_states['angle'] = df_periods_states['prop']*2*pi
        
        self.ds_geo_state = source
        self.ds_geo_state_lineage = ColumnDataSource(df_periods_states)
        
    
    def prepare_geo_county(self, state, periods=90):
        logging.info(f'Preparing datasource geographical plot in county-scale...')
        import bokeh.sampledata
        bokeh.sampledata.download()
        # import pyproj for geo coordinates
        from pyproj import Transformer, CRS
        from bokeh.transform import cumsum
        from bokeh.sampledata.us_counties import data as counties
        from bokeh.sampledata.us_states import data as states
        from datetime import datetime

        df_meta = self.data.df_meta

        # prepare data
        df_meta['location'] = df_meta['location'].str.replace(' /', '')
        df_periods = df_meta[df_meta['date'].isin(pd.date_range(pd.Timedelta(-periods, unit='d')+datetime.today().date(), periods=periods))].copy()

        lineage_list = df_periods.groupby('pango_lineage').size().reset_index().sort_values(0, ascending=False).pango_lineage.to_list()
        top10_list = lineage_list
        if len(lineage_list) > 10:
            top10_list = lineage_list[:10]
        if self.target_lineage and not self.target_lineage in top10_list:
            top10_list = top10_list+[lineage_list]

        df_periods['lineage_type'] = df_periods['pango_lineage']
        df_periods.loc[~df_periods['pango_lineage'].isin(top10_list), 'lineage_type'] = 'Others'
        
        # sequence count for each county (location)
        df_periods_cnt = df_periods.groupby(['location']).count()[['name']].rename(columns={'name': 'count'})
        df_periods_counties = df_periods.groupby(['location', 'lineage_type']).count()[['name']].rename(columns={'name': 'count'})
        df_periods_counties['prop'] = df_periods_counties.div(df_periods_cnt)

        # transformer
        transformer = Transformer.from_crs("EPSG:4326", "EPSG:3857", always_xy=True)
        
        # state info
        state_info = states[state.upper()]
        state_xs, state_ys = transformer.transform(state_info["lons"], state_info["lats"])

        # county
        counties = {
            code: county for code, county in counties.items() if county["state"] == state.lower()
        }

        # generate counties datasource
        county_lons = []
        county_lats = []
        county_names = []
        county_count = []
        county_xs = []
        county_ys = []

        for county in counties.values():
            xs, ys = transformer.transform(county["lons"], county["lats"])
            name = county['name']
            county_lons.append(xs)
            county_lats.append(ys)
            county_names.append(name)
            county_xs.append((max(xs)+min(xs))/2)
            county_ys.append((max(ys)+min(ys))/2)

            if name in df_periods_cnt.index:
                county_count.append(df_periods_cnt.loc[name, 'count'])
            else:
                county_count.append(0)

        source = ColumnDataSource(data=dict(
            lons =county_lons,
            lats =county_lats,
            x    =county_xs,
            y    =county_ys,
            name =county_names,
            count=county_count,
        ))        

        # create circle to indicate the proportion of VOC/VOI
        df_county_xy = pd.DataFrame({'location': county_names, 'x': county_xs, 'y': county_ys})

        df_periods_counties          = df_periods_counties.reset_index()
        df_periods_counties          = df_periods_counties.merge(df_county_xy, on='location')
        df_periods_counties          = df_periods_counties.rename(columns={'location': 'name'})
        df_periods_counties['angle'] = df_periods_counties['prop']*2*pi

        self.ds_geo_county = source
        self.ds_geo_county_lineage = ColumnDataSource(df_periods_counties)
    
    def prepare_protein_hotspot(self):
        logging.info(f'Preparing datasource for protein hotspots...')
        df_gene_mutation = self.data.df_mutation.groupby(['gene','pos']).count()[['acc']].reset_index()
        df_gene_mutation['pos'] = df_gene_mutation['pos'].astype(int)
        df_gene_mutation['size'] = self.get_size_bin(df_gene_mutation.acc)

        ds_protein_hotspot = ColumnDataSource(df_gene_mutation.set_index('pos'))
        self.ds_protein_hotspot = ds_protein_hotspot
        
        #display(df_gene_mutation)
        #    gene	pos	acc
        # 0	    E	1	1
        # 1	    E	10	11
        # 2	    E	11	6
        # 3	    E	12	2
        # 4	    E	13	11
        # ...	...	...	...
        # 9091	S	994	1
        # 9092	S	995	1
        # 9093	S	997	2
        # 9094	S	998	2
        # 9095	S	999	2

        
    def prepare_spike_substitution_tracking(self, freq_cutoff=0.001):
        logging.info(f'Preparing datasource for mutation tracking...')

        # mutation position on S protein
        df_mutation_S = self.data.df_mutation[self.data.df_mutation.gene=='S']
        df = df_mutation_S.groupby(['week','pos']).count()['acc'].reset_index()
        
        # circle size
        df['size'] = self.get_size_bin(df.acc, min_s=5, max_s=50)
        
        # total genomes at the position
        df_pos_total = df.groupby(['pos']).sum().reset_index()

        # week total
        df_meta_week_total = self.data.df_meta.groupby('week').count()[['acc']].rename(columns={'acc':'total'})
        self.data.df_meta_week_total = df_meta_week_total
        
        # get prop
        df['prop'] = df['acc']/df[['week','acc']].merge(df_meta_week_total.reset_index(), on='week').loc[:,'total']
        
        self.ds_spike_substitution = ColumnDataSource(df[df['prop']>=freq_cutoff])
        self.ds_week_total = ColumnDataSource(df_meta_week_total)
        self.ds_pos_total = ColumnDataSource(df_pos_total)

        
    def prepare_trends(self, **kwargs):
        logging.info(f'Preparing datasource for lineage and mutation trendings...')
        from bokeh.models import ColumnDataSource
        
        df_tl = self.data.df_meta
        df_tlmu = self.data.df_mutation
        for col, value in kwargs.items():
            df_tl = df_tl[df_tl[col]==value].copy()
            df_tlmu = df_tlmu[df_tlmu[col]==value].copy()

        # get Lineage data from -9 week to current week    
        # df_tl = pd.crosstab(df_tl.week, df_tl.lineage, normalize='index')
        df_tl = pd.crosstab(df_tl.week, df_tl.pango_lineage, normalize='index')
        df_tl = df_tl.iloc[-9:, :].copy()
        #df_tl = df_tl.iloc[-9:-1, :].copy()

        # offset the first week to 0 pct
        df_tl_offset = df_tl-df_tl.iloc[0,:]

        # select the first 10 trending lineages based on increasing freq
        # var_list = list(df_tl_offset.iloc[-4:-2,:].sum(axis=0).sort_values(ascending=False).head(10).index)
        var_list = list(df_tl_offset.apply(lambda x: np.average(x, weights=np.array(range(1,len(df_tl_offset)+1))), axis=0).sort_values(ascending=False).head(10).index)
        
        df_lineage_trend = df_tl_offset.loc[:,var_list].reset_index()

        # get Mutation data from -9 week to current week
        df_mu = df_tlmu.loc[df_tlmu.gene=='S',:].copy()
        df_mu = pd.crosstab(df_mu.week, df_mu.mutation)
        df_week_total = df_tlmu.groupby('acc').agg({'week':'first'}).reset_index().groupby('week').count()
        df_mu = df_mu.div(df_week_total['acc'], axis=0)
        df_mu = df_mu.iloc[-9:, :].copy()

        # offset the first week to 0 pct
        df_mu_offset = df_mu-df_mu.iloc[0,:]

        # select the first 10 trending lineages based on increasing freq of the past 3 weeks
        # mu_list = list(df_mu_offset.iloc[-4:-2, :].sum(axis=0).sort_values(ascending=False).head(10).index)
        mu_list = list(df_mu_offset.apply(lambda x: np.average(x, weights=np.array(range(1,len(df_mu_offset)+1))), axis=0).sort_values(ascending=False).head(10).index)

        df_position_trend = df_mu_offset.loc[:, mu_list].reset_index()

        self.ds_trend_lineage = ColumnDataSource(df_lineage_trend)
        self.ds_trend_position = ColumnDataSource(df_position_trend)
        self.trend_lineage_10 = var_list
        self.trend_position_10 = mu_list

    def prepare_mutation_tracking_per_variant(self, variant_name, gene=None):
        logging.info(f'Preparing datasource for mutations of variant: {variant_name}...')
        
        df_mutation = self.data.df_mutation
        if gene:
            # mutation in spike protein
            df_mutation = self.data.df_mutation[self.data.df_mutation.gene==gene]
        
        # mutations in a variant
        mutation_list = self.data.df_variant[self.data.df_variant.lineage==variant_name].mutation.to_list() 
        logging.debug(f'mutation_list: {mutation_list}')        

        # filter in specific variants only (heatmap source: df_lineage)
        df_mutation_lineage = df_mutation[df_mutation['mutation'].isin(mutation_list)]
        df_lineage = df_mutation_lineage.groupby(['week','mutation']).count()['acc'].reset_index()
        df_lineage = df_lineage.rename(columns={'acc':'count'})
        
        df_lineage = df_lineage.merge(self.data.df_meta_week_total, on='week', how='left')
        df_lineage['prop'] = df_lineage['count']/df_lineage['total']

        # resorting mutation by first appearance
        m_fisrt_app = list(df_lineage.sort_values(['week'])['mutation'].unique())
        mutation_list.sort(key=lambda x: m_fisrt_app.index(x) if x in m_fisrt_app else 999, reverse=True)
        logging.debug(f'df_lineage: {df_lineage}')        
        
        return df_lineage, mutation_list
