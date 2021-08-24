import os
import logging
import variant_viz.surv_viz.surv_datasource as cdata
import variant_viz.surv_viz.surv_metadata as cm
import variant_viz.surv_viz.surv_plot as cplot

import pandas as pd
import numpy as np
from datetime import date
from math import pi

from bokeh.io import output_file, show, output_notebook
from bokeh.plotting import figure, save
from bokeh.models import Div, Paragraph
from bokeh.layouts import column, row, layout, gridplot
from bokeh.models import (BoxZoomTool, Circle, HoverTool,
                          MultiLine, Plot, Range1d, ResetTool,
                          BasicTicker, ColorBar, LinearColorMapper, 
                          PrintfTickFormatter, LogColorMapper, RangeTool,
                          Panel, Tabs, Legend, ColumnDataSource,
                          DataTable, DateFormatter, TableColumn, BoxAnnotation)
from bokeh.models.annotations import Title
from bokeh.palettes import Spectral11, Blues9, d3, Category20c

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(module)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M',
)

class GISAID_stats():
    def __init__(self, gisaid_tsv=None, gisaid_pkl=None, gisaid_mutation_pkl=None, geo_type='global', country=None, state=None):
        self.country = country
        self.state = state

        surv_viz_data_path = os.path.dirname(cm.__file__)
        self.data = cm.CovidMetadata(
            filename_cnty_tsv     = f"{surv_viz_data_path}/data/country_coordinate.tsv",
            filename_v_info_tsv   = f"{surv_viz_data_path}/data/lineage_info.tsv",
            filename_v_muta_tsv   = f"{surv_viz_data_path}/data/S_protein_variant_mutation.tsv",
            filename_anno_tsv     = f"{surv_viz_data_path}/data/S_protein_anno.tsv",
            filename_meta_tsv     = gisaid_tsv.name if gisaid_tsv else None,
            filename_meta_pkl     = gisaid_pkl.name if gisaid_pkl else None,
            filename_mutation_pkl = gisaid_mutation_pkl.name if gisaid_mutation_pkl else None,
            country               = country,
            division              = state
        )

    def tsv2pkl(self, prefix):
        """
        Convert GISAID metadata.tsv to .pkl files
        """
        self.data.save_pkl(prefix)

    def generate_stats_html(self, geo_type, outfile):
        """
        Generate stats html file
        """
        region_state_map = cdata.StatesRegionsLookup()

        # init datasource
        ds = cdata.PlotDataSource(self.data)
        # init plots
        sars_cov_2_plots = cplot.CovidPlots()
        page = cplot.CovidViz(sars_cov_2_plots)

        # plot geo visualization
        ds.prepare_geo_country()

        if geo_type=='global':
            sars_cov_2_plots.plot_geo_country(ds)
        elif geo_type=='country':
            ds.prepare_geo_state()
            sars_cov_2_plots.plot_geo_state(ds)
        elif geo_type=='state':
            abbr =  region_state_map.state_abbv_lookup[self.state]
            ds.prepare_geo_county(state=abbr)
            sars_cov_2_plots.plot_geo_county(ds)
                
        # week_lineage_bar_plot
        ds.prepare_week_lineage_bar_ds()
        sars_cov_2_plots.plot_week_lineage_bar(ds)

        # protein_hotspot_plot
        ds.prepare_protein_hotspot()
        sars_cov_2_plots.plot_protein_hotspot(ds)

        # Spike_substitution_tracking
        ds.prepare_spike_substitution_tracking()
        sars_cov_2_plots.plot_spike_substitution_tracking(ds)

        # trending plot
        ds.prepare_trends()
        sars_cov_2_plots.plot_trends(ds)

        # lineage mutation plot
        sars_cov_2_plots.plot_lineage_mutation_tracking(ds)

        # layout
        page.layout()

        # save html
        page.save_html(outfile=outfile)


class EC19_Report():
    """
    Generate a visualization report for EDGE-Covid19 workflow
    """

    def __init__(self, snps, gaps, alnstats, pango, metadata, output):
        surv_viz_data_path = os.path.dirname(cm.__file__)
        self.data = cm.CovidMetadata(
            filename_cnty_tsv     = f"{surv_viz_data_path}/data/country_coordinate.tsv",
            filename_v_info_tsv   = f"{surv_viz_data_path}/data/lineage_info.tsv",
            filename_v_muta_tsv   = f"{surv_viz_data_path}/data/S_protein_variant_mutation.tsv",
            filename_anno_tsv     = f"{surv_viz_data_path}/data/S_protein_anno.tsv",
        )

        self.filename_ec19 = snps
        self.filename_alnstats = alnstats
        self.filename_pango = pango
        self.filename_metadata = metadata
        self.filename_gaps = gaps
        self.filename_output = output

    def generate_EC_repot(self):
        """
        The mehtod is used to generate following plots and save to a html file:
            p_sample
            row(column([p_var_count]+p_feature_list), column(p_pie_list))
            column([p_mut_ref, p_mut_v, p_mut_count, p_mut])
        """
        # loading EC19 SNP report
        df_ec19 = pd.read_csv(self.filename_ec19, sep="\t")
        df_ec19 = df_ec19.drop(df_ec19[df_ec19.Chromosome=="Chromosome"].index)

        # loading EC19 alnstats report 
        df_ec19_alnstats = pd.read_csv(self.filename_alnstats, sep="\t", 
                                    names=["Sample", "Ref_len", "Ref_GC_pct", "Mapped_reads", "Ref_recovery_pct", "Avg_fold_x", "Fold_std", "Num_of_Gap", "Total_Gap_bases", "Num_of_SNPs", "Num_of_INDELs"],
                                    comment='R')
        # loading EC19 pango lineage
        df_ec19_pango = pd.read_csv(self.filename_pango)
        df_ec19_pango.loc[:,'scorpio_call'] = df_ec19_pango.loc[:,'scorpio_call'].replace(regex=r' .*$', value='')
        df_ec19_pango['variant'] = df_ec19_pango['lineage']
        idx = df_ec19_pango['scorpio_call'].notnull()
        df_ec19_pango.loc[idx, 'variant'] = df_ec19_pango.loc[idx, 'scorpio_call']   

        # loading EC19 metadata
        df_ec19_meta = pd.read_csv(self.filename_metadata, sep='\t')
        df_ec19_meta.columns = [*df_ec19_meta.columns[:-1], 'sample']
        df_ec19_meta = df_ec19_meta.drop(df_ec19_meta[df_ec19_meta.virus_name=="virus_name"].index) # removing additional header lines
        df_ec19_meta = df_ec19_meta.dropna(subset=['virus_name']) # removing 'empty' lines

        # loading dropping mutations in intergenic regions
        df_ec19 = df_ec19.dropna(subset=['CDS_start'])

        # remove Synonymous changes
        df_ec19 = df_ec19.drop(df_ec19[df_ec19.Synonymous=="Yes"].index)

        # converting product names to abbreviations
        df_ec19.loc[:,'Product'] = df_ec19.loc[:,'Product'].str.lower()
        df_ec19.loc[:,'Product'] = df_ec19.loc[:,'Product'].replace(regex=[r'^\S+\d+:', ' \S*protein'], value='')
        df_ec19.loc[:,'Product'] = df_ec19.loc[:,'Product'].replace('nucleocapsid', 'N')
        df_ec19.loc[:,'Product'] = df_ec19.loc[:,'Product'].replace('surface', 'S')
        df_ec19.loc[:,'Product'] = df_ec19.loc[:,'Product'].replace('membrane', 'M')
        df_ec19.loc[:,'Product'] = df_ec19.loc[:,'Product'].replace('envelope', 'E')   
        df_ec19.loc[:,'Product'] = df_ec19.loc[:,'Product'].str.capitalize()
        df_ec19['CDS_start']    = df_ec19['CDS_start'].astype(int)
        df_ec19['CDS_end']      = df_ec19['CDS_end'].astype(int)
        df_ec19['SNP_position'] = df_ec19['SNP_position'].astype(int)

        df_ec19['AA_pos'] = np.int16((df_ec19['SNP_position']-df_ec19['CDS_start'])/3)+1
        df_ec19['AA_pos'] = df_ec19['AA_pos'].astype(str)

        df_ec19["Mutation"] = df_ec19.Product + ":" + df_ec19.aa_Ref + df_ec19.AA_pos.astype(str).str.replace(r'\.0$', '') + df_ec19.aa_Sub


        # init datasource
        ds = cdata.PlotDataSource(self.data)

        # adding lineage/variant info to SNP dataframe
        df_ec19_g = df_ec19[df_ec19.Product=='S'].merge(df_ec19_pango[['taxon', 'variant', 'lineage']], left_on='Chromosome', right_on='taxon', how='left').drop(columns=['taxon'])

        df_variant = self.data.df_variant.copy()
        df_variant['pos'] = df_variant['pos'].astype(str)
        df_variant[['aa_Ref','aa_Sub']] = df_variant['mutation'].str.extract(r'\S+:(\w)\d+(\S+)')
        df_variant[df_variant.lineage=='Alpha']

        # EC19 gaps report 
        df_ec19_gaps = pd.read_csv(self.filename_gaps, sep="\t", 
                                names=["sample","gap_start","gap_end","gap_len","Missing","cds_start","cds_end","cds_product"])

        df_ec19_gaps['gap_aa_start'] = np.int64((df_ec19_gaps['gap_start']-df_ec19_gaps['cds_start'])/3)+1
        df_ec19_gaps['gap_aa_end'] = np.int64((df_ec19_gaps['gap_end']-df_ec19_gaps['cds_start'])/3)+1
        df_ec19_gaps = df_ec19_gaps.drop(columns=['cds_start', 'cds_end', 'cds_product'])
        idx = df_ec19_gaps['gap_aa_start']<=0
        df_ec19_gaps.loc[idx, 'gap_aa_start'] = 1

        pos = list(set(df_ec19_g['AA_pos']))+list(set(df_variant['pos']))
        pos = list(set(pos))

        df_ec19_gaps['gap_pos'] = np.nan

        df_ec19_mut_missing = pd.DataFrame()

        for p in pos:
            index = (df_ec19_gaps['gap_aa_start']<=int(p)) & (df_ec19_gaps['gap_aa_end']>=int(p))
            df_ec19_gaps.loc[index,'gap_pos'] = str(p)
            df_ec19_mut_missing = pd.concat([df_ec19_mut_missing, df_ec19_gaps.loc[index,:]])

        df_ref = pd.concat([df_variant[['pos','aa_Ref']], df_ec19_g[['AA_pos','aa_Ref']].rename(columns={'AA_pos':'pos'})]).groupby('pos').agg({'aa_Ref': 'first'}).reset_index()
        df_ref['lineage']='Reference'

        # getting unique mutation for each variant
        variant_mut = df_variant.groupby('lineage')['mutation'].unique()

        df_ec19_g['novel'] = df_ec19_g.apply(lambda x: 'Y' if x.variant in variant_mut and not(x.Mutation in variant_mut[x.variant]) else 'N', axis=1)

        df_pos_aasub = pd.crosstab(df_ec19_g.AA_pos, df_ec19_g.aa_Sub, margins=True)
        ds_pos_aasub = ColumnDataSource(df_pos_aasub)

        #
        # Comparison plot for reference / named variants and EC19 samples
        #
        #
        from bokeh.transform import dodge, factor_cmap, factor_hatch, transform
        from bokeh.palettes import OrRd9, GnBu9, Category20, viridis
        from bokeh.models import CategoricalMapper, Text
        import re

        # add background color to different lineage
        df_sample_variant = df_ec19_g.groupby(['Chromosome','variant','lineage']).count()[['Mutation']].reset_index()

        df_pos_aasub = pd.crosstab(df_ec19_g.AA_pos, df_ec19_g.aa_Sub, margins=True)
        ds_pos_aasub = ColumnDataSource(df_pos_aasub)

        voc_list = df_variant.loc[df_variant.type=='VOC', 'lineage'].unique()
        voi_list = df_variant.loc[df_variant.type=='VOI', 'lineage'].unique()

        cmap = {voc_list[i]: OrRd9[-i-2] for i in range(len(voc_list))}
        cmap.update({voi_list[i]: GnBu9[-i-2] for i in range(len(voi_list))})

        def atoi(text):
            return int(text) if text.isdigit() else text

        pos = list(set(df_ec19_g['AA_pos']))+list(set(df_variant['pos']))
        pos = list(set(pos))
        pos.sort(key=atoi)

        ### Reference sequence (wildtype) plot

        TOOLTIPS = [
            ("Wildtype", "@aa_Ref"),
        ]

        p_mut_ref = figure(
            title='Spike Protein Reference',
            plot_width=len(pos)*15,
            plot_height=85,
            x_range=pos,
            y_range=['Reference'],
            toolbar_location=None, 
            x_axis_location="above",
            y_axis_location="left",
            output_backend="svg",
            tooltips=TOOLTIPS
        )

        # adding reference bases
        rv_ref = p_mut_ref.rect(
            'pos', 'lineage', 1, 1,
            color='grey', alpha=0.4, line_color='white',
            source=df_ref
        )

        # adding bases
        text_props = {"source": df_ref, "text_align": "center", "text_baseline": "middle"}
        x = dodge("pos", 0, range=p_mut_ref.x_range)
        p_mut_ref.text(x=x, y="lineage", text="aa_Ref", text_color="black", text_font_size="10px", **text_props)

        p_mut_ref.axis.axis_line_color = None
        p_mut_ref.axis.major_tick_line_color = None
        p_mut_ref.axis.major_label_standoff = 0
        p_mut_ref.axis.major_label_text_font_size = '11px'
        p_mut_ref.xaxis.major_label_orientation = pi/2
        p_mut_ref.outline_line_color = None
        p_mut_ref.grid.grid_line_color = None
        p_mut_ref.hover.renderers = [rv_ref]
        p_mut_ref.margin = [0, 20, 0, 0]

        ### Variants and type mutations plot

        lineage = df_variant.lineage.unique()

        TOOLTIPS = [
            ("Variant", "@lineage"),
            ("Pangolineage", "@lineage_pango"),
            ("Nextstrain", "@nextstrain_clade"),
            ("Type", "@type"),
            ("Mutation", "@mutation"),
            ("Optional", "@not_all"),
        ]

        p_mut_v = figure(
            plot_width=len(pos)*15,
            plot_height=len(lineage)*15,
            x_range=pos,
            y_range=lineage[::-1],
            toolbar_location=None, 
            x_axis_location="above",
            y_axis_location="left",
            output_backend="svg",
            tooltips=TOOLTIPS,
            y_axis_label='Variant Reference'
        )

        # add background color to different lineage
        p_mut_v.rect(
            [pos[0]]*len(lineage),
            lineage,  
            len(pos)*2, 
            1, 
            line_color=None, 
            alpha=0.3,
            color=[cmap[x] for x in lineage]
        )

        rv = p_mut_v.rect(
            'pos',
            'lineage',
            1,
            1,
            fill_alpha=1,
            color=factor_cmap('lineage', palette=list(cmap.values()), factors=list(cmap.keys())),
            hatch_pattern=factor_hatch('not_all', patterns=[' ', '/'], factors=['N', 'Y']),
            hatch_color='white',
            hatch_alpha=0.9,
            hatch_weight=3,
            line_color='grey',
            line_width=1,  
            source=df_variant
        )

        text_props = {"source": df_variant, "text_align": "center", "text_baseline": "middle"}
        x = dodge("pos", 0, range=p_mut_v.x_range)
        p_mut_v.text(x=x, y="lineage", text="aa_Sub", text_color="black", text_font_size="10px", **text_props)

        p_mut_v.axis.axis_line_color = None
        p_mut_v.axis.major_tick_line_color = None
        p_mut_v.axis.major_label_standoff = 0
        p_mut_v.axis.major_label_text_font_size = '11px'
        p_mut_v.xaxis.major_label_text_font_size = '0px'
        p_mut_v.xaxis.major_label_orientation = pi/2
        p_mut_v.outline_line_color = None
        p_mut_v.xgrid.grid_line_color = None
        p_mut_v.ygrid.grid_line_dash = [6, 4]
        p_mut_v.hover.renderers = [rv]
        p_mut_v.margin = [5, 20, 5, 0]


        ### This block of figure is for EC19 samples. 
        ### A summary plot will be put above the mutation plot.

        sample = list(set(df_ec19_g['Chromosome']))
        sample.sort(reverse=True)


        # raw counting plot on the top  
        pos_std = np.std(ds_pos_aasub.data['All'][:-1])
        pos_avg = np.average(ds_pos_aasub.data['All'][:-1])
        pos_max = np.max(ds_pos_aasub.data['All'][:-1])
        y_range = (0, pos_max+1)
        if (pos_max-pos_avg)/pos_std >= 3:
            y_range = (0, int(pos_avg+pos_std))


        p_mut_count = figure(
            title='Samples',
            plot_width=len(pos)*15,
            plot_height=150,
            x_range=pos,
            y_range=y_range,
            toolbar_location='above',
            tools=['hover','ywheel_zoom','ypan','reset'],
            output_backend="svg",
            tooltips="Position @AA_pos - $name: @$name"
        )

        aasub = list(ds_pos_aasub.data.keys())[1:-1]

        p_mut_count.vbar_stack(
            aasub, x='AA_pos',
            width=0.8,
            color=viridis(len(aasub)),
            alpha=0.5,
            line_width=0,
            legend_label=aasub,
            source=ds_pos_aasub
        )

        p_mut_count.ygrid.grid_line_color = None
        p_mut_count.xgrid.grid_line_color = None
        # p_mut_count.xaxis.major_tick_line_color = None
        p_mut_count.xaxis.major_label_text_font_size = '11px'
        p_mut_count.xaxis.major_label_orientation = pi/2
        p_mut_count.axis.minor_tick_line_color = None
        p_mut_count.xaxis.major_tick_line_color = None
        p_mut_count.xaxis.major_label_standoff = 0

        p_mut_count.toolbar.logo = None
        p_mut_count.legend.location = "top_left"
        p_mut_count.legend.orientation = "horizontal"
        p_mut_count.legend.label_text_font_size = '8pt'
        # Increasing the glyph's label height
        p_mut_count.legend.label_height = 10
        p_mut_count.legend.label_width = 10
        # Increasing the glyph height
        p_mut_count.legend.glyph_height = 10
        p_mut_count.legend.glyph_width = 10



        highlight_hover = HoverTool(
            names=['highlight'],
            point_policy="follow_mouse",
            tooltips=[
                ('Sample', '@Chromosome'), 
                ('Variant', '@variant'), 
                ('PANGO-lineage', '@lineage'),
                ('Total mutation(s)', '@Mutation'),
            ]
        )

        mutation_hover = HoverTool(
            names=['mutation'],
            tooltips=[
                ("Sample", "@Chromosome (@variant)"),
                ("Mutation", "@Mutation"),
                ("SNP position", "@SNP_position"),
                ("Ref codon", "@Ref_codon"),
                ("Sub codon", "@Sub_codon"),
                ("New mutation (Y/N)", "@novel")
            ]
        )

        gaps_hover = HoverTool(
            names=['gaps'],
            tooltips=[
                ("Sample", "@sample"),
                ("Position", "@gap_pos"),
                ("Gap region (genome)", "@gap_start..@gap_end (@gap_len bp)"),
                ("Gap region (protein)", "@gap_aa_start..@gap_aa_end"),
                ("Missing", "@Missing"),
            ]
        )

        p_mut = figure(
            plot_width=p_mut_count.plot_width,
            plot_height=len(sample)*15,
            x_range=pos,
            y_range=sample,
            toolbar_location=None, 
            tools=[highlight_hover, mutation_hover, gaps_hover],
            x_axis_location="above",
            y_axis_location="left",
            output_backend="svg",
            y_axis_label='Sample'
        )

        p_mut.rect(
            0,
            'Chromosome',  
            len(pos)*2, 
            1, 
            line_color=None, 
            alpha=0.3,
            hover_color='#cccccc',
            color=factor_cmap('variant', palette=list(cmap.values()), factors=list(cmap.keys()), nan_color='white'),
            name='highlight',
            source=df_sample_variant
        )

        p_mut.rect( # adding mutations falling in gaps
            'gap_pos',
            'sample',
            1,
            1,
            color='#999999',
            line_width=0,  
            name='gaps',
            source=df_ec19_mut_missing
        )

        p_mut.rect(
            'AA_pos', 'Chromosome',
            1, 1,
            color=factor_cmap('variant', palette=list(cmap.values()), factors=list(cmap.keys()), nan_color='#946D8C'),
            line_color='grey',
            alpha=0.5,
            line_width=1,
            name='mutation',
            source=df_ec19_g
        )

        text_props = {"source": df_ec19_g, "text_align": "center", "text_baseline": "middle"}
        x = dodge("AA_pos", 0, range=p_mut.x_range)

        p_mut.text(
            x=x, y="Chromosome", text="aa_Sub",
            text_color='black',
            text_font_size="10px", **text_props)


        idx=(df_ec19_g.novel=='Y')
        p_mut.circle(
            x=[pos.index(x)+1 for x in df_ec19_g.loc[idx, 'AA_pos']],
            y=[sample.index(x)+1 for x in df_ec19_g.loc[idx, 'Chromosome']],
            color='#d7301f',
            size=12
        )

        source = ColumnDataSource(dict(x=[pos.index(x)+1 for x in df_ec19_g.loc[idx, 'AA_pos']], 
                                    y=[sample.index(x)+1 for x in df_ec19_g.loc[idx, 'Chromosome']],
                                    text=['!']*len(df_ec19_g.loc[idx, 'AA_pos'])))

        text_props = {"text_align": "center", "text_baseline": "middle", "text_font_style": "bold"}

        glyph = Text(
            x='x',
            y='y',
            text='text', 
            text_color="white",
            text_font_size="10px", **text_props)

        p_mut.add_glyph(source, glyph)


        p_mut.xaxis.major_label_text_font_size = '0px'
        p_mut.xaxis.major_label_orientation = pi/2
        p_mut.toolbar.autohide = True
        p_mut.outline_line_color = None
        p_mut.axis.axis_line_color = None
        p_mut.axis.major_tick_line_color = None
        p_mut.axis.major_label_standoff = 0
        p_mut.margin = [0, 20, 0, 0]
        p_mut.xgrid.grid_line_color = None
        p_mut.ygrid.grid_line_dash = [6, 4]


        #show(column([p_mut_ref, p_mut_v, p_mut_count, p_mut]))
        df_ec19_alnstats_m = df_ec19_alnstats.merge(df_ec19_pango[['taxon', 'variant']], left_on='Sample', right_on='taxon', how='left').drop(columns=['taxon'])
        df_ec19_alnstats_m = df_ec19_alnstats_m.merge(df_ec19_meta, left_on='Sample', right_on='sample', how='left').drop(columns=['sample'])
        df_ec19_alnstats_m = df_ec19_alnstats_m.merge(df_ec19.groupby('Chromosome').count()[['Mutation']].reset_index(), left_on='Sample', right_on='Chromosome', how='left')
        df_ec19_alnstats_m.Num_of_SNPs = df_ec19_alnstats_m.Mutation
        df_ec19_alnstats_m['size'] = ds.get_size_bin(df_ec19_alnstats_m.Num_of_SNPs.fillna(0).astype(int), min_s=5, max_s=30)
        df_ec19_alnstats_m['alpha'] = [0.9 if x in cmap else 0.3 for x in df_ec19_alnstats_m['variant']]

        TOOLS = "hover,save,pan,box_zoom,reset,xwheel_zoom,ywheel_zoom"

        TOOLTIPS = [
            ("Sample", "@Sample"),
            ("Name", "@virus_name"),
            ("Collection", "@collection_date"),
            ("Location", "@location"),
            ("Sequencing tech", "@sequencing_technology"),
            ("Metadata", "@virus_passage / @gender / @age"),
            ("variant", "@variant"),
            ("Ref_len", "@Ref_len{0,0}"),
            ("Ref_GC_pct", "@Ref_GC_pct{0.00}%"),
            ("Mapped_reads", "@Mapped_reads{0,0}"),
            ("Ref_recovery_pct", "@Ref_recovery_pct{0.00}%"),
            ("Avg_fold_x", "@Avg_fold_x{0,0}"),
            ("Fold_std", "@Fold_std{0,0}"),
            ("Num_of_Gap", "@Num_of_Gap{0,0}"),
            ("Total_Gap_bases", "@Total_Gap_bases{0,0}"),
            ("Num_of_SNPs", "@Num_of_SNPs{0,0}"),
            ("Num_of_INDELs", "@Num_of_INDELs{0,0}")
        ]

        p_sample = figure(
            title="SAMPLE STATS",
            x_range=(0, 102),
            y_range=(0.1, df_ec19_alnstats_m.Avg_fold_x.max()*5),
            y_axis_type="log",
            tools=TOOLS,
            tooltips=TOOLTIPS,
            toolbar_location='below',
            output_backend="svg",
            plot_width=1000,
            plot_height=450,
            y_axis_label='Avgerage fold (x)',
            x_axis_label='Reference recovery (%)'
        )

        # add a circle renderer with a size, color, and alpha
        p_sample.circle(
            'Ref_recovery_pct', 'Fold_std', 
            size='size',
            line_width=1,
            line_color='grey',
            line_alpha=0.5,
            color=factor_cmap('variant', palette=list(cmap.values()), factors=list(cmap.keys()), nan_color='#946D8C'),
            #color='grey',
            fill_alpha='alpha',
            source=df_ec19_alnstats_m)

        p_sample.toolbar.logo = None
        p_sample.toolbar.autohide = True
        p_sample.margin = [0, 20, 0, 0]

        # show(column(p_sample))

        df_ec19_vcount = df_ec19_alnstats_m.groupby('variant').count()[['Sample']].reset_index()

        TOOLTIPS = [
            ("Lineage", "@variant"),
            ("Count", "@Sample"),
        ]

        p_var_count = figure(
            x_range=df_ec19_vcount.variant,
            y_axis_type='log',
            plot_width=800,
            plot_height=350, 
            title="LINEAGE",
            toolbar_location=None, 
            tooltips=TOOLTIPS,
            output_backend="svg",
            y_axis_label='Count')

        p_var_count.vbar(
            x='variant', 
            top='Sample',
            bottom=0.1,
            width=0.8,
            line_width=1,
            line_color='grey',
            line_alpha=0.5,
            fill_alpha=0.5,
            fill_color=factor_cmap('variant', palette=list(cmap.values()), factors=list(cmap.keys()), nan_color='#946D8C'),
            source=df_ec19_vcount)

        p_var_count.xgrid.grid_line_color = None
        p_var_count.xaxis.major_tick_line_color = None
        p_var_count.y_range.start = 0.1
        p_var_count.y_range.end = df_ec19_vcount.Sample.max()*2
        p_var_count.xaxis.major_label_orientation = pi / 2
        p_var_count.toolbar.logo = None
        p_var_count.toolbar.autohide = True
        p_var_count.margin = [0, 20, 0, 0]

        from bokeh.palettes import Set1

        df_ec19 = df_ec19_alnstats.merge(df_ec19_pango[['taxon','status','variant']], how='left', left_on='Sample', right_on='taxon')
        df_ec19 = df_ec19.merge(df_ec19_meta[['sample','virus_passage','collection_date','gender','age','sequencing_technology','location']], how='left', left_on='Sample', right_on='sample')

        def plot_meta(ds, p_var_count, title):
            p_meta = figure(
                title=title,
                x_range=p_var_count.x_range,
                plot_width=800,
                plot_height=150, 
                toolbar_location=None, 
                output_backend="svg",
                tools="hover", 
                tooltips="@variant $name: @$name")

            cates = list(ds.data.keys())[1::]    
            colors = list(Set1[9][:len(cates)])
            colors[-1] = "#AAAAAA"
            
            p_meta.vbar_stack(
                cates,
                x='variant', 
                width=0.8,
                line_width=1,
                line_color='grey',
                line_alpha=0.5,
                fill_alpha=0.5,
                color=colors,
                legend_label=cates,
                source=ds
            )
            
            p_meta.y_range.start = 0
            p_meta.grid.grid_line_color = None
            p_meta.axis.minor_tick_line_color = None
            p_meta.xaxis.major_tick_line_color = None
            p_meta.xaxis.major_label_text_font_size = '0px'
            p_meta.outline_line_color = None
            p_meta.legend.location = "top_left"
            p_meta.legend.orientation = "horizontal"
            p_meta.legend.click_policy="mute"
            p_meta.legend.label_text_font_size = '8pt'
            # Increasing the glyph's label height
            p_meta.legend.label_height = 10
            p_meta.legend.label_width = 10
            # Increasing the glyph height
            p_meta.legend.glyph_height = 10
            p_meta.legend.glyph_width = 10
            p_meta.legend.background_fill_alpha = 0.8
            p_meta.margin = [0, 20, 0, 0]

            return p_meta

        p_feature_list = []

        # age
        bins = list(range(0, 100, 10))
        df_ec19['age_range'] = 'N/A'
        idx = df_ec19.age.notnull()
        df_ec19.loc[idx, 'age_range'] = pd.cut(df_ec19.loc[idx, 'age'].astype(int), bins=bins)
        df_ec19_feature = pd.crosstab(df_ec19['age_range'], df_ec19['variant'], dropna=True, normalize='columns')

        df_ec19_feature = df_ec19_feature.reset_index()
        df_ec19_feature['age_range'] = df_ec19_feature['age_range'].astype(str)

        ds = ColumnDataSource(df_ec19_feature.set_index('age_range').T)

        p_feature_list.append(plot_meta(ds, p_var_count, 'AGE RANGE')) 

        # gender
        idx = df_ec19.gender.isnull()
        df_ec19.loc[idx, 'gender'] = 'N/A'
        df_ec19_feature = pd.crosstab(df_ec19['gender'], df_ec19['variant'], dropna=True, normalize='columns')
        ds = ColumnDataSource(df_ec19_feature.T)
        p_feature_list.append(plot_meta(ds, p_var_count, 'GENDER')) 

        # sequencing_technology
        idx = df_ec19.sequencing_technology.isnull()
        df_ec19.loc[idx, 'sequencing_technology'] = 'N/A'
        df_ec19_feature = pd.crosstab(df_ec19['sequencing_technology'], df_ec19['variant'], dropna=True, normalize='columns')
        ds = ColumnDataSource(df_ec19_feature.T)
        p_feature_list.append(plot_meta(ds, p_var_count, 'SEQUENCING TECH')) 

        # # submitting_lab
        # idx = df_ec19.submitting_lab.isnull()
        # df_ec19.loc[idx, 'submitting_lab'] = 'N/A'
        # df_ec19_feature = pd.crosstab(df_ec19['submitting_lab'], df_ec19['variant'], dropna=True, normalize='columns')
        # ds = ColumnDataSource(df_ec19_feature.T)
        # p_sub_lab = plot_meta(ds, p_var_count, 'SUBMITTING LAB')
        # p_sub_lab.xaxis.major_label_orientation = pi / 2
        # p_sub_lab.xaxis.major_label_text_font_size = '10px'
        # p_sub_lab.plot_height = p_sub_lab.plot_height+60
        # p_feature_list.append(p_sub_lab) 

        # location
        idx = df_ec19.location.isnull()
        df_ec19.loc[idx, 'location'] = 'N/A'
        df_ec19_feature = pd.crosstab(df_ec19['location'], df_ec19['variant'], dropna=True, normalize='columns')
        ds = ColumnDataSource(df_ec19_feature.T)
        p_sub_loc = plot_meta(ds, p_var_count, 'Location')
        p_sub_loc.xaxis.major_label_orientation = pi / 2
        p_sub_loc.xaxis.major_label_text_font_size = '10px'
        p_sub_loc.plot_height = p_sub_loc.plot_height+60
        p_feature_list.append(p_sub_loc)


        # show(column(p_feature_list))

        from bokeh.transform import cumsum

        p_pie_list = []

        def plot_pie_meta(df, column):
            cates = list(df.index)[:-1] 
            colors = list(Set1[9][:len(cates)])
            colors[-1] = "#AAAAAA"

            df = df.reset_index()
            df['angle'] = df['All']/df.iloc[-1,-1] * 2*pi
            df.drop(df.tail(1).index,inplace=True)
            df['color'] = colors

            p_pie_meta = figure(plot_width=280, plot_height=250, 
                                title=column.upper().replace('_',' '), 
                                toolbar_location=None,
                                tools="hover", tooltips=f"@{column}: @All", x_range=(-0.5, 1.0))

            p_pie_meta.wedge(x=0, y=1, 
                            radius=0.4,
                            alpha=0.6,
                            start_angle=cumsum('angle', include_zero=True), end_angle=cumsum('angle'),
                            line_color="white", fill_color='color', legend_field=column, source=df)

            p_pie_meta.axis.axis_label=None
            p_pie_meta.axis.visible=False
            p_pie_meta.grid.grid_line_color = None
            p_pie_meta.outline_line_color=None
            p_pie_meta.legend.label_text_font_size = '8pt'
            
            # Increasing the glyph's label height
            p_pie_meta.legend.label_height = 10
            p_pie_meta.legend.label_width = 10
            # Increasing the glyph height
            p_pie_meta.legend.glyph_height = 10
            p_pie_meta.legend.glyph_width = 10
            p_pie_meta.legend.background_fill_alpha = 0.8
            p_pie_meta.margin = [0, 20, 0, 0]
            
            return p_pie_meta


        # age
        bins = list(range(0, 100, 10))
        df_ec19['age_range'] = 'N/A'
        idx = df_ec19.age.notnull()
        df_ec19.loc[idx, 'age_range'] = pd.cut(df_ec19.loc[idx, 'age'].astype(int), bins=bins)
        df_ec19_feature = pd.crosstab(df_ec19['age_range'], df_ec19['variant'], dropna=True, margins=True, margins_name='All')
        df_ec19_feature = df_ec19_feature.reset_index()
        df_ec19_feature['age_range'] = df_ec19_feature['age_range'].astype(str)
        p_pie_list.append(plot_pie_meta(df_ec19_feature, 'age_range')) 

        # gender
        idx = df_ec19.gender.isnull()
        df_ec19.loc[idx, 'gender'] = 'N/A'
        df_ec19_feature = pd.crosstab(df_ec19['gender'], df_ec19['variant'], dropna=True, margins=True, margins_name='All')
        p_pie_list.append(plot_pie_meta(df_ec19_feature, 'gender'))

        # sequencing_technology
        idx = df_ec19.sequencing_technology.isnull()
        df_ec19.loc[idx, 'sequencing_technology'] = 'N/A'
        df_ec19_feature = pd.crosstab(df_ec19['sequencing_technology'], df_ec19['variant'], dropna=True, margins=True, margins_name='All')
        p_pie_list.append(plot_pie_meta(df_ec19_feature, 'sequencing_technology'))

        # # submitting_lab
        # idx = df_ec19.submitting_lab.isnull()
        # df_ec19.loc[idx, 'submitting_lab'] = 'N/A'
        # df_ec19_feature = pd.crosstab(df_ec19['submitting_lab'], df_ec19['variant'], dropna=True, margins=True, margins_name='All')
        # p_pie_list.append(plot_pie_meta(df_ec19_feature, 'submitting_lab'))
        
        # location
        idx = df_ec19.location.isnull()
        df_ec19.loc[idx, 'location'] = 'N/A'
        df_ec19_feature = pd.crosstab(df_ec19['location'], df_ec19['variant'], dropna=True, margins=True, margins_name='All')
        p_pie_list.append(plot_pie_meta(df_ec19_feature, 'location'))

        # show(column([p_sample, p_var_count]+p_feature_list+[p_mut_ref, p_mut_v, p_mut])) #, sizing_mode='stretch_width'

        update_date = date.today().strftime("%Y-%m-%d")

        pango_ver = df_ec19_pango.loc[0, 'version']
        pangolin_ver = df_ec19_pango.loc[0, 'pangolin_version']
        pangolearn_ver = df_ec19_pango.loc[0, 'pangoLEARN_version']
        sample_num = len(df_ec19_alnstats)
        submitted_num = len(df_ec19_meta)

        #
        div_header = Div(text=f"""
        <H2>EDGE COVID-19 Report</H2>
        <div style='margin-bottom: 1em'>
        {sample_num} samples / {submitted_num} submitted / Pango ver: {pango_ver} / Pangolin ver: {pangolin_ver} / PangoLEARN ver: {pangolearn_ver} / Report date: {update_date}
        </div>
        """, style={'width':'100%', 'margin-left':'10px', 'margin-top': '20px'})

        output_layout = gridplot([
                [div_header],
                [row(p_sample, sizing_mode='stretch_width')],
                [row(column([p_var_count]+p_feature_list, sizing_mode='stretch_width'), column(p_pie_list))],
                [column([p_mut_ref, p_mut_v, p_mut_count, p_mut])],
            ]
            , merge_tools=False
            #, sizing_mode='stretch_width'
        )


        template = """
        {% block postamble %}
        <link rel="stylesheet" href="https://stackpath.bootstrapcdn.com/bootstrap/4.3.1/css/bootstrap.min.css">
        <style type="text/css">
        @import url('https://fonts.googleapis.com/css?family=Quicksand');
        H1, H2, H3 {
            margin-top: 1em;
            font-weight: 600;
            font-family: Quicksand;
        }
        H1 {
            font-size: 2rem;
        }
        H2 {
            font-size: 1.5rem;
        }

        .bk-tabs-header .bk-tab.bk-active {
            background-color: #e6e6e6 !important;
        }
        .plot-tooltip {
            margin: 10px;
        }
        .bk-root {
            display: flex !important;
            justify-content: center !important;
        }
        </style>
        {% endblock %}
        """

        output_file(self.filename_output, title="EDGE COVID-19 Report")
        save(output_layout, template=template)