from math import pi
from bokeh.layouts import column, row, layout, gridplot
from bokeh.plotting import figure
from bokeh.models.annotations import Title
from bokeh.palettes import Spectral11, Blues9, d3, RdYlBu10
from bokeh.models import (Circle, HoverTool,
                          Plot, Range1d, ResetTool,
                          BasicTicker, ColorBar, LinearColorMapper, 
                          PrintfTickFormatter, LogColorMapper, RangeTool,
                          Panel, Tabs, Legend,
                          DataTable, DateFormatter, TableColumn, BoxAnnotation)
from bokeh.tile_providers import CARTODBPOSITRON, get_provider
from bokeh.models import Label, LabelSet, BoxAnnotation
import logging
import numpy as np

class CovidPlots(object): 
    
    def __init__(self, main_plot_width=1000):
        self.main_plot_width = main_plot_width
        self.major_label_text_font_size = "10px"

        self.geo_plot = None
        self.week_lineage_bar_plot = None
        self.trending_plot = None
        self.protein_hotspot_plot = None
        self.spike_substitution_tracking_plot = None
        self.lineage_mutation_tracking_plot = None
        self.lineage_region_tracking_plot = None
        self.lineage_states_tracking_plot = None        
        
    def plot_week_lineage_bar(self, ds):
        """ The week-lineage bar chart """
        
        logging.info(f'Plotting week lineage bar plot...')

        def lineage_plot(source):
            p = figure(
                plot_width=self.main_plot_width, 
                plot_height=300,
                x_range=ds.week_ticks,
                tools='wheel_zoom,pan,reset,hover',
                toolbar_location='below',
                output_backend="webgl",
                tooltips="$name: @$name (@week)"
            )

            p.axis.major_label_text_font_size = self.major_label_text_font_size
            p.xaxis.minor_tick_line_color = None
            p.grid.grid_line_color = None
            p.xaxis.major_label_orientation = pi/3
            p.y_range.start=0
            
            p.vbar_stack(
                ds.dis_lineage,
                x='week',
                width=0.8,
                alpha=0.7,
                color=ds.lineage_colors,
                line_color='white',
                hatch_pattern=ds.lineage_patterns,
                hatch_color='white',
                hatch_alpha=0.3,
                hatch_weight=1,
                legend_label=ds.dis_lineage,
                source=source
            )
            
            return p

        # week_lineage_bar_plot raw number
        p_tl = lineage_plot(ds.ds_week_variant)
        p_tl.legend.location = 'top_left'
        p_tl.legend.label_text_font_size = "8pt"
        p_tl.legend.orientation = "horizontal"

        t = Title()
        t.text = 'Variant'
        p_tl.title = t

        # legend = Legend(items=[(x, [p_tl[i]]) for i, x in enumerate(ds.dis_lineage)], location=(0, -30))
        # p_tl.add_layout(legend, 'left')

        # week_lineage_bar_plot raw freq
        p_tl_freq = lineage_plot(ds.ds_week_variant_freq)
        p_tl_freq.legend.visible = False
        
        # plots in grid
        p_week_lineage = gridplot([[p_tl],[p_tl_freq]], toolbar_options=dict(logo=None))

        self.week_lineage_bar_plot = p_week_lineage
        
        
    def plot_geo_country(self, ds):
        """ The country-scale geo chart """
        
        logging.info(f'Plotting geographical plot in country-scale...')

        def geo_plot(tile_source, ds):
            # set to roughly extent of points
            x_range = Range1d(start=ds.data['x'].min()-4000000, end=ds.data['x'].max()-4000000, bounds=None)
            y_range = Range1d(start=ds.data['y'].min()-5000000, end=ds.data['y'].max()+5000000, bounds=None)

            # create plot and add tools
            p = figure(title='Genome sequence distribution',
                       tools='hover,wheel_zoom,pan,reset',
                       x_range=x_range, 
                       y_range=y_range, 
                       tooltips=[
                           ("Country", "@country"), 
                           ("Count", "@acc{0,0}")
                       ],
                       toolbar_location='below',
                       output_backend="svg",
                       plot_width=self.main_plot_width, 
                       plot_height=600)
            
            p.axis.visible = False
            p.add_tile(tile_source)
            p.toolbar.logo = None
            p.toolbar.autohide = True

            colors = ["#BBBBBB"]+list(Spectral11)

            circle_color_mapper = LogColorMapper(
                palette=colors,
                low=10,
                high=ds.data['acc'].max()
            )

            # create point glyphs
            p.circle(x='x', 
                     y='y', 
                     size='size', 
                     color={'field': 'acc', 'transform': circle_color_mapper},
                     alpha=0.6,
                     line_color="white", 
                     line_width=1, 
                     source=ds)
            
            return p

        self.geo_plot = geo_plot(get_provider(CARTODBPOSITRON), ds.ds_geo_country)

    def plot_geo_state(self, ds):
        logging.info(f'Plotting geographical plot in state-scale...')

        ds_regions = ds.ds_geo_state
        ds_regions_lineage = ds.ds_geo_state_lineage
        ds_boundry = ds.ds_geo_country

        # set to roughly extent of points
        x_range = Range1d(start=min(ds_boundry.data['x'])-3000000, end=max(ds_boundry.data['x'])-3000000, bounds=None)
        y_range = Range1d(start=min(ds_boundry.data['y'])-4000000, end=max(ds_boundry.data['y'])+4000000, bounds=None)

        self.geo_plot = self.plot_geo_regions(ds_regions, ds_regions_lineage, ds_boundry, x_range, y_range)

    def plot_geo_county(self, ds):
        logging.info(f'Plotting geographical plot in county-scale...')

        ds_regions = ds.ds_geo_county
        ds_regions_lineage = ds.ds_geo_county_lineage
        ds_boundry = ds_regions
                
        # set to roughly extent of points
        x_range = Range1d(start=min(ds_boundry.data['x'])*0.95, end=max(ds_boundry.data['x'])*1.05, bounds=None)
        y_range = Range1d(start=min(ds_boundry.data['y'])*0.95, end=max(ds_boundry.data['y'])*1.05, bounds=None)

        self.geo_plot = self.plot_geo_regions(ds_regions, ds_regions_lineage, ds_regions, x_range, y_range)
    
    def plot_geo_regions(self, regions_ds, regions_lineage_ds, boundrey_ds, x_range, y_range):
        from bokeh.transform import cumsum, factor_cmap
        from bokeh.palettes import brewer, Spectral11
            
        # init variables
        tile_source = get_provider(CARTODBPOSITRON)
        
        # set color mapper
        colors = ["#BBBBBB"]+list(Spectral11)
        color_mapper = LogColorMapper(
            palette=colors,
            low=10,
            high=max(regions_ds.data['count'])
        )
        
        # create plot and add tools
        p = figure(title='Genome sequence distribution',
                   tools='wheel_zoom,pan,reset',
                   x_range=x_range, 
                   y_range=y_range, 
                   toolbar_location='below',
                   output_backend="svg",
                   plot_width=1000, 
                   plot_height=600,
        )
        p.axis.visible = False
        p.add_tile(tile_source)
        p.toolbar.logo = None

        p.patches('lons', 'lats',
                  fill_color={'field': 'count', 'transform': color_mapper},
                  fill_alpha=0.3,
                  line_color="black", 
                  line_alpha=0.1,
                  line_width=1,
                  name='region_count',
                  source=regions_ds)
        
        #all lineages
        X = regions_lineage_ds.data['lineage_type']
        Y = regions_lineage_ds.data['prop']
        Z = [x for _,x in sorted(zip(Y,X))]
        lineages = list(set(Z))

        logging.debug(f'Target lineages {lineages}...')

        palette=list(brewer["Spectral"][len(lineages)-1])

        if "Others" in lineages:
            lineages.remove("Others")
            lineages += ["Others"]
            palette += ["#BBBBBB"]

        # create annular_wedge to indicate the proportion of lineages
        p.annular_wedge(
            x='x', y='y',
            inner_radius=0,
            outer_radius=11,
            direction="anticlock",
            start_angle=cumsum('angle', include_zero=True), 
            end_angle=cumsum('angle'),
            line_color=None,
            fill_color=factor_cmap('lineage_type', palette=palette, factors=lineages),
            line_alpha=0.7,
            alpha=0.7,
            inner_radius_units='screen',
            outer_radius_units='screen',
            name='period_prop',
            legend_field='lineage_type',
            source=regions_lineage_ds)

        hover1 = HoverTool(tooltips=[("Region", "@name"),
                                     ("Total (Last 90 days)", "@count{0,0}")], names = ['region_count'])
        hover2 = HoverTool(tooltips=[("Region", "@name"),
                                     ("Trends", "Last 90 days"),
                                     ("Lineage", "@lineage_type"),
                                     ("Count (prop)", "@count{0,0} (@prop{0.1%})")], names = ['period_prop'])

        hover1.point_policy = "follow_mouse"
        hover2.point_policy = "follow_mouse"
        p.add_tools(hover1)
        p.add_tools(hover2)
        p.legend.location = "top_left"

        return p
        
        
    
    def plot_protein_hotspot(self, ds):
        logging.info(f'Plotting protein hotspots...')

        """ The protein hotspot plot """
        colors = ["#BBBBBB"]+list(Spectral11)
        mapper = LogColorMapper(
            palette = colors,
            low = ds.ds_protein_hotspot.data['acc'].min(),
            high = ds.ds_protein_hotspot.data['acc'].max()
        )

        TOOLS = "hover,save,pan,xpan,box_zoom,reset,wheel_zoom,xwheel_zoom"
        genes = list(set(ds.ds_protein_hotspot.data['gene']))
        
        p_gene = figure(
            title="Protein hotspot overview",
            x_range=(1, 650),
            y_range=genes,
            x_axis_location="above",
            plot_width=self.main_plot_width,
            plot_height=100+len(genes)*20,
            active_drag = 'xpan',
            tools=TOOLS,
            toolbar_location='above',
            output_backend="svg",
            tooltips=[('Gene', '@gene'),
                      ('Position', '@pos{0,0}'),
                      ('Count', '@acc{0,0}')]
        )

        p_gene.grid.grid_line_color = None
        p_gene.axis.axis_line_color = None
        p_gene.axis.major_tick_line_color = None
        p_gene.axis.major_label_text_font_size = self.major_label_text_font_size
        p_gene.axis.major_label_standoff = 0
        p_gene.xaxis.major_label_orientation = pi / 2
        p_gene.toolbar.logo = None
        p_gene.toolbar.autohide = True

        p_gene.circle(
            x='pos',
            y='gene',
            width=1,
            color={'field': 'acc', 'transform': mapper},
            alpha=0.4,
            line_color=None,
            size='size',
            source=ds.ds_protein_hotspot
        )

        p_gene_select = figure(
            x_range=(1, ds.ds_protein_hotspot.data['pos'].max()),
            y_range=(0, 10),
            x_axis_location="below",
            plot_width=self.main_plot_width,
            plot_height=50,
            active_drag = 'xpan',
            tools=TOOLS,
            output_backend="svg",
            toolbar_location=None
        )

        range_tool = RangeTool(x_range=p_gene.x_range)
        range_tool.overlay.fill_color = "navy"
        range_tool.overlay.fill_alpha = 0.2

        p_gene_select.grid.grid_line_color = None
        p_gene_select.ygrid.grid_line_color = None
        p_gene_select.axis.axis_line_color = None

        p_gene_select.yaxis.visible = False

        p_gene_select.add_tools(range_tool)
        p_gene_select.toolbar.active_multi = range_tool
        p_gene_select.toolbar.autohide = True

        self.protein_hotspot_plot = column(p_gene, p_gene_select)
        
    
    
    def plot_spike_substitution_tracking(self, ds):
        logging.info(f'Plotting spike substitution tracking...')

        #################################
        # Main plot
        #################################
        colors = ["#BBBBBB"]+list(Spectral11)
        mapper = LogColorMapper(
               palette=colors,
               low=ds.ds_spike_substitution.data['acc'].min(),
               high=ds.ds_spike_substitution.data['acc'].max()
        )

        TOOLS = "hover,save,pan,ypan,box_zoom,reset,wheel_zoom,ywheel_zoom"

        p = figure(
               x_range=ds.week_ticks,
               y_range=(450,650),
               x_axis_location="above",
               plot_width=self.main_plot_width,
               plot_height=700,
               active_drag = 'ypan',
               tools=TOOLS,
               toolbar_location='below',
               output_backend="svg",
               tooltips=[('Week', '@week'),
                         ('Position', '@pos{0,0}'),
                         ('Count', '@acc{0,0}'),
                         ('Proportion', '@prop{0.00%}')]
        )

        p.grid.grid_line_color = None
        p.axis.axis_line_color = None
        p.axis.major_tick_line_color = None
        p.axis.major_label_text_font_size = self.major_label_text_font_size
        p.axis.major_label_standoff = 0
        p.xaxis.major_label_orientation = pi / 2
        p.yaxis.major_label_text_font_size = '0pt'
        p.yaxis.major_tick_line_color = None # turn off y-axis major ticks
        p.yaxis.minor_tick_line_color = None # turn off y-axis minor ticks
        p.toolbar.logo = None
        p.toolbar.autohide = True
        p.yaxis.visible = False

        p.circle(x="week",
                 y="pos",
                 width=1,
                 color={'field': 'acc', 'transform': mapper},
                 alpha=0.4,
                 line_color=None,
                 size='size',
                 source=ds.ds_spike_substitution
        )


        #################################
        # TOP plot
        #################################

        p_top = figure(
               title="Spike protein mutations",
               plot_width=p.plot_width,
               plot_height=100,
               x_range=p.x_range,
               x_axis_type=None,
               tools="",
               toolbar_location=None, 
               output_backend="svg",
               tooltips=[('Week', '@week'),
                         ('Total genomes', '@total{0,0}')]
        )
        p_top.axis.major_label_text_font_size = "0px"
        p_top.axis.axis_line_color = None
        p_top.yaxis.visible = False
        p_top.grid.grid_line_color = None

        #colors = ["#75968f", "#a5bab7", "#c9d9d3", "#e2e2e2", "#dfccce", "#ddb7b1", "#cc7878", "#933b41", "#550b1d"]
        colors = list(Blues9)[1:7]
        bar_mapper = LinearColorMapper(
               palette=np.flip(colors),
               low=1,
               high=ds.ds_week_total.data['total'].max()
        )

        p_top.vbar(
               x='week',
               top='total',
               width=0.8,
               fill_color={'field': 'total', 'transform': bar_mapper},
               line_color=None,
               source=ds.ds_week_total
        )



        #################################
        # Right plot
        #################################

        p_right = figure(
               plot_width=100,
               plot_height=p.plot_height,
               y_range=(1, ds.ds_spike_substitution.data['pos'].max()),
               x_range=(10, ds.ds_pos_total.data['acc'].max()*10),
               x_axis_location="above",
               tools="",
               toolbar_location=None, 
               y_axis_location="right",
               x_axis_type="log",
               output_backend="svg",
               tooltips=[
                      ('Position', '@pos{0,0}'),
                      ('Count', '@acc{0,0}')]
        )
        p_right.axis.major_label_text_font_size = self.major_label_text_font_size
        p_right.grid.grid_line_color = None
        p_right.xaxis.minor_tick_line_color = None
        p_right.xaxis.major_label_orientation = pi / 3

        colors = list(Spectral11)
        circle_color_mapper = LogColorMapper(
            palette=colors,
            low=10,
            high=ds.ds_pos_total.data['acc'].max()
        )

        p_right.circle(
            'acc',
            'pos',
            alpha=0.6,
            size=7,
            color={'field': 'acc', 'transform': circle_color_mapper},
            source=ds.ds_pos_total
        )

        range_tool = RangeTool(y_range=p.y_range)
        range_tool.overlay.fill_color = "navy"
        range_tool.overlay.fill_alpha = 0.2

        p_right.ygrid.grid_line_color = None
        p_right.add_tools(range_tool)
        p_right.toolbar.active_multi = range_tool
        p_right.toolbar.autohide = True

        #################################
        # Left plot - lineage annotation
        #################################

        p_anno = figure(
            plot_width=150,
            plot_height=p.plot_height,
            x_range=np.sort(ds.data.df_v_info.lineage)[::-1],
            x_axis_location='above',
            # y_axis_location="right",
            y_range=p.y_range,
            output_backend="svg",
            tools="hover",
            tooltips=[('Lineage', '@lineage'),
                      ('Mutation', '@mutation'),
                      ('Type', '@type'),
                      ('Not in all genomes', '@not_all')]
        )

        p_anno.toolbar.logo = None
        p_anno.toolbar_location = None
        p_anno.toolbar.autohide = True

        p_anno.axis.major_label_standoff = 0
        p_anno.axis.axis_line_color = None

        p_anno.xaxis.major_label_text_font_size = '7pt'
        p_anno.xgrid.grid_line_alpha = 0.5
        p_anno.xgrid.grid_line_dash = [6, 4]
        p_anno.xaxis.major_label_orientation = pi/2
        p_anno.xaxis.major_tick_line_color = None

        p_anno.yaxis.axis_label = 'Spike protein position (aa)'
        p_anno.ygrid.grid_line_color = None
        p_anno.yaxis.major_label_text_font_size = '7pt'

        def add_anno(x):
            boxanno = BoxAnnotation(bottom=x.start, top=x.end, line_width=1, fill_alpha=0.05, fill_color='navy')
            p_anno.add_layout(boxanno)

            label = Label(x=0.5, y=x.end-10, text=x['name'], text_color='grey', text_font_size='12px')
            p_anno.add_layout(label)

        ds.data.df_v_anno.apply(add_anno, axis=1)

        p_anno.circle(
            y='pos',
            x='lineage',
            width=1,
            color='color',
            alpha=0.8,
            line_color='grey',
            size=7,
            source=ds.ds_variant
        )

        #############################

        s_protein_layout = gridplot([
            [None, p_top, None], 
            [p_anno, p, p_right],
        ], merge_tools=False)
        
        self.spike_substitution_tracking_plot_top = p_top
        self.spike_substitution_tracking_plot_left = p_anno
        self.spike_substitution_tracking_plot_main = p
        self.spike_substitution_tracking_plot_right = p_right
        self.spike_substitution_tracking_plot_layout = s_protein_layout
        
    def plot_trends(self, ds):
        logging.info(f'Plotting trends...')

        def plot_trend(ds, target_list, title):
            # figure of trends
            TOOLS = "save,pan,ypan,box_zoom,reset,wheel_zoom,ywheel_zoom"
            
            p_trending = figure(
                title=title,
                plot_width=int(self.main_plot_width/2),
                plot_height=600,
                x_range=ds.data['week'],
                toolbar_location='below', 
                tools=TOOLS,
                output_backend="svg",
                # y_axis_label = f'Proportion increased since {ds.data["week"][0]}'
            )

            # grey_box = BoxAnnotation(left=len(ds.data['week'])-1, fill_alpha=0.05, line_color='grey', fill_color='grey')
            # p_trending.add_layout(grey_box)
            
            colors = list(RdYlBu10)[::-1]
            
            for target, color in zip(target_list, colors):
                p_trending.circle(
                    'week',
                    target,
                    size=10, 
                    line_width=0,
                    fill_alpha=0.8,
                    fill_color=color,
                    name=target,
                    legend_label=target,
                    source=ds.data
                )

                p_trending.line(
                    'week',
                    target,
                    legend_label=target,
                    line_color=color,
                    line_width=1,
                    source=ds.data
                )

            hover = HoverTool(
                tooltips = [
                    ('Week', '@week'),
                    ('Target', '$name'),
                    ('Prop increased', '@$name{0.00%}')
                ], names = target_list)

            hover.point_policy = "follow_mouse"
            p_trending.add_tools(hover)

            p_trending.legend.location = "top_left"
            p_trending.legend.click_policy="hide"
            p_trending.toolbar.logo = None
            p_trending.toolbar.autohide = False
            p_trending.xaxis.major_label_orientation = pi / 3
            p_trending.xgrid.grid_line_color = None
            p_trending.ygrid.grid_line_color = None

            return p_trending

        week_since = ds.ds_trend_lineage.data['week'][0]
        p_trending_lin = plot_trend(ds.ds_trend_lineage, ds.trend_lineage_10, f"Variant/lineage proportion changes since {week_since}")
        p_trending_mu = plot_trend(ds.ds_trend_position, ds.trend_position_10, f"S protein position of mutation proportion changes since {week_since}")
                
        self.trending_plot = row(p_trending_lin, p_trending_mu)
    
    
    def plot_lineage_mutation_tracking(self, ds, target_variant=None, target_lineage=None):
        logging.info(f'Plotting lineage mutation tracking...')
        
        def lineage_var_plot(lineage_name, df_lineage, mut_list, lineage_list, title, display_type='count', target_lineage=None):

            ds_week = ds.ds_week_variant_freq if display_type == 'prop' else ds.ds_week_variant
            if target_lineage:
                ds_week = ds.ds_week_lineage_freq if display_type == 'prop' else ds.ds_week_lineage
                lineage_name = target_lineage

            # display raw count or proportion
            display_col = 'prop' if display_type == 'prop' else 'count'

            # prep color mapper
            colors = ["#BBBBBB"]+list(Spectral11)
            #mapper = LinearColorMapper(
            mapper = LogColorMapper(
                   palette=colors,
                   low=df_lineage[display_col].min(),
                   high=df_lineage[display_col].max()
            )

            # display VOC/VOI count on the top        
            p_count = figure(
                   title=title,
                   plot_width=1000,
                   plot_height=100,
                   x_range=ds.week_ticks,
                   x_axis_type=None,
                   output_backend="svg",
                   tooltips=[('Week', '@week'),
                             (display_type, f'@{lineage_name}')]
            )

            p_count.toolbar.logo = None
            p_count.toolbar_location = None
            p_count.toolbar.autohide = True
            p_count.ygrid.grid_line_alpha = 0.5
            p_count.ygrid.grid_line_dash = [6, 4]
            p_count.axis.axis_line_color = None
            p_count.axis.major_label_text_font_size = "10px"
            p_count.axis.major_tick_line_color = None
            p_count.axis.major_label_standoff = 0
            p_count.yaxis.minor_tick_line_color = None
            # p_count.yaxis.axis_label = 'Count (log)'
                        
            p_count.vbar(
               x='week',
               top=lineage_name,
               width=0.8,
               fill_color={'field': lineage_name, 'transform': mapper},
               line_color=None,
               source=ds_week
            )

            # week-mutation figure
            TOOLS = "hover,save,pan,xpan,box_zoom,reset,wheel_zoom,xwheel_zoom"
            
            p_aa = figure(
               x_range=p_count.x_range,
               y_range=mut_list,
               plot_width=1000,
               plot_height=90+len(mut_list)*12,
               tools=TOOLS,
               output_backend="svg",
               tooltips=[('Week', '@week'),
                         ('Mutation', '@mutation'),
                         ('Count', '@count{0,0}'),
                         ('Propotion', '@prop{0.00%}')]
            )

            p_aa.toolbar.logo = None
            p_aa.toolbar_location = 'below'
            p_aa.toolbar.autohide = True
            p_aa.grid.grid_line_color = None
            p_aa.axis.axis_line_color = None
            p_aa.axis.major_tick_line_color = None
            p_aa.axis.major_label_text_font_size = "10px"
            p_aa.axis.major_label_standoff = 0
            p_aa.xaxis.major_label_orientation = pi / 3

            p_aa.rect(
                x="week",
                y="mutation",
                width=1,
                height=1,
                alpha=0.9,
                fill_color={'field': display_col, 'transform': mapper},
                line_color=None,
                source=df_lineage,
            )

            # variant mutation figure
            p_variant = figure(
                plot_width=len(lineage_list)*12,
                plot_height=p_aa.plot_height,
                x_range=lineage_list,
                y_range=p_aa.y_range,
                output_backend="svg",
                tools="hover",
                x_axis_location="below",
                tooltips=[('Lineage', '@lineage'),
                          ('Mutation', '@mutation'),
                          ('Type', '@type'),
                          ('Not in all genomes', '@not_all')]
            )

            p_variant.circle(
                y='mutation',
                x='lineage',
                width=1,
                color='color',
                alpha=0.8,
                line_color='grey',
                size=7,
                source=ds.ds_variant
            )

            p_variant.toolbar.logo = None
            p_variant.toolbar_location = None
            p_variant.toolbar.autohide = True

            p_variant.xaxis.major_label_text_font_size = '7pt'
            p_variant.xgrid.grid_line_alpha = 0.5
            p_variant.xgrid.grid_line_dash = [6, 4]
            p_variant.xaxis.major_label_orientation = pi/2
            p_variant.xaxis.major_tick_line_color = None
            p_variant.axis.major_label_standoff = 0
            p_variant.axis.axis_line_color = None
            p_variant.ygrid.grid_line_color = None
            p_variant.yaxis.visible = False

            return (p_count, p_aa, p_variant)


        variant_list = ds.data.df_v_info.lineage
        target_variants = variant_list
        if target_variant:
            target_variants = [target_variant] # only run the target variant
        lineage_plots_raw = []
        lineage_plots_prop = []

        for variant_name in target_variants:

            df_lineage, mut_list = ds.prepare_mutation_tracking_per_variant(variant_name)

            # lineage info
            v_info = ds.data.df_v_info[ds.data.df_v_info.lineage==variant_name]
            nextstrain_clade = v_info.iloc[0,1]
            v_type = v_info.iloc[0,2]
            first_detected = v_info.iloc[0,5]
            
            title = f'{variant_name} ({nextstrain_clade}) {v_type} first detected in {first_detected}'
            
            # calculating raw counts
            (p_count, p_aa, p_voc) = lineage_var_plot(variant_name, df_lineage, mut_list, variant_list, title, 'count', target_lineage)
            lineage_plots_raw.append(
                gridplot([
                    [p_count, None],
                    [p_aa, p_voc]
                ], merge_tools=False)
            )

            # calculating propotion charts
            (p_count, p_aa, p_voc) = lineage_var_plot(variant_name, df_lineage, mut_list, variant_list, title,'prop', target_lineage)
            lineage_plots_prop.append(
                gridplot([
                    [p_count, None],
                    [p_aa, p_voc]
                ], merge_tools=False)
            )

        tab1 = Panel(child=column(lineage_plots_raw), title="Count")
        tab2 = Panel(child=column(lineage_plots_prop), title="Proportion")

        lineage_aa_tracking = Tabs(tabs=[tab1, tab2])

        self.lineage_mutation_tracking_plot = lineage_aa_tracking



class CovidViz(object):
    def __init__(self, CovidPlots, title="SARS-CoV-2 Genome Tracking"):
        self.title = title
        self.plots = CovidPlots
        self.output_layout = None

        # web components
        self.padding_div = None
        self.header = None
        self.geo_p = None
        self.week_lineage_bar_t = None
        self.week_lineage_bar_p = None
        self.trending_t = None
        self.trending_p = None
        self.spike_substitution_tracking_t = None
        self.spike_substitution_tracking_p = None
        self.lineage_mutation_tracking_t = None
        self.lineage_mutation_tracking_p = None
        self.footer = None

        self._prep_web_component()

        self.template = """
{% block postamble %}
<link rel="stylesheet" href="https://stackpath.bootstrapcdn.com/bootstrap/4.3.1/css/bootstrap.min.css">
<style type="text/css">
@import url('https://fonts.googleapis.com/css?family=Quicksand');
H1, H2, H3 {
    margin-top: 1em;
    font-family: Quicksand;
}
H1 {
    font-size: 2em;
}
H2 {
    font-size: 1.5em;
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
        

    def _prep_web_component(self):
        from bokeh.models import Div, Paragraph
        from datetime import date

        update_date = date.today().strftime("%Y-%m-%d")

        # padding_div
        self.padding_div = Div(text=' ', width=150)

        # geo_country_t
        self.header = Div(text=f"""
        <H1>{self.title}</H1>
        <div>
        Last update: {update_date} |
        Enabled by data from 
            <a href="https://www.gisaid.org/" target="_blank"><img src="http://gisaid.org/fileadmin/gisaid/img/schild.png" style="width: 71px; height: 25px; vertical-align: middle;" alt="GISAID Initiative" class="ml-2"></a>
        </div>
        """, width=self.plots.main_plot_width)

        # geo_country_p
        self.geo_p = Paragraph(text="""
        Multiple variants of SARS-CoV-2, the virus that causes COVID-19, have emerged or spread throughout different parts of the world, including the United States. Variant viruses may carry mutations that could be associated with differences in diagnostic test performance, changes in disease epidemiology, clinical outcomes, and effectiveness of certain therapeutics or vaccines. The genome information and metadata are downloaded from the Global Initiative on Sharing Avian Influenza Data (GISAID).
        """, width=self.plots.main_plot_width)

        
        # week_lineage_bar_t
        # week_lineage_bar_p

        self.week_lineage_bar_t = Div(text="""<H2>SARS-CoV-2 Variants</H2>""", width=self.plots.main_plot_width)
        self.week_lineage_bar_p = Paragraph(text="""
        Phylogenetic Assignment of Named Global Outbreak (PANGO) Lineages is software tool developed by members of the Rambaut Lab. The associated web application was developed by the Centre for Genomic Pathogen Surveillance in South Cambridgeshire and is intended to implement the dynamic nomenclature of SARS-CoV-2 lineages, known as the PANGO nomenclature. Each classification of variant includes the possible attributes of lower classes (e.g., VOC includes the possible attributes of VOI); variant status might escalate or deescalate based on scientific evidence.
        The plots below show the sequence count and the propotion for the top 10 and all VOC/VOI variants reported by CDC.
        """, width=self.plots.main_plot_width)

        # trending_t
        # trending_p

        self.trending_t = Div(text="""<H2>Trending Lineages and Mutations</H2>""", width=self.plots.main_plot_width)
        self.trending_p = Paragraph(text="""
        The trending SASRS-CoV-2 variants are estimated from weekly growth rates over the past 8 weeks.
        """, width=self.plots.main_plot_width)

        # spike_substitution_tracking_t
        # spike_substitution_tracking_p

        self.spike_substitution_tracking_t = Div(text="""<H2>Protein Substitution Tracking</H2>""", width=self.plots.main_plot_width)
        self.spike_substitution_tracking_p = Paragraph(text="""
        Notable variants of SARS-CoV-2 and notable mutations found.
        D614G: Logarithmic Prevalence of D614G in 2020 according to sequences in the GISAID database D614G is a missense mutation that affects the spike protein of SARS-CoV-2. The frequency of this mutation in the viral population has increased during the pandemic. G (glycine) has replaced D (aspartic acid) at position 614 in many countries, especially in Europe though more slowly in China and the rest of East Asia, supporting the hypothesis that G increases the transmission rate, which is consistent with higher viral titers and infectivity in vitro.
        E484K: The name of the mutation, E484K, refers to an exchange whereby the glutamic acid (E) is replaced by lysine (K) at position 484. 484K has been reported to be an escape mutation (i.e., a mutation that improves a virus's ability to evade the host's immune system[138][139]) from at least one form of monoclonal antibody against SARS-CoV-2, indicating there may be a "possible change in antigenicity".
        N501Y denotes a change from asparagine (N) to tyrosine (Y) in amino-acid position 501. Variants with N501Y include P.1 (Brazil/Japan), Variant of Concern 20DEC-01 (UK), 501.V2 (South Africa), and COH.20G/501Y (Columbus, Ohio). This last became the dominant form of the virus in Columbus in late December 2020 and January and appears to have evolved independently of other variants.
        S477G/N: A highly flexible region in the receptor binding domain (RBD) of SARS-CoV-2, starting from residue 475 and continuing up to residue 485, was identified using bioinformatics and statistical methods in several studies.
        P681H: Logarithmic Prevalence of P681H in 2020 according to sequences in the GISAID database In January 2021, scientists reported in a preprint that the mutation 'P681H', a characteristic feature of the significant novel SARS-CoV-2 variants detected in the U.K. (B.1.1.7) and Nigeria (B.1.1.207), is showing a significant exponential increase in worldwide frequency, similar to the now globally prevalent 'D614G'.
        """, width=self.plots.main_plot_width)

        # lineage_mutation_tracking_t
        # lineage_mutation_tracking_p
        
        self.lineage_mutation_tracking_t = Div(text="""<H2>Characteristics of SARS-CoV-2 Variants</H2>""", width=self.plots.main_plot_width)
        self.lineage_mutation_tracking_p = Paragraph(text="""
        Variant of Concern is a variant for which there is evidence of an increase in transmissibility, more severe disease (e.g.,  increased hospitalizations or deaths), significant reduction in neutralization by antibodies generated during previous infection or vaccination, reduced effectiveness of treatments or vaccines, or diagnostic detection failures. Variant of Interest is a variant with specific genetic markers that have been associated with changes to receptor binding, reduced neutralization by antibodies generated against previous infection or vaccination, reduced efficacy of treatments, potential diagnostic impact, or predicted increase in transmissibility or disease severity.
        """, width=self.plots.main_plot_width)

        # footer
        
        self.footer = Div(text="""
        <div>
        B-10, Bioscience Div, Los Alamos National Laboratory ©2021
        </div>
        """, width=self.plots.main_plot_width, align='center', style={'width':'100%', 'text-align':'center', 'margin-top': '1em'})

                
    def layout(self):
        logging.info(f'Layouting plots...')
        from bokeh.layouts import column, layout, gridplot

        layout = gridplot([
            [None, self.header, None],
            [None, self.geo_p, None],
            [None, self.plots.geo_plot, None],
            [None, column(self.week_lineage_bar_t, self.week_lineage_bar_p), None],
            [None, self.plots.week_lineage_bar_plot, None],
            [None, column(self.trending_t, self.trending_p), None],
            [None, self.plots.trending_plot, None],
            [None, column(self.spike_substitution_tracking_t, self.spike_substitution_tracking_p), None],
            [None, self.plots.protein_hotspot_plot, None],
            [None, self.plots.spike_substitution_tracking_plot_top, None],
            [self.plots.spike_substitution_tracking_plot_left, self.plots.spike_substitution_tracking_plot_main, self.plots.spike_substitution_tracking_plot_right],
        ], merge_tools=False)

        layout_lin = gridplot([
            [self.padding_div, column(self.lineage_mutation_tracking_t, self.lineage_mutation_tracking_p)],
            [None, self.plots.lineage_mutation_tracking_plot]
        ], merge_tools=False)
        
        self.output_layout = column(
            layout,
            layout_lin, 
            self.footer)
        
    def save_html(self, outfile='variant_stats_page.html', template=None):
        logging.info(f'Save HTML file to {outfile}...')

        from bokeh.plotting import save
        from bokeh.io import output_file, reset_output
        reset_output()
        output_file(outfile, title=self.title)
        if template is None:
            template = self.template
        save(self.output_layout, template=template)
        
class StatesRegionsLookup(object):
    """
    This is the class for conversion of US states and federal regions
    """
    def __init__(self):
        self.region_states_lookup = {
            "I": ["Connecticut", "Maine", "Massachusetts", "New Hampshire", "Rhode Island", "Vermont"],
            "II": ["New Jersey", "New York", "Puerto Rico", "US Virgin Islands"],
            "III": ["Delaware", "District of Columbia", "Maryland", "Pennsylvania", "Virginia", "West Virginia"],
            "IV": ["Alabama", "Florida", "Georgia", "Kentucky", "Mississippi", "North Carolina", "South Carolina", "Tennessee"],
            "V": ["Illinois", "Indiana", "Michigan", "Minnesota", "Ohio", "Wisconsin"],
            "VI": ["Arkansas", "Louisiana", "New Mexico", "Oklahoma", "Texas"],
            "VII": ["Iowa", "Kansas", "Missouri", "Nebraska"],
            "VIII": ["Colorado", "Montana", "North Dakota", "South Dakota", "Utah", "Wyoming"],
            "IX": ["Arizona", "California", "Hawaii", "Nevada", "Guam"],
            "X": ["Alaska", "Idaho", "Oregon", "Washington"],  
        }
        self.state_abbv_lookup = {
            'Alabama': 'AL', 'Alaska': 'AK', 'Arizona': 'AZ', 'Arkansas': 'AR', 'California': 'CA', 'Canal Zone': 'CZ',
            'Colorado': 'CO', 'Connecticut': 'CT', 'Delaware': 'DE', 'District of Columbia': 'DC', 'Florida': 'FL', 'Georgia': 'GA',
            'Guam': 'GU', 'Hawaii': 'HI', 'Idaho': 'ID', 'Illinois': 'IL', 'Indiana': 'IN', 'Iowa': 'IA', 'Kansas': 'KS',
            'Kentucky': 'KY', 'Louisiana': 'LA', 'Maine': 'ME', 'Maryland': 'MD', 'Massachusetts': 'MA', 'Michigan': 'MI',
            'Minnesota': 'MN', 'Mississippi': 'MS', 'Missouri': 'MO', 'Montana': 'MT', 'Nebraska': 'NE', 'Nevada': 'NV',
            'New Hampshire': 'NH', 'New Jersey': 'NJ', 'New Mexico': 'NM', 'New York': 'NY', 'North Carolina': 'NC', 'North Dakota': 'ND',
            'Ohio': 'OH', 'Oklahoma': 'OK', 'Oregon': 'OR', 'Pennsylvania': 'PA', 'Puerto Rico': 'PR', 'Rhode Island': 'RI',
            'South Carolina': 'SC', 'South Dakota': 'SD', 'Tennessee': 'TN', 'Texas': 'TX', 'Utah': 'UT', 'Vermont': 'VT',
            'US Virgin Islands': 'VI', 'Virginia': 'VA', 'Washington': 'WA', 'West Virginia': 'WV', 'Wisconsin': 'WI', 'Wyoming': 'WY',
        }

        self.state_region_lookup = {s:r for r in self.region_states_lookup for s in self.region_states_lookup[r]}
        self.region_abbv_lookup = {r: [self.state_abbv_lookup[s] for s in self.region_states_lookup[r]] for r in self.region_states_lookup}
        
    def region_to_states(self, region):
        try:
            return self.region_states_lookup[region]
        except:
            return None
            
    def state_to_region(self, state):
        try:
            return self.state_region_lookup[state]
        except:
            return None

    def get_all_regions(self):
        try:
            return list(self.region_states_lookup.keys())
        except:
            return None