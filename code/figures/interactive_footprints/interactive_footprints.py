import numpy as np 
import pandas as pd 
import bokeh.io
import bokeh.plotting
from bokeh.models import * 
from bokeh.themes import Theme
from bokeh.transform import linear_cmap
import glob

#############################################
# Helper functions. Credit to Griffin Chure #
#############################################

def color_palette():
    """
    Returns a dictionary of the PBOC color palette
    """
    return {'green': '#7AA974', 'light_green': '#BFD598',
              'pale_green': '#DCECCB', 'yellow': '#EAC264',
              'light_yellow': '#F3DAA9', 'pale_yellow': '#FFEDCE',
              'blue': '#738FC1', 'light_blue': '#A9BFE3',
              'pale_blue': '#C9D7EE', 'red': '#D56C55', 'light_red': '#E8B19D',
              'pale_red': '#F1D4C9', 'purple': '#AB85AC',
              'light_purple': '#D4C2D9', 'dark_green':'#7E9D90', 'dark_brown':'#905426'}

def bokeh_theme():
    """A custom bokeh theme to match PBoC 2e colors"""
    theme_json = {'attrs':
            {'figure': {
                'background_fill_color': '#E3E7E9',
                'outline_line_color': '#FFFFFF',
            },
            'Axis': {
            'major_tick_in': 4,
            'major_tick_line_width': 1,
            'axis_label_text_font': 'Lato',
            'axis_label_text_font_style': 'normal'
            },
            'Grid': {
                'grid_line_color': "white",
            },
            'Legend': {
                'background_fill_color': '#E3E7E9',
                'border_line_color': '#FFFFFF',
                'border_line_width': 1.5,
                'background_fill_alpha': 0.5
            },
            'Text': {
                'text_font_style': 'normal',
               'text_font': 'Lato'
            },
            'Title': {
                'background_fill_color': '#FFFBCE',
                'text_font_style': 'normal',
                'align': 'center',
                'text_font': 'Lato',
                'offset': 2,
            }}}

    theme = Theme(json=theme_json)
    bokeh.io.curdoc().theme = theme

    # Define the colors
    colors = color_palette()
    palette = [v for k, v in colors.items() if 'pale' not in k]
    return [colors, palette]



def load_js(fname, args):
    """
    Given external javascript file names and arguments, load a bokeh CustomJS
    object
    
    Parameters
    ----------
    fname: str or list of str
        The file name of the external javascript file. If the desired javascript
        exists in multiple external files, they can be provided as a list of
        strings.
    args: dict
        The arguments to supply to the custom JS callback. 
    
    Returns
    -------
    cb : bokeh CustomJS model object
        Returns a bokeh CustomJS model object with the supplied code and
        arguments. This can be directly assigned as callback functions.
    """
    if type(fname) == str:
        with open(fname) as f:
            js = f.read() 
    elif type(fname) == list:
        js = ''
        for _fname in fname:
            with open(_fname) as f:
                js += f.read()

    cb = CustomJS(code=js, args=args)
    return cb

bokeh_theme()

# Path to store html file
bokeh.io.output_file('interactive_footprints.html')

###############
# Data Import #
###############
collapse_df = pd.DataFrame()


# Go through all files and compact the data
for file in glob.glob("../../analysis/all_data/footprints/*"):
    df = pd.read_csv(file)
    for name, group in df.groupby(['promoter', 'replicate', 'growth_condition', 'd']):
        x = group.pos.values
        y = group.mut_info.values
        collapse_df = pd.concat([collapse_df, pd.DataFrame(data={
            'mut_info':[y], 
            'pos': [x], 
            'promoter': name[0],
            'replicate': name[1],
            'growth_condition': name[2],
            'd': name[3]
            })])


collapse_exshift_df = pd.DataFrame()
for file in glob.glob("../../analysis/all_data/expression_shifts/*"):
    df = pd.read_csv(file)
    for name, group in df.groupby(['promoter', 'replicate', 'growth_condition']):
        pos = group.pos.values 
        base = group.base.values                                   
        wt_base = group.wt_base.values
        expression_shift = group.expression_shift.values
        collapse_exshift_df = pd.concat([collapse_exshift_df, pd.DataFrame(data={
            'pos':[pos], 
            'base': [base], 
            'wt_base': [wt_base],
            'expression_shift': [expression_shift],
            'promoter': name[0],
            'replicate': name[1],
            'growth_condition': name[2],
            })])

# Import metadata for promoters
#df_meta = pd.read_csv('./20230525_footprints_meta.csv')
df_meta = pd.read_csv('./20230907_footprints_meta.csv')
df_regulonDB = pd.read_csv('./regulonDB_meta.csv')

# Transform types to strings
collapse_df['replicate'] = collapse_df['replicate'].astype(str)
collapse_df['d'] = collapse_df['d'].astype(str)
collapse_exshift_df['replicate'] = collapse_exshift_df['replicate'].astype(str)


# 
data = ColumnDataSource(collapse_df)
exshift = ColumnDataSource(collapse_exshift_df)
meta = ColumnDataSource(df_meta)
regulonDB = ColumnDataSource(df_regulonDB)
promoters = list(df['promoter'].unique())

# Set inital settings for plot
prom_ini = 'xylFp'
gc_ini = 'Glucose'
rep_ini = '1'
d_ini = '1'

wt_seq = df_meta.loc[df_meta['promoter'] == prom_ini, "promoter_seq"].values[0]

# populate datasources with initial values
_df = collapse_df.loc[(collapse_df['promoter'] == prom_ini) 
           & (collapse_df['growth_condition'] == gc_ini)
           & (collapse_df['replicate'] == rep_ini)
           & (collapse_df['d'] == d_ini), 
           ['pos', 'mut_info']]



data_display = ColumnDataSource({'pos': _df['pos'].values[0], 
                                 'mut_info': _df['mut_info'].values[0]})
                                 #'footprint': list(_df['footprint'].values), 
                                 #'footprint_test': list(_df['footprint_test'].values)})

_df_exshift = collapse_exshift_df.loc[(collapse_exshift_df['promoter'] == prom_ini) 
                                    & (collapse_exshift_df['growth_condition'] == gc_ini)
                                    & (collapse_exshift_df['replicate'] == rep_ini),
                                        ['pos', 'base', 'wt_base', 'expression_shift']]

exshift_display = ColumnDataSource({'pos': _df_exshift['pos'].values[0], 
                                   'base': _df_exshift['base'].values[0],
                                    'wt_base': _df_exshift['wt_base'].values[0],
                                    'expression_shift': _df_exshift['expression_shift'].values[0]})

# Define the selections
prom_selector = Select(options=list(np.sort(promoters)), value=prom_ini)
gc_selector = Select(options=list(np.sort(collapse_df['growth_condition'].unique())), value=gc_ini)
rep_selector = Select(options=list(collapse_df['replicate'].unique()), value=rep_ini)
d_selector = Select(options=list(collapse_df['d'].unique()), value=d_ini)


# titles for selectors
prom_title = Div(text="<b>Promoter</b>")
gc_title = Div(text="<b>Growth Condition</b>")
rep_title = Div(text="<b>Replicate</b>")
d_title = Div(text="<b>Window Width</b>")

# metadata for default choice
meta_ini = df_meta.loc[df_meta['promoter'] == prom_ini, :]

# boxes for description
prom_desc = Div(text='<div style="width:300px; overflow-wrap: break-word;"><b> Genes controlled by promoter</b>: <br/>' + meta_ini['genes'].values[0] + '<br/><b>Strand: </b><br/>' + meta_ini['direction'].values[0] + '<br/><b>5\':</b><br/>' + str(meta_ini['five_prime'].values[0]) + '<br/><b>3\':</b><br/>' + str(meta_ini['three_prime'].values[0]) + '</div>')
regulonDB_desc = Div(text="")

def update_sites(attr, old, new):
    x = '<div style="width:700px;"><b> Annotation in RegulonDB</b><br/>'
    promoter = prom_selector.value
    regulons = df_regulonDB.loc[df_regulonDB['PROMOTER_NAME'] == promoter, :]
    if len(regulons) == 0:
        x += '<br/>No Binding Sites Found'
    else:
        for index, site in regulons.iterrows():
            if -115 < site['CENTER_POSITION'] < 45:
                x += '<div style="overflow-wrap: break-word;"><br/><b>' + site['RI_FUNCTION'] + '</b><br/>Transcription Factor: ' + site['TRANSCRIPTION_FACTOR_NAME'] + '<br/>Binding Site Position Relative to TSS: ' + str(site['CENTER_POSITION']) + '</div>'#+ '<br/> Binding Site Sequence (Capital Letters): ' + site['RI_SEQUENCE'] + '<br/> Consensus Sequence: ' + site['CONSENSUS_SEQUENCE'] + '</div>';
    x += '</div>'
    regulonDB_desc.update(text=x)

update_sites("", "", "")

# initiate plot windows
TOOLS = "hover,save,pan,box_zoom,reset,wheel_zoom"

p_exshift = bokeh.plotting.figure(width=1000, height=200, 
                                  x_axis_label='sequence',
                                  title="Expression Shift upon mutation",
                                  x_range=[-0.5-115, 45 - 0.5], 
                                  y_range=[0.5, 5 - 0.5],
                                  tooltips=[('wild type base', '@wt_base')],
                                  tools=TOOLS)

p_info = bokeh.plotting.figure(width=1000, height=200, 
                               x_axis_label='position',
                               y_axis_label='mutual information [bits]',
                               title="Mutual Information from Data",)
                               #x_range=p_exshift.x_range)

#p_footprint = bokeh.plotting.figure(width=1000, height=200, 
#                                    x_axis_label='position',
#                                    y_axis_label='relative entropy [bits]',
#                                    title="Relative Entropy from MCMC",
#                                    x_range=p_exshift.x_range)



#p_footprint_test = bokeh.plotting.figure(width=1000, height=200#, 
#                                         x_axis_label='position',
 #                                        y_axis_label='relative entropy [bits]',
#                                         title="Relative Entropy from MCMC (Expression Shift as initial #condition)",
 #                                        x_range=p_exshift.x_range,
#                                         )


# update axis to show letter
p_exshift.yaxis.ticker = np.arange(1,5)
p_exshift.yaxis.major_label_overrides = {(tick+1): x_ for tick, x_ in enumerate(['A', 'C', 'G', 'T'])}

p_exshift.xaxis.major_label_overrides = {(tick-115): x_ for tick, x_ in enumerate(exshift_display.data['wt_base'][0::4])}
p_exshift.xaxis.major_label_text_font_size = "6pt"

# populate plots
p_info.vbar(x='pos', top='mut_info', source=data_display)
#p_footprint.vbar(x='pos', top='footprint', source=data_display)
#p_footprint_test.vbar(x='pos', top='footprint_test', source=data_display)


r = p_exshift.rect(x='pos', 
               y='base',
               width=1,
               height=1,
               fill_color=linear_cmap('expression_shift', 
                                      bokeh.palettes.interp_palette(["#D14241", "#FFFFFF", "#738FC1"], 100),
                                      low=-1, 
                                      high=1),
               line_color=None,
               source=exshift_display
)
p_exshift.xaxis.ticker = np.arange(-115, 45)

p_exshift.extra_x_ranges['x_above'] = Range1d(-115, 45)
#p_exshift.extra_x_ranges.update({'x_above':  p_exshift.x_range})
p_exshift.add_layout(LinearAxis(x_range_name='x_above', ticker=np.arange(-11, 5) * 10), 'above')

color_bar = r.construct_color_bar(padding=5)
p_exshift.add_layout(color_bar, "right")

#p_footprint.xaxis.ticker = np.arange(-11, 5) * 10
p_info.xaxis.ticker = np.arange(-11, 5) * 10
#p_footprint_test.xaxis.ticker = np.arange(-11, 5) * 10


# Define the callbacks
args = {
    'data_display': data_display,
    'exshift_display': exshift_display,
    'data': data,
    'exshift': exshift,
    'prom_selector': prom_selector,
    'gc_selector': gc_selector,
    'rep_selector': rep_selector,
    'd_selector': d_selector,
    'prom_desc': prom_desc,
    'meta': meta,
    'regulonDB_desc': regulonDB_desc,
    'regulonDB': regulonDB,
    'x_axis': p_exshift.xaxis[1],
    'p': p_exshift
}




prom_cb = load_js(['prom_selector.js', 'footprint_selector.js'], args=args)
prom_selector.js_on_change('value', prom_cb)

gc_cb =  load_js(['growth_condition_selector.js', 'footprint_selector.js'], args=args)
gc_selector.js_on_change('value', gc_cb)

footprint_cb = load_js('footprint_selector.js', args=args)
for s in [rep_selector, d_selector]:
    s.js_on_change('value', footprint_cb)


selector_box = bokeh.layouts.row(
    bokeh.layouts.column(
        prom_title, 
        prom_selector,
        gc_title,
        gc_selector, 
    ),
    bokeh.layouts.column( 
        rep_title, 
        rep_selector,
        d_title, 
        d_selector
    ),
    prom_desc
)
plot = bokeh.layouts.column(
    selector_box, 
    p_info, 
    p_exshift, 
    regulonDB_desc)

bokeh.io.save(plot)

# Remove first line from html document
with open(r'interactive_footprints.html', 'r+') as fp:
    # read an store all lines into list
    lines = fp.readlines()
    # move file pointer to the beginning of a file
    fp.seek(0)
    # truncate the file
    fp.truncate()

    # start writing lines except the first line
    # lines[1:] from line 2 to last line
    fp.writelines(lines[1:])


