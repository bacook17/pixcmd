# pcmd_game.py
# execute with: bokeh serve --show pcmd_first.py
import numpy as np
import sys

sys.path.append('/Users/bcook/pCMDs/pixcmd/pysrc')
import instrument as ins, isochrones as iso, galaxy as gal
import driver, fit_model, utils, gpu_utils
from scipy.misc import lena

from bokeh.io import curdoc
from bokeh.plotting import figure
from bokeh.layouts import column, row
from bokeh.palettes import Viridis256, Greys256
from bokeh.models import ColumnDataSource
from bokeh.models.mappers import LinearColorMapper, LogColorMapper
from bokeh.models.widgets import Select, TextInput, Slider, Button, Toggle, RadioButtonGroup

filters = np.array([ins.Filter.HST_F475W(1.0), ins.Filter.HST_F814W(1.0)])
iso_model = iso.Isochrone_Model(filters, path='/Users/bcook/pCMDs/pixcmd/isoc_csv/')
driv = driver.Driver(iso_model, gpu=False)

class Model():

    def __init__(self, logz, logdust, logNpix, logage, binning, N_scale, psf, psf_width):
        self.attributes = {}
        self.attributes['logz'] = logz
        self.attributes['logdust'] = logdust
        self.attributes['logNpix'] = logNpix
        self.attributes['logage'] = logage
        self.attributes['binning'] = binning
        self.attributes['N_scale'] = N_scale
        self.attributes['psf'] = psf
        self.attributes['width'] = psf_width

        self.changed = False

    def get_params(self):
        d = self.attributes
        return np.array([d['logz'], d['logdust'], d['logNpix'], d['logage']])

    def set_attr(self, attrname, value):
        self.attributes[attrname] = value
        self.changed = True

    def simulate(self):
        binning, N_scale = self.attributes['binning'], self.attributes['N_scale']
        psf = self.attributes['psf']
        if psf == 0:
            opts = dict(psf=False)
        elif psf == 1:
            opts = dict(psf=True, convolve_func='gaussian', width=self.attributes['width'])
        elif psf == 2:
            opts = dict(psf=True, convolve_func=None, multi_psf=False)
        else:
            opts = dict(psf=True, convolve_func=None, multi_psf=True)
        _, mags, _, images = driv.simulate(gal.Galaxy_SSP(self.get_params()), N_scale, fixed_seed=False, **opts)
        
        pcmd = utils.make_pcmd(mags)
        im1 = images[0] * 127. / np.mean(images[0])
        im2 = images[1] * 127. / np.mean(images[1])

        xlim = [p2.x_range.start, p2.x_range.end]
        ylim = [p2.y_range.end, p2.y_range.start]
        xnum = int(np.ceil((xlim[1] - xlim[0]) / binning))
        ynum = int(np.ceil((ylim[1] - ylim[0]) / binning))
        
        xbins = np.linspace(xlim[0], xlim[1], xnum)
        ybins = np.linspace(ylim[0], ylim[1], ynum)
        
        hess, _, _ = np.histogram2d(pcmd[0], pcmd[1], bins=np.array([xbins, ybins]))

        hess *= 255. / np.max(hess)
        hess[hess == 0.] = -1.
        return dict(hess=[hess.astype(int)[:, ::-1].T]), dict(im1=[im1.astype(int)], im2=[im2.astype(int)])

random = Model(0., 0., 0., 0., 0.1, 80, 3, 3.)
random.set_attr('logz', np.random.uniform(low=-2.5, high=1.))
random.set_attr('logdust', np.random.uniform(low=-2., high=1.))
random.set_attr('logNpix', np.random.uniform(low=-1., high=8.))
random.set_attr('logage', np.random.uniform(low=6., high=10.5))
    
current = Model(0., -6., 2., 10., 0.1, 80, 3, 3.)

plot_options = dict(width=400, plot_height=400)

p1 = figure(x_range=[0,1], y_range=[0,1], **plot_options)
p1.title.text = 'Data: %s'%(filters[0].name)
p2 = figure(x_range=[-1.5, 4.6], y_range=[15.6, -12.], **plot_options)
p2.title.text = "Data: Pixel CMD"
p2.xaxis.axis_label = '%s - %s'%(filters[0].name, filters[1].name)
p2.yaxis.axis_label = filters[1].name

p3 = figure(x_range=[0,1], y_range=[0,1], **plot_options)
p3.title.text = 'Model: %s'%(filters[0].name)
p4 = figure(x_range=p2.x_range, y_range=p2.x_range, **plot_options)
p4.title.text = "Model: Pixel CMD"
p4.xaxis.axis_label = '%s - %s'%(filters[0].name, filters[1].name)
p4.yaxis.axis_label = filters[1].name

pcmd_data, image_data = random.simulate()

pcmd_dict, image_dict = current.simulate()

pcmd_true = ColumnDataSource(pcmd_data)
images_true = ColumnDataSource(image_data)

pcmd = ColumnDataSource(pcmd_dict)
images = ColumnDataSource(image_dict)

mapper = LinearColorMapper(palette=Greys256, low=0., high=255.)
logmapper = LogColorMapper(palette=Viridis256, low=0., high=255., low_color='white')

i1 = p1.image("im1", x=0, y=0, dw=1, dh=1, source=images_true, color_mapper=mapper)
i2 = p2.image("hess", x=-1.5, y=15.6, dw=6.1, dh=27.6, source=pcmd_true, color_mapper=logmapper)

i3 = p3.image("im1", x=0, y=0, dw=1, dh=1, source=images, color_mapper=mapper)
i4 = p4.image("hess", x=-1.5, y=15.6, dw=6.1, dh=27.6, source=pcmd, color_mapper=logmapper)

slider_z = Slider(start=-2.5, end=1.0, value=current.attributes['logz'], step=0.1, title='logz', name='logz')
slider_dust = Slider(start=-2., end=1., value=current.attributes['logdust'], step=0.1, title='log E(B-V)', name='logdust')
slider_Npix = Slider(start=-1., end=8., value=current.attributes['logNpix'], step=0.1, title='log Npix', name='logNpix')
slider_age = Slider(start=5., end=10.6, value=current.attributes['logage'], step=0.05, title='log Age', name='logage')
slider_bins = Slider(start=0.01, end=0.5, value=current.attributes['binning'], step=0.01, title='Hess Binning', name='binning')
slider_Nscale = Slider(start=80, end=400, value=80, step=80, title='N_scale', name='N_scale')
slider_psf = Slider(start=0.5, end=10., value=1., step=0.5, title='PSF Width', name='width')

sliders = [slider_z, slider_dust, slider_Npix, slider_age, slider_bins, slider_Nscale, slider_psf]

button = Button(label='Resample', button_type="success")
psf_button = RadioButtonGroup(labels=['No PSF', 'Gaussian', 'HST Simple', 'HST Complex'], active=3)

def update_model(attrname, old, new):
    for i in sliders:
        current.set_attr(i.name, i.value)
    current.set_attr("psf", psf_button.active)
    #current.set_attr('N_scale', int(select_Nscale.value))

def update_width(attrname, old, new):
    if (current.attributes['psf'] == 1):
        update_model(None, None, None)
    
def update_psf(on):
    update_model(None, None, None)
    
def update_if_change():
    if current.changed:
        update()

def update():
    current.changed = False
    
    i4.glyph.x = p4.x_range.start
    i4.glyph.dw = (p4.x_range.end - p4.x_range.start)
    i4.glyph.y = p4.y_range.start
    i4.glyph.dh = (p4.y_range.start - p4.y_range.end)

    pcmd.data, images.data= current.simulate()

for i in sliders:
    if i.name=='width':
        continue
    i.on_change('value', update_model)
slider_psf.on_change('value', update_width)
#select_Nscale.on_change('value', update_interactor)

psf_button.on_click(update_psf)
button.on_click(update)
curdoc().add_periodic_callback(update_if_change, 500)

interactors = [i for i in sliders]
interactors.append(psf_button)
interactors.append(button)

layout = row([column([p1, p3]), column([p2, p4]), column(interactors)])
#layout = row([p1, p2])
curdoc().add_root(layout)
