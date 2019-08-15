import numpy as np
import astropy.wcs as wcs
import matplotlib.pyplot as pl
import astropy.coordinates as coo
import astropy.units as u



def init_coo_graph(dims,w1,fts=16):
    ra = np.array(w1.wcs_pix2world(np.arange(dims[1]),0,0))[0]
    dec = np.array(w1.wcs_pix2world(0,np.arange(dims[0]),0))[1]
    extw=[ra[0],ra[-1],dec[0],dec[-1]]      
    ax=pl.gca()
    ratx=np.array([0.2,0.5,0.8])*(ra.max()-ra.min())+ra.min()
    dectx=np.array([0.2,0.5,0.8])*(dec.max()-dec.min())+dec.min()
    axv = coo.SkyCoord(ratx*u.degree,dectx*u.degree)
    ax.set_xticks(ratx)
    ax.set_xticklabels(axv.ra.to_string(format='latex',precision=2,unit=u.h))
    ax.set_yticks(dectx)
    ax.set_yticklabels(axv.dec.to_string(format='latex',precision=1,unit=u.deg))
    pl.xlabel('RA',fontsize=fts)
    pl.ylabel('DEC',fontsize=fts)
    return extw,ra,dec



def init_coo_graph2(dims,w1,fts=16,precra=2,precd=2,delt=0.2,dval=10.0):
    ra = np.array(w1.wcs_pix2world(np.arange(dims[1]),0,0))[0]
    dec = np.array(w1.wcs_pix2world(0,np.arange(dims[0]),0))[1]
    extw=[ra[0],ra[-1],dec[0],dec[-1]]
    ax=pl.gca()
#    ratx=np.array([0.2,0.5,0.8])*(ra.max()-ra.min())+ra.min()
#    dectx=np.array([0.2,0.5,0.8])*(dec.max()-dec.min())+dec.min()
    nt = 1 + int(min(dims)/2*delt/dval)*2
    pos = np.round(np.arange(min(dims)/2-(nt-1)/2*dval/delt,min(dims)/2+(nt+1)/2*dval/delt-1.0,dval/delt)).astype('i')[0:nt]
    labels = range(int(-(nt-1)/2)*int(dval),int((nt+1)/2)*int(dval),int(dval))
    labels = [str(v) for v in labels]
    ratx=ra[pos]
    dectx=dec[pos]
    print(nt)
    print(pos)
    print(labels)
    axv = coo.SkyCoord(ratx*u.degree,dectx*u.degree)
    ax.set_xticks(ratx)
#    ax.set_xticklabels(axv.ra.to_string(format='latex',precision=precra,unit=u.h))
    ax.set_xticklabels(labels[::-1])    
    ax.set_yticks(dectx)
#    ax.set_yticklabels(axv.dec.to_string(format='latex',precision=precd,unit=u.deg))
    ax.set_yticklabels(labels)
    pl.xlabel(r'$\Delta$RA',fontsize=fts)
    pl.ylabel(r'$\Delta$DEC',fontsize=fts)
    return extw,ra,dec

def init_coo_graph3(dims,w1,fts=16,precra=2,precd=2,delt=0.2,dval=10.0):
    ra = np.array(w1.wcs_pix2world(np.arange(dims[1]),0,0))[0]
    dec = np.array(w1.wcs_pix2world(0,np.arange(dims[0]),0))[1]
    extw=[ra[0],ra[-1],dec[0],dec[-1]]      
    ax=pl.gca()
    nt = 1 + int(min(dims)/2*delt/dval)*2
    pos = np.round(np.arange(min(dims)/2-(nt-1)/2*dval/delt,min(dims)/2+(nt+1)/2*dval/delt,dval/delt)).astype('i')[0:nt]
    labels = range(-(nt-1)/2*int(dval),(nt+1)/2*int(dval),int(dval))
    labels = [str(v) for v in labels]
    ratx=ra[pos]
    dectx=dec[pos]
    print(nt)
    print(pos)
    print(labels)
    axv = coo.SkyCoord(ratx*u.degree,dectx*u.degree)
    ax.set_xticks(ratx)
    ax.set_xticklabels(labels[::-1])    
    ax.set_yticks(dectx)
    ax.set_yticklabels(['']*10)
    pl.xlabel(r'$\Delta$RA',fontsize=fts)
    return extw,ra,dec


def pophead(hh,val):
    try :
        hh.pop(val)
    except:
        pass




def make_plot(data,header_in,full_coo=False,dval=20,fts=16,*args,**kwargs) :
    from matplotlib import rc
    from copy import deepcopy
    #rc('text', usetex = True)
    #rc('axes', linewidth=2)
    header = deepcopy(header_in)

    pophead(header,'HISTORY')
    pophead(header,'COMMENT')
    pophead(header,'')

    w2=wcs.WCS(header,naxis=[1,2])
    if ('CD2_2' in header) :
        delt =  header['CD2_2']*3600.
    else :
        delt = header['CDELT2']*3600.
        
    if (full_coo) :
        extw,ra,dec = init_coo_graph(data.shape,w2,fts=fts)
    else :
        extw,ra,dec = init_coo_graph2(data.shape,w2,precra=1,precd=1,delt=delt,dval=dval,fts=fts)        
    pl.imshow(data[::-1,:],extent=extw,aspect=1./np.cos(np.mean(extw[2:])*np.pi/180.),*args,**kwargs)
    pl.xticks(fontsize = fts)
    pl.yticks(fontsize = fts)
    ax = pl.gca()
    ax.tick_params(length=20,width=4,labelsize=fts)
    ax.tick_params(length=12,width=4,which='minor')
    pl.draw()
    ax.axis(extw)
    return ax,extw,w2
