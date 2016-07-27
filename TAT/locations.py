import numpy as np
import lxml.etree as et
import datetime
import copy
import pandas


RUNPATH='/home/critechuser/DELFT3D/python/'
# function for appending to xml file

def modtxt(nitem,tag,val):
   for elem in nitem.findall(tag): elem.text=val


def locfile(points,lat,lon,plat,plon,b,ha,tsn,tend,t_ref,SAVEPATH):

  # read xml templates
  tree0=et.ElementTree(file=RUNPATH+'loc0.xml')
  root=tree0.getroot() #rss
  channel=root[0] ##channel

  channel[2][0].text=et.CDATA(channel[2][0].text)

  for elem in list(tree0.getiterator()):
   if elem.tag == 'pubDate': elem.text=datetime.datetime.strftime(t_ref,'%d %b %Y %H:%S')

  tree=et.ElementTree(file=RUNPATH+'loc.xml')
  litem=tree.getroot()[0]

  #print et.tostring(litem,pretty_print=True)


  # APPEND LOCATIONS


  for l in range(points.shape[0]):
   if (lat.min() < points['lat'][l] < lat.max()) & (lon.min() < points['long'][l] < lon.max()):
  #    print points.values[l]
       nitem=copy.deepcopy(litem)

       i=np.abs(lon-plon[l]).argmin()
       j=np.abs(lat-plat[l]).argmin()

 #     print l, i,j

       bd=b[i-1:i+2,j-1:j+2]

       bda=np.where(bd>20.)

       if np.size(bda) > 0 : 

          bloc=bd[bda].min()

          pi,pj=np.where(bd==bd[bda].min())

          ip,jp=i-1+pi.ravel()[0],j-1+pj.ravel()[0]

          tseries=ha[:,ip,jp]

          maxh=tseries.max()
        # print maxh

          if maxh-tseries[0] < 0.05 : continue

          loc=np.argwhere(tseries==maxh)
        # print loc

          dh=tseries-tseries[0]
          try: 
             arr=np.argwhere(dh>0.05)[0]
          except:
             arr=loc

          tMax=tsn[loc.ravel()[0]]/3600.
          tArr=tsn[arr.ravel()[0]]/3600.

#      print points.values[l]

          for elem in nitem.findall('title'): 
            elem.text=et.CDATA('{}: {} ({} m)'.format(points['$country'][l],points['$name'][l],maxh))

          for elem in nitem.findall('description'): 
           elem.text=et.CDATA('Country: {}\n \
Location: {}\n \
Time (hh:mm): {}\n \
Maximum Height: {} m\n \
Time Max (hh:mm): {}'.format(points['$country'][l],points['$name'][l],'{}:00'.format(tArr.astype(int)),maxh,'{}:00'.format(tMax.astype(int))))

          modtxt(nitem,'pubDate',datetime.datetime.strftime(tend,'%d %b %Y %H:00:00'))

          modtxt(nitem,'cityName',points['$name'][l])

          modtxt(nitem,'country',points['$country'][l])

          modtxt(nitem,'maxHeight',maxh.astype(str))

          modtxt(nitem,'ID',points['id'][l].astype(str))

          modtxt(nitem,'timeMaxH','{}:00'.format(tMax.astype(int)))

          modtxt(nitem,'timeMaxH_value',tMax.astype(str))

          modtxt(nitem,'timeArrival','{}:00'.format(tArr.astype(int)))

          modtxt(nitem,'timeArrival_value',tArr.astype(str))

          modtxt(nitem,'cityClass',points['cityclass'][l].astype(str))

          if not pandas.isnull(points['popest'][l]) :
              modtxt(nitem,'popEst',points['popest'][l])
          else:
              modtxt(nitem,'popEst','')

          modtxt(nitem,'{http://www.w3.org/2003/01/geo/wgs84_pos#}long',lon[ip].astype(str))
          modtxt(nitem,'{http://www.w3.org/2003/01/geo/wgs84_pos#}lat',lat[jp].astype(str))

          modtxt(nitem,'depth',bloc.astype(str))

          channel.append(nitem)


  tree0.write(SAVEPATH+'locations.xml',xml_declaration=True,encoding='utf-8',method='xml')

