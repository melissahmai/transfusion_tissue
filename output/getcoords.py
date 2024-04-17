#!/usr/bin/env python3
import sys
import json
import numpy as np

f = open(sys.argv[1]+'.json')

data = json.load(f)
f.close()

nodedat = data['nodes']

N = np.size(nodedat)
coords = np.array(np.zeros((N,3)))

for i in range(N):
   coords[i,:] = [int(nodedat[i]['key']), nodedat[i]['attributes']['x'], nodedat[i]['attributes']['y']]

coords = coords[coords[:,0].argsort()]
np.savetxt(sys.argv[1]+'.csv',coords,delimiter=',',fmt='%f',header='ID,x,y',comments='')
