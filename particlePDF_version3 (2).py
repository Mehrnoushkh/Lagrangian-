
# coding: utf-8

# In[11]:


from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation, rc
import matplotlib.colors as colors
from netCDF4 import Dataset
import seaborn as sns
from scipy.integrate import simps
from scipy.stats import skew
import scipy as sp
import pandas as pd
from scipy import interpolate


# In[12]:


from mpi4py import MPI
wRank = MPI.COMM_WORLD.Get_rank()
wSize = MPI.COMM_WORLD.Get_size()


# In[13]:


def non_nan_indices_2(array):
    non_nan_indices = np.where(np.logical_not(np.isnan(array)))[0]
    if len(non_nan_indices) > 0:
        first_non_nan_index = non_nan_indices.min()
        last_non_nan_index = non_nan_indices.max()
        return first_non_nan_index, last_non_nan_index
    else:
        return None,None


# In[14]:


fname = 'particle_fine_500km.nc'
ds = Dataset(fname)
ds.variables.keys()


# In[15]:


lon = np.array(ds.variables['longitude'])
lat = np.array(ds.variables['latitude'])
Itime = np.array(ds.variables['time'])
traj = np.array(ds.variables['trajectory'])
vel = np.array(ds.variables['vel_lon'])
u = np.array(ds.variables['vel_lon'])
v = np.array(ds.variables['vel_lat'])


# In[16]:


lon[np.where(lon < -1e3)] = np.nan
lat[np.where(lat < -1e3)] = np.nan


# In[17]:


def distance( lon1, lat1, lon2, lat2 ):
    
    acos_argument = np.sin(lat1 ) * np.sin(lat2 )  +  np.cos(lat1 ) * np.cos(lat2 ) * np.cos((lon1 - lon2))
    
    #if ( np.all(acos_argument < -1)) and (np.all(abs(acos_argument+1) < 1e-10)):
   #     delta_sigma = np.pi
   # elif (np.all(1 - acos_argument) < 1e-10):
   #     delta_sigma = 0
   # else:
    delta_sigma=np.arccos( acos_argument )
    #acos_argument = min( max( acos_argument, 0 ), np.pi )
    
    return 6371e3 * delta_sigma


# In[18]:


# *************.  MAASKIN  ********************

fname_o = '../full_AVISO_timeseries.nc'
ds2 = Dataset(fname_o)
long = np.array(ds2.variables['longitude'])
latit = np.array(ds2.variables['latitude'])



print(ds2.variables.keys())

mask = ds2['uo'][0,0,:,:].mask

print(mask.shape)


# In[19]:


print( lat[200,1].min(), lat[201,1].max() )


# In[10]:


print( lat[200,2], lat[201,2] )


# In[23]:


ii = 134
t1 = 26
t2 =t1+1 


lat1 = lat[t1,ii]
lat2 = lat[t2,ii]

lon1 = lon[t1,ii]
lon2 = lon[t2,ii]
print(distance( lon1, lat1, lon2, lat2 ))


# In[13]:


ii=0
d = distance( lon[1:,ii], lat[1:,ii], lon[:-1,ii], lat[:-1,ii])


# In[14]:


plt.plot(d[:-1] / 1e3)
plt.yscale('log')


# In[9]:


lon.shape


# In[ ]:


plt.figure(figsize=(12, 8))
plt.pcolormesh(long, latit, mask, cmap='seismic')
plt.colorbar()

#plt.scatter(lonOrig.flatten()[0:500000]*180/np.pi, latOrig.flatten()[0:500000]*180/np.pi, c = 'w', s=1 )
plt.scatter(lon.flatten()[0:500000]*180/np.pi, lat.flatten()[0:500000]*180/np.pi, c = 'g', s=3 )


# In[15]:


lon = np.rad2deg(lon)
lat = np.rad2deg(lat)

print(np.nanmin(lon), np.nanmax(lon))
print(np.nanmin(lat), np.nanmax(lat))

lon += 180.0
lat += 90
lat = abs(lat)

print(np.nanmin(lon), np.nanmax(lon))
print(np.nanmin(lat), np.nanmax(lat))


lonIndices = np.array(np.rint(lon/0.25), dtype=int)
latIndices = np.array(np.rint(lat/0.25), dtype=int)

print('lonindices', np.nanmin(lonIndices), np.nanmax(lonIndices))
print('latindices', np.nanmin(latIndices), np.nanmax(latIndices))



oneDindices = latIndices.copy() *1440 + lonIndices.copy()

print('oneDindices',oneDindices.shape, latIndices.shape, lonIndices.shape, np.nanmin(oneDindices), np.nanmax(oneDindices))


# In[27]:


#plt.pcolormesh(mask)


# In[26]:


#a = np.arange(1440*720).reshape(720,1440)
#plt.pcolormesh(a)

#a.flatten()[0:10]


# In[16]:


#flatMask = mask.T.flatten()
#flatMask = mask.T.flatten()
flatMask = mask.flatten()
#flatMask = np.flip(mask, axis=0).flatten()
print(np.sum(flatMask))
indices = np.array(np.where(flatMask == 1)[0], dtype=int)
print('indices', np.max(indices), np.min(indices))


# In[13]:


np.sum(oneDindices == indices[20000])


# In[14]:


oneDindices.shape , flatMask.shape


# In[15]:


#indices=1036800
#finalMask = np.zeros(oneDindices.shape, dtype=bool)
#for index in range(wRank, indices, wSize):

 #   finalMask[oneDindices == index] = True


# In[25]:


finalMask = np.zeros(oneDindices.shape, dtype=bool)
for index in indices[700000:1036800]:
    finalMask[oneDindices == index] = True


# In[16]:


np.max((oneDindices == indices[60000])*1)


# In[17]:


#oneDindices.dtype, indices.dtype


# In[13]:


lon = np.array(ds.variables['longitude'])
lat = np.array(ds.variables['latitude'])
lon[np.where(lon < -1e3)] = np.nan
lat[np.where(lat < -1e3)] = np.nan

#

lonOrig = np.array(ds.variables['longitude'])
latOrig = np.array(ds.variables['latitude'])

lat[finalMask] = float('nan')
lon[finalMask] = float('nan')

lat[finalMask] = float('nan')
lon[finalMask] = float('nan')


# In[26]:


lat[finalMask] = float('nan')
lon[finalMask] = float('nan')

lat[finalMask] = float('nan')
lon[finalMask] = float('nan')


# In[27]:


lonOrig = np.array(ds.variables['longitude'])
latOrig = np.array(ds.variables['latitude'])

plt.figure(figsize=(12, 8))
plt.pcolormesh(long, latit, mask, cmap='seismic')
plt.colorbar()

plt.scatter(lonOrig.flatten()[0:500000]*180/np.pi, latOrig.flatten()[0:500000]*180/np.pi, c = 'w', s=1 )
plt.scatter(lon.flatten()[0:500000]*180/np.pi, lat.flatten()[0:500000]*180/np.pi, c = 'g', s=3 )


# In[41]:


#plt.scatter(long[lonIndices.flatten()][0:100], latit[latIndices.flatten()][0:100], c = oneDindices.flatten()[0:100], cmap='nipy_spectral')


# In[9]:


## Refined version:  In the given scenario, particles are observed to move within a range of -180 to 180 degrees, oscillating back and forth


ii = 0
epsilon = 0.01
epsilon_2 = 0.01

##
## First, for-loop to split the particle records into individual particles (i.e. split up where they are recycled)
##

for ii in range(1600):
   
    c_1 = np.abs( np.sin(lon[1:, ii]) - np.sin(lon[:-1, ii]) ) > epsilon
    c_2 = np.abs( np.cos(lon[1:, ii]) - np.cos(lon[:-1, ii]) ) > epsilon
    c_3 = np.abs(        lon[1:, ii]  -        lon[:-1, ii]  ) > epsilon_2
              

    recycle_indices = np.where( np.logical_or(c_1,np.logical_or(c_2,c_3)))[0]  #times at which we have extremum in data

    if len(recycle_indices > 0):
                 
        Ndeath = len(recycle_indices)   # number of times a particle disappears and re-appears
        
        # Add in the appropriate number of new particles
        padding = np.tile( lon[:,ii].reshape((len(Itime),1)), (1,Ndeath) ) 
        lon = np.concatenate( (lon, padding) , axis = 1)  #copy the column of each traj 
        
        padding = np.tile( lat[:,ii].reshape((len(Itime),1)), (1,Ndeath) ) 
        lat = np.concatenate( (lat, padding) , axis = 1) 
        # Fill in with nan values where the particle doesn't exist 
        j = recycle_indices[0]
        lon[j:,ii] = np.nan         # put nan value after the recycle event
        lat[j:,ii] = np.nan  
        for Ij in range(Ndeath-1):
            lon[:recycle_indices[Ij]+1,-(1+Ij)] = np.nan
            lon[recycle_indices[Ij+1]:,-(1+Ij)] = np.nan
            
            lat[:recycle_indices[Ij]+1,-(1+Ij)] = np.nan
            lat[recycle_indices[Ij+1]:,-(1+Ij)] = np.nan
            
        lon[:recycle_indices[-1],-Ndeath] = np.nan
        lat[:recycle_indices[-1],-Ndeath] = np.nan
    print(ii, len(recycle_indices))
   # print(ii, len(recycle_indices))

print('done')


# In[37]:


with Dataset('No_Land_particle_file.nc', 'r', format='NETCDF4') as ds3:
    lon = np.array(ds3.variables['longitude'])
    lat = np.array(ds3.variables['latitude'])
    print(ds3.variables.keys())
    print(lon.shape, lat.shape)


# In[39]:


plt.figure(figsize=(18, 10))
#plt.pcolormesh(long, latit, mask,cmap= 'YlGn')

#plt.colorbar()
plt.scatter(lon[:-1, :] * 180 / np.pi, lat[:-1,:] * 180 / np.pi, s=0.1, c='blue', linewidths=0.001)


# In[52]:


plt.figure(figsize=(18, 10))
plt.pcolormesh(long, latit, mask,cmap= 'YlGn')

plt.colorbar()
plt.scatter(lon[:-1, :] * 180 / np.pi, lat[:-1,:] * 180 / np.pi, s=0.1, c='yellow', linewidths=0.001)


# In[9]:


for ii in range(1600):
   
   # c_1 = np.abs( np.sin(lon[1:, ii]) - np.sin(lon[:-1, ii]) ) > epsilon
   # c_2 = np.abs( np.cos(lon[1:, ii]) - np.cos(lon[:-1, ii]) ) > epsilon
    #c_3 = np.abs(        lon[1:, ii]  -        lon[:-1, ii]  ) > epsilon_2
    c_1 = distance( lon[1:,ii], lat[1:,ii], lon[:-1,ii], lat[:-1,ii] ) > 20000
              

    recycle_indices = np.where( c_1)[0]  #times at which we have extremum in data

    if len(recycle_indices )> 0:
                 
        Ndeath = len(recycle_indices)   # number of times a particle disappears and re-appears
        
        # Add in the appropriate number of new particles
        padding = np.tile( lon[:,ii].reshape((len(Itime),1)), (1,Ndeath) ) 
        lon = np.concatenate( (lon, padding) , axis = 1)  #copy the column of each traj 
        
        padding = np.tile( lat[:,ii].reshape((len(Itime),1)), (1,Ndeath) ) 
        lat = np.concatenate( (lat, padding) , axis = 1) 
        # Fill in with nan values where the particle doesn't exist 
        j = recycle_indices[0]
        lon[j:,ii] = np.nan         # put nan value after the recycle event
        lat[j:,ii] = np.nan  
        for Ij in range(Ndeath-1):
            lon[:recycle_indices[Ij]+1,-(1+Ij)] = np.nan
            lon[recycle_indices[Ij+1]:,-(1+Ij)] = np.nan
            
            lat[:recycle_indices[Ij]+1,-(1+Ij)] = np.nan
            lat[recycle_indices[Ij+1]:,-(1+Ij)] = np.nan
            
        lon[:recycle_indices[-1],-Ndeath] = np.nan
        lat[:recycle_indices[-1],-Ndeath] = np.nan
    print(ii, len(recycle_indices))
   # print(ii, len(recycle_indices))

print('done')


# In[10]:


lon.shape


# In[10]:


lon.shape


# In[11]:


plt.figure(figsize=(12, 8))
plt.pcolormesh(long, latit, mask)
plt.colorbar()
ii= range(2)

plt.scatter(lon[:-1, ii] * 180 / np.pi, lat[:-1, ii] * 180 / np.pi, s=4, c='red', linewidths=0.1)


# In[21]:


lon.shape


# In[27]:


8889+38


# In[30]:


particle_in_land = my_interp(((lon[10, 6500]*180/np.pi)), ((lat[10, 6500]*180/np.pi))) 
print(particle_in_land)


# In[20]:


my_interp = sp.interpolate.interp2d(long, latit, mask, kind='linear')
for ii in range(3):
    lon[:,ii] = np.logical_not(np.isnan(lon[:,ii]))
    lat[:,ii] = np.logical_not(np.isnan(lat[:,ii]))
particle_in_land = my_interp(((lon[10, 2]*180/np.pi)), ((lat[10, 2]*180/np.pi))) ==1
print(particle_in_land)
if np.all(particle_in_land):
        lon[:, ii] = np.nan
        lat[:, ii] = np.nan
        indices_to_remove.append(ii)
        lon = np.delete( lon, indices_to_remove, axis = 1 )
        lat = np.delete( lat, indices_to_remove, axis = 1 )
print(ii)


# In[19]:


indices_to_remove = []
my_interp = sp.interpolate.interp2d(long, latit, mask, kind='linear')
for ii in range(lon.shape[1]):
    #lon[:,ii] = np.logical_not(np.isnan(lon[:,ii]))
    #lat[:,ii] = np.logical_not(np.isnan(lat[:,ii]))
    particle_in_land = my_interp(((lon[:, ii]*180/np.pi)), ((lat[:, ii]*180/np.pi))) == True
    particle_in_land[np.isnan(lon[:, ii]) | np.isnan(lat[:, ii])] = True

    if np.all(particle_in_land):
        lon[:, ii] = np.nan
        lat[:, ii] = np.nan
        indices_to_remove.append(ii)
    print(ii)


# In[20]:


indices_to_remove


# In[18]:


lonOrig = np.array(ds.variables['longitude'])
latOrig = np.array(ds.variables['latitude'])

plt.figure(figsize=(12, 8))
plt.pcolormesh(long, latit, mask, cmap='seismic')
plt.colorbar()

#plt.scatter(lonOrig[:-1, :] * 180 / np.pi, latOrig[:-1,:] * 180 / np.pi, s=0.1, c='white', linewidths=0.001)

plt.scatter(lon[:-1, 3] * 180 / np.pi, lat[:-1,3] * 180 / np.pi, s=3, c='green', linewidths=0.001)


# In[12]:


indices_to_remove = []
my_interp = sp.interpolate.interp2d(long, latit, mask, kind='linear')
for ii in range(lon.shape[1]):
    particle_in_land = my_interp(((lon[:, ii]*180/np.pi)), ((lat[:, ii]*180/np.pi))) == 1
    particle_in_land[np.isnan(lon[:, ii]) | np.isnan(lat[:, ii])] = True

    if np.all(particle_in_land):
        lon[:, ii] = np.nan
        lat[:, ii] = np.nan
        indices_to_remove.append(ii)
    print(ii)


# In[24]:


# The 'land particles' are just nans, so remove them from the list of particles since they're meaningless
lon = np.delete( lon, indices_to_remove, axis = 1 )
lat = np.delete( lat, indices_to_remove, axis = 1 )


# In[21]:


#df = pd.DataFrame( mean_zonal_speed)
#df2 = pd.DataFrame(negative_edges) 





#for jj in range(len(lon[1])):
   # if np.any(np.abs(mean_zonal_speed[jj]) <1e-3):
   #     mean_zonal_mask[jj] = mean_zonal_speed[jj]
   # else:
   #     continue
        
#q = np.where(np.logical_not(np.isnan(mean_zonal_mask)))
#q = np.array(q)[0]        
#print(q)


#print(mean_zonal_speed[144])
#print(len(q))


# In[22]:


#plt.figure(figsize=(12, 8))
#plt.pcolormesh(long, latit, mask)
#plt.colorbar()
#lon = np.rad2deg(lon)
#lat = np.rad2deg(lat)
#for ii in range(len(q)): 
  

   # plt.scatter(lon[:-1, q[ii]] * 180 / np.pi, lat[:-1, q[ii]] * 180 / np.pi, s=4, c='red', linewidths=0.01)

    #plt.show()


# In[23]:


# remove the particles over land
#lon = np.delete( lon, q, axis = 1 )
#lat = np.delete( lat, q, axis = 1 )
#print(lon.shape)


# In[ ]:


#np.sum(lon[:,6500] == lat[:,6500]), np.sum(np.logical_not(np.isnan(lat[:,6500])))


# In[45]:


plt.figure(figsize=(12, 8))
plt.pcolormesh(long, latit, mask)

plt.colorbar()
plt.scatter(lon[:-1, :] * 180 / np.pi, lat[:-1,:] * 180 / np.pi, s=10, c='red', linewidths=0.001)

plt.savefig('scatter.png')


# In[51]:


##
## Second, shift each particle so that we're looking at deviations from the initial position
##
Ntime, num_particles = lon.shape
print(lon.shape)
for ii in range(num_particles):

    first_index, last_index = non_nan_indices_2(lon[:, ii])
    if first_index is not None:
        begin = first_index
    else:
        begin = 0 

    # Re-centre the starting longitude to 0
    lon_tmp = lon[:,ii] - lon[begin,ii]

    # Periodically roll any values that are outside of [-pi,pi]
    lon_tmp[ lon_tmp >  np.pi ] = lon_tmp[ lon_tmp >  np.pi ] - 2 * np.pi
    lon_tmp[ lon_tmp < -np.pi ] = lon_tmp[ lon_tmp < -np.pi ] + 2 * np.pi

    # Put the new trajectory back in
    lon[:,ii] = lon_tmp

    # If you also want to centre latitudes, uncomment this
    lat[:,ii] = lat[:,ii] - lat[begin,ii]

##
## If you now plot the trajectories, that should all start at 0 longitude. If you uncommented to centre the latitude as well, then all trajectories should start at the origin (might be an interesting figure)
##
 #   print(ii, len(recycle_indices))
print('done')


# In[5]:


#ofile = 'new_particle_file.nc'

# Assuming you have defined the two-dimensional arrays: lat and lon

#Nlat, Nlon = lat.shape  # Get the shape of the lat and lon arrays

# Define the data types for dimensions and variables
dtype_dim = np.float64
dtype_var = np.float32

# Create a new NetCDF file
with Dataset(ofile, 'w', format='NETCDF4') as fp:
    # Create dimensions
    fp.createDimension('latitude', Nlat)
    fp.createDimension('longitude', Nlon)
    #fp.createDimension('longitude', Nlon)

    # Create variables
    lat_var = fp.createVariable('latitude', dtype_dim, ('latitude', 'longitude'))
    lon_var = fp.createVariable('longitude', dtype_dim, ('latitude', 'longitude'))

    # Write data to variables
    lat_var[:] = lat
    lon_var[:] = lon


# In[35]:


with Dataset('new_particle_file.nc', 'r', format='NETCDF4') as ds3:
    lon = np.array(ds3.variables['longitude'])
    lat = np.array(ds3.variables['latitude'])
    print(ds3.variables.keys())
    


# In[10]:


print(lon.shape)
print(np.nanmin(lon))


# In[11]:


4380*3908


# In[12]:


plt.figure(figsize=(12, 8))
#plt.pcolormesh(long, latit, mask)

#plt.colorbar()
plt.scatter(lon[:-1, :] * 180 / np.pi, lat[:-1,:] * 180 / np.pi, s=10, c='red', linewidths=0.001)

#plt.savefig('scatter.png')


# In[13]:


# find the particles on land : step one : finding the noises in data

def non_nan_indices_2(array):
    non_nan_indices = np.where(np.logical_not(np.isnan(array)))[0]
    if len(non_nan_indices) > 0:
        first_non_nan_index = non_nan_indices.min()
        last_non_nan_index = non_nan_indices.max()
        return first_non_nan_index, last_non_nan_index
    else:
        return None,None
num = len(lon[1])
delta_v = np.zeros(num, dtype=float)
mean_zonal_speed = np.zeros( lon.shape[1] )
    


for ii in range( lon.shape[1] ):
    first_index, last_index = non_nan_indices_2(lon[:, ii])
    if first_index is not None and last_index is not None:
        
        
        particle_lon_rad = lon[:,ii]
        particle_lon_m = particle_lon_rad * 6371e3 * np.cos( lat[:,ii] ) 
        particle_del_lon_m = particle_lon_m[last_index] - particle_lon_m[first_index]
        particle_lifespan = Itime[last_index] - Itime[first_index]

        particle_zonal_speed = particle_del_lon_m / particle_lifespan

        mean_zonal_speed[ii] = particle_zonal_speed
    ii += 1
print(mean_zonal_speed)


plt.figure(figsize=(16, 8))

plt.grid(True, which='major')

negative_bins = np.pi * np.concatenate( ( -np.geomspace( 1e-9, 20, 200 )[::-1], np.array([0,]) ) )
positive_bins = np.pi * np.concatenate( ( np.array([0,]), np.geomspace( 1e-9, 20, 200 ) ) )


mean_zonal_mask = np.full(8644, np.nan)  

#for jj in range(8644):
#    if np.any(np.abs(mean_zonal_speed[jj]) > 1e-8):
 #       mean_zonal_mask[jj] = mean_zonal_speed[jj]
 #   else:
 #       continue
    
    
negative_hist, negative_edges = np.histogram( mean_zonal_speed, bins = negative_bins )
positive_hist, positive_edges = np.histogram( mean_zonal_speed, bins = positive_bins )

plt.plot(  0.5 * ( negative_edges[1:] + negative_edges[:-1] ), negative_hist, label = 'Westward' )
plt.plot(   0.5 * ( positive_edges[1:] + positive_edges[:-1] ), positive_hist, label = 'Eastward' )
plt.scatter(0.5 * ( negative_edges[1:] + negative_edges[:-1] ), negative_hist, label = 'Westward')
plt.scatter(   0.5 * ( positive_edges[1:] + positive_edges[:-1] ), positive_hist, label = 'Eastward' )
plt.legend()
#plt.xscale('symlog', linthresh = 1e-4)
#plt.xlim( -60, 60 )

plt.xlabel('velocity (m/s)')
plt.ylabel('Particle density')

#plt.ylim(0,3)

q = np.where(np.logical_not(np.isnan(mean_zonal_speed)))
q = np.array(q)[0]

print(skew(q, axis=0, bias=True))
#plt.xlim(-0.01,0.01)



skewness = skew(q, axis=0, bias=True)
plt.text(0.5, 0.92, f'Skewness: {skewness:.5f} m/s', transform=plt.gca().transAxes, ha='right', va='top',weight='bold', fontsize = '14', color = 'black',fontname='paratype-pt-sans')
#plt.savefig('vel_PDf_noLand.png')


# In[33]:


lonD = lon[:,:]* 180 / np.pi
latD = lat[:,:]* 180 / np.pi

lonD[np.where(lonD > 0)] = np.nan
lonD[np.where(lonD < -50)] = np.nan
lonD.shape

latD[np.where(lonD > 0)] = np.nan
latD[np.where(lonD <-50)] = np.nan


# In[36]:


lonD = lon[:,:]* 180 / np.pi
latD = lat[:,:]* 180 / np.pi


# In[12]:


lonD.shape


# In[43]:


plt.figure(figsize=(12, 8))
#plt.pcolormesh(long, latit, mask)

#plt.colorbar()
plt.scatter(lonD[:-1, :] , latD[:-1,:], s=1, c='black', linewidths=0.001)


# In[44]:


lonD.shape


# In[42]:



for jj in range(3908):
    if np.all(lonD[:,jj]<0) :
        lonD[:,jj] = lonD[:,jj]
    else:
        lonD[:,jj] = np.nan


# In[1]:


mask = np.logical_not(np.isnan(lonD))


lonD = lonD[mask].reshape(-1, np.sum(mask, axis=0))


print(lonD.shape)


# In[17]:


LONd = lonD[np.where(np.logical_not(np.isnan(lonD[:,:])))]
LONd.shape


# In[ ]:


# find the particles on land : step one : finding the noises in data

def non_nan_indices_2(array):
    non_nan_indices = np.where(np.logical_not(np.isnan(array)))[0]
    if len(non_nan_indices) > 0:
        first_non_nan_index = non_nan_indices.min()
        last_non_nan_index = non_nan_indices.max()
        return first_non_nan_index, last_non_nan_index
    else:
        return None,None

mean_zonal_speed = np.empty( lonD.shape[1] )
    


for ii in range( lonD.shape[1] ):
    if np.all(lonD[:,ii] )is not None:
        first_index, last_index = non_nan_indices_2(lonD[:, ii])
        if first_index is not None and last_index is not None:
        
            particle_lon_rad = lonD[:,ii]*np.pi/180
            particle_lon_m = particle_lon_rad * 6371e3 * np.cos( latD[:,ii] ) 
            particle_del_lon_m = particle_lon_m[last_index] - particle_lon_m[first_index]
            particle_lifespan = Itime[last_index] - Itime[first_index]

            particle_zonal_speed = particle_del_lon_m /particle_lifespan

            mean_zonal_speed[ii] = particle_zonal_speed
        ii += 1
print(mean_zonal_speed)


plt.figure(figsize=(16, 8))

plt.grid(True, which='major')

negative_bins = np.pi * np.concatenate( ( -np.geomspace( 1e-9, 20, 400 )[::-1], np.array([0,]) ) )
positive_bins = np.pi * np.concatenate( ( np.array([0,]), np.geomspace( 1e-9, 20, 400 ) ) )


#mean_zonal_mask = np.full(8644, np.nan)  

#for jj in range(8644):
#    if np.any(np.abs(mean_zonal_speed[jj]) > 1e-8):
 #       mean_zonal_mask[jj] = mean_zonal_speed[jj]
 #   else:
 #       continue
    
    
negative_hist, negative_edges = np.histogram( mean_zonal_speed, bins = negative_bins )
positive_hist, positive_edges = np.histogram( mean_zonal_speed, bins = positive_bins )

plt.plot(  0.5 * ( negative_edges[1:] + negative_edges[:-1] ), negative_hist, label = 'Westward' )
plt.plot(   0.5 * ( positive_edges[1:] + positive_edges[:-1] ), positive_hist, label = 'Eastward' )
plt.scatter(0.5 * ( negative_edges[1:] + negative_edges[:-1] ), negative_hist, label = 'Westward')
plt.scatter(   0.5 * ( positive_edges[1:] + positive_edges[:-1] ), positive_hist, label = 'Eastward' )
plt.legend()
#plt.xscale('symlog', linthresh = 1e-4)
#plt.xlim( -60, 60 )

plt.xlabel('velocity (m/s)')
plt.ylabel('Particle density')

#plt.ylim(0,300)

q = np.where(np.logical_not(np.isnan(mean_zonal_speed)))
q = np.array(q)[0]

print(skew(q, axis=0, bias=True))
plt.xlim(-0.05,0.05)



skewness = skew(q, axis=0, bias=True)
plt.text(0.5, 0.92, f'Skewness: {skewness:.5f} m/s', transform=plt.gca().transAxes, ha='right', va='top',weight='bold', fontsize = '14', color = 'black',fontname='paratype-pt-sans')
#plt.savefig('vel_PDf_noLand.png')


# In[31]:


#lonD = lon[:,:]*180/np.pi
#print(lonD.shape)
#for ii in range(5485):
#    lonD[:,ii]=lonD[np.where(lonD[:,ii] < 5)]
#print(lonD.shape)
#print(np.nanmin(lonD[0:,10]),np.nanmax(lonD[0:,10]))


# In[24]:


np.nanmax(particle_lifespan)


# In[73]:


#plt.figure(figsize=(12, 8))
#plt.pcolormesh(long, latit, mask)

#plt.colorbar()
#plt.scatter(lonD[:, :] , lat[:,:] * 180 / np.pi, s=10, c='red', linewidths=0.001)

    #plt.show()

