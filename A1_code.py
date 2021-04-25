# -*- coding: utf-8 -*-
"""
Created on Wed Mar 10 13:32:29 2021

@author: Nikos
"""

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

hdulist=fits.open("C:/Users/Nikos/Desktop/mosaic.fits")
values=hdulist[0].data
zp=2.530
zp_er=2.000E-02
plt.hist(values, range=(3350,3500),bins=50)
#plt.yscale('log')
#plt.grid()
#plt.plot()
plt.ylim(0,1000)
#%%
#removing bleeding stars
def remove_circle(xc,yc,rad,values):
    X=np.arange(0,4610,1)
    Y=np.arange(0,2569,1)
    for x in X:
        for y in Y:
            d2 = (x - xc)**2 + (y - yc)**2
            if d2<rad**2:
                values[x][y]=0
    return values

def remove_box(x1,x2,y1,y2,values):
    X=np.arange(0,4610,1)
    Y=np.arange(0,2569,1)
    for y in Y:
        for x in X:
            if (y>y1) and (y<y2):
                if (x>x1) and (x<x2):
                    values[x][y]=0
    return values

def remove_edges(thickness,values):
    X=np.arange(0,4610,1)
    Y=np.arange(0,2569,1)
    for x in X:
        if x<thickness or x>4610-thickness:
            for y in Y:
                values[x,y]=0
    for y in Y:
        if y<thickness or y>2569-thickness:
            for x in X:
                values[x,y]=0
    return values
values=remove_edges(250,values)
values=remove_box(0,4610,1382,1494,values)
values=remove_box(0,545,966,1783,values)

values=remove_circle(3214,1434,365,values)
values=remove_circle(3316,788,106,values)
values=remove_circle(3267,2255,89,values)
values=remove_circle(3755,2147,91,values)
values=remove_circle(4097,565,73,values)
values=remove_circle(2886,862,86,values)
values=remove_circle(2759,989,96,values)
values=remove_circle(2291,915,93,values)
values=remove_circle(2291,2141,99,values)
values=remove_circle(1439,2088,80,values)


#%%
#fits.writeto('newtable7.fits',values)
#%%
std=40
background=3420+5*std
def detect_source(values):
    a=background
    for i in range(0,4611):
        for j in range (0,2570):
            if values[i][j]>a:
                a=values[i][j]
                location=(i,j)
    return a,location

def remove_galaxy(location,values):
    x=location[0]
    y=location[1]
    values1=np.copy(values)
    values1[x][y]=0
    return values1

def remove_adjacent(values,x,y):
    Adj=[]
    removed=0
    if x+1 <= 4610:
        right=values[x+1][y]
        if right> background:
            values[x+1][y]=0
            removed+=1
            Adj.append([x+1,y])
    if x-1 >= 0:
        left=values[x-1][y]
        if left>background:
            values[x-1][y]=0
            removed+=1
            Adj.append([x-1,y])
    if y+1 <= 2569:
        up=values[x][y+1]
        if up>background:
            values[x][y+1]=0
            removed+=1
            Adj.append([x,y+1])
    if y-1>= 0:
        down=values[x][y-1]
        if down>background:
            values[x][y-1]=0
            removed+=1
            Adj.append([x,y-1])
    return (values,Adj,removed)
        
def source_count(loc,values,aperture_rad):
    #X=np.arange(0,4610,1)
    #Y=np.arange(0,2569,1)
    er=0
    S=0
    for x in np.arange(loc[0]-aperture_rad,loc[0]+aperture_rad,1):
        for y in np.arange(loc[1]-aperture_rad,loc[1]+aperture_rad,1):
            d2 = (x - loc[0])**2 + (y - loc[1])**2
            if d2<aperture_rad**2:
                if values[x][y] !=0:
                    S+=values[x][y]-3420
                    er+=std
    return S,er


def magnitude(S,erS):
    m=[]
    er=[]
    for i in range (len(S)):
        m.append(zp-2.5*np.log10(S[i]))
        er.append(zp_er+2.5*(erS[i]/S[i]))
    return m,er

"""    
def count_galaxies(values,minimum_pixels):
    sky=np.copy(values)
    V=[]
    L=[]
    while sky.max()>background:
        det=detect_source(sky)
        #val=det[0]
        loc=det[1]
        sky=remove_galaxy(loc,sky)
        sky,Adj,rem=remove_adjacent(sky,loc[0],loc[1])
        while len(Adj)>0:
            Adj1=[]
            for i in Adj:
                #print(Adj)
                #print(i)
                sky,Adj2,rem1=remove_adjacent(sky,i[0],i[1])
                rem+=rem1
                Adj1.append(Adj2)
            Adj=Adj1[0]
        if rem>minimum_pixels:
            print(loc)
            L.append(loc)
            V.append(source_count(loc,values,3))
    return sky,V,L   
"""
#%%
def catalogue_fast(values,minimum_pixels,removing_rad):
    sky=np.copy(values)
    V=[]
    L=[]
    erV=[]
    while sky.max()>background:
        det=detect_source(sky)
        #val=det[0]
        loc=det[1]
        sky=remove_galaxy(loc,sky)
        v,erv=(source_count(loc,values,3))
        rem=0
        for x in np.arange(loc[0]-removing_rad,loc[0]+removing_rad,1):
            for y in np.arange(loc[1]-removing_rad,loc[1]+removing_rad,1):
                    if sky[x][y]>background:
                        rem+=1
                        sky[x][y]=0
        if rem>minimum_pixels:
            V.append(v)
            erV.append(erv)
            L.append(loc)
            #print(magnitude(v,erv))
    return sky,V,L,erV

def N_smaller_than_m(m,erM,Marray):
    Y=[]
    E=[]
    for i in range(len(m)):
        s=0
        er=0
        for j in Marray:
            if j<m[i]:
                s+=1
                if m[i]+erM[i]<j:
                    er+=1
        Y.append(s)
        E.append(er)
    return Y,E
#%%

sky,V,L,erV=catalogue_fast(values,9,30)
#%%
M,erM=magnitude(V,erV)
m=np.arange(-6,-12,-0.3)
Y,erY=N_smaller_than_m(m,erM,M)
a=np.polyfit(m,np.log(Y),1)

def bestfit(m,a):
    return np.exp(a[0]*m+a[1])

#plt.plot(m,Y,"x",label="data")
plt.errorbar(m,Y,yerr=np.sqrt(Y),fmt="x",label="N(<m)")
plt.plot(m,bestfit(m,a),":",label="lnY=am+c")
ax=plt.gca()
plt.grid()
plt.legend()
plt.ylabel("N(<m)")
plt.xlabel("m")
ax.set_yscale("log")
plt.show()

print("a=",a[0])
print("c=",a[1])
#%%
#export catalogue
import csv

with open('catalogueA1_ck3818.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["location","Magnitude","Magnitude_Error"])
    for i in range(len(L)):
        writer.writerow([L[i],M[i],erM[i]])

