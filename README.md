# Coddingme (calculate geosaptial coordinate using Guass method)
"""The Gauss method, also known as the Gauss-Krüger method, is used to calculate geospatial coordinates. Here's a simplified summary of the steps involved:

Determine the central meridian for the Gauss-Krüger zone.
Find the difference between the central meridian and the longitude of the point you want to calculate coordinates for.
Convert the delta longitude to radians.
Calculate the northing by multiplying the latitude, scale factor, and Earth's radius.
Calculate the easting using a formula involving arc length, delta longitude, and false easting.
Add the false easting and false northing to obtain the final coordinates.
Note that this summary provides a basic overview, and actual implementations may involve more complexity and specialized tools"""
import numpy as np

def nptan(a) :
    nptan = (np.sin(a)/np.cos(a))
    return nptan
def npsec(a):
    npsec = 1 /(np.cos(a))
    return npsec

def derece(der,dak,san):
    derece = der + (dak/60) + (san/3600)
    return derece

def derdaksan(der):
    derece = int(der)
    dakika= int((der-derece)*60)
    saniye= (der-derece-dakika/60)*3600
    
    return derece,dakika,saniye

#Direct problem with Gauss 
def gauss_direct(lat,lon,alpha12,s12,a,b):
    lat = np.deg2rad(lat)
    lon = np.deg2rad(lon)
    alpha12 = np.deg2rad(alpha12)
    e1 = np.sqrt((a**2 - b**2) / (a**2))
    M1 = a * (1 - e1**2) / (1 - e1**2 * np.sin(lat)**2)**(3/2)
    N1 = a / (1 - e1**2 * np.sin(lat)**2)**(1/2)
    lat2 = lat
    lon2 = lon
    alpha2 = alpha12
    dalpha =0
    M2 = M1
    N2 = N1
    i = 0
    while True :
        if i == 10:
            break
        i +=1
        alpha_m = alpha12 + (dalpha) / 2
        Mm = (M1 + M2) / 2
        dlat = (s12 * np.cos(alpha_m))/ Mm
        lat2 =  lat + dlat
        M2 = a * (1 - e1**2) / (1 - e1**2 * np.sin(lat2)**2)**(3/2)
        N2 = a / (1 - e1**2 * np.sin(lat2)**2)**(1/2)
        Mm = (M1 + M2) / 2
        Nm = (N1 + N2) /2
        lat_m = (lat + lat2)/2
        dlamda = (s12 * np.sin(alpha_m)) / (Nm*np.cos(lat_m))
        lon2 = lon + dlamda
        dalpha = dlamda * np.sin(lat_m)
        alpha2 = alpha12 + dalpha
        print(derdaksan(np.rad2deg(lat2)), derdaksan(np.rad2deg(lon2)))
    return(lat2,lon2)

#For example    
lat = derece(41,15,56.25)
lon = derece(26,0,37.11)
alpha12 =derece(37,26,43.1348)
s12 = 40000
a = 6378137.000
b = 6356752.3141
gauss_direct(lat,lon,alpha12,s12,a,b)
