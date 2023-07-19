import cdsapi

c = cdsapi.Client()

# years to access
years_to_access = [str(i) for i in range(1950, 2021)]
months_to_access = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']

# viet geographical extent (ymax, xmin, ymin, xmax)
geo_extent = [24, 100, 7.5, 110]

# mean total preciptiation 1997_2019
c.retrieve(
    'reanalysis-era5-land-monthly-means',
    {
        'format': 'netcdf',
        'product_type': 'monthly_averaged_reanalysis',
        'variable': '2m_temperature',
        'year': years_to_access,
        'month': months_to_access,
        'time': '00:00',
        'area': geo_extent,
    },
    'tmean_month_viet.nc')
