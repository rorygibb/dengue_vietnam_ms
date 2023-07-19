import cdsapi

c = cdsapi.Client()

# years to access
#years_to_access = [str(i) for i in range(2020, 2021)]
years_to_access = 2021
months_to_access = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']

# viet geographical extent (ymax, xmin, ymin, xmax)
geo_extent = [24, 100, 7.5, 110]

# mean total preciptiation 1997_2019
c.retrieve(
    'reanalysis-era5-land-monthly-means',
    {
        'format': 'netcdf',
        'product_type': 'monthly_averaged_reanalysis',
        'variable': 'total_precipitation',
        'year': years_to_access,
        'month': months_to_access,
        'time': '00:00',
        'area': geo_extent,
    },
    'meanprecip_month_viet_20202021.nc')

# mean wind u 1997_2019
c.retrieve(
    'reanalysis-era5-land-monthly-means',
    {
        'format': 'netcdf',
        'product_type': 'monthly_averaged_reanalysis',
        'variable': '10m_u_component_of_wind',
        'year': years_to_access,
        'month': months_to_access,
        'time': '00:00',
        'area': geo_extent,
    },
    'windu_month_viet_20202021.nc')

# mean wind v 1997_2019
c.retrieve(
    'reanalysis-era5-land-monthly-means',
    {
        'format': 'netcdf',
        'product_type': 'monthly_averaged_reanalysis',
        'variable': '10m_v_component_of_wind',
        'year': years_to_access,
        'month': months_to_access,
        'time': '00:00',
        'area': geo_extent,
    },
    'windv_month_viet_20202021.nc')

# temperature and PEV 1981-present
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
    'meantemp_month_viet_20202021.nc')

c.retrieve(
    'reanalysis-era5-land-monthly-means',
    {
        'format': 'netcdf',
        'product_type': 'monthly_averaged_reanalysis',
        'variable': 'potential_evaporation',
        'year': years_to_access,
        'month': months_to_access,
        'time': '00:00',
        'area': geo_extent,
    },
    'meanpotentialevap_month_viet_20202021.nc')
