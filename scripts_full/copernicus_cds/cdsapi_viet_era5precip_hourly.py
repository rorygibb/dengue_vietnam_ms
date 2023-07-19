import cdsapi

c = cdsapi.Client()

years_to_access = [str(i) for i in range(2017, 2020)]
months_to_access = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']

# geographical extent (ymax, xmin, ymin, xmax)
geo_extent = [24, 100, 7.5, 110]

c.retrieve(
    'reanalysis-era5-land',
    {
        'format': 'netcdf',
        'variable': 'total_precipitation',
        'year': years_to_access,
        'month': months_to_access,
        'day': [
            '01', '02', '03',
            '04', '05', '06',
            '07', '08', '09',
            '10', '11', '12',
            '13', '14', '15',
            '16', '17', '18',
            '19', '20', '21',
            '22', '23', '24',
            '25', '26', '27',
            '28', '29', '30',
            '31',
        ],
        'time': [
            '00:00', '01:00', '02:00',
            '03:00', '04:00', '05:00',
            '06:00', '07:00', '08:00',
            '09:00', '10:00', '11:00',
            '12:00', '13:00', '14:00',
            '15:00', '16:00', '17:00',
            '18:00', '19:00', '20:00',
            '21:00', '22:00', '23:00',
        ],
        'area': geo_extent,
    },
    'viet_preciphourly_viet_era5land_20172020.nc')