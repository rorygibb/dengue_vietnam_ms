
# =============== CDS api query to access WFDE5 bias-corrected precip data =============

# n.b. https://cds.climate.copernicus.eu/cdsapp#!/dataset/derived-near-surface-meteorological-variables?tab=overview
# n.b. units = kg/m2/sec
# n.b. only runs until 2018 but can cross-ref w/ ERA5-land data

import cdsapi

c = cdsapi.Client()

# years to access
years = [str(i) for i in range(1979, 2020)]
months = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']

# viet geographical extent (ymax, xmin, ymin, xmax)
geo_extent = [24, 100, 7.5, 110]

# iterate through years 
for yx in years:

    # reporting
    print(" ".join(["Accessing data for", yx]))
    file_name = "".join(['wfde5_precip_hourly_', yx, ".tar.gz"])

    c.retrieve(
    'derived-near-surface-meteorological-variables',
    {
        'format': 'tgz',
        'variable': 'rainfall_flux',
        'reference_dataset': 'cru_and_gpcc',
        'year': yx,
        'month': months,
        'version': '2.1',
        'area': geo_extent,
    },
    file_name)

