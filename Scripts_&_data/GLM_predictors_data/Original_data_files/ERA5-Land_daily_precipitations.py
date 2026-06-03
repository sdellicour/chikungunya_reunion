import cdsapi

c = cdsapi.Client()

c.retrieve(
    'reanalysis-era5-land',
    {
        'variable': [
            'total_precipitation',
        ],
        'year': ['2024', '2025'],

        'month': [
            '01','02','03','04','05','06',
            '07','08','09','10','11','12'
        ],
        'day': [
            '01','02','03','04','05','06','07','08','09','10',
            '11','12','13','14','15','16','17','18','19','20',
            '21','22','23','24','25','26','27','28','29','30','31'
        ],
        'time': [
            '12:00'
        ],
        'data_format': 'netcdf',
        'download_format': 'unarchived',
        'area': [-20.8, 55.1, -21.5, 55.9],
    },
    'reunion_era5land.nc'
)