import cdsapi
import os

c = cdsapi.Client()

dataset = "reanalysis-era5-single-levels"

var_dict = {
    "prec": "total_precipitation",
    "evap": "evaporation",
    "tcwv": "total_column_water_vapour", 
    "viwve": "vertical_integral_of_eastward_water_vapour_flux",
    "viwvn": "vertical_integral_of_northward_water_vapour_flux"
}

out_dir = '../era5/test_cam'
os.makedirs(out_dir, exist_ok=True)

year = '2025'
month = '07'
day = ['01', '02', '03', '04', '05', '06', '07']
area = [40, -140, -10, -80]
time = [f'{h:02d}:00' for h in range(24)]
grid = [1.0, 1.0]

for vshort, vlong in var_dict.items():
    outfile = os.path.join(out_dir, f"{vshort}.nc")
    
    request = {
        "product_type": "reanalysis",
        "variable": vlong,
        "year": year,
        "month": month,
        "day": day,
        "time": time,
        "area": area,
        "grid": grid,
        "data_format": "netcdf"
    }
    print("Downloading", outfile)

    c.retrieve(
        dataset,
        request,
        outfile
    )
    print(f"DONE: {vlong}")
