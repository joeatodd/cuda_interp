# Example 1

This is a sample 3D dataset representing temperature of Store Glacier, West Greenland.

eg1.ini is a configuration file which specifies the input filename & the required output grid.

## Running the example

First fetch the dataset (400Mb):
```
wget "https://www.dropbox.com/s/s04mfwth9yx4sn9/Seasonal_Spinup_8_body_temp.nc"
```
then run it:
```
/your/install/prefix/CudaInterp eg1.ini
```
