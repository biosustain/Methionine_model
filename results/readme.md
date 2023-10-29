Results go here.

Note that, even though Maud's stores `InferenceData` output as a netcdf file
`idata.nc`, here these are stored as [zarr](https://zarr.readthedocs.io/en/stable/#) 
stores. This is to avoid large files that cannot be stored on GitHub.

To load an `InferenceData`, use the arviz method
[`InferenceData.from_zarr`](https://python.arviz.org/en/stable/api/generated/arviz.InferenceData.from_zarr.html),
e.g.:

```python
import arviz as az

idata = az.InferenceData.from_zarr("methionine/idata")
```

To save a new output as zarr, use the arviz function
[`to_zarray`](https://python.arviz.org/en/stable/api/generated/arviz.to_zarr.html),
e.g.:


```python
import arviz as az

idata = az.from_netcdf("my_new_run/idata.nc")
az.to_zarr(idata, "idata")
```
