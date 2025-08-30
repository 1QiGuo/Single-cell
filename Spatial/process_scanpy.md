# Load

## h5

```{py}
ad = sc.read_visium(path=datadir, count_file="filtered_feature_bc_matrix.h5")
```

## matrix triple

```{py}
ad = sc.read_10x_mtx(dir_path, var_names="gene_symbols")
```

### add tissue position to obs and obsm

If there is a matrix triple, we need to add the spatial coordinate to obs (Observation Annotations), and obsm (Observation Multi-dimensional Annotations) manually

```{py}
#read file
coords = pd.read_csv(
                os.path.join(spatial_dir, "tissue_positions_list.csv"),
                header=None
)
#rename
coords.columns = [
                "barcode", "in_tissue", "array_row", "array_col",
                "pxl_col_in_fullres", "pxl_row_in_fullres"
]

coords.set_index("barcode", inplace=True)

# merge into ad.obs
ad.obs = ad.obs.join(coords, how="left")

# ensure coords matches ad.obs.index exactly
coords = coords.reindex(ad.obs.index)

# force float type for spatial coords
coords[["pxl_row_in_fullres", "pxl_col_in_fullres"]] = coords[
   ["pxl_row_in_fullres", "pxl_col_in_fullres"]
].astype(float)

# assign spatial coordinates (for plotting)
ad.obsm["spatial"] = coords[["pxl_row_in_fullres", "pxl_col_in_fullres"]].to_numpy()

```

