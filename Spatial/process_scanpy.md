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

### add image

uns (unstructured annotation slot ) is used to save images, including high-resolution and low-resolution images.
```{py}
images = {}
        image_files = {
            "hires": "tissue_hires_image.png",
            "lowres": "tissue_lowres_image.png",
            "detected": "detected_tissue_image.jpg",
            "fiducial": "aligned_fiducials.jpg"
        }
        for key, fname in image_files.items():
            fpath = os.path.join(spatial_dir, fname)
            if os.path.exists(fpath):
                if load_images:
                    try:
                        images[key] = np.array(Image.open(fpath))
                        print(f"[INFO] Loaded image: {fname} -> {key}")
                    except Exception as e:
                        print(f"[WARN] Failed to load image {fname}: {e}")
                else:
                    images[key] = fpath
        ad.uns["spatial"]["sample"]["images"] = images

    else:
        raise FileNotFoundError("No filtered_feature_bc_matrix.h5 or filtered_feature_bc_matrix/ found")

```


## visualize image

? how to make sure 

```{py}
sc.pl.spatial(ad)
```

## check image resolution

```{py}
library_id = list(ad.uns["spatial"].keys())[0]  
hires_img = ad.uns["spatial"][library_id]["images"]["hires"]

hires_shape = hires_img.shape  # (height, width, channels)
```
