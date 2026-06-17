# Data Directory

Place downloaded public input files here before running the scripts.

The repository does not include large raw sequencing files, processed expression matrices from third-party repositories, or atlas objects. Download these files from GEO or the original public resource named in the manuscript.

Suggested local layout:

```text
data/
  public_datasets/
    GSE216376/
    GSE183841/
    GSE233078/
    GSE175759/
    GSE30122/
```

Some scripts expect the original project folder structure. If your downloaded files are placed elsewhere, update the path variables near the top of the corresponding script.


