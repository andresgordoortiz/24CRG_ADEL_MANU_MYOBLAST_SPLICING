# betAS Local Installation Instructions

This project uses a custom fork of betAS with the `getDataset()` function.

## Setup (first time or new environment)

1. **Initialize renv** (if not already done):
   ```r
   renv::init()
   ```

2. **Install betAS from local submodule**:
   ```r
   renv::install("./betAS")
   ```

3. **Restore other dependencies**:
   ```r
   renv::restore()
   ```

## Daily usage

Simply load the library as normal:
```r
library(betAS)
getDataset(pathTables = "/path/to/your/data.tab", tool = "vast-tools")
```

## Updating betAS

To update the betAS submodule to the latest version:
```bash
git submodule update --remote betAS
git add betAS
git commit -m "Update betAS submodule"
```

Then reinstall:
```r
renv::install("./betAS", force = TRUE)
```

## Technical details

- betAS is included as a git submodule pointing to: `andresgordoortiz/betAS@dev-add-getDataset`
- The fork includes the missing `getDataset()`, `getEvents()`, and related functions
- Local installation avoids .rdb corruption issues seen with GitHub installs