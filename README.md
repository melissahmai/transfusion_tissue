# Transfusion Tissue Network Model

## Codebase

The following scripts and functions are included within the codebase for this project. Please see file headers for more detailed information about inputs and outputs.

### Network Construction

`build_network.m` is the main script used to iterate through images in the `data/` folder to manually construct networks. It depends on the following scripts/functions:`editpoints.m`, `editchi.m`, `get_pairs.m`, `shownetwork.m`.

For _de novo_ construction of simple, hypothetical models, `bottleneck.m` can be used.

### Main model

The main function is `tt_model.m`. Full documentation of its inputs and outputs are found in the file header.

The following scripts are used prior to `tt_model.m` to prepare inputs or manage configurations:
- `get_configs.m` compiles all available configurations
- `read_connections.m` extracts the identity vector and connectivity matrix from the outputs of `build_network.m` for input into `tt_model.m`.

### Analysis

The following scripts and functions are used in analyzing the model output: `needle_analysis.m`, `summarize.m`.

These files depend on subscripts `toSI.m`, `findpath.m`, `findpath_edge.m`, `get_thetas.m`.

### Visualization

The following files are used for data visualizations: `get_pathfigs.m`, `removeSeries.m`.

### Miscellaneous Utility

The following functions are included for utility or aesthetic convenience:`renamefield.m`, `unpackStruct.m`, `mstart.m`.

The colormap used for visualization is included in `cmap.mat`.