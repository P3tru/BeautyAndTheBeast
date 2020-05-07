# BeautyAndTheBeast

Reconstruct offline the events produced by rat-pac-*
Works with ROOT 6, and prior one should source both ROOT and RATROOT libs.
Then, an additional sourcing is necessary (ROOT_INCLUDE_PATH and LD_LIBRARY_PATH):
```bash
source SetROOTIncludePath.sh
```

At the moment it's a ~~3~~ 4 steps process:

## Step 1 : Flatten the hits
```bash
./FlattenHits -i output_from_rat.root -o FLAT_output.npz
```

## Step 2 : Create PDFs
```bash
python OffPDFs.py FLAT_output.npz -v -o PREFIX
```

## Step 3 : PROFIT!
```bash
python OffRecon.py FLAT_output.npz --TResPDF PREFIX_TRes_PDF.npy --CTPDF PREFIX_CTheta_PDF.npy -r PREFIX_RIndex_PDF.npy -o OUTPUT
```

## Step 4 : Rewrite new EV object inside original RAT file
```bash
./FillBackRATMC -i--RAT output_from_rat.root -i--FLAT OUTPUT.root
```
