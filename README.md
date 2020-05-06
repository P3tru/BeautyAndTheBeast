# BeautyAndTheBeast

Reconstruct offline the events produced by rat-pac-*
At the moment it's a 3 steps process:

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
python OffRecon.py --TResPDF PREFIX_TRes_PDF.npy --CTPDF PREFIX_CTheta_PDF.npy -r  PREFIX_RIndex_PDF.npy -o FLAT_output.npz
```
