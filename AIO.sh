#!/bin/bash

./FlattenHits -i $1 -o $1_FLAT.npz
python OffPDFs.py $1_FLAT.npz -o $1
python OffRecon.py $1_FLAT.npz --TResPDF $1_TRes_PDF.npy --CTPDF $1_CTheta_PDF.npy -r $1_RIndex_PDF.npy -o $1
./FillBackRATMC -i--RAT $1 -i--FLAT $1_RECON.root
