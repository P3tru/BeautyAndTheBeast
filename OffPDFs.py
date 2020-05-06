import OffAnaFunct

import argparse
import numpy as np
import os
from progress.bar import IncrementalBar

''' Start by parsing input file name, and true MC informations '''
parser = argparse.ArgumentParser()
parser.add_argument("input", help="FLAT npz file to create PDF")
parser.add_argument("-x", "--x_true", help="X true position", type=float, default=0.)
parser.add_argument("-y", "--y_true", help="Y true position", type=float, default=0.)
parser.add_argument("-z", "--z_true", help="Z true position", type=float, default=0.)
parser.add_argument("-t", "--t_true", help="T true time", type=float, default=0.)
parser.add_argument("-th", "--theta_true", help="Theta true", type=float, default=0.)
parser.add_argument("-p", "--phi_true", help="Phi true", type=float, default=0.)
parser.add_argument("-o", "--output", help="PDFs output path", type=str, default="")
parser.add_argument("-v", "--verbose", help="increase output verbosity", action="store_true")
args = parser.parse_args()

evts = np.load(args.input)

dirname, filename = os.path.split(args.input)
output_basename = os.path.splitext(filename)

output = args.output if len(args.output) > 0 else output_basename[0]

if args.verbose:
    print("PDFs will be saved to\n {0}\n {1}".format(output + "_TRes_PDF", output + "_CTheta_PDF"))
    print(" rindex will be saved to {0}".format(output + "_RIndex_PDF"))

bar = IncrementalBar("", max=len(evts))

X_TRUE = args.x_true
Y_TRUE = args.y_true
Z_TRUE = args.z_true
T_TRUE = args.t_true
THETA_TRUE = args.theta_true
PHI_TRUE = args.phi_true

TrueOrigin = [X_TRUE, Y_TRUE, Z_TRUE]
TrueDir = [np.cos(PHI_TRUE) * np.sin(THETA_TRUE), np.sin(PHI_TRUE) * np.sin(THETA_TRUE), np.cos(THETA_TRUE)]

rindex = OffAnaFunct.GetAverageRindex(evts, x0_true=[X_TRUE, Y_TRUE, Z_TRUE])
SoL = 299.792 / rindex  # mm/ns

vTRes = []
vCT = []
for evt in evts:
    TRes, CT = OffAnaFunct.CreateEvtTResAndCTheta(evts[evt], SoL, TrueOrigin, T_TRUE, TrueDir)
    if args.verbose:
        bar.next()
    vTRes.append(TRes)
    vCT.append(CT)

bar.finish()

tresids = np.concatenate(vTRes)
costhetas = np.concatenate(vCT)

np.save(output + "_TRes_PDF", tresids)
np.save(output + "_CTheta_PDF", costhetas)
np.save(output + "_RIndex_PDF", np.asarray(rindex))
