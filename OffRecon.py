import OffAnaFunct

import argparse
import numpy as np
import os

''' Parsing input arguments '''

parser = argparse.ArgumentParser()
parser.add_argument("input",
                    help="FLAT npz file to reconstruct",
                    type=str)
parser.add_argument("-o", "--output",
                    help="RECON output path",
                    type=str, default="")
parser.add_argument("--TResPDF",
                    help="FLAT npy file with TRes PDF",
                    type=str, default="")
parser.add_argument("--CTPDF",
                    help="FLAT npy file with cos theta PDF",
                    type=str, default="")
parser.add_argument("-r", "--rindex",
                    help="rindex for PDF material",
                    type=str, default="")
parser.add_argument("-c", "--cores",
                    help="nb of cores to run parallel recon",
                    type=int, default=8)
parser.add_argument("-m", "--max",
                    help="max nb of fits",
                    type=int, default=None)
parser.add_argument("-pc", "--prompt-cut", dest="pc",
                    help="prompt cut value for dir fit",
                    type=float, default=None)
parser.add_argument("-x", "--x_true",
                    help="X true position",
                    type=float, default=None)
parser.add_argument("-y", "--y_true",
                    help="Y true position",
                    type=float, default=None)
parser.add_argument("-z", "--z_true",
                    help="Z true position",
                    type=float, default=None)
parser.add_argument("-t", "--t_true",
                    help="T true time",
                    type=float, default=None)
parser.add_argument("-th", "--theta_true",
                    help="Theta true",
                    type=float, default=None)
parser.add_argument("-p", "--phi_true",
                    help="Phi true",
                    type=float, default=None)

args = parser.parse_args()

evts = np.load(args.input)
dirname, filename = os.path.split(args.input)
output_basename = os.path.splitext(filename)

output = args.output if len(args.output) > 0 else output_basename[0]
tag = output + '_RECON'

pdfs = [np.load(args.TResPDF), np.load(args.CTPDF)]

X_TRUE = args.x_true
Y_TRUE = args.y_true
Z_TRUE = args.z_true
T_TRUE = args.t_true
THETA_TRUE = args.theta_true
PHI_TRUE = args.phi_true

cores = args.cores

max_fits = args.max

prompt_cut = args.pc

if args.rindex is not None:
    rindex = np.load(args.rindex)
else:
    rindex = None

res = OffAnaFunct.reconstruct_from_ratpac(tag=tag, evts=evts,
                                          x0_true=None,
                                          cores=cores,
                                          max_fits=max_fits,
                                          prompt_cut=prompt_cut,
                                          pdfs=pdfs,
                                          rindex=rindex)
OffAnaFunct.Export2TTree(res[0] + ".root", res[2])
