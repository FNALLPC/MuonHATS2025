import ROOT
import correctionlib

correctionlib.register_pyroot_binding()

ROOT.EnableImplicitMT()

ROOT.gROOT.ProcessLine(
    'auto cset = correction::CorrectionSet::from_file("schemaV2.json");'
)
ROOT.gROOT.ProcessLine('#include "MuonScaRe.cc"')

df_data = ROOT.RDataFrame("Events", "data.root")
df_mc = ROOT.RDataFrame("Events", "mc.root")

df_data = df_data.Define(
    'pt_1_corr',
    'pt_scale(1, pt_1, eta_1, phi_1, charge_1, "nom")'
)

# MC apply scale shift and resolution smearing
df_mc = df_mc.Define(
    'pt_1_scale_corr',
    'pt_scale(0, pt_1, eta_1, phi_1, charge_1, "nom")'
)
df_mc = df_mc.Define(
    "pt_1_corr",
    'pt_resol(pt_1_scale_corr, eta_1, float(nTrkLayers_1), "nom")'
)