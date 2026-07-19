#root -b -q -l TrackStateAnalysis.C'("all_acts.root","hist_acts.root")'
#root -b -q -l TrackStateAnalysis.C'("all_genfit.root","hist_genfit.root")'

root -b -q -l TrackStateAnalysis.C'("all_acts_includesecondaries_fullphi_minpt0p2.root","hist_acts_includesecondaries_fullphi_minpt0p2.root",200000)'
root -b -q -l TrackStateAnalysis.C'("all_genfit_includesecondaries_fullphi_minpt0p2.root","hist_genfit_includesecondaries_fullphi_minpt0p2.root",200000)'
