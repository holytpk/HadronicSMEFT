python plots.py --plot_directory delphes-v4                                                          --sample TT01j1l_HT800
python plots.py --plot_directory delphes-v4  --selection singlelep-njet4p-btag1p                     --sample TT01j1l_HT800
python plots.py --plot_directory delphes-v4  --selection singlelep-AK8merged-njet4p-btag1p           --sample TT01j1l_HT800
python plots.py --plot_directory delphes-v4  --selection singlelep-AK8pt500-AK8merged-njet4p-btag1p  --sample TT01j1l_HT800
python plots.py --plot_directory delphes-v4                                                          --sample TT01j1l_HT800 --reweight
python plots.py --plot_directory delphes-v4  --selection singlelep-njet4p-btag1p                     --sample TT01j1l_HT800 --reweight
python plots.py --plot_directory delphes-v4  --selection singlelep-AK8merged-njet4p-btag1p           --sample TT01j1l_HT800 --reweight
python plots.py --plot_directory delphes-v4  --selection singlelep-AK8pt500-AK8merged-njet4p-btag1p  --sample TT01j1l_HT800 --reweight

python plots.py --plot_directory delphes-v4                                      --sample TT01j1l
python plots.py --plot_directory delphes-v4  --selection singlelep               --sample TT01j1l
python plots.py --plot_directory delphes-v4  --selection singlelep-njet4p        --sample TT01j1l
python plots.py --plot_directory delphes-v4  --selection singlelep-njet4p-btag1p --sample TT01j1l

python plots.py --plot_directory delphes-v4                                      --sample TT01j1l --reweight
python plots.py --plot_directory delphes-v4  --selection singlelep               --sample TT01j1l --reweight
python plots.py --plot_directory delphes-v4  --selection singlelep-njet4p        --sample TT01j1l --reweight
python plots.py --plot_directory delphes-v4  --selection singlelep-njet4p-btag1p --sample TT01j1l --reweight
