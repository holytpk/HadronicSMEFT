python plots.py                                                         --sample TT01j1l_HT800
python plots.py --selection singlelep-njet4p-btag1p                     --sample TT01j1l_HT800
python plots.py --selection singlelep-AK8merged-njet4p-btag1p           --sample TT01j1l_HT800
python plots.py --selection singlelep-AK8pt500-AK8merged-njet4p-btag1p  --sample TT01j1l_HT800
python plots.py                                                         --sample TT01j1l_HT800 --reweight
python plots.py --selection singlelep-njet4p-btag1p                     --sample TT01j1l_HT800 --reweight
python plots.py --selection singlelep-AK8merged-njet4p-btag1p           --sample TT01j1l_HT800 --reweight
python plots.py --selection singlelep-AK8pt500-AK8merged-njet4p-btag1p  --sample TT01j1l_HT800 --reweight

python plots.py                                     --sample TT01j1l
python plots.py --selection singlelep               --sample TT01j1l
python plots.py --selection singlelep-njet4p        --sample TT01j1l
python plots.py --selection singlelep-njet4p-btag1p --sample TT01j1l

python plots.py                                     --sample TT01j1l --reweight
python plots.py --selection singlelep               --sample TT01j1l --reweight
python plots.py --selection singlelep-njet4p        --sample TT01j1l --reweight
python plots.py --selection singlelep-njet4p-btag1p --sample TT01j1l --reweight
