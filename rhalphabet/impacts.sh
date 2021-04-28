text2workspace.py ttHbb_r-20to20.txt
combineTool.py -M Impacts -d ttHbb_r-20to20.root -m 125 --doInitialFit --robustFit 1 --rMin -20 --rMax 20
combineTool.py -M Impacts -d ttHbb_r-20to20.root -m 125 --doFits --robustFit 1 --rMin -20 --rMax 20
combineTool.py -M Impacts -d ttHbb_r-20to20.root -m 125 -o impacts.json
plotImpacts.py -i impacts.json -o impacts --blind
