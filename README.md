## Predict chemical shieldings of membranes

Example file of computing chemical shieldings for a membrane model system.

Install the required packages, ideally in a fresh virtual environment:
```bash
pip install -r requirements.txt
```

Run the `get_data.py` script from within the `./data/` directory to download the example trajectory:

```bash
cd data
python get_data.py
```

This will download a small example trajectory from Zenodo.
Then from within the `./analysis/` directory, run the `analysis.py`script to compute the chemical shieldings:   

```bash
cd analysis
python analysis.py
```


