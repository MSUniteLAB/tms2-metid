# tms2-metid

Code for annotating biotransformation products from **targeted MS2 (tMS2)** acquisitions driven
by an inclusion list of predicted metabolite masses.

The workflow takes the mzML files of a tMS2 run, sorts every MS2 scan by the precursor it was
fragmented for, reconstructs a pseudo-extracted ion chromatogram for each precursor, picks the
chromatographic peaks on that trace, and builds one consensus MS2 spectrum per peak. Nothing in
it is specific to a compound: the inputs are a set of predicted masses and the mzML files of the
runs that fragmented them.

It is the code behind *Pendimethalin biotransformation in zebrafish larvae* (see **Citation**
below), where it was used on zebrafish larvae exposed to pendimethalin.

## What it does, step by step

1. **Parse.** Every MS2 scan is read from the mzML, with its precursor target, retention time,
   collision energy and centroided peak list.
2. **Group.** Scans are grouped by the precursor mass the instrument was told to isolate, not by
   what it actually found.
3. **Pseudo-XIC.** For each precursor, a trace is built over retention time. Each scan
   contributes its intensity weighted by how similar its fragments are to the group's own
   consensus, and the contribution is doubled when the intact precursor survives the isolation
   within 5 ppm.
4. **Peak picking.** Peaks on that trace are found with a prominence and a width criterion, so a
   precursor that was fragmented throughout the run but never eluted produces no peak.
5. **Consensus spectrum.** The scans across one peak are averaged on a 0.01 Da grid, with the
   scans in which the precursor survived counted three times, and the precursor region excluded
   from the fragment list.
6. **Comparison.** Each consensus spectrum is compared with the spectrum of the parent compound
   as a mirror plot, and candidates that co-elute with a heavier conjugate are examined as
   possible in-source fragments rather than metabolites.

## Requirements

Python 3.10. Tested with the versions in `requirements.txt`:

```
pip install -r requirements.txt
```

## Running it

Put the mzML files in a folder named `converted/` next to the script, and run:

```
python met_id_pipeline.py
```

Two conventions are built in: a file is processed as a sample when its name starts with `Zeb`,
and a file named `parrent_spectrum.mzML` is read as the spectrum of the parent compound.

The script writes to `output/` one MGF with every consensus spectrum, a summary CSV with one
row per feature, and a retention-time plot per precursor. It prints as it goes how many scans
it loaded, how many precursors it grouped them into, and how many peaks it found on each.

## Input data

mzML converted from the instrument raw files with ThermoRawFileParser. The runs are targeted
MS2 acquisitions in positive electrospray: each mass on the inclusion list is fragmented
whenever it appears, whatever its intensity, which is what makes a pseudo-XIC possible.

No acquisition is included here. One run of the published study is 160 MB, and the complete
files are deposited separately; see the data availability statement of the paper.

## Not in this repository

- The scripts written to produce the tables and figures of the paper, and to check its
  numbers. They are bookkeeping around one study, not the method.
- The scripts that assemble the manuscript and its supplementary document. They work on Word
  files, not on data, and one of them carries the authors' contact details.
- The mass spectrometry files, for their size.

## Citation

If you use this code, please cite the paper and the archived release of this repository.

<!-- DOI and full reference to be added when the release is archived. -->

## License

MIT. See `LICENSE`.
