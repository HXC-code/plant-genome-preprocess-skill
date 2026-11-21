# Third-Party Notices

The repository does not redistribute the WGDI source tree or the external
bioinformatics executables. They are installed separately by the user or by
the bootstrap scripts.

## Python runtime

- WGDI 0.75, installed through `requirements-wgdi-validated.txt`.
- Matplotlib, NumPy, SciPy, pandas, Biopython and tabulate, installed from the
  same pinned requirements file.

Consult each upstream distribution for its copyright and license terms. Python
package metadata and installed license files remain authoritative.

## Bundled third-party script: PAL2NAL v14

`code/tools/pal2nal.v14/` redistributes the PAL2NAL v14 Perl script and its
README, used by WGDI's `-ks` stage to convert protein alignments into codon
alignments. The files are bundled unmodified.

- Author: Mikita Suyama (Kyoto University), with David Torrents and Peer Bork.
- Citation: Suyama, Torrents & Bork (2006), PAL2NAL: robust conversion of
  protein sequence alignments into the corresponding codon alignments.
  *Nucleic Acids Research* 34:W609–W612.
- Upstream: <http://www.bork.embl.de/pal2nal/> (the v14 tarball is still served
  there; newer releases exist).
- License caveat: the v14 distribution ships **no explicit license file**, and
  the upstream project page does not state license terms. Internal use with
  attribution is straightforward; **before any public redistribution of this
  archive, confirm the redistribution terms with the upstream author**.

## External executables

- DIAMOND: <https://github.com/bbuchfink/diamond>
- NCBI BLAST+: <https://blast.ncbi.nlm.nih.gov/>
- MUSCLE: <https://drive5.com/muscle/>
- PAML `yn00`: <https://github.com/abacus-gene/paml>
- Perl: <https://www.perl.org/>

These tools are invoked at runtime and are not bundled in the release archive.
Their own licenses govern their use and redistribution.

The workflow source in this repository is covered by the root `LICENSE`. That
license does not extend to PAL2NAL or to any external tool listed above.
