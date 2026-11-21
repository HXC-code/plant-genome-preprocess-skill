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

## PAL2NAL v14 (installed separately)

PAL2NAL is used by WGDI's `-ks` stage to convert protein alignments into codon
alignments. The workflow release does **not** bundle PAL2NAL. Install it from
the upstream distribution when Ks analysis is needed:

```bash
python code/tools/install_pal2nal.py
```

The installer downloads the official v14 tarball, verifies SHA-256
`0e458d38571ab6d61078047a9e80d843bb2e733414eeff8b1144d1c6bc13848f`, and
stores it under the ignored `external/pal2nal.v14/` directory.

- Author: Mikita Suyama (Kyoto University), with David Torrents and Peer Bork.
- Citation: Suyama, Torrents & Bork (2006), PAL2NAL: robust conversion of
  protein sequence alignments into the corresponding codon alignments.
  *Nucleic Acids Research* 34:W609–W612.
- Upstream: <http://www.bork.embl.de/pal2nal/> (the v14 tarball is still served
  there; newer releases exist).
- License caveat: the historical v14 distribution ships **no explicit license
  file**, and the upstream project page does not state license terms for that
  tarball. The workflow does not redistribute it. Users obtain it directly from
  upstream and are responsible for complying with the applicable terms.
- The newer GitHub project (<https://github.com/grp-bork/pal2nal>) declares an
  MIT license for its current code. That declaration is not being assumed to
  cover the historical v14 tarball.

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
