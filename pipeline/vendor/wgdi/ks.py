import os
import sys
import shutil
import numpy as np
import pandas as pd
from Bio import SeqIO
import subprocess
from Bio.Phylo.PAML import yn00
import wgdi.base as base


class ks:
    def __init__(self, options):
        base_conf = base.config()
        self.pair_pep_file = 'pair.pep'
        self.pair_cds_file = 'pair.cds'
        self.prot_align_file = 'prot.aln'
        self.mrtrans = 'pair.mrtrans'
        self.pair_yn = 'pair.yn'

        for k, v in base_conf:
            setattr(self, str(k), v)
        for k, v in options:
            setattr(self, str(k), v)
            print(f'{str(k)} = {v}')

    def auto_file(self):
        pairs = []
        with open(self.pairs_file) as f:
            p = ' '.join(f.readlines()[:30])

        # Detect file format and process accordingly
        if 'path length' in p or 'MAXIMUM GAP' in p:
            collinearity = base.read_colinearscan(self.pairs_file)
            pairs = [[v[0], v[2]] for k in collinearity for v in k[1]]
        elif 'MATCH_SIZE' in p or '## Alignment' in p:
            collinearity = base.read_mcscanx(self.pairs_file)
            pairs = [[v[0], v[2]] for k in collinearity for v in k[1]]
        elif '# Alignment' in p:
            collinearity = base.read_collinearity(self.pairs_file)
            pairs = [[v[0], v[2]] for k in collinearity for v in k[1]]
        elif '###' in p:
            collinearity = base.read_jcvi(self.pairs_file)
            pairs = [[v[0], v[2]] for k in collinearity for v in k[1]]
        elif ',' in p:
            collinearity = pd.read_csv(self.pairs_file, header=None)
            pairs = collinearity.values.tolist()
        else:
            collinearity = pd.read_csv(self.pairs_file, header=None, sep='\t')
            pairs = collinearity.values.tolist()

        df = pd.DataFrame(pairs).drop_duplicates()
        df[0] = df[0].astype(str)
        df[1] = df[1].astype(str)
        df.index = df[0] + ',' + df[1]
        return df

    def run(self):
        # Load sequence data
        cds = SeqIO.to_dict(SeqIO.parse(self.cds_file, "fasta"))
        pep = SeqIO.to_dict(SeqIO.parse(self.pep_file, "fasta"))
        df_pairs = self.auto_file()

        # Keep only completed KS rows. Older WGDI versions wrote ID-only rows
        # after a failed alignment, which made later runs skip those pairs.
        ks_columns = ['id1', 'id2', 'ka_NG86', 'ks_NG86', 'ka_YN00', 'ks_YN00']
        completed = pd.DataFrame(columns=ks_columns)
        if os.path.exists(self.ks_file):
            existing = pd.read_csv(self.ks_file, sep='\t').drop_duplicates()
            required = ks_columns
            if all(column in existing.columns for column in required):
                for column in required[2:]:
                    existing[column] = pd.to_numeric(existing[column], errors='coerce')
                completed = existing[
                    existing['ks_NG86'].notna() | existing['ks_YN00'].notna()
                ][required].copy()
                if len(completed):
                    completed['id'] = completed['id1'].astype(str) + ',' + completed['id2'].astype(str)
                    reverse = completed.copy()
                    reverse['id'] = completed['id2'].astype(str) + ',' + completed['id1'].astype(str)
                    valid_ids = pd.concat([completed['id'], reverse['id']])
                    df_pairs.drop(np.intersect1d(df_pairs.index, valid_ids.to_numpy()), inplace=True)

        ks_file = open(self.ks_file, 'w')
        ks_file.write('\t'.join(['id1', 'id2', 'ka_NG86', 'ks_NG86', 'ka_YN00', 'ks_YN00']) + '\n')
        for row in completed[ks_columns].itertuples(index=False, name=None):
            ks_file.write('\t'.join([str(i) for i in row]) + '\n')
        ks_file.flush()

        # Filter valid pairs based on sequence data
        df_pairs = df_pairs[
            (df_pairs[0].isin(cds.keys())) & (df_pairs[1].isin(cds.keys())) &
            (df_pairs[0].isin(pep.keys())) & (df_pairs[1].isin(pep.keys()))
        ]

        pairs = df_pairs[[0, 1]].to_numpy()

        if len(pairs) > 0 and pairs[0][0][:3] == pairs[0][1][:3]:
            allpairs = []
            pair_hash = {}
            for k in pairs:
                if k[0] + ',' + k[1] in pair_hash or k[1] + ',' + k[0] in pair_hash:
                    continue
                else:
                    pair_hash[k[0] + ',' + k[1]] = 1
                    pair_hash[k[1] + ',' + k[0]] = 1
                    allpairs.append(k)
            pairs = allpairs

        calculated = 0
        for k in pairs:
            cds_gene1, cds_gene2 = cds[k[0]], cds[k[1]]
            cds_gene1.id, cds_gene2.id = 'gene1', 'gene2'
            pep_gene1, pep_gene2 = pep[k[0]], pep[k[1]]
            pep_gene1.id, pep_gene2.id = 'gene1', 'gene2'

            # Write sequences to files
            SeqIO.write([cds[k[0]], cds[k[1]]], self.pair_cds_file, "fasta")
            SeqIO.write([pep[k[0]], pep[k[1]]], self.pair_pep_file, "fasta")

            # Compute Ka/Ks values
            kaks = self.pair_kaks(['gene1', 'gene2'])
            if not kaks:
                print(f'Warning: KS calculation skipped {k[0]} {k[1]}')
                continue

            ks_file.write('\t'.join([str(i) for i in list(k) + list(kaks)]) + '\n')
            ks_file.flush()
            calculated += 1

        ks_file.close()

        if len(pairs) > 0 and calculated == 0 and len(completed) == 0:
            raise RuntimeError(
                f"No valid KS values were calculated. Check Muscle, PAL2NAL and YN00 paths. "
                f"Configured paths: muscle={self.muscle_path}, "
                f"pal2nal={self.pal2nal_path}, yn00={self.yn00_path}"
            )

        # Clean up temporary files
        for file in [
            self.pair_pep_file, self.pair_cds_file, self.mrtrans, self.pair_yn,
            self.prot_align_file, '2YN.dN', '2YN.dS', '2YN.t', 'rst', 'rst1', 'yn00.ctl', 'rub'
        ]:
            try:
                os.remove(file)
            except OSError:
                pass

    def pair_kaks(self, k):
        if not self.align():
            return None
        pal = self.pal2nal()
        if not pal:
            return None

        kaks = self.run_yn00()
        if kaks is None:
            return None

        kaks_new = [
            kaks[k[0]][k[1]]['NG86']['dN'], kaks[k[0]][k[1]]['NG86']['dS'],
            kaks[k[0]][k[1]]['YN00']['dN'], kaks[k[0]][k[1]]['YN00']['dS']
        ]
        return kaks_new

    def align(self):
        if self.align_software == 'mafft':
            try:
                command = [self.mafft_path, '--quiet', self.pair_pep_file, '>', self.prot_align_file]
                subprocess.run(" ".join(command), shell=True, check=True)
            except subprocess.CalledProcessError as e:
                print(f"Error while running MAFFT: {e}")
                return False

        elif self.align_software == 'muscle':
            try:
                command = [self.muscle_path, '-align', self.pair_pep_file, '-output', self.prot_align_file, '-quiet']
                subprocess.run(" ".join(command), shell=True, check=True)
            except subprocess.CalledProcessError as e:
                print(f"Error while running Muscle: {e}")
                return False

        return os.path.exists(self.prot_align_file) and os.path.getsize(self.prot_align_file) > 0

    def pal2nal(self):
        perl_executable = shutil.which('perl') or 'perl'
        command = [
            perl_executable,
            self.pal2nal_path,
            self.prot_align_file,
            self.pair_cds_file,
            '-output',
            'paml',
            '-nogap'
        ]
        try:
            with open(self.mrtrans, 'w') as output:
                subprocess.run(command, stdout=output, stderr=subprocess.PIPE, check=True)
        except Exception as e:
            print(f"Error while running PAL2NAL: {e}")
            return False
        return os.path.exists(self.mrtrans) and os.path.getsize(self.mrtrans) > 0

    def run_yn00(self):
        yn = yn00.Yn00()
        yn.alignment = self.mrtrans
        yn.out_file = self.pair_yn
        yn.set_options(icode=0, commonf3x4=0, weighting=0, verbose=1)

        try:
            run_result = yn.run(command=self.yn00_path)
        except Exception as e:
            print(f"Error while running YN00: {e}")
            run_result = None
        return run_result
