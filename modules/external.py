"""Adapters and explicit result schemas for specialist transcriptome tools.
No shell execution, simulated biological calls, or inferred PSI from DEJU logFC.
"""
import csv
import gzip
import math
import re
from pathlib import Path

EVENTS = {'SE': 'exon_skipping', 'A5SS': 'alternative_donor',
          'A3SS': 'alternative_acceptor', 'MXE': 'mutually_exclusive_exons', 'RI': 'intron_retention'}


def write_tsv(path, rows, fields):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path, required=()):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        if not set(required) <= set(reader.fieldnames or []):
            raise ValueError(f'Missing columns in {path}: {set(required)-set(reader.fieldnames or [])}')
        return list(reader)


def number(value):
    try:
        result = float(value)
        return result if math.isfinite(result) else None
    except (ValueError, TypeError):
        return None


def bh(values):
    out = [None] * len(values)
    valid = sorted((p, i) for i, p in enumerate(values) if p is not None)
    running = 1.0
    for rank in range(len(valid), 0, -1):
        p, i = valid[rank - 1]
        running = min(running, p * len(valid) / rank)
        out[i] = running
    return out


def event_junctions(row, kind):
    """rMATS exon starts are 0-based, exon ends are 1-based exclusive/closed equivalent."""
    def i(key):
        return int(row[key])
    if kind == 'SE':
        return [(i('upstreamEE'), i('exonStart_0base') + 1, 'inclusion'),
                (i('exonEnd'), i('downstreamES') + 1, 'inclusion'),
                (i('upstreamEE'), i('downstreamES') + 1, 'skipping')]
    if kind == 'MXE':
        return [(i('upstreamEE'), i('1stExonStart_0base') + 1, 'first_exon'),
                (i('1stExonEnd'), i('downstreamES') + 1, 'first_exon'),
                (i('upstreamEE'), i('2ndExonStart_0base') + 1, 'second_exon'),
                (i('2ndExonEnd'), i('downstreamES') + 1, 'second_exon')]
    if kind == 'RI':
        return [(i('upstreamEE'), i('downstreamES') + 1, 'spliced_exclusion')]
    if kind in ('A5SS', 'A3SS'):
        if i('flankingEE') <= i('longExonStart_0base'):
            return [(i('flankingEE'), i('longExonStart_0base') + 1, 'long_exon'),
                    (i('flankingEE'), i('shortES') + 1, 'short_exon')]
        return [(i('longExonEnd'), i('flankingES') + 1, 'long_exon'),
                (i('shortEE'), i('flankingES') + 1, 'short_exon')]
    raise ValueError(kind)


def summarize_events(directory, out, treatment_ids, reference_ids, fdr=.05, delta_filter=.1):
    """b1 must be treatment: rMATS IncLevelDifference is b1 minus b2."""
    all_rows, psi, links = [], [], []
    for mode in ('JC', 'JCEC'):
        mode_rows = []
        for kind, label in EVENTS.items():
            path = Path(directory) / f'{kind}.MATS.{mode}.txt'
            if not path.is_file():
                raise ValueError(f'Missing rMATS result: {path}')
            for raw in read_tsv(path, ('ID','GeneID','chr','strand','PValue','FDR','IncLevel1','IncLevel2','IncLevelDifference',
                      'IncFormLen','SkipFormLen','IJC_SAMPLE_1','IJC_SAMPLE_2','SJC_SAMPLE_1','SJC_SAMPLE_2')):
                gene = raw['GeneID'].strip('"')
                event_id = f'{kind}:{raw["ID"]}'
                p = number(raw['PValue'])
                if p is not None and not 0 <= p <= 1:
                    raise ValueError(f'Invalid event p-value: {event_id}')
                group_values = []
                for group, ids in ((1, treatment_ids), (2, reference_ids)):
                    vals = raw[f'IncLevel{group}'].split(',')
                    inc = raw[f'IJC_SAMPLE_{group}'].split(',')
                    skip = raw[f'SJC_SAMPLE_{group}'].split(',')
                    if not len(vals) == len(inc) == len(skip) == len(ids):
                        raise ValueError(f'rMATS replicate order/number mismatch: {event_id}')
                    parsed = []
                    li, ls = float(raw['IncFormLen']), float(raw['SkipFormLen'])
                    if li <= 0 or ls <= 0:
                        raise ValueError('Invalid effective event lengths')
                    for sid, val, ic, sc in zip(ids, vals, inc, skip):
                        value, ic, sc = number(val), float(ic), float(sc)
                        if ic < 0 or sc < 0 or not ic.is_integer() or not sc.is_integer():
                            raise ValueError('Invalid event counts')
                        expected = (ic / li) / (ic / li + sc / ls) if ic + sc else None
                        if value is not None and (not 0 <= value <= 1 or expected is None or abs(value - expected) > .002):
                            raise ValueError(f'PSI disagrees with length-normalized event counts: {event_id}')
                        parsed.append(value)
                        psi.append(dict(event_id=event_id, event_type=label, count_type=mode, gene_id=gene,
                                        sample_id=sid, group='treatment' if group == 1 else 'reference',
                                        PSI=value, inclusion_count=ic, skipping_count=sc,
                                        inclusion_effective_length=li, skipping_effective_length=ls))
                    group_values.append([v for v in parsed if v is not None])
                delta = number(raw['IncLevelDifference'])
                if all(group_values):
                    expected = sum(group_values[0])/len(group_values[0]) - sum(group_values[1])/len(group_values[1])
                    if delta is None or abs(delta - expected) > .003:
                        raise ValueError(f'Incorrect delta PSI: {event_id}')
                row = dict(event_id=event_id, event_type=label, count_type=mode, gene_id=gene,
                           chr=raw['chr'], strand=raw['strand'], delta_PSI=delta, PValue=p,
                           rMATS_type_FDR=number(raw['FDR']))
                mode_rows.append(row)
                if mode == 'JC':
                    for left, right, role in event_junctions(raw, kind):
                        links.append(dict(event_id=event_id,event_type=label,gene_id=gene,chr=raw['chr'],
                                          strand=raw['strand'],left=left,right=right,role=role))
        for row, q in zip(mode_rows, bh([r['PValue'] for r in mode_rows])):
            row['all_event_types_FDR'] = q
            row['significant'] = q is not None and q <= fdr
            row['passes_delta_filter'] = row['significant'] and row['delta_PSI'] is not None and abs(row['delta_PSI']) >= delta_filter
        all_rows.extend(mode_rows)
    fields = ['event_id','event_type','count_type','gene_id','chr','strand','delta_PSI','PValue','rMATS_type_FDR',
              'all_event_types_FDR','significant','passes_delta_filter']
    write_tsv(Path(out)/'events.tsv', all_rows, fields)
    write_tsv(Path(out)/'PSI_by_sample.tsv', psi, ['event_id','event_type','count_type','gene_id','sample_id','group','PSI',
              'inclusion_count','skipping_count','inclusion_effective_length','skipping_effective_length'])
    write_tsv(Path(out)/'event_junctions.tsv', links, ['event_id','event_type','gene_id','chr','strand','left','right','role'])
    return all_rows, links


def link_dju_events(dju_file, links, output):
    by_key = {}
    for row in links:
        key = tuple(str(row[k]) for k in ('gene_id','chr','strand','left','right'))
        by_key.setdefault(key, []).append(row)
    rows = []
    for j in read_tsv(dju_file):
        if not all(k in j for k in ('chr','strand','left','right')):
            continue
        key = tuple(j[k] for k in ('gene_id','chr','strand','left','right'))
        matches = by_key.get(key, [dict(event_id='',event_type='unclassified',role='')])
        for event in matches:
            rows.append(dict(feature_id=j['feature_id'],gene_id=j['gene_id'],junction_FDR=j['junction_FDR'],
                event_id=event['event_id'],event_type=event['event_type'],role=event['role']))
    write_tsv(output, rows, ['feature_id','gene_id','junction_FDR','event_id','event_type','role'])


def gtf_transcripts(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) != 9 or f[2] != 'transcript':
                continue
            attrs = dict(re.findall(r'(\S+)\s+"([^"]*)"', f[8]))
            if not attrs.get('transcript_id') or not attrs.get('gene_id'):
                raise ValueError('Transcript GTF lacks IDs')
            rows.append(dict(attrs, chr=f[0], start=int(f[3]), end=int(f[4]), strand=f[6]))
    return rows


def summarize_isoforms(gtfs, annotated_gtf, out):
    matrices = {}
    for sid, path in gtfs.items():
        transcripts = gtf_transcripts(path)
        seen = set()
        for row in transcripts:
            tid = row['transcript_id']
            if tid in seen:
                raise ValueError(f'Duplicate quantified transcript {tid}')
            seen.add(tid)
            val = number(row.get('TPM'))
            if val is None or val < 0:
                raise ValueError(f'Missing/invalid StringTie TPM for {tid}')
            matrices.setdefault(tid, {})[sid] = val
    # Missing transcript records remain missing; do not invent zero abundance.
    matrix = [dict(transcript_id=tid,**{sid: vals.get(sid, 'NA') for sid in gtfs}) for tid, vals in sorted(matrices.items())]
    write_tsv(Path(out)/'transcript_TPM.tsv',matrix,['transcript_id',*gtfs])
    candidates=[]
    for row in gtf_transcripts(annotated_gtf):
        code=row.get('class_code','unknown')
        row['candidate_type']='intergenic_transcript_novel_gene_candidate' if code=='u' else ('annotated_intron_chain' if code=='=' else 'other_transcript_candidate')
        candidates.append(row)
    write_tsv(Path(out)/'transcript_classification.tsv',candidates,
              ['transcript_id','gene_id','chr','start','end','strand','class_code','ref_gene_id','cmp_ref','candidate_type'])
    write_tsv(Path(out)/'novel_gene_candidates.tsv',[r for r in candidates if r.get('class_code')=='u'],
              ['transcript_id','gene_id','chr','start','end','strand','class_code','candidate_type'])


def summarize_fusions(files, output):
    rows=[]
    for sid,path in files.items():
        for row in read_tsv(path, ('#gene1','gene2','breakpoint1','breakpoint2','confidence')):
            if '#gene1' in row:
                row['gene1']=row.pop('#gene1')
            if not all(k in row for k in ('gene1','gene2','breakpoint1','breakpoint2','confidence')):
                raise ValueError('Unrecognized Arriba fusion schema')
            rows.append(dict(row,sample_id=sid,interpretation='fusion_candidate_requires_review'))
    fields=['sample_id','gene1','gene2','breakpoint1','breakpoint2','confidence','type','split_reads1','split_reads2',
            'discordant_mates','filters','interpretation']
    write_tsv(output,rows,fields)
    return rows
