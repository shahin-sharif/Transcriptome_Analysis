#!/usr/bin/env python3
"""Deterministic synthetic Salmon-like quantifications, GTF and spliced sequences.
Not real Salmon alignment output or a biological benchmark.
"""
import csv, json, math, random
from pathlib import Path


def make_fixture(out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True);rng=random.Random(2048)
    annotation=[];sequences=[];transcripts=[]
    def row(feature,start,end,gene,tx=None,extra='',phase='.'):
        attr=f'gene_id "{gene}"; gene_name "{gene}"; gene_type "protein_coding";'
        if tx:attr+=f' transcript_id "{tx}"; transcript_type "protein_coding";'
        annotation.append(f'chrSynthetic\tsynthetic\t{feature}\t{start}\t{end}\t.\t+\t{phase}\t{attr} {extra}'.rstrip())
    for g in range(60):
        gene=f'SYN_G{g+1:03d}';start=1000*g+1;exons=[(start,start+101),(start+300,start+401),(start+600,start+701)]
        row('gene',start,start+701,gene)
        for label,indices in [('A',[0,1,2]),('B',[0,2])]:
            tx=f'SYN_T{g+1:03d}_{label}.1';seq=''
            row('transcript',start,start+701,gene,tx)
            for n,j in enumerate(indices,1):
                a,b=exons[j];piece='GCT'*34
                if j==0:piece='ATG'+piece[3:]
                if j==2:piece=piece[:-3]+'TAA'
                seq+=piece
                row('exon',a,b,gene,tx,f'exon_number "{n}";')
                row('CDS',a,b-3 if j==2 else b,gene,tx,phase='0')
            row('start_codon',start,start+2,gene,tx,phase='0')
            row('stop_codon',start+699,start+701,gene,tx,phase='0')
            sequences.append(f'>{tx}\n{seq}\n');transcripts.append((tx,len(seq),g,label))
    (out/'annotation.gtf').write_text('# Synthetic validation annotation; not a human reference.\n'+'\n'.join(annotation)+'\n')
    (out/'transcripts.fa').write_text(''.join(sequences))
    samples=[]
    for condition in ('control','treated'):
        for replicate in range(1,4):
            sid=f'{condition}{replicate}';dest=out/sid;dest.mkdir(exist_ok=True)
            counts=[]
            for g in range(60):
                total=max(100,round((600+g*15)*math.exp(rng.gauss(0,.22))))
                prop=.85 if g<5 and condition=='control' else .15 if g<5 else .45
                prop=max(.03,min(.97,prop+rng.gauss(0,.025)))
                a=round(total*prop);counts.extend([a,total-a])
            eff=[length-80 for _,length,_,_ in transcripts]
            rates=[c/l for c,l in zip(counts,eff)];norm=sum(rates)
            with (dest/'quant.sf').open('w') as fh:
                writer=csv.writer(fh,delimiter='\t',lineterminator='\n');writer.writerow(['Name','Length','EffectiveLength','TPM','NumReads'])
                for (tx,length,_,_),l,rate,c in zip(transcripts,eff,rates,counts):writer.writerow([tx,length,l,rate/norm*1e6,c])
            samples.append([sid,condition,f'{sid}/quant.sf'])
    with (out/'samples.tsv').open('w') as fh:
        w=csv.writer(fh,delimiter='\t',lineterminator='\n');w.writerow(['sample_id','condition','quant']);w.writerows(samples)
    config={'samples':'samples.tsv','out':'../../../example_isoform_results','gtf':'annotation.gtf','transcript_fasta':'transcripts.fa',
      'comparisons':[{'name':'treated_vs_control','reference':'control','treatment':'treated'}],
      'covariates':{},'alpha':.05,'delta_if':.1,'filter':{'iso_count':10,'count_proportion':.7,'if_cutoff':.01,'if_proportion':.5},
      'strip_pipe':False,'consequences':False,'predict_novel_orfs':False,'plots':0}
    (out/'config.json').write_text(json.dumps(config,indent=2)+'\n')
    return out

if __name__=='__main__':make_fixture(Path(__file__).resolve().parents[1]/'examples'/'isoform_switch')
