#!/usr/bin/env python3
"""Integrated RNA-seq workflow. Python standard library + R/Bioconductor.
Specialist modules use explicit external programs; missing modules never silently pass.
"""
import argparse
import csv
import hashlib
import html
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import tempfile
import time
from modules.external import (summarize_events, link_dju_events, summarize_isoforms,
                              summarize_fusions, write_tsv)

VERSION = '0.2.0'
ROOT = Path(__file__).resolve().parent
MODULES = {'dge','dju','events','isoforms','fusions','loci','isoform_switch'}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def safe_name(value):
    import re
    return isinstance(value, str) and re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', value) is not None


def load_config(path):
    path = Path(path).resolve()
    with path.open() as fh:
        cfg = json.load(fh)
    allowed = {'samples','out','modules','gene_counts','dju_counts','dju_features','gtf','fasta','design','contrasts',
               'bam','dge','dju','events','isoforms','isoform_switch','fusions','loci','tools','threads','fdr'}
    require(not set(cfg)-allowed, f'Unknown configuration keys: {set(cfg)-allowed}')
    schemas = {
        'bam': {'paired','strand','min_mapq','remove_duplicates'},
        'dge': {'min_count','min_samples','fit_type'},
        'dju': {'engine','min_anchor','min_intron','max_intron','min_count','min_total_count'},
        'events': {'read_length','variable_read_length','min_anchor','delta_filter','comparisons'},
        'isoforms': set(),
        'isoform_switch': {'config'},
        'fusions': {'star_chimeric_bam','blacklist','known_fusions','protein_domains'},
        'design': {'formula','factors','numeric'},
        'tools': {'Rscript','rmats','stringtie','gffcompare','arriba'}
    }
    for section, keys in schemas.items():
        if section in cfg:
            require(isinstance(cfg[section],dict) and not set(cfg[section])-keys, f'Unknown/invalid {section} options')
    require(cfg.get('dge',{}).get('fit_type','parametric') in ('parametric','local','mean'), 'Invalid DESeq2 fit_type')
    require(cfg.get('dju',{}).get('engine','modern') in ('modern','legacy'), 'Invalid DEJU engine')
    require(type(cfg.get('events',{}).get('variable_read_length',False)) is bool, 'variable_read_length must be boolean')
    def absolute(p):
        p=Path(p).expanduser()
        return str((path.parent/p).resolve() if not p.is_absolute() else p.resolve())
    for k in ('samples','out','gene_counts','dju_counts','dju_features','gtf','fasta'):
        if k in cfg:
            cfg[k]=absolute(cfg[k])
            if k!='out':
                require(Path(cfg[k]).is_file(), f'Missing {k}: {cfg[k]}')
    require('samples' in cfg and 'out' in cfg, 'samples and out are required')
    mods=cfg.get('modules',[])
    require(mods and len(mods)==len(set(mods)) and set(mods)<=MODULES, f'modules must be a nonempty unique list drawn from {sorted(MODULES)}')
    require(isinstance(cfg.get('threads',1),int) and 1<=cfg.get('threads',1)<=64, 'threads must be an integer 1..64')
    require(0<cfg.get('fdr',.05)<1,'fdr must be between zero and one')
    with open(cfg['samples'],newline='') as fh:
        samples=list(csv.DictReader(fh,delimiter='\t'))
    ids=[s.get('sample_id') for s in samples]
    require(samples and all(safe_name(s) for s in ids) and len(ids)==len(set(ids)), 'Invalid or duplicate sample_id')
    need_bams = bool(set(mods)&{'events','isoforms','fusions','loci'}) or ('dge' in mods and 'gene_counts' not in cfg) or ('dju' in mods and 'dju_counts' not in cfg)
    if 'dju_counts' in cfg:
        require('dju_features' in cfg,'dju_features is required with dju_counts')
    if need_bams:
        require('gtf' in cfg, 'BAM modules require gtf')
        require(not cfg['gtf'].endswith('.gz'), 'For cross-tool compatibility supply an uncompressed GTF')
        b=cfg.get('bam',{})
        require(type(b.get('paired')) is bool and type(b.get('strand')) is int and b['strand'] in (0,1,2), 'bam.paired and bam.strand must be explicit')
        for s in samples:
            require(s.get('bam'),'BAM path missing from sample sheet')
            p=Path(s['bam']).expanduser()
            s['bam']=str((Path(cfg['samples']).parent/p).resolve() if not p.is_absolute() else p.resolve())
            require(Path(s['bam']).is_file(),f'Missing BAM: {s["bam"]}')
            require(',' not in s['bam'] and '\n' not in s['bam'], 'BAM paths must not contain comma/newline for external tool lists')
        require(len(set(s['bam'] for s in samples))==len(samples),'Duplicate BAM paths cannot be biological replicates')
    b=cfg.get('bam',{})
    require(0<=b.get('min_mapq',10)<=255 and type(b.get('remove_duplicates',False)) is bool,'Invalid BAM filters')
    for k in ('min_anchor','min_intron','max_intron','min_count','min_total_count'):
        if k in cfg.get('dju',{}):
            require(isinstance(cfg['dju'][k],(int,float)) and cfg['dju'][k]>=0,f'Invalid dju.{k}')
    for k in ('min_count','min_samples'):
        if k in cfg.get('dge',{}):
            require(isinstance(cfg['dge'][k],int) and cfg['dge'][k]>=1,f'Invalid dge.{k}')
    if 'fusions' in mods:
        f=cfg.get('fusions',{})
        require('fasta' in cfg and f.get('star_chimeric_bam') is True, 'Arriba requires fasta and fusions.star_chimeric_bam=true for STAR BAMs generated with chimeric evidence')
        require('blacklist' in f,'Supply the Arriba blacklist matching the reference assembly')
        for k in ('blacklist','known_fusions','protein_domains'):
            if k in f:
                f[k]=absolute(f[k]);require(Path(f[k]).is_file(),f'Missing fusion resource: {k}')
    if 'events' in mods:
        e=cfg.get('events',{})
        require(isinstance(e.get('read_length'),int) and e['read_length']>0,'events.read_length is required for effective-length PSI')
        require(e.get('comparisons'),'Specify pairwise events.comparisons explicitly')
        require(0<=e.get('delta_filter',.1)<=1 and isinstance(e.get('min_anchor',8),int) and e.get('min_anchor',8)>=1,'Invalid event delta filter/anchor')
        names=[]
        for cmp in e['comparisons']:
            require(safe_name(cmp.get('name')),'Invalid event comparison name');names.append(cmp['name'])
            group=cmp.get('column','condition')
            require(cmp.get('treatment')!=cmp.get('reference'),'Event treatment and reference must differ')
            for label in ('treatment','reference'):
                require(sum(s.get(group)==cmp[label] for s in samples)>=2,'Event comparisons need >=2 biological replicates in each group')
            # Event models never inherit the general DESeq2/DEJU design silently.
            require(cmp.get('model') in ('unpaired','paired'),'Set event model explicitly to unpaired or paired')
            if cmp['model']=='paired':
                require(cmp.get('subject'),'Paired event comparison requires a subject column')
        require(len(names)==len(set(names)),'Duplicate event comparison names')
    if 'loci' in mods:
        require(cfg.get('loci'),'Supply regions for the loci module')
        require(len({l.get('name') for l in cfg['loci']})==len(cfg['loci']),'Duplicate locus names')
        for l in cfg['loci']:
            require(safe_name(l.get('name')) and l.get('chr') and isinstance(l.get('start'),int) and isinstance(l.get('end'),int)
                    and 1<=l['start']<l['end'] and l['end']-l['start']<=1000000,'Invalid locus coordinates/name')
    if 'isoform_switch' in mods:
        require(cfg.get('isoform_switch',{}).get('config'),'isoform_switch.config is required')
        cfg['isoform_switch']['config']=absolute(cfg['isoform_switch']['config'])
        require(Path(cfg['isoform_switch']['config']).is_file(),'Missing isoform-switch configuration')
    cfg['_samples']=samples
    return cfg


def tool(cfg,name,default):
    value=cfg.get('tools',{}).get(name,default)
    require(isinstance(value,list) and value and all(isinstance(x,str) and x for x in value),f'tools.{name} must be an argument array')
    return value


def event_groups(samples,comparison):
    col=comparison.get('column','condition')
    groups=[[s for s in samples if s.get(col)==comparison[k]] for k in ('treatment','reference')]
    if comparison['model']=='paired':
        subject=comparison['subject']
        subject_groups=[[s.get(subject) for s in g] for g in groups]
        require(all(all(x for x in g) and len(g)==len(set(g)) for g in subject_groups),'Missing/duplicated event subject IDs')
        require(set(subject_groups[0])==set(subject_groups[1]),'Paired event subject sets differ')
        groups=[sorted(g,key=lambda s:s[subject]) for g in groups]
    return groups


class Workflow:
    def __init__(self,cfg):
        self.cfg=cfg;self.out=Path(cfg['out']);self.samples=cfg['_samples'];self.commands=[];self.states={m:'not_requested' for m in sorted(MODULES)}
        for m in cfg['modules']: self.states[m]='pending'

    def command(self,args,label,expected=()):
        self.commands.append({'label':label,'argv':list(map(str,args))})
        log=self.out/'logs'/f'{label}.log';log.parent.mkdir(parents=True,exist_ok=True)
        print(f'Running {label}',flush=True)
        with log.open('w') as fh:
            fh.write(shlex.join(list(map(str,args)))+'\n');fh.flush()
            subprocess.run(list(map(str,args)),stdout=fh,stderr=subprocess.STDOUT,check=True)
        for p in expected:
            require(Path(p).is_file(),f'{label} did not produce {p}; see {log}')

    def report(self):
        self.out.mkdir(parents=True,exist_ok=True)
        payload={'version':VERSION,'modules':self.states,'commands':self.commands,'updated_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())}
        (self.out/'run_manifest.json').write_text(json.dumps(payload,indent=2))
        rows=''.join(f'<tr><td>{html.escape(m)}</td><td>{html.escape(s)}</td></tr>' for m,s in self.states.items())
        summaries=[]
        for title, relative in [('Differential gene expression','dge/summary.tsv'),('Alignment QC','native/alignment_qc.tsv')]:
            source=self.out/relative
            if source.is_file():
                with source.open(newline='') as fh:
                    reader=csv.DictReader(fh,delimiter='\t'); fields=reader.fieldnames or []; records=list(reader)
                cells=''.join('<tr>'+''.join('<td>'+html.escape(str(r.get(k,'')))+'</td>' for k in fields)+'</tr>' for r in records[:30])
                summaries.append('<h2>'+title+'</h2><table><tr>'+''.join('<th>'+html.escape(k)+'</th>' for k in fields)+'</tr>'+cells+'</table>')
        links=[]
        for p in sorted(self.out.rglob('*')):
            if p.is_file() and (p.suffix in ('.tsv','.pdf','.json') or p.name.endswith('.tsv.gz')) and p.name!='report.html':
                relative=p.relative_to(self.out).as_posix()
                links.append(f'<li><a href="{html.escape(relative,quote=True)}">{html.escape(relative)}</a></li>')
        (self.out/'report.html').write_text('''<!doctype html><meta charset="utf-8"><title>Transcriptome analysis report</title>
<style>body{font:16px system-ui;max-width:1000px;margin:40px auto;padding:16px}td,th{padding:8px;border:1px solid #ddd}table{border-collapse:collapse}</style>
<h1>Transcriptome analysis</h1><p>Version '''+VERSION+'''</p><table><tr><th>Module</th><th>Status</th></tr>'''+rows+'''</table>
<p>Read each module's QC before interpreting results. General designs apply to DESeq2/DEJU only. rMATS comparisons use their explicitly selected two-group model.
PSI and relative junction log-fold change are different quantities. Transcript assemblies, novel loci and fusions are candidates requiring review.</p>
<p>DGE: independent gene counts; paired-end fragments. DEJU: aligned-read exon/junction features. External tools use their own documented filtering policies on original BAMs.
Locus plots show raw coverage and junction evidence, not normalized abundance or a PSI estimate.</p>'''+''.join(summaries)+'''<h2>Outputs</h2><ul>'''+''.join(links)+'</ul>')

    def native(self,resolved):
        requested=set(self.cfg['modules'])&{'dge','dju','loci'}
        if not requested: return
        for m in requested:self.states[m]='running'
        self.report()
        try:
            self.command(tool(self.cfg,'Rscript',['Rscript'])+[str(ROOT/'TranscriptomeAnalysis.R'),str(resolved)],'native',
                         [self.out/'native'/'COMPLETE.txt'])
        except Exception:
            for m in requested:self.states[m]='failed_or_incomplete'
            raise
        for m in requested:self.states[m]='completed'

    def events(self):
        cfg=self.cfg;opt=cfg['events']
        for cmp in opt['comparisons']:
            directory=self.out/'events'/cmp['name'];directory.mkdir(parents=True)
            groups=event_groups(self.samples,cmp)
            files=[]
            for n,group in enumerate(groups,1):
                file=directory/f'b{n}.txt';file.write_text(','.join(s['bam'] for s in group)+'\n');files.append(file)
            raw=directory/'raw';tmp=directory/'tmp'
            cmd=tool(cfg,'rmats',['rmats.py'])+['--b1',str(files[0]),'--b2',str(files[1]),'--gtf',cfg['gtf'],
                '-t','paired' if cfg['bam']['paired'] else 'single','--readLength',str(opt['read_length']),
                '--libType',{0:'fr-unstranded',1:'fr-secondstrand',2:'fr-firststrand'}[cfg['bam']['strand']],
                '--nthread',str(cfg.get('threads',1)),'--od',str(raw),'--tmp',str(tmp),'--anchorLength',str(opt.get('min_anchor',8))]
            if opt.get('variable_read_length',False):cmd+=['--variable-read-length']
            if cmp['model']=='paired':cmd+=['--paired-stats']
            self.command(cmd,'rmats_'+cmp['name'])
            _,links=summarize_events(raw,directory,[s['sample_id'] for s in groups[0]],[s['sample_id'] for s in groups[1]],
                cfg.get('fdr',.05),opt.get('delta_filter',.1))
            # Explicit mapping prevents associating unrelated statistical contrasts by name.
            if cmp.get('dju_contrast'):
                source=self.out/'dju'/cmp['dju_contrast']/'junctions.tsv.gz'
                require(source.is_file(),'Event dju_contrast does not have a completed junction result')
                link_dju_events(source,links,directory/'dju_event_links.tsv')
            write_tsv(directory/'sample_order.tsv',
                [dict(sample_id=s['sample_id'],group=g,position=i+1) for g,grp in zip(('treatment','reference'),groups) for i,s in enumerate(grp)],
                ['sample_id','group','position'])

    def isoforms(self):
        cfg=self.cfg;out=self.out/'isoforms';out.mkdir(parents=True)
        stringtie=tool(cfg,'stringtie',['stringtie']);strand={0:[],1:['--fr'],2:['--rf']}[cfg['bam']['strand']]
        assembled=[]
        for s in self.samples:
            p=out/f'{s["sample_id"]}.assembled.gtf';assembled.append(p)
            self.command(stringtie+[s['bam'],'-G',cfg['gtf'],'-p',str(cfg.get('threads',1)),'-o',str(p)]+strand,
                         'assemble_'+s['sample_id'],[p])
        lst=out/'assemblies.txt';lst.write_text('\n'.join(map(str,assembled))+'\n')
        merged=out/'merged.gtf'
        self.command(stringtie+['--merge','-G',cfg['gtf'],'-p',str(cfg.get('threads',1)),'-o',str(merged),str(lst)],'merge_isoforms',[merged])
        prefix=out/'comparison';annotated=out/'comparison.annotated.gtf'
        self.command(tool(cfg,'gffcompare',['gffcompare'])+['-r',cfg['gtf'],'-o',str(prefix),str(merged)],'classify_isoforms',[annotated])
        quantified={}
        for s in self.samples:
            p=out/f'{s["sample_id"]}.quantified.gtf';quantified[s['sample_id']]=p
            self.command(stringtie+[s['bam'],'-e','-G',str(merged),'-p',str(cfg.get('threads',1)),'-o',str(p),
                         '-A',str(out/f'{s["sample_id"]}.gene_abundance.tsv')]+strand,'quantify_'+s['sample_id'],[p])
        summarize_isoforms(quantified,annotated,out)

    def isoform_switch(self):
        config=prepare_switch_config(self.cfg,self.out/'isoform_switch')
        path=self.out/'isoform_switch.config.json';path.write_text(json.dumps(config,indent=2))
        self.command(tool(self.cfg,'Rscript',['Rscript'])+[str(ROOT/'IsoformSwitchAnalysis.R'),str(path)],
                     'isoform_switch',[self.out/'isoform_switch'/'COMPLETE.txt'])

    def fusions(self):
        cfg=self.cfg;out=self.out/'fusions';out.mkdir(parents=True);files={}
        opt=cfg['fusions']
        for s in self.samples:
            p=out/f'{s["sample_id"]}.fusions.tsv';files[s['sample_id']]=p
            cmd=tool(cfg,'arriba',['arriba'])+['-x',s['bam'],'-a',cfg['fasta'],'-g',cfg['gtf'],'-b',opt['blacklist'],
                 '-o',str(p),'-O',str(out/f'{s["sample_id"]}.discarded.tsv')]
            for k,flag in (('known_fusions','-k'),('protein_domains','-p')):
                if k in opt:cmd += [flag,opt[k]]
            self.command(cmd,'fusions_'+s['sample_id'],[p])
        summarize_fusions(files,out/'fusion_candidates.tsv')

    def run(self):
        require(not self.out.exists(),'Output directory exists; choose a new path')
        self.out.mkdir(parents=True)
        (self.out/'INCOMPLETE.txt').write_text('Requested analysis is incomplete. Inspect logs and run_manifest.json.\n')
        cfg={k:v for k,v in self.cfg.items() if not k.startswith('_')}
        resolved=self.out/'resolved_config.json';resolved.write_text(json.dumps(cfg,indent=2))
        shutil.copyfile(cfg['samples'],self.out/'samples.input.tsv')
        hashes=[]
        for key in ('samples','gtf','fasta','gene_counts','dju_counts','dju_features'):
            if key in cfg:
                p=Path(cfg[key]);digest=hashlib.sha256()
                with p.open('rb') as fh:
                    for block in iter(lambda:fh.read(1024*1024),b''):digest.update(block)
                hashes.append(dict(input=key,path=str(p),sha256=digest.hexdigest(),size=p.stat().st_size))
        write_tsv(self.out/'input_checksums.tsv',hashes,['input','path','sha256','size'])
        write_tsv(self.out/'bam_manifest.tsv',[dict(sample_id=s['sample_id'],path=s['bam'],size=Path(s['bam']).stat().st_size,
                  modified=Path(s['bam']).stat().st_mtime) for s in self.samples if s.get('bam') and Path(s['bam']).is_file()],
                  ['sample_id','path','size','modified'])
        try:
            self.native(resolved)
            for m in ('events','isoforms','isoform_switch','fusions'):
                if m not in self.cfg['modules']:continue
                self.states[m]='running';self.report()
                try:getattr(self,m)()
                except Exception:self.states[m]='failed';raise
                self.states[m]='completed';self.report()
            self.report()
            (self.out/'COMPLETE.txt').write_text('All requested modules completed. Review report.html and module QC.\n')
            (self.out/'INCOMPLETE.txt').unlink()
        finally:self.report()


def prepare_switch_config(cfg,out):
    source=Path(cfg['isoform_switch']['config'])
    config=json.loads(source.read_text())
    def absolute(p):
        p=Path(p).expanduser()
        return str((source.parent/p).resolve() if not p.is_absolute() else p.resolve())
    for key in ('samples','gtf','transcript_fasta'):
        if key in config:config[key]=absolute(config[key])
    if config.get('go',{}).get('mapping'):config['go']['mapping']=absolute(config['go']['mapping'])
    require(config.get('samples')==cfg['samples'],'Integrated isoform-switch module must use the same sample sheet')
    if cfg.get('gtf'):require(config.get('gtf')==cfg['gtf'],'Integrated modules must use the same reference GTF')
    config['out']=str(out)
    return config


def preflight(cfg):
    wanted={'Rscript':['Rscript']} if set(cfg['modules'])&{'dge','dju','loci','isoform_switch'} else {}
    if 'events' in cfg['modules']:wanted['rmats']=['rmats.py']
    if 'isoforms' in cfg['modules']:wanted.update(stringtie=['stringtie'],gffcompare=['gffcompare'])
    if 'fusions' in cfg['modules']:wanted['arriba']=['arriba']
    for name,default in wanted.items():
        cmd=tool(cfg,name,default)
        require(shutil.which(cmd[0]) is not None,f'Missing executable for {name}: {cmd[0]}; install it or set tools.{name}')
        for arg in cmd[1:]:
            if arg.endswith('.py'):require(Path(arg).is_file(),f'Missing tool script: {arg}')
    if 'events' in cfg['modules']:
        for c in cfg['events']['comparisons']:event_groups(cfg['_samples'],c)
    if 'isoform_switch' in cfg['modules']:
        with tempfile.TemporaryDirectory(prefix='isoform-switch-check-') as tmp:
            path=Path(tmp)/'config.json'
            path.write_text(json.dumps(prepare_switch_config(cfg,Path(cfg['out'])/'isoform_switch')))
            subprocess.run(tool(cfg,'Rscript',['Rscript'])+[str(ROOT/'IsoformSwitchAnalysis.R'),str(path),'--check'],check=True)
    if set(cfg['modules'])&{'dge','dju','loci'}:
        with tempfile.TemporaryDirectory(prefix='transcriptome-check-') as tmp:
            path=Path(tmp)/'config.json';path.write_text(json.dumps({k:v for k,v in cfg.items() if not k.startswith('_')}))
            subprocess.run(tool(cfg,'Rscript',['Rscript'])+[str(ROOT/'TranscriptomeAnalysis.R'),str(path),'--check'],check=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config',required=True)
    mode=p.add_mutually_exclusive_group();mode.add_argument('--check',action='store_true');mode.add_argument('--plan',action='store_true')
    args=p.parse_args();cfg=load_config(args.config)
    if args.plan:
        print(json.dumps({k:v for k,v in cfg.items() if not k.startswith('_')},indent=2))
        print('Plan only. Executables, R packages and biological results have not been validated.')
        return
    preflight(cfg)
    if args.check:print('Preflight passed; no analysis outputs created.');return
    Workflow(cfg).run();print(f'Completed. Report: {cfg["out"]}/report.html')


if __name__=='__main__':
    try:main()
    except (ValueError,OSError,subprocess.CalledProcessError) as e:
        print(f'ERROR: {e}',file=sys.stderr);sys.exit(1)
