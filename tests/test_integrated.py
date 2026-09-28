"""Offline contract/regression tests. These do not substitute for external caller validation."""
import csv
import json
from pathlib import Path
import sys
import tempfile
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from modules.external import (EVENTS,write_tsv,read_tsv,bh,summarize_events,event_junctions,
                              summarize_isoforms,summarize_fusions,link_dju_events)
from TranscriptomePipeline import load_config,event_groups,Workflow,prepare_switch_config


def fixture_events(path):
    row=dict(ID='1',GeneID='G1',chr='chr1',strand='+',PValue='.01',FDR='.05',
             IncLevel1='.8,.8',IncLevel2='.2,.2',IncLevelDifference='.6',
             IncFormLen='100',SkipFormLen='100',IJC_SAMPLE_1='80,80',IJC_SAMPLE_2='20,20',
             SJC_SAMPLE_1='20,20',SJC_SAMPLE_2='80,80',upstreamEE='100',downstreamES='400',
             exonStart_0base='200',exonEnd='250',**{'1stExonStart_0base':'200','1stExonEnd':'250',
             '2ndExonStart_0base':'300','2ndExonEnd':'350','longExonStart_0base':'200','longExonEnd':'270',
             'shortES':'200','shortEE':'250','flankingES':'400','flankingEE':'450'})
    for kind in EVENTS:
        for mode in ('JC','JCEC'):
            write_tsv(path/f'{kind}.MATS.{mode}.txt',[row],list(row))


class Contracts(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.p=Path(self.temp.name)
    def tearDown(self):
        self.temp.cleanup()
    def test_isoform_switch_integration_paths_and_sample_identity(self):
        sheet=self.p.resolve()/'samples.tsv'
        write_tsv(sheet,[dict(sample_id='A',condition='control',quant='A/quant.sf')],['sample_id','condition','quant'])
        switch=self.p/'switch.json'
        switch.write_text(json.dumps(dict(samples='samples.tsv',gtf='reference.gtf',out='unused')))
        cfg=dict(samples=str(sheet),out=str(self.p/'out'),isoform_switch=dict(config=str(switch)),modules=['isoform_switch'])
        prepared=prepare_switch_config(cfg,self.p/'out/isoform_switch')
        self.assertEqual(prepared['samples'],str(sheet))
        self.assertEqual(prepared['out'],str(self.p/'out/isoform_switch'))
        self.assertEqual(prepared['gtf'],str(self.p.resolve()/'reference.gtf'))
        cfg['samples']=str(self.p/'different.tsv')
        with self.assertRaisesRegex(ValueError,'same sample sheet'):prepare_switch_config(cfg,self.p/'out')

    def test_bh(self):
        self.assertEqual(bh([.01,.04,.03,None]),[.03,.04,.04,None])
    def test_events_psi_direction_and_coordinates(self):
        fixture_events(self.p)
        rows,links=summarize_events(self.p,self.p/'out',['K1','K2'],['W1','W2'])
        self.assertEqual(len(rows),10)
        self.assertTrue(all(abs(r['delta_PSI']-.6)<1e-8 for r in rows))
        self.assertEqual(set(r['event_type'] for r in rows),set(EVENTS.values()))
        self.assertEqual(len(read_tsv(self.p/'out/PSI_by_sample.tsv')),40)
        se=[r for r in links if r['event_id']=='SE:1']
        self.assertEqual([(r['left'],r['right']) for r in se],[(100,201),(250,401),(100,401)])
        self.assertEqual([r['role'] for r in links if r['event_type']=='intron_retention'],['spliced_exclusion'])
        dju=[dict(feature_id='J1',gene_id='G1',chr='chr1',strand='+',left=100,right=401,junction_FDR=.01),
             dict(feature_id='J2',gene_id='G1',chr='chr1',strand='+',left=500,right=600,junction_FDR=.2)]
        write_tsv(self.p/'dju.tsv',dju,list(dju[0]));link_dju_events(self.p/'dju.tsv',links,self.p/'linked.tsv')
        self.assertTrue(any(r['event_type']=='exon_skipping' for r in read_tsv(self.p/'linked.tsv')))
        self.assertTrue(any(r['event_type']=='unclassified' for r in read_tsv(self.p/'linked.tsv')))
    def test_reject_invalid_psi_or_missing_results(self):
        fixture_events(self.p)
        file=self.p/'SE.MATS.JC.txt';file.write_text(file.read_text().replace('.8,.8','.1,.1'))
        with self.assertRaisesRegex(ValueError,'PSI disagrees'):summarize_events(self.p,self.p/'out',['K1','K2'],['W1','W2'])
        file.unlink()
        with self.assertRaisesRegex(ValueError,'Missing rMATS'):summarize_events(self.p,self.p/'out',['K1','K2'],['W1','W2'])
    def test_paired_subject_order(self):
        s=[dict(sample_id='K2',condition='KO',subject='2'),dict(sample_id='W1',condition='WT',subject='1'),
           dict(sample_id='K1',condition='KO',subject='1'),dict(sample_id='W2',condition='WT',subject='2')]
        c=dict(treatment='KO',reference='WT',model='paired',subject='subject')
        g=event_groups(s,c)
        self.assertEqual([[r['sample_id'] for r in a] for a in g],[['K1','K2'],['W1','W2']])
        s[0]['subject']='3'
        with self.assertRaisesRegex(ValueError,'sets differ'):event_groups(s,c)
    def test_assembly_candidates_and_missing_quantification(self):
        template='chr1\tStringTie\ttranscript\t1\t100\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; TPM "12"; class_code "u";\n'
        a=self.p/'a.gtf';a.write_text(template);b=self.p/'b.gtf';b.write_text('')
        summarize_isoforms({'A':a,'B':b},a,self.p/'out')
        self.assertEqual(read_tsv(self.p/'out/transcript_TPM.tsv')[0]['B'],'NA')
        self.assertEqual(read_tsv(self.p/'out/novel_gene_candidates.tsv')[0]['candidate_type'],'intergenic_transcript_novel_gene_candidate')
    def test_fusion_schema(self):
        row={'#gene1':'A','gene2':'B','breakpoint1':'chr1:100','breakpoint2':'chr2:200','confidence':'high'}
        p=self.p/'f.tsv';write_tsv(p,[row],list(row))
        rows=summarize_fusions({'S':p},self.p/'summary.tsv')
        self.assertEqual(rows[0]['interpretation'],'fusion_candidate_requires_review')
        p.write_text('invalid\n')
        with self.assertRaises(ValueError):summarize_fusions({'S':p},self.p/'summary.tsv')
    def test_html_report(self):
        cfg={'out':str(self.p/'out'),'modules':['dge'],'_samples':[]}
        w=Workflow(cfg);w.states['dge']='completed'
        write_tsv(w.out/'dge/summary.tsv',[{'contrast':'A_vs_B','significant':5}],['contrast','significant'])
        (w.out/'dge/results.tsv.gz').touch()
        w.report();report=(w.out/'report.html').read_text()
        self.assertIn('A_vs_B',report);self.assertIn('dge/results.tsv.gz',report)
        self.assertEqual(json.loads((w.out/'run_manifest.json').read_text())['modules']['dge'],'completed')

    def test_config_and_tool_contracts(self):
        samples=[]
        for sid,condition in [('K1','KO'),('K2','KO'),('W1','WT'),('W2','WT')]:
            p=self.p/(sid+'.bam');p.touch();samples.append(dict(sample_id=sid,condition=condition,bam=str(p)))
        write_tsv(self.p/'s.tsv',samples,['sample_id','condition','bam']);(self.p/'a.gtf').touch()
        cfg=dict(samples='s.tsv',gtf='a.gtf',out='out',modules=['events'],bam=dict(paired=True,strand=2),
                 events=dict(read_length=100,comparisons=[dict(name='KO_WT',treatment='KO',reference='WT',model='unpaired')]))
        path=self.p/'config.json';path.write_text(json.dumps(cfg));resolved=load_config(path)
        workflow=Workflow(resolved);workflow.out.mkdir()
        seen=[]
        def stub(args,label,expected=()):
            # Fake only the executable boundary; validates actual adapter arguments/parsing.
            seen.append(args);fixture_events(Path(args[args.index('--od')+1]))
        workflow.command=stub;workflow.events()
        cmd=seen[0]
        self.assertEqual(cmd[cmd.index('--libType')+1],'fr-firststrand')
        self.assertIn('K1.bam',Path(cmd[cmd.index('--b1')+1]).read_text())
        self.assertIn('W1.bam',Path(cmd[cmd.index('--b2')+1]).read_text())
        self.assertTrue((workflow.out/'events/KO_WT/events.tsv').is_file())
        cfg['modules']=['fusions'];cfg['fusions']={};path.write_text(json.dumps(cfg))
        with self.assertRaisesRegex(ValueError,'STAR BAMs'):load_config(path)


if __name__=='__main__':unittest.main()
