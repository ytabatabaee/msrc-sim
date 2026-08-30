from __future__ import annotations
from pathlib import Path
import csv,json,yaml
from dataclasses import asdict

from .frequency_history_plot import A0_COLOR, A1_COLOR, plot_frequency_history_tree

def write_outputs(config,tree,rearrangement,history,sampled,results):
    out=Path(config['output']['directory']); out.mkdir(parents=True,exist_ok=True)
    frequency_history_path = out / 'frequency_history.csv'
    sampled_arrangements_path = out / 'sampled_arrangements.csv'
    if config['output']['record_resolved_config']:
        with open(out/'config.resolved.yaml','w') as f: yaml.safe_dump(config,f,sort_keys=False)
    if config['output']['record_frequency_history']:
        rows=[asdict(r) for r in history.records]
        with open(frequency_history_path,'w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)
    if config['output']['record_sampled_arrangements']:
        with open(sampled_arrangements_path,'w',newline='') as f:
            w=csv.writer(f); w.writerow(['taxon','arrangement','terminal_frequency_A1']);
            for t,z in sampled.items(): w.writerow([t,z,history.terminal_frequency(t)])
    if config['output'].get('make_history_plot', False):
        if not config['output'].get('record_frequency_history', True):
            raise ValueError('output.make_history_plot requires output.record_frequency_history')
        plot_cfg = config['output'].get('history_plot', {}) or {}
        filename = plot_cfg.get('filename', 'wright_fisher_history.png')
        figure_path = out / filename
        fmt = plot_cfg.get('format') or figure_path.suffix.lower().lstrip('.') or 'png'
        if fmt not in {'png', 'pdf', 'svg'}:
            raise ValueError("output.history_plot.format must be one of 'png', 'pdf', or 'svg'")
        if figure_path.suffix.lower().lstrip('.') not in {'png', 'pdf', 'svg'}:
            figure_path = figure_path.with_suffix(f'.{fmt}')
        figure_path.parent.mkdir(parents=True, exist_ok=True)
        records = [asdict(r) for r in history.records]
        sampled_rows = {str(t): {'arrangement': str(z), 'terminal_frequency_A1': str(history.terminal_frequency(t))} for t, z in sampled.items()}
        tip_order = plot_cfg.get('tip_order')
        fig, _ax = plot_frequency_history_tree(
            records,
            sampled_arrangements=sampled_rows,
            tip_order=tip_order,
            glyphs_per_row=int(plot_cfg.get('glyphs_per_row', 12)),
            max_rows_per_branch=None if plot_cfg.get('all_generations', False) else int(plot_cfg.get('max_rows_per_branch', 30)),
            width_mode=plot_cfg.get('width_mode', 'constant'),
            show_internal_labels=bool(plot_cfg.get('show_internal_labels', True)),
            show_frequency_trace=bool(plot_cfg.get('show_frequency_trace', False)),
            title=plot_cfg.get('title'),
        )
        fig.savefig(figure_path, format=fmt, dpi=int(plot_cfg.get('dpi', 300)), bbox_inches='tight')
        metadata = getattr(fig, '_msrc_history_metadata', {})
        import matplotlib.pyplot as plt
        plt.close(fig)
        origin = next((r for r in records if r.get('is_origin')), None)
        manifest = {
            'frequency_history': str(frequency_history_path),
            'config': str(out / 'config.resolved.yaml') if config['output'].get('record_resolved_config', True) else None,
            'sampled_arrangements': str(sampled_arrangements_path) if config['output'].get('record_sampled_arrangements', True) else None,
            'figure': str(figure_path),
            'num_records_total': len(records),
            'num_records_displayed': metadata.get('num_records_displayed'),
            'glyphs_per_row': int(plot_cfg.get('glyphs_per_row', 12)),
            'max_rows_per_branch': None if plot_cfg.get('all_generations', False) else int(plot_cfg.get('max_rows_per_branch', 30)),
            'width_mode': plot_cfg.get('width_mode', 'constant'),
            'tip_order': metadata.get('tip_order', tip_order),
            'origin_branch': str(origin['branch_id']) if origin else None,
            'origin_age': float(origin['absolute_age']) if origin else None,
            'arrangement_colors': {'A0': A0_COLOR, 'A1': A1_COLOR},
        }
        figure_path.with_name(f'{figure_path.stem}.figure.json').write_text(json.dumps(manifest, indent=2))
    if config['output']['record_gene_trees']:
        with open(out/'true_gene_trees.nwk','w') as f:
            for r in results: f.write(r.newick+'\n')
    with open(out/'coalescence_times.csv','w',newline='') as f:
        w=csv.writer(f); w.writerow(['locus_id','coalescence_index','time'])
        for r in results:
            for i,t in enumerate(r.coalescence_times): w.writerow([r.locus_id,i,t])
    with open(out/'coalescence_events.csv','w',newline='') as f:
        rows=[asdict(x) for r in results for x in r.coalescences]
        if rows: w=csv.DictWriter(f,fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)
    if config['output']['record_backward_events']:
        rows=[asdict(x) for r in results for x in r.events]
        if rows:
            with open(out/'genealogy_events.csv','w',newline='') as f: w=csv.DictWriter(f,fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)
    counts=[sum(r.topology_index==i for r in results) for i in range(3)]
    with open(out/'summary.json','w') as f: json.dump({'version':'0.3.0','taxa':tree.taxa,'sampled_arrangements':sampled,'topology_counts':counts,'topology_frequencies':[c/len(results) for c in counts]},f,indent=2)
    return out
