"""Check the miniprot-predicted models as a separate, non-independent evidence tier."""

import json
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from analyze_models import WORKERS, analyze, attach_junctions
from exact_position import OUT, read_json_gz, sha, write_json_gz


def main():
    # inputs: focal junctions, predicted models, the annotation review
    evidence = read_json_gz('evidence.json.gz')
    predictions = read_json_gz('predictions.json.gz')
    certificate = json.loads((OUT / 'prediction_extraction_complete.json').read_text())
    assert sha(OUT / 'predictions.json.gz') == certificate['output_sha256']
    annotation_review = {row['key']: row for row in read_json_gz('model_review.json.gz')['rows']}

    # the annotation analysis, run on the predicted models
    attach_junctions(predictions['rows'], evidence['loci'])
    with ProcessPoolExecutor(max_workers=WORKERS) as pool:
        results = list(pool.map(analyze, predictions['rows'], chunksize=4))

    # relabel phase sources as miniprot predictions
    for result in results:
        result['evidence_type'] = 'miniprot_prediction_not_independent_annotation'
        result['original_annotation_status'] = annotation_review[result['key']]['status']
        result['original_annotated_phase'] = False
        for model in result['models']:
            if model.get('phase_source') == 'annotation':
                model['phase_source'] = 'miniprot_GFF_prediction'
    write_json_gz('prediction_review.json.gz', dict(rows=results))

    # prediction status cross-tabulated against the annotation status
    cross = Counter((r['original_annotation_status'], r['status']) for r in results)
    inputs = [
        OUT / 'predictions.json.gz',
        OUT / 'prediction_extraction_complete.json',
        Path(__file__).resolve().parent / 'analyze_models.py',
        OUT / 'model_review.json.gz',
    ]
    summary = dict(
        complete=True,
        pairs=len(results),
        statuses=dict(Counter(r['status'] for r in results)),
        cross_tab=[
            dict(annotation=annotation, prediction=prediction, n=n)
            for (annotation, prediction), n in sorted(cross.items())
        ],
        source_sha256={str(path): sha(path) for path in inputs},
        output_sha256=sha(OUT / 'prediction_review.json.gz'),
        script_sha256=sha(__file__),
    )
    (OUT / 'prediction_review_complete.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
