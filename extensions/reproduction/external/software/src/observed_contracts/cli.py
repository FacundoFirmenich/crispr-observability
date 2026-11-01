import argparse, gzip, json, sys
from pathlib import Path
from .core import factorization, compile_partition, require_scope, PROPERTIES, TARGETS

def main():
    if len(sys.argv)>1 and sys.argv[1]=='compile':
        return compile_main(sys.argv[2:])
    p = argparse.ArgumentParser(description='Audit whether a declared target factors through a representation on an actual finite coding catalog.')
    p.add_argument('catalog', type=Path)
    p.add_argument('--operator', required=True)
    p.add_argument('--run', default='ALL')
    p.add_argument('--representation', choices=PROPERTIES, required=True)
    p.add_argument('--target', choices=TARGETS, required=True)
    p.add_argument('--exact-coding', required=True)
    p.add_argument('--require-scope', choices=['finite_observed_catalog','full_allele','biological_function','future_sufficiency'])
    a = p.parse_args()
    data = json.loads(a.catalog.read_text())
    rows = data['catalogs'][a.operator][a.run]
    result = factorization(rows, a.representation, a.target, a.exact_coding)
    if a.require_scope:
        try: require_scope(a.require_scope, result)
        except ValueError as exc:
            result['claim_authorization'] = 'REJECTED'
            result['reason'] = str(exc)
            print(json.dumps(result, ensure_ascii=False, indent=2))
            raise SystemExit(2)
        result['claim_authorization'] = 'ACCEPTED_WITHIN_FINITE_CATALOG'
    print(json.dumps(result, ensure_ascii=False, indent=2))

def compile_main(argv):
    p = argparse.ArgumentParser(description='Compile the coarsest observed partition preserving a family of declared targets.')
    p.add_argument('catalog',type=Path)
    p.add_argument('--operator',required=True)
    p.add_argument('--run',default='ALL')
    p.add_argument('--targets',nargs='+',choices=TARGETS,required=True)
    p.add_argument('--exact-coding',required=True)
    a = p.parse_args(argv)
    data = json.loads(a.catalog.read_text())
    result = compile_partition(data['catalogs'][a.operator][a.run],a.targets,a.exact_coding)
    print(json.dumps(result,ensure_ascii=False,indent=2))

if __name__ == '__main__': main()
