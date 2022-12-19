#!/usr/bin/env python
from __future__ import print_function
import os
import sys
import argparse
import getpass
try:
    from tqdm.auto import tqdm
    hasTQDM = True
except ImportError:
    hasTQDM = False

from dataset_allowed_definitions import (get_allowed_resonances, 
                                         get_allowed_eras,
                                         get_allowed_sub_eras)
from config import Configuration


# parallel processing
def _futures_handler(futures_set, status=True, unit='items', desc='Processing',
                     add_fn=None, output=None):
    try:
        from tqdm.auto import tqdm
        hasTQDM = True
    except ImportError:
        hasTQDM = False
    import time

    def _handle(pbar=None):
        while len(futures_set) > 0:
            finished = set(job for job in futures_set if job.done())
            futures_set.difference_update(finished)
            while finished:
                res = finished.pop().result()
                if add_fn:
                    add_fn(output, res)
                if pbar is not None:
                    pbar.update(1)
            time.sleep(0.5)

    try:
        if hasTQDM:
            with tqdm(disable=not status, unit=unit, total=len(futures_set),
                      desc=desc, ncols=80) as pbar:
                _handle(pbar)
        else:
            _handle()
    except KeyboardInterrupt:
        for job in futures_set:
            job.cancel()
    except Exception:
        for job in futures_set:
            job.cancel()
        raise


# argparse common functions
def add_common_multi(parser):
    parser.add_argument('--workers', '-j', type=int, default=1,
                        help='Number of cores')
    parser.add_argument('--dryrun', action='store_true',
                        help='Don\'t run, just print number of jobs')
    parser.add_argument('--condor', action='store_true',
                        help='Prepare condor submit script')
    parser.add_argument('--condorTag', type=str, default='',
                        help='Tag for condor submit script')
    parser.add_argument('--jobsPerSubmit', '-nj', type=int, default=1,
                        help='Number of jobs to run per submit')


def add_common_flatten(parser):
    parser.add_argument('--numerator', nargs='*',
                        help='Filter by numerator')
    parser.add_argument('--denominator', nargs='*',
                        help='Filter by denominator')
    parser.add_argument('--shiftType', nargs='*',
                        help='Filter by shift type')
    parser.add_argument('--dataOnly', action='store_true',
                        help='Only flatten data')
    parser.add_argument('--bySubEraAlso', action='store_true',
                        help='Flatten each sub-era')


def add_common_fit(parser):
    parser.add_argument('--numerator', nargs='*',
                        help='Filter by numerator')
    parser.add_argument('--denominator', nargs='*',
                        help='Filter by denominator')
    parser.add_argument('--fitType', nargs='*',
                        help='Filter by fit type')
    parser.add_argument('--shiftType', nargs='*',
                        help='Filter by shift type')
    parser.add_argument('--sampleType', nargs='*',
                        help='Filter by sample type (data, mc)')
    parser.add_argument('--efficiencyBin', nargs='*',
                        help='Filter by efficiency bin')
    parser.add_argument('--recover', action='store_true',
                        help='Auto recover failed fits')


def add_common_prepare(parser):
    parser.add_argument('--numerator', nargs='*',
                        help='Filter by numerator')
    parser.add_argument('--denominator', nargs='*',
                        help='Filter by denominator')
    parser.add_argument('--skipPlots', action='store_true',
                        help='Skip efficiency plots')
    parser.add_argument('--cutAndCount', action='store_true',
                        help='Use cut and count rather than fits')

def add_common_prepare_corrected(parser):
    parser.add_argument('--numerator', nargs='*',
                        help='Filter by numerator')
    parser.add_argument('--denominator', nargs='*',
                        help='Filter by denominator')
    parser.add_argument('--skipPlots', action='store_true',
                        help='Skip efficiency plots')
    parser.add_argument('--cutAndCount', action='store_true',
                        help='Use cut and count rather than fits')

def add_common_compare(parser):
    parser.add_argument('--subera1', default='',
                        help='Introduce first subera')
    parser.add_argument('--subera2', default='',
                        help='Introduce second subera')
    parser.add_argument('--era2', default='',
                        help='Introduce de second era, if needed')
    parser.add_argument('--condor_submit', action='store_true',
                        help='Prepare condor submit script')

def add_common_particle(parser):
    parser.add_argument('particle', choices=['muon', 'electron'],
                        help='Particle for scalefactors')


def add_common_probe(parser):
    parser.add_argument('probe', choices=['generalTracks', 'standAloneMuons', 'dSAMuons'],
                        help='Probe for scalefactors')


def add_common_resonance(parser):
    allowed = sorted(get_allowed_resonances())
    parser.add_argument('resonance', choices=allowed,
                        help='Resonance for scalefactors')


def add_common_era(parser):
    a = get_allowed_resonances()
    allowed = []
    for r in a:
        allowed += get_allowed_eras(r)
    allowed = sorted(set(allowed))
    parser.add_argument('era', choices=allowed,
                        help='Scale factor set to produce')


def add_common_data_tier(parser, isPositional=True):
    allowed = ['AOD', 'MINIAOD']
    if isPositional:
        parser.add_argument('dataTier', choices=allowed,
                            help='Data-tier of underlying ntuples')
    else:
        parser.add_argument('--dataTier', choices=allowed,
                            help='Data-tier of underlying ntuples')


def add_common_sub_era(parser):
    a_rs = get_allowed_resonances()
    a_sub = []
    for r in a_rs:
        a_es = get_allowed_eras(r)
        for e in a_es:
            a_sub += get_allowed_sub_eras(r, e)
    a_sub = sorted(set(a_sub))
    
    parser.add_argument('subEra', nargs='?', choices=a_sub,
                        help='Sub-era if desired')


def add_common_config(parser):
    parser.add_argument('config',
                        help='Efficiency configuration file')


def add_common_options(parser):
    parser.add_argument('--baseDir', default='',
                        help='Working directory')
    parser.add_argument('--usePOGSpace', action='store_true',
                        help='Use POG central space instead of user')
    parser.add_argument('--useLocalSpark', action='store_true',
                        help='Use local spark master for tests')
    parser.add_argument('--registry', default=None,
                        help='Specify data registry')
    parser.add_argument('--ntupleVer', default=None,
                        help='Specify ntuple version defined in data registry')


def parse_command_line(argv):
    parser = argparse.ArgumentParser(description='TnP Fitter')

    subparsers = parser.add_subparsers(help='Fitting step', dest='command')
    subparsers.required = True

    parser_convert = subparsers.add_parser(
        'convert',
        help='Convert ROOT to parquet',
    )
    add_common_particle(parser_convert)
    add_common_resonance(parser_convert)
    add_common_era(parser_convert)
    if argv[0] == 'convert':
        add_common_data_tier(parser_convert, isPositional=True)
    add_common_sub_era(parser_convert)
    add_common_options(parser_convert)

    parser_flatten = subparsers.add_parser(
        'flatten',
        help='Flatten to histograms',
    )
    add_common_particle(parser_flatten)
    add_common_probe(parser_flatten)
    add_common_resonance(parser_flatten)
    add_common_era(parser_flatten)
    add_common_config(parser_flatten)
    if argv[0] != 'convert':
        add_common_data_tier(parser_flatten, isPositional=False)
    add_common_options(parser_flatten)
    add_common_flatten(parser_flatten)

    parser_fit = subparsers.add_parser(
        'fit',
        help='Fit pass/fail histograms',
    )
    add_common_particle(parser_fit)
    add_common_probe(parser_fit)
    add_common_resonance(parser_fit)
    add_common_era(parser_fit)
    add_common_config(parser_fit)
    add_common_options(parser_fit)
    add_common_multi(parser_fit)
    add_common_fit(parser_fit)

    parser_prepare = subparsers.add_parser(
        'prepare',
        help='Prepare efficiencies',
    )
    add_common_particle(parser_prepare)
    add_common_probe(parser_prepare)
    add_common_resonance(parser_prepare)
    add_common_era(parser_prepare)
    add_common_config(parser_prepare)
    if argv[0] != 'convert':
        add_common_data_tier(parser_prepare, isPositional=False)
    add_common_options(parser_prepare)
    add_common_multi(parser_prepare)
    add_common_prepare(parser_prepare)

    parser_prepare_corrected = subparsers.add_parser(
        'prepare_corrected',
        help='Prepare efficiencies',
    )
    add_common_particle(parser_prepare_corrected)
    add_common_probe(parser_prepare_corrected)
    add_common_resonance(parser_prepare_corrected)
    add_common_era(parser_prepare_corrected)
    add_common_config(parser_prepare_corrected)
    if argv[0] != 'convert':
        add_common_data_tier(parser_prepare_corrected, isPositional=False)
    add_common_options(parser_prepare_corrected)
    add_common_multi(parser_prepare_corrected)
    add_common_prepare_corrected(parser_prepare_corrected)

    parser_compare = subparsers.add_parser(
        'compare',
        help='Compare distributions',
    )
    add_common_particle(parser_compare)
    add_common_probe(parser_compare)
    add_common_resonance(parser_compare)
    add_common_era(parser_compare)
    add_common_config(parser_compare)
    add_common_options(parser_compare)
    add_common_multi(parser_compare)
    add_common_compare(parser_compare)

    return parser.parse_args(argv)


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    job_fn = None
    unit = 'unit'
    desc = 'Processing'
    add_fn = None
    output = None

    if args.baseDir:
        baseDir = args.baseDir
    elif args.particle == 'muon':
        baseDir = os.path.join(
            '/eos/cms/store/group/phys_muon',
            f'{getpass.getuser()}/TagAndProbe',
        )
    else:
        baseDir = os.path.join(
            '/eos/cms/store/user',
            f'{getpass.getuser()}/TagAndProbe/{args.particle}',
        )

    if args.command == 'convert':
        from converter import run_all
        run_all(args.particle, args.resonance, args.era, args.dataTier, args.subEra, '', args.baseDir, args.usePOGSpace, args.useLocalSpark)
        return 0
    elif args.command == 'flatten':
        from flattener import run_spark
        run_spark(
            args.particle, args.probe, args.resonance, args.era,
            Configuration(args.config),
            numerator=args.numerator,
            denominator=args.denominator,
            shiftType=args.shiftType,
            baseDir=baseDir,
            dataOnly=args.dataOnly,
            bySubEraAlso=args.bySubEraAlso,
            useLocalSpark=args.useLocalSpark,
            registry=args.registry,
            dataTier=args.dataTier,
            ntupleVer=args.ntupleVer,
        )
        return 0
    elif args.command == 'fit':
        from fitter import run_single_fit, build_fit_jobs, build_condor_submit
        job_fn = run_single_fit
        jobs = build_fit_jobs(
            args.particle, args.probe, args.resonance, args.era,
            Configuration(args.config),
            baseDir=baseDir,
            numerator=args.numerator,
            denominator=args.denominator,
            fitType=args.fitType,
            sampleType=args.sampleType,
            shiftType=args.shiftType,
            efficiencyBin=args.efficiencyBin,
            recover=args.recover,
        )
        unit = 'fit'
        desc = 'Fitting'
    elif args.command == 'prepare':
        from prepare import prepare, build_prepare_jobs
        job_fn = prepare
        jobs = build_prepare_jobs(
            args.particle,
            args.probe,
            args.resonance,
            args.era,
            Configuration(args.config),
            numerator=args.numerator,
            denominator=args.denominator,
            baseDir=baseDir,
            registry=args.registry,
            dataTier=args.dataTier,
            ntupleVer=args.ntupleVer,
        )
        jobs = [job + [args.skipPlots, args.cutAndCount] for job in jobs]
        unit = 'efficiency'
        desc = 'Preparing'
    elif args.command == 'prepare_corrected':
        from prepare_corrected import prepare_corrected, build_prepare_corrected_jobs
        job_fn = prepare_corrected
        jobs = build_prepare_corrected_jobs(
            args.particle,
            args.probe,
            args.resonance,
            args.era,
            Configuration(args.config),
            numerator=args.numerator,
            denominator=args.denominator,
            baseDir=baseDir,
            registry=args.registry,
            dataTier=args.dataTier,
            ntupleVer=args.ntupleVer,
        )
        jobs = [job + [args.skipPlots, args.cutAndCount] for job in jobs]
        unit = 'efficiency'
        desc = 'Preparing'

    elif args.command == 'compare':
        from compare import compare
        if args.condor_submit == True:
            from run_multiple_compare import compare_multiple
            compare_multiple(
                args.particle,
                args.probe,
                args.resonance,
                args.era,
                args.config,
                subera1=args.subera1,
                subera2=args.subera2,
                era2=args.era2,
                baseDir=baseDir,
            )
        else:
            compare(
                args.particle,
                args.probe,
                args.resonance,
                args.era,
                Configuration(args.config),
                subera1=args.subera1,
                subera2=args.subera2,
                era2=args.era2,
                baseDir=baseDir,
            )

        return 0
        
    if args.dryrun:
        print('Will run {} {} jobs'.format(len(jobs), args.command))
    elif args.condor:
        test = False
        submit_dir = ''
        joblist = os.path.join(
            submit_dir,
            '{}joblist_{}_{}_{}_{}{}.txt'.format(
                'test_' if test else '',
                args.particle,
                args.probe,
                args.resonance,
                args.era,
                '_'+args.condorTag if args.condorTag != '' else ''
            )
        )
        config = build_condor_submit(joblist,
                                     test=test,
                                     jobsPerSubmit=args.jobsPerSubmit,
                                     njobs=len(jobs))
        if test:
            os.makedirs('condor', exist_ok=True)
        configpath = os.path.join(
            submit_dir,
            '{}condor_{}_{}_{}_{}{}.sub'.format(
                'test_' if test else '',
                args.particle,
                args.probe,
                args.resonance,
                args.era,
                '_'+args.condorTag if args.condorTag != '' else ''
            )
        )
        with open(configpath, 'w') as f:
            f.write(config)
        with open(joblist, 'w') as f:
            for job in jobs:
                f.write(','.join([str(j) for j in job])+'\n')
        print('Condor submit script written to {}'.format(configpath))
        print('To submit:')
        print('    condor_submit {}'.format(configpath))
    elif args.workers > 1:
        import concurrent.futures
        with concurrent.futures.ProcessPoolExecutor(args.workers) as executor:
            futures = set(executor.submit(job_fn, *job) for job in jobs)
            _futures_handler(futures, status=True, unit=unit, desc=desc,
                             add_fn=add_fn, output=output)
    else:
        if hasTQDM:
            for job in tqdm(jobs, ncols=80, unit=unit, desc=desc):
                result = job_fn(*job)
                if add_fn is not None:
                    add_fn(output, result)
        else:
            for job in jobs:
                print(f'Submitting job: {job}')
                result = job_fn(*job)
                if add_fn is not None:
                    add_fn(output, result)


if __name__ == "__main__":
    status = main()
    sys.exit(status)
