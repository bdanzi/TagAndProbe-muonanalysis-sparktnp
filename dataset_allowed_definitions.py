
# allowed choices
def get_allowed_resonances():
    
    resonances = [
        'Z',
        'JPsi'
    ]
    return resonances


def get_allowed_eras(resonance):

    eras = {
        'Z': [
            # ultra legacy
            'Run2016_UL_HIPM',
            'Run2016_UL',
            'Run2017_UL',
            'Run2018_UL',
            # Double muon PD
            'Run2016_UL_HIPM_DM',
            'Run2016_UL_DM',
            'Run2017_UL_DM',
            'Run2018_UL_DM',
            # rereco (i.e. legacy)
            'Run2016',
            'Run2017',
            'Run2018',
            'Run2022'
        ],
        'JPsi': [
            # heavy ion
            'Run2016_HI_pPb_8TeV',
            # ultra legacy
            'Run2016_UL_HIPM',
            'Run2016_UL',
            'Run2017_UL',
            'Run2018_UL',
            # rereco (i.e. legacy)
            'Run2016',
            'Run2017',
            'Run2018',
            'Run2022'
        ],
    }
    return eras.get(resonance, [])


def get_allowed_sub_eras(resonance, era):
    
    subEras = {
        'Z': {
            # ultra legacy
            'Run2016_UL_HIPM': ['Run2016_UL_HIPM'] + [
                f'Run2016{b}' for b in 'BCDEF']+['DY_madgraph'],
            'Run2016_UL': ['Run2016_UL'] + [
                f'Run2016{b}' for b in 'FGH']+['DY_madgraph'],
            'Run2017_UL': ['Run2017_UL'] + [
                f'Run2017{b}' for b in 'BCDEF']+['DY_madgraph'],
            'Run2018_UL': ['Run2018_UL'] + [
                #f'Run2018{b}' for b in 'ABCD']+['DY_madgraph'],
                f'Run2018{b}' for b in 'ABCD']+['DY_madgraph'],
                #f'Run2018{b}' for b in 'B'],#+['DY_madgraph', 'DY_powheg'],
            # Double muon PD
            'Run2016_UL_HIPM_DM': ['Run2016_UL_HIPM_DM'] + [
                f'Run2016{b}' for b in 'BCDEF']+['DY_madgraph'],
            'Run2016_UL_DM': ['Run2016_UL_DM'] + [
                f'Run2016{b}' for b in 'FGH']+['DY_madgraph'],
            'Run2017_UL_DM': ['Run2017_UL_DM'] + [
                f'Run2017{b}' for b in 'BCDEF']+['DY_madgraph'],
            'Run2018_UL_DM': ['Run2018_UL_DM'] + [
                f'Run2018{b}' for b in 'ABCD']+['DY_madgraph'],
            # rereco (i.e. legacy)
            'Run2016': ['Run2016'] + [
               f'Run2016{b}' for b in 'BCDEFGH']+['DY_madgraph'],
            'Run2017': ['Run2017'] + [
               f'Run2017{b}' for b in 'BCDEF']+['DY_madgraph'],
            'Run2018': ['Run2018'] + [
               f'Run2018{b}' for b in 'ABCD']+['DY_madgraph'],
            'Run2022': ['Run2022'] + [
               f'Run2022{b}' for b in 'BCD']+['DY_madgraph']
               #f'Run2022{b}' for b in 'D']+['DY_madgraph']
        },
        'JPsi': {
            # ultra legacy
            'Run2016_UL_HIPM': ['Run2016_UL_HIPM'] + [
                f'Run2016{b}' for b in 'BCDEF']+['JPsi_pythia8'],
            'Run2016_UL': ['Run2016_UL'] + [
                f'Run2016{b}' for b in 'FGH']+['JPsi_pythia8'],
            'Run2017_UL': ['Run2017_UL'] + [
                f'Run2017{b}' for b in 'BCDEF']+['JPsi_pythia8'],
            'Run2018_UL': ['Run2018_UL'] + [
                f'Run2018{b}' for b in 'ABCD']+['JPsi_pythia8'],
            # heavy ion
            'Run2016_HI_pPb_8TeV': ['Run2016'],
            # rereco (i.e. legacy)
            'Run2016': ['Run2016'] + [
               f'Run2016{b}' for b in 'BCDEFGH']+['JPsi_pythia8'],
            'Run2017': ['Run2017'] + [
               f'Run2017{b}' for b in 'BCDEF']+['JPsi_pythia8'],
            'Run2018': ['Run2018'] + [
               f'Run2018{b}' for b in 'ABCD']+['JPsi_pythia8'],
            'Run2022': ['Run2022'] + [
               f'Run2022{b}' for b in 'BC']+['DY_madgraph']
        },
    }
    return subEras.get(resonance, {}).get(era, [])


def get_data_mc_sub_eras(resonance, era):
    eraMap = {
        'Z': {
            # TODO: decide how to handle alternate generators
            # ultra legacy
            'Run2016_UL_HIPM': ['Run2016_UL_HIPM', 'DY_madgraph'],
            'Run2016_UL': ['Run2016_UL', 'DY_madgraph'],
            'Run2017_UL': ['Run2017_UL', 'DY_madgraph'],
            'Run2018_UL': ['Run2018_UL', 'DY_madgraph'],
            # Double muon PD
            'Run2016_UL_HIPM_DM': ['Run2016_UL_HIPM_DM', 'DY_madgraph'],
            'Run2016_UL_DM': ['Run2016_UL_DM', 'DY_madgraph'],
            'Run2017_UL_DM': ['Run2017_UL_DM', 'DY_madgraph'],
            'Run2018_UL_DM': ['Run2018_UL_DM', 'DY_madgraph'],
            # rereco (i.e. legacy)
            'Run2016': ['Run2016', 'DY_madgraph'],
            'Run2017': ['Run2017', 'DY_madgraph'],
            'Run2018': ['Run2018', 'DY_madgraph'],
            'Run2022': ['Run2022', 'DY_madgraph']
        },
        'JPsi': {
            # ultra legacy
            'Run2016_UL_HIPM': ['Run2016_UL_HIPM', 'JPsi_pythia8'],
            'Run2016_UL': ['Run2016_UL', 'JPsi_pythia8'],
            'Run2017_UL': ['Run2017_UL', 'JPsi_pythia8'],
            'Run2018_UL': ['Run2018_UL', 'JPsi_pythia8'],
            # heavy ion
            'Run2016_HI_pPb_8TeV': ['Run2016', None],
            # rereco (i.e. legacy)
            'Run2016': ['Run2016', 'JPsi_pythia8'],
            'Run2017': ['Run2017', 'JPsi_pythia8'],
            'Run2018': ['Run2018', 'JPsi_pythia8'],
        },
    }
    return eraMap.get(resonance, {}).get(era, [None, None])

