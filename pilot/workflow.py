import os
from gwf import Workflow, AnonymousTarget

gwf = Workflow(defaults={'account': 'chrXh2'})


def modpath(p, parent=None, base=None, suffix=None, relpath=None):
    """
    Modifies dir, basename, or suffix of path.
    """
    assert not (parent is not None and relpath is not None), "Use either parent or relpath"
    par, name = os.path.split(p)
    name_no_suffix, suf = os.path.splitext(name)
    if type(suffix) is str:
        suf = suffix
    if parent is not None:
        par = parent
    if relpath is not None:
        os.path.dirname(os.path.relpath(p, relpath))
    if base is not None:
        name_no_suffix = base
    new_path = os.path.join(par, name_no_suffix + suf)
    if type(suffix) is tuple:
        assert len(suffix) == 2
        new_path, nsubs = re.subn(r'{}$'.format(suffix[0]), suffix[1], new_path)
        assert nsubs == 1, nsubs
    return new_path



def ld_scores(input_base_name, output_base_name):
    """
    Segment based LD score
    """
    inputs = [modpath(input_base_name, suffix=x) for x in ['.bed', '.bim', '.fam']]
    outputs = {'score.ld': output_base_name + '.score.ld'}
    cores = 1
    options = {'memory': f'{8*cores}g',
               'walltime': '6-00:00:00',
               'cores': cores,
              } 
    spec = f"""
    {gcta} \
        --bfile {input_base_name} \
        --ld-score-region 200 \
        --ld-wind 10000 \
        --ld-rsq-cutoff 0 \
        --thread-num {cores} \
        --out {output_base_name}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)
    
def segment_ld_scores(ls_scores_file_name, output_base_name):
    """

    """
    inputs = [ls_scores_file_name]
    outputs = [output_base_name + f'_snp_group_{i}.txt' for i in range(1,5)]
    cores = 1
    options = {'memory': f'{8*cores}g',
               'walltime': '6-00:00:00',
               'cores': cores,
              } 
    spec = f"""
    Rscript --vanilla scripts/segment_ld_scores.r {ls_scores_file_name} {output_base_name}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)
    
def grm(input_base_name, snp_subset_file_name, output_base_name):
    """

    """
    inputs = [modpath(input_base_name, suffix=x) for x in ['.bed', '.bim', '.fam']] + [snp_subset_file_name]
    outputs = {'grm.bin': output_base_name + '.grm.bin',
               'grm.N.bin': output_base_name + '.grm.N.bin',
               'grm.id': output_base_name + '.grm.id'
               }
    cores = 1
    options = {'memory': f'{8*cores}g',
               'walltime': '6-00:00:00',
               'cores': cores,
              } 
    spec = f"""
    {gcta} \
        --bfile {input_base_name} \
        --extract {snp_subset_file_name} \
        --chr {chrom} \
        --make-grm \
        --out {output_base_name} \
        --thread-num {cores}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)
    

def reml(input_base_name, grm_file_list, pheno_file_name, output_base_name):
    """

    """
    inputs = [modpath(input_base_name, suffix=x) for x in ['.bed', '.bim', '.fam']]
    inputs += [pheno_file_name]
    for f in grm_file_list:
        inputs += [modpath(f, suffix=x) for x in ['.grm.bin', '.grm.N.bin', '.grm.id']]
    outputs = []
    cores = 1
    options = {'memory': f'{8*cores}g',
               'walltime': '6-00:00:00',
               'cores': cores,
              } 
    multi_grms_file_name = output_base_name + 'multi_grms.txt'
    spec = f"""
    rm -f {multi_grms_file_name}
    for N in {' '.join(grm_file_list)}; do echo $N >> {multi_grms_file_name}; done

    {gcta} \
        --reml --mgrm {multi_grms_file_name} \
        --pheno {pheno_file_name} \
        --out {output_base_name} \
        --prevalence {prevalence} \
        --thread-num {cores}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)
    
gcta = 'gcta_v1.94.0Beta_macOS/gcta_v1.94.0Beta_macOS'
#gcta = 'gcta_v1.94.1_linux_kernel_3_x86_64/gcta_v1.94.1_linux_kernel_3_x86_64_static'

#chromosomes = list(map(str, range(1,23))) + ['X']
chromosomes = ['1']

prevalence = 0.01

output_base_name = 'steps/gwa'
pheno_file_name = modpath(output_base_name, suffix='.phen')

target = gwf.target_from_template(f'ld_scores',
    ld_scores(output_base_name, output_base_name))
target = gwf.target_from_template(f'segment_ld_scores',
    segment_ld_scores(target.outputs['score.ld'], output_base_name))
segment_ld_files = target.outputs
for chrom in chromosomes:
    grm_base_name_list = []
    for i, sgm_ld_file_name in enumerate(segment_ld_files):
        grm_base_name = output_base_name + f'_{chrom}_{i+1}'
        target = gwf.target_from_template(f'grm_{chrom}_snp_group_{i+1}',
            grm(output_base_name, sgm_ld_file_name, grm_base_name))
        grm_base_name_list.append(grm_base_name)

    target = gwf.target_from_template(f'reml',
        reml(output_base_name, grm_base_name_list, pheno_file_name, output_base_name))



# test.hsq (rows are
# header line;
# name of genetic variance, estimate and standard error (SE);
# residual variance, estimate and SE;
# phenotypic variance, estimate and SE;
# ratio of genetic variance to phenotypic variance, estimate and SE;
# log-likelihood;
# sample size). If there are multiple GRMs included in the REML analysis, there will be multiple rows for the genetic variance (as well as their ratios to phenotypic variance) with the names of V(1), V(2), …

# Source  Variance        SE
# V(1)    0.389350        0.161719
# V(e)    0.582633        0.160044
# Vp      0.971984        0.031341
# V(1)/Vp 0.400573        0.164937
# The estimate of variance explained on the observed scale is transformed to that on the underlying scale:
# (Proportion of cases in the sample = 0.5; User-specified disease prevalence = 0.1)
# V(1)/Vp_L       0.657621      0.189123
# logL    -945.65
# logL0    -940.12
# LRT    11.06
# Pval    4.41e-4
# n       2000