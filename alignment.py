#!/usr/bin/env python2.7
# We abstract the analysis module of alignment for NGS reads,
# the main step including mapping to reference genome, indel realignment,
# mark the duplication and quality control!
# Written By Schaudge King
# Date : 2017-03-21
from __future__ import print_function
import argparse
import os, sys
import subprocess
from toil.common import Toil
from toil.job import Job
from toil_lib import UserError

def align_construct(job, config, parameters, cores=4):
    """
    Define the whole alignment steps, including realign local indel region and duplication analysis

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param Namespace config: Argparse Namespace object containing argument inputs
    :param parameters: the parameter used to bwa align software.
    """
    config = argparse.Namespace(**vars(config))
    work_dir = job.fileStore.getLocalTempDir()
    output_bam = os.path.join(work_dir, config.uuid + '_sorted.bam')
    
    if config.toolVendor and config.toolVendor == "sentieon":
        #map_sort_command = [config.alignTool, "mem", "-M", "-t", str(cores), "-K", "10000000", "-R" ]+ parameters \
        #+["||", "echo", "-n", 'error', "|", config.bamSortTool, "util", "sort", "-r", config.bwa_ref, "-o", output_bam, "-t", str(cores),"--sam2bam", "-i", "-"]
        map_command = [config.alignTool, "mem", "-M", "-t", str(cores), "-K", "10000000"]+ parameters
        sort_command = [config.bamSortTool, "util", "sort", "-r", config.bwa_ref, "-o", output_bam, "-t", str(cores),"--sam2bam", "-i", "-"]
        #subprocess.check_call(map_sort_command)
        align_ps = subprocess.Popen(map_command, stdout=subprocess.PIPE)
        subprocess.check_output(sort_command, stdin=align_ps.stdout)
        align_ps.wait()
    else:
        map_command = [config.alignTool, "mem", "-t", str(cores), "-Y", "-M"] + parameters
        convert_command = [config.sam2bamTool, "view", "-u"]
        sort_command = [config.bamSortTool, "sort", "-@", str(cores), "-o", output_bam]
        align_ps = subprocess.Popen(map_command, stdout=subprocess.PIPE)
        bam_create = subprocess.Popen(convert_command, stdin=align_ps.stdout, stdout=subprocess.PIPE)
        subprocess.check_output(sort_command, stdin=bam_create.stdout)
        align_ps.wait()
        bam_create.wait()

    if align_ps.returncode != 0 or bam_create.returncode != 0:
        raise UserError('E20030: The mapping subprocess (bwa) by some unkown reasons, mybe system overhead !')
    subprocess.check_call([config.sam2bamTool, "index", output_bam])
    
    sorted_bam = job.fileStore.writeGlobalFile(output_bam)
    bam_index = job.fileStore.writeGlobalFile(output_bam + ".bai")
    dedupID = job.addChildJobFn(duplication_bam, config, sorted_bam, bam_index)
    if config.realign and config.baseRecalibration:
        realignID = job.addFollowOnJobFn(realign_indel, config, dedupID.rv(0), dedupID.rv(1))
        recalID = job.addFollowOnJobFn(base_recalibration, config, realignID.rv(0), realignID.rv(1))
        return { "seq_sorted_bam_file-" + config.uuid + "_recal.bam": recalID.rv(0)}, { "seq_bam_index_file-" + config.uuid + "_recal.bam.bai": recalID.rv(1)}, \
               { "qc_dup_metrics-" + config.uuid + "_dup.stat": dedupID.rv(2)}
    elif config.realign:
        realignID = job.addFollowOnJobFn(realign_indel, config, dedupID.rv(0), dedupID.rv(1))
        return { "seq_sorted_bam_file-" + config.uuid + "_realigned.bam": realignID.rv(0)}, { "seq_bam_index_file-" + config.uuid + "_realigned.bam.bai": realignID.rv(1)}, \
               { "qc_dup_metrics-" + config.uuid + "_dup.stat": dedupID.rv(2)}
    elif config.baseRecalibration:
        recalID = job.addFollowOnJobFn(base_recalibration, config, dedupID.rv(0), dedupID.rv(1))
        return { "seq_sorted_bam_file-" + config.uuid + "_recal.bam": recalID.rv(0)}, { "seq_bam_index_file-" + config.uuid + "_recal.bam.bai": recalID.rv(1)}, \
               { "qc_dup_metrics-" + config.uuid + "_dup.stat": dedupID.rv(2)}
    else:
        return { "seq_sorted_bam_file-" + config.uuid + "_dupmarked.bam": dedupID.rv(0)}, { "seq_bam_index_file-" + config.uuid + "_dupmarked.bam.bai": dedupID.rv(1)}, \
               { "qc_dup_metrics-" + config.uuid + "_dup.stat": dedupID.rv(2)}

def realign_indel(job, config, bamID=None, indexID=None):
    """
    Implement local realign for indel et al complex genome region!

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param Namespace config: Argparse Namespace object containing argument inputs
    :param bamID: input deduped or dupmarked bam fileID
    :param indexID: the fileID of input deduped or dupmarked bam index file
    """
    config = argparse.Namespace(**vars(config))
    work_dir = job.fileStore.getLocalTempDir()
    input_bam = os.path.join(work_dir, config.uuid + '_deduped.bam')
    input_bam_idx = os.path.join(work_dir, config.uuid + '_deduped.bam.bai')
    if bamID and indexID:
        job.fileStore.readGlobalFile(bamID, input_bam)
        job.fileStore.readGlobalFile(indexID, input_bam_idx)
    else :
        raise UserError('E20005: The realign need sorted bam. Here is right file ID !')
    realign_bam = os.path.join(work_dir, config.uuid + '_realign.bam')
    job.fileStore.logToMaster('### ------ local realign for sample {}.'.format(config.uuid))
    vcf_args = []

    if config.realignTool.endswith(".jar"):
        for vcf in config.realign_vcf.split(','):
            vcf_args += ["-known", vcf]
        target_intervals = os.path.join(work_dir, config.uuid + '_realign.intervals')
        target_maker = ["java", "-Xmx2g", "-jar", config.realignTool, "-T", "RealignerTargetCreator", "-nt", "2", "-I", input_bam] \
                          + ["-R", config.bwa_ref, "-L", config.bedfile] + vcf_args + ["-o", target_intervals]
        subprocess.check_call(target_maker)
        realign_process = ["java", "-Xmx2g", "-jar", config.realignTool, "-T", "IndelRealigner", "--filter_bases_not_stored", "-I", input_bam, "-R", config.bwa_ref] \
                        + ["-targetIntervals", target_intervals] + vcf_args + ["-o", realign_bam]
        subprocess.check_call(realign_process)
        subprocess.check_call([config.bamTool, "index", realign_bam])
    elif config.realignTool.endswith("sentieon"):
        for vcf in config.realign_vcf.split(','):
            vcf_args += ["-k", vcf]
        realign_process = [config.realignTool, "driver", "-r", config.bwa_ref, "-t", "4", "-i", input_bam, "--algo", "Realigner" ] \
                        + vcf_args + ["--interval_list", config.bedfile, realign_bam]
        subprocess.check_call(realign_process)        
    else:
        realign_process = [config.dedupTool, "I=" + input_bam, "O=" + realign_bam]
        subprocess.check_call(realign_process)
        subprocess.check_call([config.bamTool, "index", realign_bam])
    
    job.fileStore.logToMaster('------ ### local realign is OK!')
    bamRealignID = job.fileStore.writeGlobalFile(realign_bam)
    idxRealignID = job.fileStore.writeGlobalFile(realign_bam + ".bai")
    return bamRealignID, idxRealignID

def duplication_bam(job, config, bamID=None, indexID=None, testFlag=False):
    """
    Implement the (PCR et al) duplication analysis, to gain the library complex for the reads!!

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param Namespace config: Argparse Namespace object containing argument inputs
    :param bamID: input sorted bam fileID
    :param indexID: the fileID of input sorted bam index file
    """
    config = argparse.Namespace(**vars(config))
    work_dir = job.fileStore.getLocalTempDir()
    input_bam = os.path.join(work_dir, config.uuid + '_sorted.bam')
    input_bam_idx = os.path.join(work_dir, config.uuid + '_sorted.bam.bai')
    if bamID and indexID:
        job.fileStore.readGlobalFile(bamID, input_bam)
        job.fileStore.readGlobalFile(indexID, input_bam_idx)
    elif not testFlag :
        raise UserError('E20004: The duplication statistic need sorted bam and bam index. Here is right file ID !')
    else :
        input_bam = config.bam
        input_bam_idx = config.bai
    dedup_bam = os.path.join(work_dir, config.uuid + '_dedup.bam')
    dep_metric = os.path.join(work_dir, config.uuid + '_deplication.metrics')
    job.fileStore.logToMaster('### ------ PCR duplication statistic for sample {}.'.format(config.uuid))
    if config.dedupTool.endswith(".jar"):
        duplication_analysis = ["java", "-Xmx2g", "-jar", config.dedupTool, "MarkDuplicates", "VALIDATION_STRINGENCY=LENIENT", "I="+input_bam, "O="+dedup_bam, "M="+dep_metric]
        subprocess.check_call(duplication_analysis)
        subprocess.check_call([config.bamTool, "index", dedup_bam])
    elif config.dedupTool.endswith("sentieon"):
        score_txt = os.path.join(work_dir, config.uuid + '_score.txt')
        score_maker = [config.dedupTool, "driver", "-t", "8", "-i", input_bam, "--algo", "LocusCollector", "--fun", "score_info", score_txt]
        subprocess.check_call(score_maker)
        #duplication_analysis = [config.dedupTool, "driver", "-t", "8", "-i", input_bam, "--algo", "Dedup", "--rmdup", "--score_info", score_txt, "--metrics", dep_metric, dedup_bam]
        duplication_analysis = [config.dedupTool, "driver", "-t", "8", "-i", input_bam, "--algo", "Dedup", "--score_info", score_txt, "--metrics", dep_metric, dedup_bam]
        subprocess.check_call(duplication_analysis)
    else:
        # This was used for duplication tool --- sambamba
        duplication_analysis = [config.dedupTool, "markdup", "-t", "2", input_bam, dedup_bam]
        subprocess.check_call(duplication_analysis, stderr=open(dep_metric, 'w'))
        # The following snippet was used for dupliacation tool ---- elPrep
        #extra_params = ["--filter-unmapped-reads", "--mark-duplicates", "--sorting-order", "keep", "--nr-of-threads", "2"]
        #duplication_analysis = [config.dedupTool, input_bam, dedup_bam] + extra_params
        #subprocess.check_call([config.bamTool, "stat", dedup_bam], stdout=open(dep_metric,'w'))
        subprocess.check_call([config.bamTool, "index", dedup_bam])
    
    job.fileStore.logToMaster('------ ### duplication information analysis is OK!')
    bamDedupID = job.fileStore.writeGlobalFile(dedup_bam)
    idxDedupID = job.fileStore.writeGlobalFile(dedup_bam + ".bai")
    dupQcID = job.fileStore.writeGlobalFile(dep_metric)
    return bamDedupID, idxDedupID, dupQcID

def base_recalibration(job, config, bamID=None, indexID=None):
    """
    Implement the base recalibration for improving SNP&INDEL calling accurately !

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param Namespace config: Argparse Namespace object containing argument inputs
    :param bam: input sorted bam file path (the whole path)
    :param idx: the input sorted bam's index file (whole path also)
    """
    config = argparse.Namespace(**vars(config))
    work_dir = job.fileStore.getLocalTempDir()
    if config.realign:
        input_bam = os.path.join(work_dir, config.uuid + '_realign.bam')
        input_bam_idx = os.path.join(work_dir, config.uuid + '_realign.bam.bai')
    else:
        input_bam = os.path.join(work_dir, config.uuid + '_dedup.bam')
        input_bam_idx = os.path.join(work_dir, config.uuid + '_dedup.bam.bai')
    if bamID and indexID:
        job.fileStore.readGlobalFile(bamID, input_bam)
        job.fileStore.readGlobalFile(indexID, input_bam_idx)
    else:
        raise UserError('E20007: The base recalibration need sorted bam. Here is the right file ID !')
    recal_bam = os.path.join(work_dir, config.uuid + '_recal.bam')
    job.fileStore.logToMaster('### ------ base recalibration for sample {}.'.format(config.uuid))
    vcf_args = []
    recal_grp = os.path.join(work_dir, config.uuid + '_recal.grp')
    for vcf in config.baseRecalibrationVcf.split(','):
        vcf_args += ["-knownSites", vcf]
    if config.baseRecalibrationTool.endswith(".jar"):

        recal_grp_maker = ["java", "-Xmx2g", "-jar", config.baseRecalibrationTool, "-T", "BaseRecalibrator", "-I", input_bam] \
                          + ["-R", config.fa_ref, "-nct", "2"] + vcf_args + ["-L", config.bedfile, "-o", recal_grp]
        subprocess.check_call(recal_grp_maker)
        recal_process = ["java", "-Xmx2g", "-jar", config.baseRecalibrationTool, "-T", "PrintReads", "-I", input_bam, "-R", config.fa_ref] \
                        + ["-BQSR", recal_grp, "-nct", "2", "-L", config.bedfile, "-o", recal_bam]
        subprocess.check_call(recal_process)
    elif config.baseRecalibrationTool.endswith("sentieon"):
        recal_table_maker = [config.baseRecalibrationTool, "driver", "-r", config.fa_ref, "-t", "4"] \
                            + ["--interval", config.bed_target, "-i", input_bam, "--algo", "QualCal"] + vcf_args + [recal_grp]
        subprocess.check_call(recal_table_maker)
        recal_process = [config.baseRecalibrationTool, "driver", "-r", config.fa_ref, "-t", "2"] \
                        + ["-i", input_bam, "-q", recal_grp, "--algo", "ReadWriter", recal_bam]
        subprocess.check_call(recal_process)
        subprocess.check_call([config.bamTool, "index", recal_bam])
    else:
        recal_process = [config.baseRecalibrationTool, "I=" + input_bam, "O=" + recal_bam]
        subprocess.check_call(recal_process)
    subprocess.check_call([config.bamTool, "index", recal_bam])
    job.fileStore.logToMaster('------ ### base recalibration is OK!')
    bamID = job.fileStore.writeGlobalFile(recal_bam)
    idxID = job.fileStore.writeGlobalFile(recal_bam + ".bai")
    return bamID, idxID


def main():

    parser = argparse.ArgumentParser(description=main.__doc__, formatter_class=argparse.RawTextHelpFormatter)
    subparsers = parser.add_subparsers(dest='command')
    parser_run = subparsers.add_parser('run', help='Runs the analysis module ...')
    Job.Runner.addToilOptions(parser_run)
    # If no arguments provided, print full help menu
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    parsed_config = {'uuid' : 'demo', 'output_dir': 'file:///gcbi/storage/toil',\
                     'bamTool' : '/home/gcbi/eclipse/workspace/GCSAS_svn_2/trunk/software/samtools',\
                     'toolVendor' : 'sentieon',\
                     'dedupTool' : '/home/gcbi/eclipse/workspace/GCSAS_svn_2/software/sentieon/bin/sentieon',\
                     'bwa_ref' : '/gcbi/ref/genome/hg19/hg19.ucsc.fa',\
                     'realignTool' : '/home/gcbi/eclipse/workspace/GCSAS_svn_2/software/sentieon/bin/sentieon',\
                     'realign_vcf' : '/gcbi/ref/genome/hg19/1000G_phase1.indels.hg19.sites.vcf,/gcbi/ref/genome/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf',\
                     'bedfile' : '/home/gcbi/llq/storage/Test_samples/1600215/bed_target_hg19.bed',\
                     'realign' : True,\
                     'sam2bamTool': '/home/gcbi/eclipse/workspace/GCSAS_svn_2/software/samtools' ,\
                     'bamSortTool' : '/home/gcbi/eclipse/workspace/GCSAS_svn_2/software/sentieon/bin/sentieon',\
                     'alignTool' : '/opt/GATK/software/bwa'
                     #'alignTool' : '/home/gcbi/eclipse/workspace/GCSAS_svn_2/software/sentieon/bin/bwa'
                     }
    config = argparse.Namespace(**parsed_config)
    if not config.output_dir.endswith('/'):
        config.output_dir += '/'
    # Start the workflow
    os.environ["SENTIEON_LICENSE"] = "192.168.2.72:8990"
    readsGroup = "@RG\\tID:UQID\\tSM:" + config.uuid + "\\tPL:illumina\\tPU:1\\tLB:PATHO"
    fastq_path = ["/home/gcbi/llq/storage/Test_samples/1600215/1600215_L005_R1.fastq.gz",\
    "/home/gcbi/llq/storage/Test_samples/1600215/1600215_L005_R2.fastq.gz"]
    parameters = ["-R", readsGroup, config.bwa_ref] + fastq_path
    with Toil(args) as toil:
        toil.start(Job.wrapJobFn(align_construct, config, parameters))

def main1():

    parser = argparse.ArgumentParser(description=main.__doc__, formatter_class=argparse.RawTextHelpFormatter)
    subparsers = parser.add_subparsers(dest='command')
    parser_run = subparsers.add_parser('run', help='Runs the analysis module ...')
    Job.Runner.addToilOptions(parser_run)
    # If no arguments provided, print full help menu
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    parsed_config = {'uuid' : 'demo', 'output_dir': 'file:///gcbi/storage/toil',\
                     'bamTool' : '/home/gcbi/eclipse/workspace/GCSAS_svn_2/trunk/software/samtools',\
                     'toolVendor' : 'sentieon',\
                     'dedupTool' : '/home/gcbi/eclipse/workspace/GCSAS_svn_2/software/sentieon/bin/sentieon',\
                     'bwa_ref' : '/gcbi/ref/genome/hg19/hg19.ucsc.fa',\
                     'realignTool' : '/home/gcbi/eclipse/workspace/GCSAS_svn_2/software/sentieon/bin/sentieon',\
                     'realign_vcf' : '/gcbi/ref/genome/hg19/1000G_phase1.indels.hg19.sites.vcf,/gcbi/ref/genome/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf',\
                     'bedfile' : '/home/gcbi/llq/storage/Test_samples/1600215/bed_target_hg19.bed',\
                     'realign' : True,\
                     'sam2bamTool': '/home/gcbi/eclipse/workspace/GCSAS_svn_2/software/samtools' ,\
                     'bamSortTool' : '/home/gcbi/eclipse/workspace/GCSAS_svn_2/software/sentieon/bin/sentieon',\
                     'alignTool' : '/opt/GATK/software/bwa',\
                     'bam':'/home/gcbi/llq/storage/Test_samples/1600215/sorted.bam',\
                     'bai':'/home/gcbi/llq/storage/Test_samples/1600215/sorted.bai'
                     }
    config = argparse.Namespace(**parsed_config)
    if not config.output_dir.endswith('/'):
        config.output_dir += '/'
    # Start the workflow
    os.environ["SENTIEON_LICENSE"] = "192.168.2.72:8990"
    readsGroup = "@RG\\tID:UQID\\tSM:" + config.uuid + "\\tPL:illumina\\tPU:1\\tLB:PATHO"
    fastq_path = ["/home/gcbi/llq/storage/Test_samples/1600215/1600215_L005_R1.fastq.gz",\
    "/home/gcbi/llq/storage/Test_samples/1600215/1600215_L005_R2.fastq.gz"]
    parameters = ["-R", readsGroup, config.bwa_ref] + fastq_path
    with Toil(args) as toil:
        toil.start(Job.wrapJobFn(duplication_bam, config, testFlag=True))

if __name__ == '__main__':
    try:
        main()
    except UserError as e:
        print(e.message, file=sys.stderr)
        sys.exit(1)
