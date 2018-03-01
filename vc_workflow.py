#!/usr/bin/env python2.7
from __future__ import print_function

import argparse
import sys, os
import subprocess
import yaml
from utils import parameter,message
from json import loads
from shutil import copyfile
from urlparse import urlparse
from toil.common import Toil
from toil.job import Job
from toil_lib import require, UserError
from toil_lib import partitions
from utils.output import consolidate_tuple_output
from variation import variation_construct, joint_calling


def construct_job(job, func, inputs, *args):
    """
    Spawns a tree of jobs to avoid overloading the number of jobs spawned by a single parent.
    This function is appropriate to use when batching samples greater than 100.

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param function func: Function to spawn dynamically, passes one sample as first argument
    :param list inputs: Array of bam samples to be batched
    :param list args: any arguments to be passed to the function
    """
    num_partitions = 100
    partition_size = len(inputs)/num_partitions
    if partition_size > 1:
        # For larger sample number (100), we only process the first 100 samples for joint calling! 
        for partition in partitions(inputs, partition_size):
            return job.addChildJobFn(construct_job, func, partition, *args).rv()
    else:
        total_gvcf = [None] * len(inputs)
        for i in range(len(inputs)):
            total_gvcf[i] = job.addChildJobFn(func, inputs[i], *args).rv()
    total_vars_out = job.addFollowOnJobFn(joint_calling, config, upstream_output=total_gvcf, errMarker='30901').rv()
    return job.addFollowOnJobFn(consolidate_tuple_output, config, upstream_output=total_vars_out, errMarker='30902').rv()

def run_variation(job, sample, config):
    config = argparse.Namespace(**vars(config))
    config.bam = sample
    disk = 2 * os.path.getsize(config.bam)
    job.fileStore.logToMaster('********** -------------- ************* Variation Calling Module for sample: {}'.format(config.uuid))
    ### Now, we would construct the whole workflow to complete the analysis.
    return job.addChildJobFn(variation_construct, config, cores=2, disk=disk).rv()    

def bamfile_check(inputFileList,checkidname="StoreId",bampath=None, baipath=None):
    for flag_file in inputFileList:
        if flag_file.has_key(checkidname):
            if urlparse(flag_file[checkidname]).scheme == "file":
                require(os.path.exists(urlparse(flag_file[checkidname]).path), 'E30002: The input file {} does not existed for panelOfNormal creating process!'.format(flag_file[checkidname]))
            else:
                return 1
            if flag_file["resultType"] == "seq_sorted_bam_file":
                bampath = urlparse(flag_file[checkidname]).path
                baipath = bampath + ".bai"
            elif flag_file["resultType"] == "seq_bam_index_file":
                copyfile(urlparse(flag_file[checkidname]).path, baipath)
    return bampath, baipath

def main():
    """
    The module WES/TRS module analysis pipeline, ***** variation calling module *****

    ----------------------------------------------
    Bioinformatics Center, Shanghai Geniminix Co. Ltd.
    Author :  Schaudge King
    Date   :  2017-03-15
    """
    
    os.environ["SENTIEON_LICENSE"] = "192.168.2.72:8990"
    
    parser = argparse.ArgumentParser(description=main.__doc__, formatter_class=argparse.RawTextHelpFormatter)
    subparsers = parser.add_subparsers(dest='command')
    parser_run = subparsers.add_parser('run', help='Runs the variation analysis module ...')
    parser_run.add_argument('--config', default='config-toil.yaml', type=str, help='Path to the (filled in) config file.')
    parser_run.add_argument('--message', default='', type=str, help='The input message string.')
    Job.Runner.addToilOptions(parser_run)
    # If no arguments provided, print full help menu
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    cwd = os.path.dirname(sys.argv[0])
    if args.command == 'run':
        args.config = cwd + "/" + args.config
        require(os.path.exists(args.config), '{} not found.'.format(args.config))
        run_env = "private" if os.getenv('GCBI_PROFILES') is None else os.getenv('GCBI_PROFILES')
        # Parse config
        tool_required = ['bamTool', 'baseRecalibrationTool', 'variantCallTool']
        message_keys = ["processId", "processType", "accession", "resultPath",
                        "genomeVersion", "bedTarget", "baseRecalibration", "includeDuplication"]
        file_required = ['fa_ref', 'dbsnp', 'baseRecalibrationVcf']
        whole_config = { x : y for x, y in yaml.load(open(args.config).read()).iteritems()}
        if whole_config.has_key(run_env) and whole_config[run_env] is not None:
            whole_config['default'].update(whole_config[run_env])
        parsed_config = whole_config['default']
        if urlparse(parsed_config["mq_conf"]).scheme == "file":
            mq_conf = urlparse(parsed_config["mq_conf"]).path
            require(os.path.exists(urlparse(parsed_config["mq_conf"]).path),
                    'E30001: The message configure file {} does not existed for variation calling process!'.format("mq_conf"))
        else:
            return 1
        message_dict = loads(args.message)

        parsed_config = parameter.tool_check(parsed_config, tool_required, subToolNeeded=True, toolName=message_dict["variantCallingTool"], parsConfFile=message_dict["callingToolArguFile"], errMarker=30101)
        parsed_config.update(parameter.arguments_modify(message_dict, message_keys, check='bedTarget', errMarker=30301))
        parsed_config.update(parameter.reference_check(parsed_config[parsed_config["genomeVersion"]], file_required, errMarker=30201))
        config = argparse.Namespace(**parsed_config)
        config.maxCores = int(args.maxCores) if args.maxCores else sys.maxint
        config.mq_conf = mq_conf
        config.uuid = config.accession
        config.bam, config.idx = bamfile_check(message_dict["inputFileList"],checkidname="storeId")

        require(config.resultPath, 'E30004: No output location specified: {}'.format(config.resultPath))
        if not config.resultPath.endswith('/'):
            config.resultPath += '/'

        # Start the workflow
        with Toil(args) as toil:
            return toil.start(Job.wrapJobFn(construct_job, run_variation, [config.bam], config))



if __name__ == '__main__':
    try:
        uuid = (sys.argv[2]).split('/')[-1]
        resultFileList = main()
        message_context = {"succ": True, "processId": uuid, "processType": "variationCalling","resultFileList": resultFileList}
        send_state = message.send_message(os.path.dirname(sys.argv[0])+"/rabbitmq.conf", message_context)
        if send_state != 0:
            connect_addr = message.config.main(os.path.dirname(sys.argv[0]) + "/rabbitmq.conf")['address']
            print("The connection to \"{}\" for rabbitmq does not successfully!".format(connect_addr))
        else:
            print("The task " + (sys.argv[2]).split(':')[1] + " send message successfully! --- end")
    except UserError as e:
        print(e.message, file=sys.stderr)
        uuid = (sys.argv[2]).split('/')[-1]
        errMarker = e.message.split(":")[0]
        message_context = {"succ": False, "processId": uuid, "errCode": errMarker, "processType": "variationCalling",
                           "outputFile": "", "resultFileList": []}
        message.send_message(os.path.dirname(sys.argv[0])+"/rabbitmq.conf", message_context)
    except Exception as ue:
        print("--------- Baddly, we meet the ******** RUNNING TIME ERROR ********\n", file=sys.stderr)
        print(str(ue) + "\n", file=sys.stderr)
        uuid = (sys.argv[2]).split('/')[-1]
        message_context = {"succ": False, "processId": uuid, "errCode": "E30000", "processType": "variationCalling",
                           "outputFile": "", "resultFileList": []}
        message.send_message(os.path.dirname(sys.argv[0])+"/rabbitmq.conf", message_context)
    finally:
        uuid = (sys.argv[2]).split(':')[1]
        if os.path.exists(uuid):
            subprocess.check_call(["toil", "clean", uuid])
        sys.exit(0)
