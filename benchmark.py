import os
import pathlib
import sys
import logging
import click
import subprocess
from typing import *


logging.basicConfig(stream=sys.stdout, format='%(asctime)s - %(levelname)s - %(message)s', level=logging.DEBUG)

class Benchmark():
    def __init__(self, genomefa, tefa, tegtf, bwagenome, bowtiete, outputdir, threads, sConfig, tefamily, coverage, repeatsmask):
        # Comment
        self.genomefa  = genomefa
        self.tefa      = tefa
        self.outputdir = outputdir
        self.threads   = threads

        # Alignment
        self.bwagenome = bwagenome
        self.bowtiete  = bowtiete

        # Simulation
        self.tegtf    = tegtf
        self.sConfig  = sConfig
        self.tefamily = tefamily
        self.coverate = coverage

        # relocate2
        self.repeatsmask = repeatsmask

        # intermediate files and results

        # output fastq
        # default output
        #deffastq = pathlib.Path(self.outputdir, "data", "forward")
        dirfastq = pathlib.Path(self.outputdir, "fastq")
        insert   = pathlib.Path(self.outputdir, "insert")

        if not dirfastq.exists():
            os.makedirs(str(dirfastq))
        if not insert.exists():
            os.makedirs(insert)
        self.dirfastq = dirfastq
        self.insert   = insert

        # outbam
        outbam = pathlib.Path(self.outputdir, "bam")
        if not outbam.exists():
            os.makedirs(str(outbam))
        self.outbam = outbam
        self.bwabam = pathlib.Path(self.outbam, "bwa_alignment.bam")

        # tei output
        teiout = pathlib.Path(self.outputdir, "tei")
        if not teiout.exists():
            os.makedirs(str(teiout))
        self.teiout = teiout

        # ngs_te_mapper2
        ngs_te_mapper2 = pathlib.Path(self.outputdir, "ngstemapper2")
        if not ngs_te_mapper2.exists():
            os.makedirs(str(ngs_te_mapper2))
        self.ngs_te_mapper2 = ngs_te_mapper2

        # RelocaTE2
        temp = pathlib.Path(self.outputdir, "temp")
        if not temp.exists():
            os.makedirs(str(temp))
        self.temp = temp
        relocate2 = pathlib.Path(self.outputdir, "relocate2")
        if not relocate2.exists():
            os.makedirs(str(relocate2))
        self.relocate2 = relocate2

        # Retroseq
        retroseq = pathlib.Path(self.outputdir, "retroseq")
        if not retroseq.exists():
            os.mkdir(str(retroseq))
        self.retroseq = retroseq

        # TEMP2
        temp2 = pathlib.Path(self.outputdir, "temp2")
        if not temp2.exists():
            os.makedirs(str(temp2))
        self.temp2 = temp2

    def run(self, cmd):
        logging.info(cmd)
        status = subprocess.run(cmd, stderr=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)
        return (status)
    def simulation(self):
        commands = f'python3 /mnt/Storage2/home/xiaoyihan/1129/development/simulation.py -r {self.genomefa} -c {self.tefa} -g {self.tegtf} -t {self.tefamily} -j {self.sConfig} -p {self.threads} --end 1 -C {self.coverate} -o {self.outputdir}'
        stauts = self.run(commands)
        if stauts.returncode == 0:
            self.run("mv " + str(pathlib.Path(self.outputdir, "data", "forward", "*fastq")) + " " + str(self.dirfastq))
            self.run("mv " + str(pathlib.Path(self.outputdir, "data", "forward", "*bed")) + " " + str(self.insert))
            self.run("rm -rf " + str(pathlib.Path(self.outputdir, "data")))
    def bwa_mapping(self):
        fastq = pathlib.Path(self.dirfastq, "*fastq")
        commands = f'bwa mem {self.genomefa} {str(fastq)}  -t {self.threads} | samtools view -bS - | samtools sort - > {self.bwabam}'
        self.run(commands)
        self.run(f'samtools index {self.bwabam}')
    def TEi(self):
        commands = f'Rscript /mnt/Storage2/home/xiaoyihan/1129/development/tei_pipeline.R -i {self.bwabam} -r {self.bowtiete} -o {self.teiout}'
        self.run(commands)
        file = str(self.teiout) + "/" + "soft_clip*"
        self.run(f'rm {file}')

    def neg_te_mapper2(self):
        fastqs = os.listdir(str(self.dirfastq))
        fastq1 = str(pathlib.Path(self.dirfastq, fastqs[0]))
        fastq2 = str(pathlib.Path(self.dirfastq, fastqs[1]))
        commands = f'ngs_te_mapper2 -f {fastq1},{fastq2} -r {self.genomefa} -l {self.tefa} -t {self.threads} -o {self.ngs_te_mapper2}'
        status = self.run(commands)

    def RelocaTE2(self):
        commands = f'python2 /mnt/Storage2/home/xiaoyihan/miniconda3/envs/RelocaTE2/scripts/relocaTE2.py --genome_fasta {self.genomefa} --fq_dir {str(self.dirfastq)} --te_fasta {self.tefa} -c {self.threads} --reference_ins {self.repeatsmask} --outdir {str(self.temp)} --run'
        self.run(commands)
        targetfile = pathlib.Path(self.temp, "repeat/results/ALL.all_nonref_insert.txt")
        self.run(f'mv {str(targetfile)} {str(self.relocate2)}')
        self.run(f'rm -rf {str(self.temp)}')

    def RetroSeq(self, te_fasta_fastafile):
        # first step
        output = str(pathlib.Path(self.retroseq, "ouptput.txt"))
        commands = f'perl /mnt/Storage2/home/xiaoyihan/Teddy/benchmark/RetroSeq/bin/retroseq.pl -discover -bam {self.bwabam} -eref {te_fasta_fastafile} -output {output}'
        self.run(commands)

        # second step
        output2 = str(pathlib.Path(self.retroseq, "retroseq"))
        commands = f'perl /mnt/Storage2/home/xiaoyihan/Teddy/benchmark/RetroSeq/bin/retroseq.pl -call -bam {self.bwabam} -input {output} -ref {str(self.genomefa)} -output {output2}'
        self.run(commands)
    def TEMP2(self):
        fastqs = os.listdir(str(self.dirfastq))
        fastq1 = str(pathlib.Path(self.dirfastq, fastqs[0]))
        fastq2 = str(pathlib.Path(self.dirfastq, fastqs[1]))
        commands = f'/mnt/Storage2/home/xiaoyihan/Teddy/benchmark/TEMP2/TEMP2 insertion -l {fastq1} -r {fastq2} -i {self.bwabam} -I {self.bwagenome} -g {self.genomefa} -R {self.tefa} -t {self.tegtf} -o {self.temp2} -p temp -c {self.threads}'
        self.run(commands)

@click.command(short_help="Benchmark of transposon new insertion site detection")
@click.argument("genomefa",
                type=click.Path(exists=True,
                                file_okay=True,
                                dir_okay=False,
                                readable=False,
                                resolve_path=True))
@click.argument("tefa",
                type=click.Path(exists=True,
                                file_okay=True,
                                dir_okay=False,

                                readable=True,
                                resolve_path=True))
@click.argument("tegtf",
                type=click.Path(exists=True,
                                file_okay=True,
                                dir_okay=False,
                                readable=True,
                                resolve_path=True))
@click.argument("bwagenome",
                type=click.Path(exists=True,
                                file_okay=True,
                                dir_okay=False,
                                readable=False,
                                resolve_path=True))
@click.argument("bowtiete",
                type=click.Path(exists=False,
                                file_okay=False,
                                dir_okay=False,
                                readable=False,
                                resolve_path=True))
@click.argument("outputdir",
                type=click.Path(exists=False,
                                file_okay=False,
                                dir_okay=True,
                                readable=False,
                                resolve_path=True))
@click.option("--threads", "-t",
              help = "Number of cpus to use",
              default=10)
@click.option("--sconfig", "-s",
              help="Use to generate simulation data",
              default="/mnt/Storage2/home/xiaoyihan/Teddy/reference/simulation/yeast/yeast_config.json",
              )
@click.option("--tefamily", "-f",
              help="File have two clomuns, contain teid and te family",
              default="/mnt/Storage2/home/xiaoyihan/Teddy/reference/sac_cer_te_families.tsv")
@click.option("--coverage", "-c",
              help="simulation data coverage",
              default=10)
@click.option("--repeatsmask", "-r",
              help="Repeats marsk output using in relocated",
              default="/mnt/Storage2/home/xiaoyihan/Teddy/reference/sacCer2.fasta.out")

def benchmark(genomefa: str, tefa: str, tegtf: str, bwagenome: str, bowtiete: str, outputdir: str, threads: int, sconfig: str, tefamily: str, coverage: int, repeatsmask: str) -> None:
    pip = Benchmark(genomefa=genomefa, tefa=tefa, tegtf=tegtf, bwagenome=bwagenome, bowtiete=bowtiete, outputdir=outputdir, threads=threads, sConfig=sconfig, tefamily=tefamily, coverage=coverage, repeatsmask=repeatsmask)
    pip.simulation()
    pip.bwa_mapping()
    pip.TEi()
    pip.neg_te_mapper2()
    #pip.RelocaTE2()
    pip.RetroSeq("/mnt/Storage2/home/xiaoyihan/Teddy/reference/te_fasta.tab")
    pip.TEMP2()

if __name__ == '__main__':
    benchmark()
