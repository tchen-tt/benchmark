import logging
import click
import os
import subprocess
import pathlib
import sys


logging.basicConfig(stream=sys.stdout, format='%(asctime)s - %(levelname)s - %(message)s', level=logging.DEBUG)
class Pipelines():
    def __init__(self, threads, gtffile, output, filename):
        self.threads = threads
        self.ouput = output
        self.gtffile = gtffile
        self.filename = filename

    def run(self, cmd):
        logging.info(cmd)
        status = subprocess.run(cmd.split(sep=" "), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return(status)

    def star_alignment(self, fastq1, fastq2, genomeindex):
        fastq = str(fastq1) if (len(fastq2.suffix) == 0) else str(fastq1) + " " + str(fastq2)
        compress = "UncompressionCommand" if (fastq1.suffix != ".gz") else "zcat"
        if compress == "UncompressionCommand":
            compress = ""
        else:
            compress = "--readFilesCommand " + "zcat"

        outbam = pathlib.Path(self.ouput, "bam")
        if not outbam.exists():
            os.makedirs(str(outbam))
        outbamfile = str(pathlib.Path(outbam, self.filename))

        commands = f'STAR --runThreadN {self.threads} --genomeDir {genomeindex} --readFilesIn {fastq} {compress} --outFileNamePrefix {outbamfile} --outSAMtype BAM SortedByCoordinate'
        self.run(commands)
        # Aligned.sortedByCoord.out.bam
    def rm_duplicated(self):
        inputbam  = str(pathlib.Path(self.ouput, "bam", self.filename)) + "Aligned.sortedByCoord.out.bam"
        outputbam = str(pathlib.Path(self.ouput, "bam", self.filename)) + ".rmdup.bam"
        commands = f'samtools markdup -r -@ {self.threads} -O BAM {inputbam} {outputbam}'
        status = self.run(commands)
        if status.returncode == 0:
            self.run('rm ' + inputbam)
        else:
            logging.ERROR(status.stderr)


    def htseq_count(self):
        bamfile = str(pathlib.Path(self.ouput, "bam", self.filename)) + ".rmdup.bam"

        outcount = pathlib.Path(self.ouput, "count")
        if not outcount.exists():
            os.makedirs(str(outcount))
        outcount = pathlib.Path(outcount, self.filename).with_suffix(".txt")

        commands = f'htseq-count -f bam -s no -i gene_id {bamfile} {self.gtffile}'
        status = self.run(commands)
        print(status.stdout.decode('utf8'), file=open(outcount, "w"))
    def stringtie_denovo(self):
        bamfile = str(pathlib.Path(self.ouput, "bam", self.filename)) + ".rmdup.bam"
        outgtf  = pathlib.Path(self.ouput, "gtf")
        if not outgtf.exists():
            os.makedirs(str(outgtf))
        outgtf = pathlib.Path(outgtf, self.filename).with_suffix(".gtf")
        commands = f'stringtie -p {self.threads} -G {self.gtffile} -o {outgtf} {bamfile}'
        self.run(commands)



@click.command(short_help='Sequence alignment, gene quantification  and transcript assembly')
@click.argument("genomeindex",
                type=click.Path(exists=True,
                                file_okay=False,
                                dir_okay=True,
                                readable=True,
                                resolve_path=False))
@click.argument("gtffile",
                type=click.Path(exists=True,
                                file_okay=True,
                                dir_okay=False,
                                readable=True,
                                resolve_path=True))
@click.argument("output",
                type=click.Path(exists=False,
                                file_okay=False,
                                dir_okay=True,
                                readable=False,
                                resolve_path=True))


@click.argument("fastq1",
                type=click.Path(exists=True,
                                file_okay=True,
                                dir_okay=False,
                                readable=True,
                                resolve_path=True))
@click.argument("fastq2",
              type=click.Path(exists=False,
                              file_okay=True,
                              dir_okay=False,
                              readable=True,
                              resolve_path=True),
              default="")
@click.option("--threads", "-t",
              help="Number of cpus use for pipelines",
              default=10)
def star_alignment_denovo(genomeindex: str, gtffile: str, output: str, fastq1: str, fastq2: str, threads: int):
    file = pathlib.Path(fastq1)
    filename = file.stem
    output = pathlib.Path(output)

    fastq1 = pathlib.Path(fastq1)
    fastq2 = pathlib.Path(fastq2)

    pipe = Pipelines(threads=threads, gtffile=gtffile, output=output, filename=filename)
    pipe.star_alignment(fastq1, fastq2, genomeindex=genomeindex)
    pipe.rm_duplicated()
    pipe.htseq_count()
    pipe.stringtie_denovo()


if __name__ == '__main__':
    star_alignment_denovo()
