# Mapping-simulations

Repo for scripts and files for read mapping simulations.

## Getting started

Make sure [Anaconda](https://www.anaconda.com/products/individual#download-section) is installed and activate it if it hasn't been already:

```bash
    source /path/to/anaconda3/bin/activate
```

Then, all dependencies (except one, see below), can be installed using the pre-set environment, `mapping-sims`. Create that environment to install this software:

```bash
    conda env create --file mapping-sims.yml
```

And then activate the environment:

```bash
    conda activate mapping-sims
```

The only software not installed via the conda environment is the read simulation software, [NEAT](https://github.com/ncsa/NEAT). NEAT is just a python script, so you can install it just by cloning the github repo, though this must be cloned in the correct spot so the snakemake pipeline can find it. First, while in the root project directory, make the directory `pkgs` and navigate to it:

```bash
    mkdir pkgs
    cd pkgs
```

Then, clone the NEAT repo:

```bash
    git clone https://github.com/ncsa/NEAT.git
```
That should be everything!

## Running the pipeline

To make things easier, let's just switch to the `scripts` directory.

```bash
    cd scripts
```
If you'd rather not switch to `scripts`, just be sure to add that to all paths listed below.

### Important note: use `--dryrun` on any snakemake commands before running them.

`--dryrun` will show you the commands that snakemake will run without actually running them. This can let you address errors before you start the whole pipeline.

### Running snakemake with a config file

The commands for the pipeline are compiled in a snakemake workflow, `mapping_sims.smk`

To run a given set of simulations, you will have to create a <b>config file</b> with the parameters and inputs of that simulation. An example config file for the mouse genome is provided here: `mm39-config.yaml`

Take a look at that config file and make any changes, including paths to data directories and reference genomes, as well as simulation parameters.

After you have edited the config file, you can run the pipeline directly:

```
snakemake -j <number of jobs to run in parallel> \
    -s mapping_sims.smk                          \
    -p                                           \
    --configfile <your config file>
```

I like to include the `-p` option that prints out the exact shell command run for each rule, which I find helpful to understand exactly what's going on.

For example, the exact command for the provided mm39 config file to launch 20 jobs would be:

```bash
snakemake -j 20 -s mapping_sims.smk -p --configfile mm30-config.yaml
```

### Using a SLURM profile to launch snakemake jobs on the cluster

Snakemake also allows jobs to be submitted to a cluster rather than run directly on the terminal. This is advantageous as it can allow the pipeline to access many more resources than using a single node with the command above. To do this, you must set-up a <b>profile</b> for your cluster. An example profile for the Harvard SLURM cluster is provided in `scripts/profiles/slurm_profile/config.yaml`. Once your profile is set-up, you can launch the pipeline as follows:

```
snakemake -s mapping_sims.smk                    \
    -p                                           \
    --configfile <your config file>              \
    --profile <path to directory containing config.yaml profile>
```

Similar, to the previous command, but without specifying the number of jobs with `-j`, which is now done with the `--profile` option. Note that the provided path must be to a directory containing your profile, which is saved as `config.yaml`.

So, the command to run the pipeline with the provided config and profile would be:

```bash
snakemake -s mapping_sims.smk -p --configfile mm30-config.yaml --profile profiles/slurm_profile/
```