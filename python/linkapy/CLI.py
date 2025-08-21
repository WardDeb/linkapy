from importlib import metadata
import click
from rich import print

@click.group()
@click.help_option("-h", "--help")
@click.version_option(version=metadata.version('linkapy'), prog_name='linkapy')
def linkapy() -> None:
    '''
    Linkapy CLI - A command line interface to process and analyze single-cell multiome data.
    '''
    pass

@linkapy.command(context_settings={"show_default": True})
@click.help_option("-h", "--help")
@click.option('--methylation_path', '-m', type=click.Path(exists=True), help='Path to the directory containing methylation data. Will be searched recursively to match pattern.')
@click.option('--transcriptome_path', '-t', type=click.Path(exists=True), help='Path to the directory containing transcriptome data. Will be searched recursively to match pattern. Note that these should be featureCounts files.')
@click.option('--output' ,'-o', type=click.Path(), default='linkapy_output', help='Output directory for the results. Default is "linkapy_output". RNA matrices will be written in arrow format, methylation derived matrices will be written in mtx format. Additionaly, if mudata is set, a MuData object is created as well.)')
@click.option('--methylation_pattern', type=str, default=("*CG*.tsv",), multiple=True, help='Pattern to match methylation files. Can be specified multiple times. Note that every pattern yields a separate matrix.')
@click.option('--transcriptome_pattern', type=str, default=("*.tsv",), multiple=True, help='Pattern to match transcriptome files. Can be specified multiple times. Note that every pattern yields a separate matrix.')
@click.option('--methylation_pattern_names', type=str, default=(), multiple=True, help='Labels for every methylation pattern provided. Can be specified multiple times. The name will be used to name the output files. If not provided. The asterisks will be removed from the pattern to yield labels.')
@click.option('--transcriptome_pattern_names', type=str, default=(), multiple=True, help='Labels for every transcriptome pattern provided. Can be specified multiple times. The name will be used to name the output files. If not provided. The asterisks will be removed from the pattern to yield labels.')
@click.option('--NOMe', is_flag=True, help='Assumes data under methylation_path is NOMe data. Setting this flag is the same as using "--methylation_pattern *GCHN* --methylation_pattern *WCGN*"')
@click.option('--threads', '-j', type=int, default=1, help='Number of threads to use for processing. Default is 1.')
@click.option('--chromsizes', '-c', type=click.Path(exists=True), help='Path to the chromsizes file for genome reference. Only needed if no regions are provided.')
@click.option('--regions', '-r', type=click.Path(exists=True), multiple=True, default=(), help='Path to regions file (bed format) to aggregate methylation data over. Can be specified multiple times.')
@click.option('--blacklist', type=click.Path(exists=True), multiple=True, default=(), help='Path of regions (bed format) to exclude from aggregation. Can be specified multiple times. Note that these are only relevant for methylation data.')
@click.option('--binsize', '-b', type=int, default=10000, help='Size of bins for aggregating methylation data over. Only used if chromsizes are provided.')
@click.option('--project', '-p', type=str, default='linkapy', help='Project name. Effectively used as a prefix for the output files.')
@click.option('--verbose', '-v', is_flag=True, help='Enable debugging output.')
@click.pass_context
def parsing(ctx, **kwargs) -> None:
    '''
    Parse single-cell scmethylation - / scNOMe -  and/or scRNA data.
    Either methylation_path or transcriptome_path must be provided.
    '''

    if not any((kwargs.get('methylation_path'), kwargs.get('transcriptome_path'), )):
        click.echo(ctx.get_help())
        print("Provide either a methylation path and/or a transcriptome path.")
        return
    

    if kwargs.get('methylation_path') and not any((kwargs.get('chromsizes'), kwargs.get('regions'))):
        click.echo(ctx.get_help())
        print("Methylation data requires either a chromsizes file or at least one regions file.")
        return
    
    if kwargs.get('nome'):
        # methylation pattern and names are set in Linkapy_Parser
        kwargs['methylation_pattern'] = ()
        kwargs['methylation_pattern_names'] = ()

    try:
        from linkapy.parsing import Linkapy_Parser
        lp = Linkapy_Parser(
            methylation_path=kwargs.get('methylation_path'),
            transcriptome_path=kwargs.get('transcriptome_path'),
            output=kwargs.get('output'),
            methylation_pattern=kwargs.get('methylation_pattern'),
            methylation_pattern_names=kwargs.get('methylation_pattern_names'),
            transcriptome_pattern=kwargs.get('transcriptome_pattern'),
            transcriptome_pattern_names=kwargs.get('transcriptome_pattern_names'),
            NOMe=kwargs.get('nome'),
            threads=kwargs.get('threads'),
            chromsizes=kwargs.get('chromsizes'),
            regions=kwargs.get('regions'),
            blacklist=kwargs.get('blacklist'),
            binsize=kwargs.get('binsize'),
            project=kwargs.get('project'),
            verbose=kwargs.get('verbose')
        )
        lp.parse()
    except ValueError as e:
        raise click.ClickException(f"ERROR: {str(e)}")
    except FileNotFoundError as e:
        raise click.ClickException(f"ERROR: {str(e)}. Please check the provided paths.")
    except AssertionError as e:
        raise click.ClickException(f"ERROR: {str(e)}. Please check the provided parameters")

    return