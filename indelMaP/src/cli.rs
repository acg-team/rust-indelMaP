use clap::Parser;
use std::path::PathBuf;

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
pub(super) struct Cli {
    /// Sequence file in fasta format
    #[arg(short, long, value_name = "SEQ_FILE")]
    pub(super) seq_file: PathBuf,

    /// Tree file in newick format
    #[arg(short, long, value_name = "TREE_FILE")]
    pub(super) tree_file: PathBuf,

    /// Tree file in newick format
    #[arg(short, long, value_name = "OUTPUT_MSA_FILE")]
    pub(super) output_msa_file: Option<PathBuf>,

    /// Sequence evolution model
    #[arg(short, long, value_name = "MODEL", rename_all = "UPPER")]
    pub(super) model: String,

    /// Sequence evolution model parameters, e.g. alpha and beta for k80 and
    /// f_t f_c f_a f_g r_tc r_ta r_tg r_ca r_cg r_ag for GTR (in this specific order)
    #[arg(short = 'p', long, value_name = "MODEL_PARAMS", num_args = 0..)]
    pub(super) model_params: Vec<f64>,

    /// Gap opening penalty
    #[arg(short = 'g', long, default_value_t = 2.5)]
    pub(super) go: f64,

    /// Gap extension penalty
    #[arg(short = 'e', long, default_value_t = 0.5)]
    pub(super) ge: f64,

    /// Number of percentile categories to use for branch length approximation
    #[arg(short, long, default_value_t = 4)]
    pub(super) categories: u32,
}
