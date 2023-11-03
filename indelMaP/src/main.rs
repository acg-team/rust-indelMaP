#![allow(non_snake_case)]
use crate::cli::Cli;
use anyhow::Error;
use clap::Parser;
use log::{error, info, LevelFilter};
use parsimony::parsimony_alignment::pars_align_on_tree;
use parsimony::parsimony_alignment::parsimony_costs::parsimony_costs_model::GapMultipliers;
use parsimony::parsimony_alignment::parsimony_costs::parsimony_costs_model::{
    DNAParsCosts, ProteinParsCosts,
};
use phylo::alignment::{compile_alignment_representation, Alignment};
use phylo::io;
use phylo::phylo_info::{setup_phylogenetic_info, PhyloInfo};
use phylo::sequences::{get_sequence_type, SequenceType};
use phylo::tree::{get_percentiles_rounded, NodeIdx};
use phylo::Rounding;
use pretty_env_logger::env_logger::Builder;
use std::path::PathBuf;
use std::result::Result::Ok;

mod cli;

type Result<T> = std::result::Result<T, Error>;

pub fn indel_map_align_dna(
    info: &PhyloInfo,
    model_name: String,
    model_params: Vec<f64>,
    gap_mult: &GapMultipliers,
    categories: u32,
) -> Result<(Vec<Alignment>, Vec<f64>)> {
    indel_map_align_dna_rounded(
        info,
        model_name,
        model_params,
        gap_mult,
        categories,
        &Rounding::none(),
    )
}

pub(crate) fn indel_map_align_dna_rounded(
    info: &PhyloInfo,
    model_name: String,
    model_params: Vec<f64>,
    gap_mult: &GapMultipliers,
    categories: u32,
    rounding: &Rounding,
) -> Result<(Vec<Alignment>, Vec<f64>)> {
    let times = get_percentiles_rounded(&info.tree.get_all_branch_lengths(), categories, rounding);
    let scoring = DNAParsCosts::new(
        &model_name,
        &model_params,
        gap_mult,
        &times,
        false,
        rounding,
    )?;
    Ok(pars_align_on_tree(&scoring, info))
}

pub fn indel_map_align_protein(
    info: &PhyloInfo,
    model_name: String,
    _: Vec<f64>,
    gap_mult: &GapMultipliers,
    categories: u32,
) -> Result<(Vec<Alignment>, Vec<f64>)> {
    indel_map_align_protein_rounded(
        info,
        model_name,
        vec![],
        gap_mult,
        categories,
        &Rounding::none(),
    )
}

pub(crate) fn indel_map_align_protein_rounded(
    info: &PhyloInfo,
    model_name: String,
    _: Vec<f64>,
    gap_mult: &GapMultipliers,
    categories: u32,
    rounding: &Rounding,
) -> Result<(Vec<Alignment>, Vec<f64>)> {
    let times = get_percentiles_rounded(&info.tree.get_all_branch_lengths(), categories, rounding);
    let scoring = ProteinParsCosts::new(&model_name, gap_mult, &times, false, rounding)?;
    Ok(pars_align_on_tree(&scoring, info))
}

fn main() -> Result<()> {
    Builder::new()
        .filter_level(LevelFilter::Info)
        .format_timestamp_secs()
        .format_module_path(false)
        .init();
    info!("IndelMaP run started");
    let cli = Cli::try_parse()?;
    info!("Successfully parsed the command line parameters");
    let info = setup_phylogenetic_info(cli.seq_file, cli.tree_file);
    match info {
        Ok(info) => {
            let (alignment, scores) = match get_sequence_type(&info.sequences) {
                SequenceType::DNA => {
                    info!("Working on DNA data -- please ensure that data type is inferred correctly.");
                    indel_map_align_dna(
                        &info,
                        cli.model,
                        cli.model_params,
                        &GapMultipliers::new(cli.go, cli.ge),
                        cli.categories,
                    )?
                }
                SequenceType::Protein => {
                    info!("Working on protein data -- please ensure that data type is inferred correctly.");
                    indel_map_align_protein(
                        &info,
                        cli.model,
                        cli.model_params,
                        &GapMultipliers::new(cli.go, cli.ge),
                        cli.categories,
                    )?
                }
            };
            info!(
                "Final alignment score is: \n{:?}",
                scores.iter().sum::<f64>()
            );
            let out_msa_path = match cli.output_msa_file {
                Some(path) => path,
                None => {
                    let path = PathBuf::from("msa.fasta");
                    info!(
                        "No output file name provided, writing MSA to default file {}.",
                        path.display()
                    );
                    path
                }
            };
            io::write_sequences_to_file(
                &compile_alignment_representation(&info, &alignment, None::<NodeIdx>),
                out_msa_path,
            )?;
            info!("IndelMAP alignment done, quitting.");
        }
        Err(error) => {
            error!("Error in input data: {}", error);
            return Ok(());
        }
    }
    Ok(())
}

#[cfg(test)]
mod real_data_tests;
