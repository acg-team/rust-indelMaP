use crate::indel_map_align_protein_rounded;
use parsimony::parsimony_alignment::parsimony_costs::parsimony_costs_model::GapMultipliers;
use phylo::phylo_info::phyloinfo_from_files;
use phylo::Rounding;
use std::path::PathBuf;

#[test]
fn align_HIV_example_wag() {
    let info = phyloinfo_from_files(
        PathBuf::from("./data/HIV_subset.fas"),
        PathBuf::from("./data/HIV_subset.nwk"),
    )
    .unwrap();
    let (_, scores) = indel_map_align_protein_rounded(
        &info,
        "WAG".to_string(),
        vec![],
        &GapMultipliers::new(2.5, 0.5),
        4,
        &Rounding::four(),
    )
    .unwrap();
    assert_eq!(scores.iter().sum::<f64>(), 350.64988524999995);
}
