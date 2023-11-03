use self::parsimony_costs::{BranchParsimonyCosts, ParsimonyCosts};
use self::parsimony_info::ParsimonySiteInfo;
use self::parsimony_matrices::ParsimonyAlignmentMatrices;
use self::parsimony_sets::get_parsimony_sets;
use log::{debug, info};
use phylo::alignment::Alignment;
use phylo::phylo_info::PhyloInfo;
use phylo::sequences::get_sequence_type;
use phylo::tree::{NodeIdx::Internal as Int, NodeIdx::Leaf};
use rand::prelude::*;

pub mod parsimony_costs;
pub mod parsimony_info;
pub mod parsimony_matrices;
pub(crate) mod parsimony_sets;

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) enum Direction {
    Matc,
    GapInY,
    GapInX,
}

fn rng_len(l: usize) -> usize {
    random::<usize>() % l
}

fn pars_align_w_rng(
    x_info: &[ParsimonySiteInfo],
    x_scoring: &dyn BranchParsimonyCosts,
    y_info: &[ParsimonySiteInfo],
    y_scoring: &dyn BranchParsimonyCosts,
    rng: fn(usize) -> usize,
) -> (Vec<ParsimonySiteInfo>, Alignment, f64) {
    let mut pars_mats = ParsimonyAlignmentMatrices::new(x_info.len() + 1, y_info.len() + 1, rng);
    debug!(
        "x_scoring: {} {} {}",
        x_scoring.avg_cost(),
        x_scoring.gap_open_cost(),
        x_scoring.gap_ext_cost()
    );
    debug!(
        "y_scoring: {} {} {}",
        y_scoring.avg_cost(),
        y_scoring.gap_open_cost(),
        y_scoring.gap_ext_cost()
    );
    pars_mats.fill_matrices(x_info, x_scoring, y_info, y_scoring);
    pars_mats.traceback(x_info, y_info)
}

fn pars_align(
    x_info: &[ParsimonySiteInfo],
    x_scoring: &dyn BranchParsimonyCosts,
    y_info: &[ParsimonySiteInfo],
    y_scoring: &dyn BranchParsimonyCosts,
) -> (Vec<ParsimonySiteInfo>, Alignment, f64) {
    pars_align_w_rng(x_info, x_scoring, y_info, y_scoring, rng_len)
}

pub fn pars_align_on_tree(
    scoring: &dyn ParsimonyCosts,
    info: &PhyloInfo,
) -> (Vec<Alignment>, Vec<f64>) {
    info!("Starting the IndelMAP alignment.");

    let tree = &info.tree;
    let sequences = &info.sequences;
    let sequence_type = &get_sequence_type(&info.sequences);
    let order = &tree.postorder;

    debug_assert_eq!(tree.internals.len() + tree.leaves.len(), order.len());

    let mut internal_info = vec![Vec::<ParsimonySiteInfo>::new(); tree.internals.len()];
    let mut leaf_info = vec![Vec::<ParsimonySiteInfo>::new(); tree.leaves.len()];
    let mut alignments = vec![Alignment::empty(); tree.internals.len()];
    let mut scores = vec![0.0; tree.internals.len()];

    for &node_idx in order {
        info!(
            "Processing {}{}.",
            node_idx,
            tree.get_node_id_string(&node_idx)
        );
        match node_idx {
            Int(idx) => {
                let (x_info, x_branch) = match tree.internals[idx].children[0] {
                    Int(idx) => (&internal_info[idx], tree.internals[idx].blen),
                    Leaf(idx) => (&leaf_info[idx], tree.leaves[idx].blen),
                };
                debug!("x_info: {:?}", x_info);
                let (y_info, y_branch) = match tree.internals[idx].children[1] {
                    Int(idx) => (&internal_info[idx], tree.internals[idx].blen),
                    Leaf(idx) => (&leaf_info[idx], tree.leaves[idx].blen),
                };
                debug!("y_info: {:?}", y_info);
                info!(
                    "Aligning sequences at nodes: \n1. {}{} with branch length {} \n2. {}{} with branch length {}",
                    tree.internals[idx].children[0],
                    tree.get_node_id_string(&tree.internals[idx].children[0]),
                    x_branch,
                    tree.internals[idx].children[1],
                    tree.get_node_id_string(&tree.internals[idx].children[1]),
                    y_branch
                );
                let (info, alignment, score) = pars_align(
                    x_info,
                    scoring.get_branch_costs(x_branch),
                    y_info,
                    scoring.get_branch_costs(y_branch),
                );
                internal_info[idx] = info;
                alignments[idx] = alignment;
                scores[idx] = score;
                info!("Alignment complete with score {}.\n", score);
            }
            Leaf(idx) => {
                let pars_sets = get_parsimony_sets(&sequences[idx], sequence_type);
                leaf_info[idx] = pars_sets
                    .into_iter()
                    .map(ParsimonySiteInfo::new_leaf)
                    .collect();
                info!("Processed leaf node.\n");
            }
        }
    }
    info!("Finished IndelMAP alignment.");
    (alignments, scores)
}

#[cfg(test)]
mod parsimony_alignment_tests;
