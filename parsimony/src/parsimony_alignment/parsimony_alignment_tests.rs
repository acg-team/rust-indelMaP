use super::parsimony_info::SiteFlag::{self, GapOpen, NoGap};
use crate::parsimony_alignment::parsimony_costs::{
    parsimony_costs_simple::ParsimonyCostsSimple, ParsimonyCosts,
};
use crate::parsimony_alignment::{
    pars_align_on_tree, pars_align_w_rng, parsimony_info::ParsimonySiteInfo,
    parsimony_sets::get_parsimony_sets,
};
use bio::io::fasta::Record;
use phylo::phylo_info::PhyloInfo;
use phylo::sequences::SequenceType;
use phylo::tree::{NodeIdx::Internal as I, NodeIdx::Leaf as L, Tree};

macro_rules! align {
    (@collect -) => { None };
    (@collect $l:tt) => { Some($l) };
    ( $( $e:tt )* ) => {vec![ $( align!(@collect $e), )* ]};
}

#[test]
pub(crate) fn align_two_first_outcome() {
    let mismatch_cost = 1.0;
    let gap_open_cost = 2.0;
    let gap_ext_cost = 0.5;
    let scoring = ParsimonyCostsSimple::new(mismatch_cost, gap_open_cost, gap_ext_cost);

    let sequences = [
        Record::with_attrs("A", None, b"AACT"),
        Record::with_attrs("B", None, b"AC"),
    ];
    let leaf_info1: Vec<ParsimonySiteInfo> = get_parsimony_sets(&sequences[0], &SequenceType::DNA)
        .into_iter()
        .map(ParsimonySiteInfo::new_leaf)
        .collect();
    let leaf_info2: Vec<ParsimonySiteInfo> = get_parsimony_sets(&sequences[1], &SequenceType::DNA)
        .into_iter()
        .map(ParsimonySiteInfo::new_leaf)
        .collect();
    let (_info, alignment, score) = pars_align_w_rng(
        &leaf_info1,
        &scoring.get_branch_costs(1.0),
        &leaf_info2,
        &scoring.get_branch_costs(1.0),
        |l| l - 1,
    );
    assert_eq!(score, 3.5);
    assert_eq!(alignment.map_x().len(), 4);
    assert_eq!(alignment.map_y().len(), 4);
    assert_eq!(alignment.map_x(), align!(0 1 2 3));
    assert_eq!(alignment.map_y(), align!(0 1 - -));
}

#[test]
pub(crate) fn align_two_second_outcome() {
    let mismatch_cost = 1.0;
    let gap_open_cost = 2.0;
    let gap_ext_cost = 0.5;
    let scoring = ParsimonyCostsSimple::new(mismatch_cost, gap_open_cost, gap_ext_cost);
    let sequences = [
        Record::with_attrs("A", None, b"AACT"),
        Record::with_attrs("B", None, b"AC"),
    ];
    let leaf_info1: Vec<ParsimonySiteInfo> = get_parsimony_sets(&sequences[0], &SequenceType::DNA)
        .into_iter()
        .map(ParsimonySiteInfo::new_leaf)
        .collect();
    let leaf_info2: Vec<ParsimonySiteInfo> = get_parsimony_sets(&sequences[1], &SequenceType::DNA)
        .into_iter()
        .map(ParsimonySiteInfo::new_leaf)
        .collect();
    let (_info, alignment, score) = pars_align_w_rng(
        &leaf_info1,
        &scoring.get_branch_costs(1.0),
        &leaf_info2,
        &scoring.get_branch_costs(1.0),
        |_| 0,
    );
    assert_eq!(score, 3.5);
    assert_eq!(alignment.map_x().len(), 4);
    assert_eq!(alignment.map_y().len(), 4);
    assert_eq!(alignment.map_x(), align!(0 1 2 3));
    assert_eq!(alignment.map_y(), align!(0 - -1));
}

#[test]
pub(crate) fn align_two_on_tree() {
    let mismatch_cost = 1.0;
    let gap_open_cost = 2.0;
    let gap_ext_cost = 0.5;

    let sequences = [
        Record::with_attrs("A", None, b"AACT"),
        Record::with_attrs("A", None, b"AC"),
    ];
    let mut tree = Tree::new(2, 0);
    tree.add_parent(0, L(0), L(1), 1.0, 1.0);
    tree.create_postorder();
    let info = PhyloInfo::new(tree, sequences.to_vec());

    let scoring = ParsimonyCostsSimple::new(mismatch_cost, gap_open_cost, gap_ext_cost);

    let (alignment_vec, score) = pars_align_on_tree(&Box::new(&scoring), &info);
    assert_eq!(score[Into::<usize>::into(info.tree.root)], 3.5);
    let alignment = &alignment_vec[Into::<usize>::into(info.tree.root)];
    assert_eq!(alignment.map_x().len(), 4);
    assert_eq!(alignment.map_y().len(), 4);
}

#[test]
pub(crate) fn internal_alignment_first_outcome() {
    let mismatch_cost = 1.0;
    let gap_open_cost = 2.0;
    let gap_ext_cost = 0.5;
    let scoring = ParsimonyCostsSimple::new(mismatch_cost, gap_open_cost, gap_ext_cost);

    let leaf_info1 = [
        (vec![b'A'], NoGap),
        (vec![b'C', b'A'], NoGap),
        (vec![b'C'], GapOpen),
        (vec![b'T'], GapOpen),
    ]
    .map(create_site_info);

    let leaf_info2 = [([b'G'], GapOpen), ([b'A'], NoGap)].map(create_site_info);

    let (_, alignment, score) = pars_align_w_rng(
        &leaf_info1,
        &scoring.get_branch_costs(1.0),
        &leaf_info2,
        &scoring.get_branch_costs(1.0),
        |_| 0,
    );
    assert_eq!(score, 1.0);
    assert_eq!(alignment.map_x(), align!(0 1 2 3));
    assert_eq!(alignment.map_y(), align!(0 1 - -));
}

#[allow(dead_code)]
pub(crate) fn create_site_info(
    args: (impl IntoIterator<Item = u8>, SiteFlag),
) -> ParsimonySiteInfo {
    ParsimonySiteInfo::new(args.0, args.1)
}

#[test]
pub(crate) fn internal_alignment_second_outcome() {
    let mismatch_cost = 1.0;
    let gap_open_cost = 2.0;
    let gap_ext_cost = 0.5;
    let scoring = ParsimonyCostsSimple::new(mismatch_cost, gap_open_cost, gap_ext_cost);

    let leaf_info1 = [
        (vec![b'A'], NoGap),
        (vec![b'A'], GapOpen),
        (vec![b'C'], GapOpen),
        (vec![b'T', b'C'], NoGap),
    ]
    .map(create_site_info);

    let leaf_info2 = [([b'G'], GapOpen), ([b'A'], NoGap)].map(create_site_info);

    let (_info, alignment, score) = pars_align_w_rng(
        &leaf_info1,
        &scoring.get_branch_costs(1.0),
        &leaf_info2,
        &scoring.get_branch_costs(1.0),
        |_| 0,
    );
    assert_eq!(score, 2.0);
    assert_eq!(alignment.map_x(), align!(0 1 2 3));
    assert_eq!(alignment.map_y(), align!(0 - -1));
}

#[test]
pub(crate) fn internal_alignment_third_outcome() {
    let mismatch_cost = 1.0;
    let gap_open_cost = 2.0;
    let gap_ext_cost = 0.5;
    let scoring = ParsimonyCostsSimple::new(mismatch_cost, gap_open_cost, gap_ext_cost);

    let leaf_info1 = [
        (vec![b'A'], NoGap),
        (vec![b'A'], GapOpen),
        (vec![b'C'], GapOpen),
        (vec![b'C', b'T'], NoGap),
    ]
    .map(create_site_info);

    let leaf_info2 = [(vec![b'G'], GapOpen), (vec![b'A'], NoGap)].map(create_site_info);

    let (_info, alignment, score) = pars_align_w_rng(
        &leaf_info1,
        &scoring.get_branch_costs(1.0),
        &leaf_info2,
        &scoring.get_branch_costs(1.0),
        |l| l - 1,
    );
    assert_eq!(score, 2.0);
    assert_eq!(alignment.map_x(), align!(- 0 1 2 3));
    assert_eq!(alignment.map_y(), align!(0 1 - - -));
}

#[test]
pub(crate) fn align_four_on_tree() {
    let a = 2.0;
    let b = 0.5;
    let c = 1.0;

    let sequences = [
        Record::with_attrs("A", None, b"AACT"),
        Record::with_attrs("B", None, b"AC"),
        Record::with_attrs("C", None, b"A"),
        Record::with_attrs("D", None, b"GA"),
    ];

    let mut tree = Tree::new(4, 2);
    tree.add_parent(0, L(0), L(1), 1.0, 1.0);
    tree.add_parent(1, L(2), L(3), 1.0, 1.0);
    tree.add_parent(2, I(0), I(1), 1.0, 1.0);
    tree.create_postorder();

    let info = PhyloInfo::new(tree, sequences.to_vec());

    let scoring = ParsimonyCostsSimple::new(c, a, b);

    let (alignment_vec, score) = pars_align_on_tree(&Box::new(&scoring), &info);
    // first cherry
    assert_eq!(score[0], 3.5);
    assert_eq!(alignment_vec[0].map_x().len(), 4);
    // second cherry
    assert_eq!(score[1], 2.0);
    assert_eq!(alignment_vec[1].map_x().len(), 2);
    // root, three possible alignments
    assert!(score[2] == 1.0 || score[2] == 2.0);
    if score[2] == 1.0 {
        assert_eq!(alignment_vec[2].map_x().len(), 4);
    } else {
        assert!(alignment_vec[2].map_x().len() == 4 || alignment_vec[2].map_x().len() == 5);
    }
}
