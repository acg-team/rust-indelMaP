pub trait ParsimonyCosts {
    fn get_branch_costs(&self, branch_length: f64) -> Box<&dyn BranchParsimonyCosts>;
}

pub trait BranchParsimonyCosts {
    fn match_cost(&self, i: u8, j: u8) -> f64;
    fn gap_open_cost(&self) -> f64;
    fn gap_ext_cost(&self) -> f64;
    fn avg_cost(&self) -> f64;
}

pub mod parsimony_costs_simple;
pub mod parsimony_costs_model;
