use anyhow::Error;

pub mod parsimony_alignment;

type Result<T> = std::result::Result<T, Error>;

#[allow(non_camel_case_types)]
type f64_h = ordered_float::OrderedFloat<f64>;

pub(crate) fn cmp_f64() -> impl Fn(&f64, &f64) -> std::cmp::Ordering {
    |a, b| a.partial_cmp(b).unwrap()
}

#[cfg(test)]
pub(crate) fn assert_float_relative_slice_eq(actual: &[f64], expected: &[f64], epsilon: f64) {
    use approx::relative_eq;
    assert_eq!(
        actual.len(),
        expected.len(),
        "Must have the same number of entries."
    );
    for (i, (&act, &exp)) in actual.iter().zip(expected.iter()).enumerate() {
        assert!(
            relative_eq!(exp, act, epsilon = epsilon),
            "Entries at position {} do not match",
            i,
        );
    }
}
