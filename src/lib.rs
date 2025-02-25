use pyo3::prelude::*;
mod keep_cool;

#[pymodule]
fn linkars(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(keep_cool::parse_cools, m)?)?;
    Ok(())
}