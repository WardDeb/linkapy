use pyo3::prelude::*;
mod keep_cool;

#[pymodule]
fn linkapy(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(keep_cool::parse_cools, m)?)?;
    Ok(())
}