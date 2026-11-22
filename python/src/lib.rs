//! The Python binding of LOSAT's library API: `pylosat._pylosat`, which `pylosat/_api.py` wraps.
//!
//! This crate hands `losat::api` its inputs and hands Python what comes back. It computes no
//! BLAST-defined value and writes no report (`PD-LOSAT-PYTHON-BINDING-BOUNDARY`); what NCBI does
//! not have, Arrow and pickle, stays here and out of `LOSAT/`.

mod search;
mod table;

use losat::api::{self, Program, ValueKind};
use pyo3::prelude::*;
use pyo3::types::PyBytes;

pyo3::create_exception!(
    pylosat,
    LosatError,
    pyo3::exceptions::PyException,
    "A failure that LOSAT reports. The message is the command line's text for the same failure, \
     with one U+FFFD per maximal ill-formed subpart of its bytes; `message_bytes` holds the bytes."
);
pyo3::create_exception!(
    pylosat,
    ArgumentError,
    LosatError,
    "An option that the command line's parser refuses, or one the library does not take."
);
pyo3::create_exception!(
    pylosat,
    InputError,
    LosatError,
    "A FASTA input that the command line's reader refuses or cannot open."
);
pyo3::create_exception!(
    pylosat,
    SearchError,
    LosatError,
    "A search that fails, with the command line's message."
);

/// The Python exception of an API error: its message is the error's text (`Display`: one U+FFFD
/// per maximal ill-formed subpart), and its attribute `message_bytes` the error's bytes
/// (`Error::message`), what the command line writes.
pub(crate) fn exception(error: api::Error) -> PyErr {
    let text = error.to_string();
    let raised = match &error {
        api::Error::Argv(_) => ArgumentError::new_err(text),
        api::Error::Input(_) => InputError::new_err(text),
        api::Error::Search(_) => SearchError::new_err(text),
        _ => LosatError::new_err(text),
    };
    Python::attach(|py| {
        let bytes = PyBytes::new(py, error.message());
        match raised.value(py).setattr("message_bytes", bytes) {
            Ok(()) => raised,
            Err(failed) => failed,
        }
    })
}

/// One option of a program: flag, help, whether it takes a value, the kind of the value, the
/// default of the parser, the listed values.
type ParameterRow = (
    String,
    String,
    bool,
    &'static str,
    Option<String>,
    Vec<String>,
);

/// The name that `pylosat` gives a `ValueKind`.
fn kind_name(kind: ValueKind) -> &'static str {
    match kind {
        ValueKind::Flag => "flag",
        ValueKind::Boolean => "boolean",
        ValueKind::Integer => "integer",
        ValueKind::Real => "real",
        ValueKind::Path => "path",
        _ => "text",
    }
}

/// The options of `program` from the engine's own table, and the genetic codes it accepts.
#[pyfunction]
fn describe(program: &str) -> PyResult<(Vec<ParameterRow>, Vec<u32>, Vec<u32>)> {
    let description = api::describe(Program::parse(program).map_err(exception)?);
    let parameters = description
        .parameters
        .into_iter()
        .map(|parameter| {
            (
                parameter.flag,
                parameter.help,
                parameter.takes_value,
                kind_name(parameter.value_kind),
                parameter.default,
                parameter.choices,
            )
        })
        .collect();
    Ok((
        parameters,
        description.query_gencodes,
        description.subject_gencodes,
    ))
}

/// Checks an argv with the command line's parser.
#[pyfunction]
fn validate(program: &str, argv: Vec<String>) -> PyResult<()> {
    search::checked(program, &argv).map_err(exception)?;
    Ok(())
}

/// Whether a checked argv has `-dryrun` (`api::Search::dry_run`).
#[pyfunction]
fn dry_run(program: &str, argv: Vec<String>) -> PyResult<bool> {
    Ok(search::checked(program, &argv)
        .map_err(exception)?
        .dry_run())
}

#[pymodule]
fn _pylosat(module: &Bound<'_, PyModule>) -> PyResult<()> {
    let py = module.py();
    module.add("LosatError", py.get_type::<LosatError>())?;
    module.add("ArgumentError", py.get_type::<ArgumentError>())?;
    module.add("InputError", py.get_type::<InputError>())?;
    module.add("SearchError", py.get_type::<SearchError>())?;
    module.add("ENGINE_VERSION", api::VERSION)?;
    module.add_class::<table::Hsp>()?;
    module.add_class::<table::HspTable>()?;
    module.add_class::<search::Batches>()?;
    module.add_function(wrap_pyfunction!(describe, module)?)?;
    module.add_function(wrap_pyfunction!(validate, module)?)?;
    module.add_function(wrap_pyfunction!(dry_run, module)?)?;
    module.add_function(wrap_pyfunction!(search::run, module)?)?;
    module.add_function(wrap_pyfunction!(search::run_to_file, module)?)?;
    module.add_function(wrap_pyfunction!(search::makeblastdb, module)?)?;
    module.add_function(wrap_pyfunction!(search::blastdbcmd, module)?)?;
    module.add_function(wrap_pyfunction!(search::dustmasker, module)?)?;
    module.add_function(wrap_pyfunction!(search::segmasker, module)?)?;
    module.add_function(wrap_pyfunction!(search::windowmasker, module)?)?;
    module.add_function(wrap_pyfunction!(search::convert2blastmask, module)?)?;
    module.add_function(wrap_pyfunction!(search::makembindex, module)?)?;
    module.add_function(wrap_pyfunction!(search::blast_formatter, module)?)?;
    Ok(())
}
