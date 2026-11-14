//! The HSPs of a search as one Arrow record batch.
//!
//! The batch is the only store of the results on this side: `HspTable` hands it to Arrow's
//! consumers without a copy, an `Hsp` is a row of it, and pickle writes it as an Arrow IPC stream.
//! Its columns are named as `-outfmt` names them, the default report's first (the `HspRecord`
//! fields that have no column keep their names); an `Hsp` has `HspRecord`'s field names. Each
//! value is the `HspRecord`'s, not rounded and never formatted here.

use std::io::Cursor;
use std::sync::Arc;

use arrow_array::cast::AsArray;
use arrow_array::ffi_stream::FFI_ArrowArrayStream;
use arrow_array::types::{Float64Type, Int64Type};
use arrow_array::{
    Array, ArrayRef, Float64Array, Int64Array, LargeStringArray, RecordBatch, RecordBatchIterator,
};
use arrow_ipc::reader::StreamReader;
use arrow_ipc::writer::StreamWriter;
use arrow_schema::ArrowError;
use losat::api::{HspRecord, Strand};
use pyo3::exceptions::{PyIndexError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyCapsule, PySlice, PyType};

fn ints(hsps: &[HspRecord], value: impl Fn(&HspRecord) -> usize) -> ArrayRef {
    Arc::new(Int64Array::from_iter_values(
        hsps.iter().map(|hsp| value(hsp) as i64),
    ))
}

fn optional_ints(hsps: &[HspRecord], value: impl Fn(&HspRecord) -> Option<i64>) -> ArrayRef {
    Arc::new(hsps.iter().map(value).collect::<Int64Array>())
}

fn floats(hsps: &[HspRecord], value: impl Fn(&HspRecord) -> f64) -> ArrayRef {
    Arc::new(Float64Array::from_iter_values(hsps.iter().map(value)))
}

fn texts<'a>(hsps: &'a [HspRecord], value: impl Fn(&'a HspRecord) -> Option<&'a str>) -> ArrayRef {
    Arc::new(hsps.iter().map(value).collect::<LargeStringArray>())
}

fn strand_name(strand: Option<Strand>) -> Option<&'static str> {
    strand.map(|strand| match strand {
        Strand::Plus => "plus",
        Strand::Minus => "minus",
    })
}

/// The table of `hsps`.
pub fn of_records(hsps: &[HspRecord]) -> RecordBatch {
    RecordBatch::try_from_iter_with_nullable([
        // The columns of `-outfmt 6`, in its order.
        ("qaccver", texts(hsps, |hsp| Some(&hsp.query_id)), false),
        ("saccver", texts(hsps, |hsp| Some(&hsp.subject_id)), false),
        ("pident", floats(hsps, HspRecord::percent_identity), false),
        ("length", ints(hsps, |hsp| hsp.length), false),
        (
            "mismatch",
            optional_ints(hsps, |hsp| Some(hsp.mismatch)),
            false,
        ),
        ("gapopen", ints(hsps, |hsp| hsp.gap_open), false),
        ("qstart", ints(hsps, |hsp| hsp.q_start), false),
        ("qend", ints(hsps, |hsp| hsp.q_end), false),
        ("sstart", ints(hsps, |hsp| hsp.s_start), false),
        ("send", ints(hsps, |hsp| hsp.s_end), false),
        ("evalue", floats(hsps, |hsp| hsp.e_value), false),
        ("bitscore", floats(hsps, |hsp| hsp.bit_score), false),
        // Other columns of `-outfmt`.
        (
            "score",
            optional_ints(hsps, |hsp| Some(i64::from(hsp.raw_score))),
            false,
        ),
        ("nident", ints(hsps, |hsp| hsp.num_ident), false),
        ("positive", ints(hsps, |hsp| hsp.num_positives), false),
        ("gaps", ints(hsps, |hsp| hsp.gaps), false),
        ("qlen", ints(hsps, |hsp| hsp.query_length), false),
        ("slen", ints(hsps, |hsp| hsp.subject_length), false),
        (
            "qframe",
            optional_ints(hsps, |hsp| hsp.query_frame.map(i64::from)),
            true,
        ),
        (
            "sframe",
            optional_ints(hsps, |hsp| hsp.subject_frame.map(i64::from)),
            true,
        ),
        (
            "sstrand",
            texts(hsps, |hsp| strand_name(hsp.subject_strand)),
            true,
        ),
        (
            "qseq",
            texts(hsps, |hsp| hsp.query_aligned.as_deref()),
            true,
        ),
        (
            "sseq",
            texts(hsps, |hsp| hsp.subject_aligned.as_deref()),
            true,
        ),
        // The fields of `HspRecord` without a column.
        (
            "query_strand",
            texts(hsps, |hsp| strand_name(hsp.query_strand)),
            true,
        ),
        (
            "sum_n",
            optional_ints(hsps, |hsp| hsp.sum_n.map(i64::from)),
            true,
        ),
        (
            "comp_adjust_method",
            optional_ints(hsps, |hsp| hsp.comp_adjust_method.map(i64::from)),
            true,
        ),
        ("index", ints(hsps, |hsp| hsp.index), false),
        ("q_idx", ints(hsps, |hsp| hsp.q_idx), false),
        ("s_idx", ints(hsps, |hsp| hsp.s_idx), false),
        ("rank", ints(hsps, |hsp| hsp.rank), false),
    ])
    .expect("the columns have one length")
}

/// The batch as the bytes of an Arrow IPC stream.
fn ipc_bytes(batch: &RecordBatch) -> Result<Vec<u8>, ArrowError> {
    let mut bytes = Vec::new();
    let mut writer = StreamWriter::try_new(&mut bytes, batch.schema_ref())?;
    writer.write(batch)?;
    writer.finish()?;
    drop(writer);
    Ok(bytes)
}

/// The batch that `ipc_bytes` wrote.
fn of_ipc_bytes(bytes: &[u8]) -> Result<RecordBatch, ArrowError> {
    let mut reader = StreamReader::try_new(Cursor::new(bytes), None)?;
    reader
        .next()
        .ok_or_else(|| ArrowError::IpcError("no batch in the stream".to_string()))?
}

fn value_error(error: ArrowError) -> PyErr {
    PyValueError::new_err(error.to_string())
}

/// The sorted positions where the values of a sorted column change, as `(value, start, end)`.
fn runs(values: &[i64]) -> impl Iterator<Item = (i64, usize, usize)> + '_ {
    let mut start = 0;
    values
        .chunk_by(|left, right| left == right)
        .map(move |run| {
            let this = (run[0], start, start + run.len());
            start += run.len();
            this
        })
}

/// The HSPs of a search: one row per HSP, in the order of the final hit list.
///
/// A table has `len`, indexing (an `Hsp`), slicing (a table) and iteration, and any Arrow
/// consumer reads it without a copy through `__arrow_c_stream__`: `pyarrow.table(hsps)`,
/// `polars.DataFrame(hsps)`, `pandas.DataFrame.from_arrow(hsps)`.
#[pyclass(frozen, module = "pylosat._pylosat")]
pub struct HspTable {
    batch: Arc<RecordBatch>,
}

impl HspTable {
    pub fn new(batch: RecordBatch) -> Self {
        Self {
            batch: Arc::new(batch),
        }
    }

    fn rows(&self, start: usize, end: usize) -> Self {
        Self::new(self.batch.slice(start, end - start))
    }

    fn sorted_ints(&self, name: &str) -> &[i64] {
        self.batch
            .column_by_name(name)
            .expect("a column of the table")
            .as_primitive::<Int64Type>()
            .values()
    }
}

#[pymethods]
impl HspTable {
    /// The table of an Arrow IPC stream of a table (pickle).
    #[new]
    fn from_ipc(bytes: &[u8]) -> PyResult<Self> {
        let batch = of_ipc_bytes(bytes).map_err(value_error)?;
        if batch.schema() != of_records(&[]).schema() {
            return Err(PyValueError::new_err(
                "the stream does not hold a table of this version of pylosat",
            ));
        }
        Ok(Self::new(batch))
    }

    /// The column names.
    #[getter]
    fn columns(&self) -> Vec<String> {
        let schema = self.batch.schema();
        schema
            .fields()
            .iter()
            .map(|field| field.name().clone())
            .collect()
    }

    fn __len__(&self) -> usize {
        self.batch.num_rows()
    }

    fn __getitem__(&self, py: Python<'_>, key: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
        if let Ok(slice) = key.cast::<PySlice>() {
            let indices = slice.indices(self.__len__() as isize)?;
            if indices.step != 1 {
                return Err(PyValueError::new_err("a table slices with step 1"));
            }
            let start = indices.start as usize;
            let table = self.rows(start, start + indices.slicelength);
            return Ok(table.into_pyobject(py)?.into_any().unbind());
        }
        let index: isize = key.extract()?;
        let length = self.__len__() as isize;
        let row = if index < 0 { index + length } else { index };
        if !(0..length).contains(&row) {
            return Err(PyIndexError::new_err("the table has no such row"));
        }
        let hsp = Hsp {
            batch: self.batch.clone(),
            row: row as usize,
        };
        Ok(hsp.into_pyobject(py)?.into_any().unbind())
    }

    fn __eq__(&self, other: &Self) -> bool {
        self.batch == other.batch
    }

    fn __repr__(&self) -> String {
        format!("<HspTable: {} HSPs>", self.__len__())
    }

    /// The Arrow PyCapsule interface. A requested schema may be ignored (best effort, as the
    /// interface says): the columns have one type each.
    #[pyo3(signature = (requested_schema=None))]
    fn __arrow_c_stream__<'py>(
        &self,
        py: Python<'py>,
        requested_schema: Option<Bound<'py, PyAny>>,
    ) -> PyResult<Bound<'py, PyCapsule>> {
        let _ = requested_schema;
        let schema = self.batch.schema();
        let batches = [Ok(self.batch.as_ref().clone())];
        let reader = RecordBatchIterator::new(batches, schema);
        let stream = FFI_ArrowArrayStream::new(Box::new(reader));
        PyCapsule::new_with_value(py, stream, c"arrow_array_stream")
    }

    fn __reduce__<'py>(
        &self,
        py: Python<'py>,
    ) -> PyResult<(Bound<'py, PyType>, (Bound<'py, PyBytes>,))> {
        let bytes = py.detach(|| ipc_bytes(&self.batch)).map_err(value_error)?;
        Ok((py.get_type::<Self>(), (PyBytes::new(py, &bytes),)))
    }

    /// The rows of query `q_idx`: the rows of a query lie together, in the order of the queries.
    fn _query(&self, q_idx: i64) -> Self {
        let values = self.sorted_ints("q_idx");
        let start = values.partition_point(|&value| value < q_idx);
        let end = values.partition_point(|&value| value <= q_idx);
        self.rows(start, end)
    }

    /// `(s_idx, rows)` of each subject of a query's rows, in their order.
    fn _subjects(&self) -> Vec<(i64, Self)> {
        runs(self.sorted_ints("s_idx"))
            .map(|(s_idx, start, end)| (s_idx, self.rows(start, end)))
            .collect()
    }
}

/// What pickle calls to make an `Hsp` again: the class, its table of one row, the row.
type HspReduction<'py> = (Bound<'py, PyType>, (Py<HspTable>, usize));

/// One HSP: a row of its table, with the field names of LOSAT's `HspRecord`.
#[pyclass(frozen, module = "pylosat._pylosat")]
pub struct Hsp {
    batch: Arc<RecordBatch>,
    row: usize,
}

impl Hsp {
    fn column(&self, name: &str) -> &dyn Array {
        self.batch
            .column_by_name(name)
            .expect("a column of the table")
            .as_ref()
    }

    fn int(&self, name: &str) -> i64 {
        self.column(name)
            .as_primitive::<Int64Type>()
            .value(self.row)
    }

    fn float(&self, name: &str) -> f64 {
        self.column(name)
            .as_primitive::<Float64Type>()
            .value(self.row)
    }

    fn optional_int(&self, name: &str) -> Option<i64> {
        let column = self.column(name).as_primitive::<Int64Type>();
        column.is_valid(self.row).then(|| column.value(self.row))
    }

    fn text(&self, name: &str) -> Option<String> {
        let column = self.column(name).as_string::<i64>();
        column
            .is_valid(self.row)
            .then(|| column.value(self.row).to_string())
    }

    fn required_text(&self, name: &str) -> String {
        self.column(name)
            .as_string::<i64>()
            .value(self.row)
            .to_string()
    }

    fn one_row(&self) -> RecordBatch {
        self.batch.slice(self.row, 1)
    }

    /// The 0-based half-open interval of a row's `start` and `end` on the plus strand, as
    /// `HspRecord::query_range` and `subject_range` give it.
    fn plus_interval<'py>(&self, py: Python<'py>, start: i64, end: i64) -> Bound<'py, PySlice> {
        PySlice::new(
            py,
            (start.min(end) - 1) as isize,
            start.max(end) as isize,
            1,
        )
    }
}

#[pymethods]
impl Hsp {
    /// Row `row` of `table` (pickle).
    #[new]
    fn new(table: &HspTable, row: usize) -> PyResult<Self> {
        if row >= table.batch.num_rows() {
            return Err(PyIndexError::new_err("the table has no such row"));
        }
        Ok(Self {
            batch: table.batch.clone(),
            row,
        })
    }

    #[getter]
    fn index(&self) -> i64 {
        self.int("index")
    }

    #[getter]
    fn q_idx(&self) -> i64 {
        self.int("q_idx")
    }

    #[getter]
    fn s_idx(&self) -> i64 {
        self.int("s_idx")
    }

    #[getter]
    fn rank(&self) -> i64 {
        self.int("rank")
    }

    #[getter]
    fn query_id(&self) -> String {
        self.required_text("qaccver")
    }

    #[getter]
    fn subject_id(&self) -> String {
        self.required_text("saccver")
    }

    #[getter]
    fn raw_score(&self) -> i64 {
        self.int("score")
    }

    #[getter]
    fn bit_score(&self) -> f64 {
        self.float("bitscore")
    }

    #[getter]
    fn e_value(&self) -> f64 {
        self.float("evalue")
    }

    #[getter]
    fn percent_identity(&self) -> f64 {
        self.float("pident")
    }

    #[getter]
    fn q_start(&self) -> i64 {
        self.int("qstart")
    }

    #[getter]
    fn q_end(&self) -> i64 {
        self.int("qend")
    }

    #[getter]
    fn s_start(&self) -> i64 {
        self.int("sstart")
    }

    #[getter]
    fn s_end(&self) -> i64 {
        self.int("send")
    }

    #[getter]
    fn query_strand(&self) -> Option<String> {
        self.text("query_strand")
    }

    #[getter]
    fn subject_strand(&self) -> Option<String> {
        self.text("sstrand")
    }

    #[getter]
    fn query_frame(&self) -> Option<i64> {
        self.optional_int("qframe")
    }

    #[getter]
    fn subject_frame(&self) -> Option<i64> {
        self.optional_int("sframe")
    }

    #[getter]
    fn query_length(&self) -> i64 {
        self.int("qlen")
    }

    #[getter]
    fn subject_length(&self) -> i64 {
        self.int("slen")
    }

    #[getter]
    fn length(&self) -> i64 {
        self.int("length")
    }

    #[getter]
    fn num_ident(&self) -> i64 {
        self.int("nident")
    }

    #[getter]
    fn num_positives(&self) -> i64 {
        self.int("positive")
    }

    #[getter]
    fn mismatch(&self) -> i64 {
        self.int("mismatch")
    }

    #[getter]
    fn gap_open(&self) -> i64 {
        self.int("gapopen")
    }

    #[getter]
    fn gaps(&self) -> i64 {
        self.int("gaps")
    }

    #[getter]
    fn sum_n(&self) -> Option<i64> {
        self.optional_int("sum_n")
    }

    #[getter]
    fn comp_adjust_method(&self) -> Option<i64> {
        self.optional_int("comp_adjust_method")
    }

    #[getter]
    fn query_aligned(&self) -> Option<String> {
        self.text("qseq")
    }

    #[getter]
    fn subject_aligned(&self) -> Option<String> {
        self.text("sseq")
    }

    /// The query interval as a 0-based slice on the record's plus strand.
    #[getter]
    fn query_range<'py>(&self, py: Python<'py>) -> Bound<'py, PySlice> {
        self.plus_interval(py, self.q_start(), self.q_end())
    }

    /// The subject interval as a 0-based slice on the record's plus strand.
    #[getter]
    fn subject_range<'py>(&self, py: Python<'py>) -> Bound<'py, PySlice> {
        self.plus_interval(py, self.s_start(), self.s_end())
    }

    fn __eq__(&self, other: &Self) -> bool {
        self.one_row() == other.one_row()
    }

    fn __repr__(&self) -> String {
        format!(
            "Hsp(query_id={:?}, subject_id={:?}, q_start={}, q_end={}, s_start={}, s_end={}, \
             raw_score={}, e_value={:?})",
            self.query_id(),
            self.subject_id(),
            self.q_start(),
            self.q_end(),
            self.s_start(),
            self.s_end(),
            self.raw_score(),
            self.e_value(),
        )
    }

    fn __reduce__<'py>(&self, py: Python<'py>) -> PyResult<HspReduction<'py>> {
        let table = Py::new(py, HspTable::new(self.one_row()))?;
        Ok((py.get_type::<Self>(), (table, 0)))
    }
}
