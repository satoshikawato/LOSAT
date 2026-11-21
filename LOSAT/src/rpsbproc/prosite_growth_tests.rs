//! The growth model against the observation tables of
//! `LOSAT/tests/fixtures/rpsbproc_growth/` (`README.md` names the programs, the compiler, the
//! library and the dates): every row of the string tables (`string_ops.tsv`,
//! `string_ops_large.tsv`), the rows that hold a gibibyte or more at once (`string_ops_host.tsv`),
//! and every row of the vector table (`vector_ops.tsv`). The allocations of NCBI's rpsbproc 0.5.1
//! binary for edited motifs are compared in `sites_tests.rs`.

use super::super::RECORD_BYTES;
use super::*;

const STRING_TABLE: &str = include_str!("../../tests/fixtures/rpsbproc_growth/string_ops.tsv");
const STRING_LARGE: &str =
    include_str!("../../tests/fixtures/rpsbproc_growth/string_ops_large.tsv");
const STRING_HOST: &str = include_str!("../../tests/fixtures/rpsbproc_growth/string_ops_host.tsv");
const VECTOR_TABLE: &str = include_str!("../../tests/fixtures/rpsbproc_growth/vector_ops.tsv");

// The numbers of rows and operations of the tables.
const STRING_COUNTS: (usize, usize) = (7620, 45502);
const STRING_HOST_ROWS: usize = 12;
const VECTOR_COUNTS: (usize, usize) = (11237, 98501);

fn rows(table: &str) -> impl Iterator<Item = (&str, &str)> {
    table
        .lines()
        .filter(|line| !line.starts_with('#') && !line.is_empty())
        .map(|line| line.split_once('\t').expect(line))
}

/// Applies the row's operations to a string and compares each result with the observed token
/// `len/capacity/peak` (the most heap bytes live at once during the operation) or
/// `!what|len/capacity/peak` (an exception, the string unchanged). A success that held
/// `HOST_DEPENDENT_BYTES` or more is `MapEnd::HostDependent` of those bytes, and ends the row.
/// Returns the number of operations compared.
fn check_string_row(ops: &str, observed: &str) -> usize {
    let mut string = CxxString::new();
    let mut count = 0;
    for (op, token) in ops.split(' ').zip(observed.split(' ')) {
        let number: u64 = op[1..].parse().unwrap_or(0);
        let (before, before_capacity) = (string.heap_bytes(), string.capacity);
        let result = match &op[..1] {
            "R" => string.reserve(number),
            "A" => string.append(number),
            "P" => string.push_back(),
            "C" => {
                string.clear();
                Ok(())
            }
            other => panic!("{other}"),
        };
        let (exception, numbers) = match token.strip_prefix('!') {
            Some(rest) => {
                let (what, numbers) = rest.split_once('|').unwrap();
                (Some(what), numbers)
            }
            None => (None, token),
        };
        let numbers: Vec<u64> = numbers.split('/').map(|n| n.parse().unwrap()).collect();
        let [len, capacity, peak] = numbers[..] else {
            panic!("{token}");
        };
        count += 1;
        let context = format!("{ops} | {op} -> {token}");
        match exception {
            Some(what) => {
                let expected = match what {
                    "basic_string::_M_replace_aux" => LENGTH_ERROR,
                    "basic_string::_M_create" => RESERVE_LENGTH_ERROR,
                    "std::bad_alloc" => BAD_ALLOC,
                    other => panic!("{other}"),
                };
                assert_eq!(result, Err(MapEnd::Throw(expected)), "{context}");
            }
            None if peak >= HOST_DEPENDENT_BYTES => {
                assert_eq!(result, Err(MapEnd::HostDependent(peak)), "{context}");
                return count;
            }
            None => assert_eq!(result, Ok(()), "{context}"),
        }
        assert_eq!((string.len, string.capacity), (len, capacity), "{context}");
        // The peak: the buffer before, and the new one beside it when the capacity changed.
        let grown = if string.capacity != before_capacity {
            string.heap_bytes()
        } else {
            0
        };
        assert_eq!(before + grown, peak, "{context}");
    }
    count
}

#[test]
fn the_string_follows_every_observed_row() {
    let (mut row_count, mut op_count) = (0, 0);
    for table in [STRING_TABLE, STRING_LARGE] {
        for (ops, observed) in rows(table) {
            row_count += 1;
            op_count += check_string_row(ops, observed);
        }
    }
    assert_eq!((row_count, op_count), STRING_COUNTS);
}

#[test]
fn the_string_refuses_what_needs_a_gibibyte_held_at_once() {
    let mut rows_checked = 0;
    for (ops, observed) in rows(STRING_HOST) {
        check_string_row(ops, observed);
        rows_checked += 1;
    }
    assert_eq!(rows_checked, STRING_HOST_ROWS);
}

/// The vector rows: `G n` (clear, push n), `K`, `Z`, `S` (push a copy onto the stack), `O` (take
/// the top back), `N n` (the vector replaced by one of n records). Each token is
/// `capacity/size/stack entries/stack records/peak` (peak: the most records of heap live at once
/// during the operation).
fn check_vector_row(ops: &str, observed: &str) -> usize {
    let mut live = Live::new(0);
    let (mut size, mut stack) = (0u64, Vec::<u64>::new());
    let mut count = 0;
    for (op, token) in ops.split(' ').zip(observed.split(' ')) {
        let number: u64 = op[1..].parse().unwrap_or(0);
        let before = live;
        // The text on the heap that makes the allocations with the peak fail by one byte.
        let run =
            |live: &mut Live, size: &mut u64, stack: &mut Vec<u64>| -> Result<(), MatchError> {
                match &op[..1] {
                    "G" => {
                        live.grow(number)?;
                        *size = number;
                    }
                    "K" => *size = 0,
                    "Z" => *size -= 1,
                    "S" => {
                        live.push(*size)?;
                        stack.push(*size);
                    }
                    "O" => {
                        let top = *stack.last().unwrap();
                        live.pop(top)?;
                        stack.pop();
                        *size = top;
                    }
                    "N" => {
                        live.capacity = number;
                        *size = number;
                    }
                    other => panic!("{other}"),
                }
                Ok(())
            };
        let numbers: Vec<u64> = token.split('/').map(|n| n.parse().unwrap()).collect();
        let [capacity, observed_size, entries, records, peak] = numbers[..] else {
            panic!("{token}");
        };
        let context = format!("{ops} | {op} -> {token}");
        count += 1;
        // The peak is exact: with `peak` records and `text` bytes of text held, the operation
        // fails exactly when `peak * 8 + text` reaches `HOST_DEPENDENT_BYTES`.
        let limit = HOST_DEPENDENT_BYTES;
        if &op[..1] == "N" {
            // Not an operation of the model: the harness sets the capacity of a vector that
            // another path of NCBI's code made (the peak is not the model's).
            run(&mut live, &mut size, &mut stack).unwrap();
        } else if peak * RECORD_BYTES >= limit {
            let mut trial = before;
            let mut trial_size = size;
            let mut trial_stack = stack.clone();
            assert_eq!(
                run(&mut trial, &mut trial_size, &mut trial_stack),
                Err(MatchError::Exhausted),
                "{context}"
            );
            return count;
        } else {
            let allocates = peak > before.capacity + before.stacked;
            if allocates {
                let mut trial = before;
                trial.text = limit - peak * RECORD_BYTES;
                let (mut trial_size, mut trial_stack) = (size, stack.clone());
                assert_eq!(
                    run(&mut trial, &mut trial_size, &mut trial_stack),
                    Err(MatchError::Exhausted),
                    "{context}: fails at the peak"
                );
            }
            let mut trial = before;
            trial.text = limit - peak * RECORD_BYTES - 1;
            let (mut trial_size, mut trial_stack) = (size, stack.clone());
            assert_eq!(
                run(&mut trial, &mut trial_size, &mut trial_stack),
                Ok(()),
                "{context}: fits below the peak"
            );
            run(&mut live, &mut size, &mut stack).unwrap();
        }
        assert_eq!(
            (live.capacity, size, stack.len() as u64, live.stacked),
            (capacity, observed_size, entries, records),
            "{context}"
        );
    }
    count
}

#[test]
fn the_vector_follows_every_observed_row() {
    let (mut row_count, mut op_count) = (0, 0);
    for (ops, observed) in rows(VECTOR_TABLE) {
        row_count += 1;
        op_count += check_vector_row(ops, observed);
    }
    assert_eq!((row_count, op_count), VECTOR_COUNTS);
}
