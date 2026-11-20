//! Tests of `lower_bound_by`: the replay of the observed calls of `std::lower_bound`
//! (`LOSAT/tests/unit/lower_bound_observed/README.md`).

use super::lower_bound_by;

/// A tree of `probe_trees.txt`: the index of a call, then the tree after the answer false, then
/// the tree after the answer true; or the position returned.
enum ProbeTree {
    Call {
        index: usize,
        if_false: Box<ProbeTree>,
        if_true: Box<ProbeTree>,
    },
    Returned(usize),
}

fn parse_probe_tree<'a>(tokens: &mut impl Iterator<Item = &'a str>) -> ProbeTree {
    let token = tokens.next().expect("a token of the tree");
    match token.strip_prefix('p') {
        Some(index) => {
            let index = index.parse().unwrap();
            let if_false = Box::new(parse_probe_tree(tokens));
            let if_true = Box::new(parse_probe_tree(tokens));
            ProbeTree::Call {
                index,
                if_false,
                if_true,
            }
        }
        None => ProbeTree::Returned(token.strip_prefix('r').unwrap().parse().unwrap()),
    }
}

/// Every path of a tree: its answers, the indices of its calls, the position returned.
type ProbePath = (Vec<bool>, Vec<usize>, usize);

fn probe_paths(
    tree: &ProbeTree,
    answers: &mut Vec<bool>,
    calls: &mut Vec<usize>,
    out: &mut Vec<ProbePath>,
) {
    match tree {
        ProbeTree::Returned(position) => out.push((answers.clone(), calls.clone(), *position)),
        ProbeTree::Call {
            index,
            if_false,
            if_true,
        } => {
            calls.push(*index);
            for (answer, subtree) in [(false, if_false), (true, if_true)] {
                answers.push(answer);
                probe_paths(subtree, answers, calls, out);
                answers.pop();
            }
            calls.pop();
        }
    }
}

// For every sequence of answers of the recorded trees, `lower_bound_by` asks about the indices
// that `std::lower_bound` gave its comparator, in the same order, and returns the same position.
#[test]
fn lower_bound_asks_in_the_observed_order() {
    let trees = include_str!("../../tests/unit/lower_bound_observed/probe_trees.txt");
    let (mut tree_count, mut path_count) = (0, 0);
    for line in trees.lines() {
        let mut tokens = line.split(' ');
        assert_eq!(tokens.next(), Some("tree"), "{line}");
        let n: usize = tokens.next().unwrap().parse().unwrap();
        let start: usize = tokens.next().unwrap().parse().unwrap();
        let total: usize = tokens.next().unwrap().parse().unwrap();
        assert!(start + n <= total, "{line}");
        let tree = parse_probe_tree(&mut tokens);
        assert_eq!(tokens.next(), None, "{line}");
        let mut paths = Vec::new();
        probe_paths(&tree, &mut Vec::new(), &mut Vec::new(), &mut paths);
        assert_eq!(paths.len(), n + 1, "{line}");
        for (answers, calls, returned) in paths {
            let mut asked = Vec::new();
            let position = lower_bound_by(start, start + n, |i| {
                asked.push(i);
                answers.get(asked.len() - 1).copied().unwrap_or(false)
            });
            assert_eq!(
                (asked, position),
                (calls, returned),
                "{n} from {start}, answers {answers:?}"
            );
            path_count += 1;
        }
        tree_count += 1;
    }
    assert_eq!((tree_count, path_count), (137, 5535));
}
