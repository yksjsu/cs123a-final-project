use pyo3::prelude::*;
use std::alloc::{alloc, dealloc, Layout};


const NUM_DIRECTIONS: usize = 3;


/// Formats the sum of two numbers as string.
#[pyfunction]
fn check_import() {
    println!("Rust module imported!");
}


/// Match-Mismatch-Linear Gap alignment
#[pyfunction]
fn linear_gap_alignment(seq1: &str, seq2: &str, match_score: i32, mismatch_penalty: i32, gap_penalty: i32, is_local: bool) -> Vec<Vec<(i32, [bool; NUM_DIRECTIONS])>> {
    // prepare individual ASCII-character sequences
    let seq1_bytes= seq1.as_bytes();
    let seq2_bytes= seq2.as_bytes();

    scoring_matrix_linear_gap(seq1_bytes, seq2_bytes, score_comparison, match_score, mismatch_penalty, gap_penalty, is_local)
}

/// Match-Mismatch-Affine Gap Alignment
#[pyfunction]
fn affine_gap_alignment(seq1: &str, seq2: &str, match_score: i32, mismatch_penalty: i32, gap_open_penalty: i32, gap_extend_penalty: i32, is_local: bool) -> Vec<Vec<(i32, [bool; NUM_DIRECTIONS])>> {
    // prepare individual ASCII-character sequences
    let seq1_bytes= seq1.as_bytes();
    let seq2_bytes= seq2.as_bytes();

    scoring_matrix_affine_gap(seq1_bytes, seq2_bytes, score_comparison, match_score, mismatch_penalty, gap_open_penalty, gap_extend_penalty, is_local)
}


/// Comparison function for linear gap alignment. It
/// expects a 3-element array, such that the first element
/// is the "left" score, the second is the "diagonal" score, and
/// the third is the "top" score.
/// 
/// Returns a 3-element boolean array, with each element
/// being True if its respective direction (i.e. first element
/// corresponds to "left", second element corresponds to "diagonal", etc.) is
/// the direction that is "maximum."
fn score_comparison(list: &[i32; NUM_DIRECTIONS], is_local: bool) -> (i32, [bool; NUM_DIRECTIONS]) {
    // first find maximum
    let mut max;

    max = if list[0] > list[1] {list[0]} else {list[1]};
    max = if max > list[2] {max} else {list[2]};

    // if local alignment is specified, compare maximum against "0"
    if is_local {
        max = if max > 0 {max} else {0};
    }

    // now, determine which directions the alignments will go
    let mut directions = [false; NUM_DIRECTIONS];

    if max == 0 {
        return (max, directions)
    }

    for i in 0..list.len() {
        if list[i] == max {directions[i] = true};
    }

    (max, directions)
}


/// Build scoring matrix for either global or local alignment
/// with an affine gap scoring system.
fn scoring_matrix_affine_gap(seq1: &[u8], seq2: &[u8], max_function: fn(&[i32; NUM_DIRECTIONS], bool) -> (i32, [bool; NUM_DIRECTIONS]),
match_score: i32, mismatch_penalty: i32, gap_open_penalty: i32, gap_extend_penalty: i32, is_local: bool) -> Vec<Vec<(i32, [bool; NUM_DIRECTIONS])>> {
    // prepare sizes
    let n = seq1.len() + 1;
    let m = seq2.len() + 1;

    // the layouts for each of the "three matrices":
    // one for when a match/mismatch is forced between seq1 and seq2,
    // one for when seq2 opens (or extends) a gap for seq1, and one for
    // when seq1 opens (or extends) a gap for seq2
    let score_matrix_layout = Layout::array::<i32>(m).unwrap();
    let force_match_layout = Layout::array::<i32>(m).unwrap();
    let seq1_to_seq2_gap_layout = Layout::array::<i32>(m).unwrap();
    let seq2_to_seq1_gap_layout = Layout::array::<i32>(m).unwrap();

    // allocate arrays
    let score_matrix_ptr = unsafe{alloc(score_matrix_layout)} as *mut i32;
    let force_match_ptr = unsafe {alloc(force_match_layout)} as *mut i32;
    let seq1_to_seq2_gap_ptr = unsafe {alloc(seq1_to_seq2_gap_layout)} as *mut i32;
    let seq2_to_seq1_gap_ptr = unsafe {alloc(seq2_to_seq1_gap_layout)} as *mut i32;

    // create vector of vectors
    let mut vectors = Vec::with_capacity(n);

    // initialize first row of penalties
    unsafe {
        // create first vector row
        let mut vec_to_add = Vec::with_capacity(m);

        // create very first cell (0)
        vec_to_add.push((0, [false, false, false]));
        *score_matrix_ptr.add(0) = 0;

        // set the memory pointer values to the first row of penalties
        for j in 1..m {
            if is_local {
                // affine local alignment means we still want the initial rows
                // to be 0
                *score_matrix_ptr.add(j) = 0;
                vec_to_add.push((0, [false, false, false]));
            } else {
                // calculate the current score for the first row
                let value = gap_open_penalty + (j as i32 - 1) * gap_extend_penalty;

                // set the value to the ptr
                *score_matrix_ptr.add(j) = value;
                vec_to_add.push((value, [true, false, false]));
            }
        }
        
        vectors.push(vec_to_add);
    }

    // create mutable array of force_mismatch, seq1_to_seq2_gap, and
    // seq2_to_seq1_gap scores, respectively
    let mut scores = [0, 0, 0];

    // temp variables for current max and the prev
    // score[0] value in the force_mismatch matrix
    let mut curr;
    let mut prev_score_matrix;

    // obtain the optimal alignment score
    // calculate row-by-row dynamic programming matrix
    for i in 1..n {
        prev_score_matrix = if is_local {0} else {gap_open_penalty + (i as i32 - 1) * gap_extend_penalty};
        
        // prepare a vector to add to vectors
        let mut vec_to_add = Vec::with_capacity(m);

        vec_to_add.push(if is_local {(0, [false, false, false])} else {(prev_score_matrix, [false, false, true])});

        for j in 1..m {
            // set the pointer at current position in force match
            unsafe {
                *seq1_to_seq2_gap_ptr.add(j) = if i == 1 {*score_matrix_ptr.add(j)} else {
                    if *seq1_to_seq2_gap_ptr.add(j) > *score_matrix_ptr.add(j) {
                        *seq1_to_seq2_gap_ptr.add(j) + gap_extend_penalty
                    } else {
                        *score_matrix_ptr.add(j) + gap_open_penalty
                    }
                };

                *force_match_ptr.add(j) = *score_matrix_ptr.add(j - 1) + if seq1[i - 1] == seq2[j - 1] {match_score} else {mismatch_penalty};

                *seq2_to_seq1_gap_ptr.add(j) = if j == 1 {prev_score_matrix} else {
                    if *seq2_to_seq1_gap_ptr.add(j - 1) > prev_score_matrix {
                        *seq2_to_seq1_gap_ptr.add(j - 1) + gap_extend_penalty
                    } else {
                        prev_score_matrix + gap_open_penalty
                    }
                };

                scores[0] = *seq1_to_seq2_gap_ptr.add(j);
                scores[1] = *force_match_ptr.add(j);
                scores[2] = *seq2_to_seq1_gap_ptr.add(j);
            }


            // get max score and directions
            let output_tup = max_function(&scores, is_local);
            curr = output_tup.0;

            // push current score and directions to the current row
            vec_to_add.push(output_tup);

            // write down the "prev" score
            unsafe {
                *score_matrix_ptr.add(j) = curr;
                *score_matrix_ptr.add(j - 1) = prev_score_matrix;
            }

            // set new "prev" score
            prev_score_matrix = curr;
        }

        // add new vector to overall vector of vectors
        vectors.push(vec_to_add);

        // update the last cell with the "left" score in the "previous row"
        unsafe {
            *score_matrix_ptr.add(m - 1) = prev_score_matrix;
        }
    }

    // deallocate all unsafe memory
    unsafe {
        dealloc(score_matrix_ptr as *mut u8, score_matrix_layout);
        dealloc(force_match_ptr as *mut u8, force_match_layout);
        dealloc(seq1_to_seq2_gap_ptr as *mut u8, seq1_to_seq2_gap_layout);
        dealloc(seq2_to_seq1_gap_ptr as *mut u8, seq2_to_seq1_gap_layout);
    }

    vectors
}


/// Build scoring matrix for either global or local alignment
/// with a linear gap scoring system.
fn scoring_matrix_linear_gap(seq1: &[u8], seq2: &[u8], max_function: fn(&[i32; NUM_DIRECTIONS], bool) -> (i32, [bool; NUM_DIRECTIONS]),
match_score: i32, mismatch_penalty: i32, gap_penalty: i32, is_local: bool) -> Vec<Vec<(i32, [bool; NUM_DIRECTIONS])>> {
    // prepare sizes
    let n = seq1.len() + 1;
    let m = seq2.len() + 1;

    // prev and curr layouts and arrays
    let prev_layout = Layout::array::<i32>(m).unwrap();

    // allocate arrays
    let prev_ptr = unsafe {alloc(prev_layout)} as *mut i32;

    // create vector of vectors
    let mut vectors = Vec::with_capacity(n);

    // initialize first row of penalties
    unsafe {
        // create first vector row
        let mut vec_to_add = Vec::with_capacity(m);

        // create very first cell (0)
        vec_to_add.push((0, [false, false, false]));
        *prev_ptr.add(0) = 0;

        // set the memory pointer values to the first row of penalties
        for j in 1..m {
            if is_local {
                *prev_ptr.add(j) = 0;
                vec_to_add.push((0, [false, false, false]));
            } else {
                // calculate the current score for the first row
                let value = j as i32 * gap_penalty;

                // set the value to the ptr
                *prev_ptr.add(j) = value;
                vec_to_add.push((value, [true, false, false]));
            }
        }

        vectors.push(vec_to_add);
    }

    // create mutable array of left, diag, and top
    // scores, respectively
    let mut scores = [0, 0, 0];


    // temp variables for current max
    // and the prev score[0] value
    let mut curr;
    let mut prev;
    
    // obtain the optimal alignment score
    // calculate row-by-row dynamic programming matrix
    for i in 1..n {
        // "left" value
        scores[0] = if is_local {0} else {i as i32 * gap_penalty};
        
        // prepare a vector to add to vectors
        let mut vec_to_add = Vec::with_capacity(m);
        
        // prepare gap penalty at current row, first column
        let obj_to_add = if is_local {(scores[0], [false, false, false])} else {(scores[0], [false, false, true])};
        vec_to_add.push(obj_to_add);

        prev = scores[0];
        scores[0] += gap_penalty;

        // calculate scores
        for j in 1..m {
            unsafe {
                // get diagonal score
                scores[1] = *prev_ptr.add(j - 1) + (if seq1[i - 1] == seq2[j - 1] {match_score} else {mismatch_penalty});
    
                // get top score
                scores[2] = *prev_ptr.add(j) + gap_penalty;
            }

            // get max score and directions
            let output_tup = max_function(&scores, is_local);
            curr = output_tup.0;

            // push current score and directions to the current row
            vec_to_add.push(output_tup);

            // write down the "prev" score
            unsafe {
                *prev_ptr.add(j - 1) = prev;
            }

            // set new "prev" score
            prev = curr;

            // update the "left" score to the current score
            scores[0] = curr + gap_penalty;
        }

        // add new vector to overall vector of vectors
        vectors.push(vec_to_add);

        // update the last cell with the "left" score in the "previous row"
        unsafe {
            *prev_ptr.add(m - 1) = prev;
        }
    }

    // deallocate all unsafe memory
    unsafe {
        dealloc(prev_ptr as *mut u8, prev_layout);
    }

    vectors
}

/// A Python module implemented in Rust.
#[pymodule]
fn biorust(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(check_import, m)?)?;
    m.add_function(wrap_pyfunction!(linear_gap_alignment, m)?)?;
    m.add_function(wrap_pyfunction!(affine_gap_alignment, m)?)?;
    Ok(())
}
