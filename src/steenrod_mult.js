let Vector = require("./vector.js");
let combinatoricsjs = require("./combinatorics.js");
let multinomial = combinatoricsjs.multinomial;
let binomial = combinatoricsjs.binomial;

function find_next_milnor_matrix(r, s, M, p){
    let found;
    let rows = r.length + 1;
    let cols = s.length + 1;
    for(let i = 1; !found && i < rows; i++){
        let sum = M[i][0];
        for(let j = 1; !found && j < cols; j++){
            let p_to_the_j = Math.pow(p,j);
            // check to see if column index j is small enough
            if(sum >= p_to_the_j){
                // now check to see if there's anything above this entry
                // to add to it
                let temp_col_sum = 0;
                for(let k = 0; k < i; k++){
                    temp_col_sum += M[k][j]
                }
                if(temp_col_sum != 0){
                    found = true;
                    for(let row = 1; row < i; row ++ ){
                        M[row][0] = r[row-1];
                        for(let col = 1; col < cols; col++){
                            M[0][col] = M[0][col] + M[row][col];
                            M[row][col] = 0;
                        }
                    }
                    for(let col = 1; col < j; col ++){
                        M[0][col] = M[0][col] + M[i][col];
                        M[i][col] = 0;
                    }
                    M[0][j] = M[0][j] - 1;
                    M[i][j] = M[i][j] + 1;
                    M[i][0] = sum - p_to_the_j;
                } else {
                    sum = sum + M[i][j] * p_to_the_j;
                }
            } else {
                sum = sum + M[i][j] * p_to_the_j;
            }
        }
    }
    return found;
}


function milnor_multiplication(r, s, p){
    let result = new Vector(p);
    let rows = r.length + 1;
    let cols = s.length + 1;
    let diags = r.length + s.length;
    // initialize matrix
    let M = Array(rows).fill().map(()=> Array(cols).fill(0));
    for(let j = 1; j < cols; j++){
        M[0][j] = s[j-1];
    }
    for(let i = 1; i < rows; i++){
        M[i][0] = r[i-1]
    }
    let found = true;
    do {
        // check diagonals
        let coeff = 1;
        let diagonal = Array(diags).fill(0);
        for(let n = 1; n <= diags; n++){
            let nth_diagonal = [];
            for(let i = Math.max(0,n-cols+1); i < Math.min(1+n,rows); i++){
                nth_diagonal.push(M[i][n-i]);
            }
            coeff *= multinomial(nth_diagonal,p);
            if(coeff === 0){
                break;
            }
            diagonal[n-1] = nth_diagonal.reduce((a, b) => a + b);
        }

        if(coeff !== 0){
            let i = diags - 1;
            while(i >= 0 && diagonal[i] == 0){
                i--;
            }
            let monomial = diagonal.slice(0,i+1);
            result.addTerm(monomial, coeff);
        }
    } while(find_next_milnor_matrix(r, s, M, p));
    return result
}


function milnor_multiplication_odd(m1, m2, p){
    let f = m2[0];
    let s = m2[1];
    // First compute Q_e0 Q_e1 ... P(r1, r2, ...) Q_f0 Q_f1 ...
    // Store results (as dictionary of pairs of tuples) in 'answer'.
    let answer = new Vector(p);
    answer.set(m1,1);
    for(let k of f){
        let old_answer = answer;
        answer = new Vector(p);
        for(let mono of old_answer.keys()){
            for(let i = 0; i < 1 + mono[1].length; i++){
                if((!mono[0].includes(k+i)) && ( i == 0 || Math.pow(p,k) <= mono[1][i-1])){
                    let q_mono = new Set(mono[0]);
                    q_mono = [...q_mono];
                    let ind;
                    if(q_mono.length > 0){
                        ind = q_mono.filter((x) => k+i<= x && x <= Math.max(...q_mono)).length;
                    } else {
                        ind = 0
                    }
                    let coeff = (2*(ind%2) - 1) * old_answer.get(mono); // (-1)^ind * old_answer[mono]
                    let lst = mono[0];
                    if(ind === 0){
                        lst.push(k+i)
                    } else {
                        lst.splice(lst.length -ind, 0, k+i);
                    }
                    q_mono = lst;
                    let p_mono = mono[1];
                    if(i>0){
                        p_mono[i-1] = p_mono[i-1] - Math.pow(p,k);
                    }

                    // The next two lines were added so that p_mono won't
                    // have trailing zeros. This makes p_mono uniquely
                    // determined by P(*p_mono).

                    while(p_mono.length>0 && p_mono[p_mono.length - 1] === 0){
                        p_mono.pop()
                    }

                    answer.set([q_mono, p_mono], (coeff % p + p)%p);
                }
            }
        }
    }
    // Now for the Milnor matrices.  For each entry '(e,r): coeff' in answer,
    // multiply r with s.  Record coefficient for matrix and multiply by coeff.
    // Store in 'result'.
    let result;
    if(s.length == 0){
        result = answer;
    } else {
        result = new Vector(p);
        for(const [e,r] of answer.keys()){
            let coeff = answer.get([e,r]);
            let prod = milnor_multiplication(r, s, p);
            for(let k of prod.keys()){
                result.addTerm([e,k], coeff*prod.get(k));
            }
        }
    }
    return result;
}


