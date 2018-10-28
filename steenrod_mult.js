let StringifyingMap = require("./StringifyingMap");
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

function adem(a, b, c, p=2, generic=undefined){
    if(generic===undefined){
        generic = p!==2;
    }
    result = new Vector(p);
    if(!generic){
        if(b === 0){
            result.set([a],1);
            return result;
        } else if(a === 0){
            result.set([b],1);
            return result;
        } else if(a >= 2*b) {
            result.set([a,b],1);
            return result;
        }
        for(let c=0; c < 1+a/2; c++){
            if(binomial(b-c-1, a-2*c) === 1){
                if(c === 0){
                    result.set([a+b],1);
                } else {
                    result.set([a+b-c,c],1);
                }
            }
        }
        return result
    }
    // p odd
    if(a === 0 && b === 0){
        result.set([c],1);
        return result;
    }
    let A;
    let B;
    let bockstein;
    if(c === 0){
        bockstein = 0;
        A = a;
        B = b;
    } else {
        A = a;
        B = c;
        bockstein = b; // should be 0 or 1
    }
    if(A === 0){
        result.set([bockstein, B, 0],1);
        return result;
    }
    if(B === 0){
        result.set([0, A, bockstein], 1);
        return result;
    }
    if(A >= p*B+bockstein){ // admissible
        result.set([0,A,bockstein,B,0],1);
        return result;
    }
    for(let j = 0; j < 1 + a/p; j++){
        let coeff = -((A+j)%2*2-1) * binomial((B-j) * (p-1) - 1 + bockstein, A - p*j, p);
        coeff = (coeff + p) % p;
        if(coeff % p !== 0){
            if(j === 0){
                result.set([bockstein,A+B,0], coeff);
            } else {
                result.set([bockstein,A+B-j,0,j,0], coeff);
            }
        }
    }
    if(bockstein !== 0){
        for(let j = 0; j < 1 + (a - 1)/p; j++){
            let coeff = -((A+j-1)%2*2-1) * binomial((B-j) * (p-1) - 1, A - p*j - 1, p);
            coeff = (coeff % p + p) % p;
            if(coeff !== 0){
                if(j === 0){
                    result.set([0,A+B,1], coeff);
                } else {
                    result.set([0,A+B-j,1,j,0], coeff);
                }
            }
        }
    }
    return result
}

function make_mono_admissible(mono, p=2, generic=undefined){
    if(generic===undefined){
        generic = p!==2;
    }
    let result = new Vector(p);
    if(mono.length === 1){
        result.set(mono,1);
        return result;
    }
    if(!generic && mono.length === 2){
        return adem(...mono)
    }
    if(!generic){
        // check to see if admissible:
        let admissible = true;
        let j;
        for(j = 0; j < mono.length - 1; j++){
            if(mono[j] < 2*mono[j+1]){
                admissible = false;
                break
            }
        }
        if(admissible){
            result.set(mono,1);
            return result;
        }
        // else j is the first index where admissibility fails
        let y = adem(mono[j], mono[j+1]);
        for(let x of y.keys()){
            let new_vec = mono.slice(0,j).concat(x,mono.slice(j+2));
            new_vec = make_mono_admissible(new_vec);
            result.addVector(new_vec, y.get(x));
        }
        return result;
    }
    // p odd
    // check to see if admissible:
    let admissible = true;
    let j;
    for(j = 1; j < mono.length - 2; j += 2){
        if(mono[j] < mono[j+1] + p*mono[j+2]){
            admissible = false;
            break;
        }
    }
    if(admissible){
        result.set(mono,1);
        return result;
    }
    // else j is the first index where admissibility fails
    let y = adem(...mono.slice(j,j+3), p, true);
    for(let x of y.keys()){
        let new_vec_x = x.slice();
        new_vec_x[0] = mono[j-1] + x[0];
        if(mono.length >= j+3){
            new_vec_x[new_vec_x.length - 1] = mono[j + 3] + x[x.length-1];
        }
        if(new_vec_x[0] <= 1 && new_vec_x[new_vec_x.length-1] <= 1){
            let new_vec = mono.slice(0,j-1).concat(new_vec_x,mono.slice(j+4));
            new_vec = make_mono_admissible(new_vec, p, true);
            result.addVector(new_vec, y.get(x));
        }
    }
    return result
}

console.log(adem(2,2).m);
console.log( make_mono_admissible([2,2]).m);
console.log( make_mono_admissible([2,2,2]).m);
console.log( make_mono_admissible([0,2,0,1,0],7).m);
console.log("____")

console.log(adem(2,2).m);
// console.log(adem(4,2).m);
// console.log(adem(2,4).m);

console.log(adem(3,0,1, 3).m); //{(0, 3, 0, 1, 0): 1}
console.log(adem(1,0,1, 7).m); // {(0, 2, 0): 2}
console.log(adem(1,1,1, 5).m); // {(0, 2, 1): 1, (1, 2, 0): 1}
console.log(adem(1,1,2, 5).m); // {(0, 3, 1): 1, (1, 3, 0): 2}

console.log(milnor_multiplication([6],[2],3).m);


console.log("Tests");

console.log(milnor_multiplication_odd([[0,2],[5]], [[1],[1]], 5).m);   // {((0, 1, 2), (0, 1)): 4, ((0, 1, 2), (6,)): 4}
console.log(milnor_multiplication_odd([[0,2,4],[]], [[1,3],[]], 7).m); //  {((0, 1, 2, 3, 4), ()): 6}
console.log(milnor_multiplication_odd([[0,2,4],[]], [[1,5],[]], 7).m); // {((0, 1, 2, 4, 5), ()): 1}
console.log(milnor_multiplication_odd([[],[6]], [[],[2]], 3).m);     //  {((), (0, 2)): 1, ((), (4, 1)): 1, ((), (8,)): 1}


console.log("2: ")
console.log(milnor_multiplication([2],[1],2).m); // {(0, 1): 1, (3,): 1}
console.log(milnor_multiplication([4],[2,1],2).m); // {(0, 3): 1, (2, 0, 1): 1, (6, 1): 1}
console.log(milnor_multiplication([2,4],[0,1],2).m); //  {(2, 0, 0, 1): 1, (2, 5): 1}




console.log(multinomial([1,2,4], 2)) // 1
console.log(multinomial([1,2,4], 7)) //0
console.log(multinomial([1,2,4], 11)) //6
//console.log(multinomial_odd([1,2,4], 101)) //4
//console.log(multinomial_odd([1,2,4], 107)) // 105
