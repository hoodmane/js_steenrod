let Vector = require("./vector.js");
let error = require("./errors.js");
let combinatorics = require("./combinatorics.js");
let multinomial = combinatorics.multinomial;


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
                if(temp_col_sum !== 0){
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


function milnor_multiplication_even(r, s, p){
    let result = new MilnorBasis(p);
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


function milnor_multiplication_full(m1, m2, p){
    let f = m2[0];
    let s = m2[1];
    // First compute Q_e0 Q_e1 ... P(r1, r2, ...) Q_f0 Q_f1 ...
    // Store results (as dictionary of pairs of tuples) in 'answer'.
    let answer = new MilnorBasis(p);
    answer.set(m1,1);
    for(let k of f){
        let old_answer = answer;
        answer = new MilnorBasis(p);
        for(let mono of old_answer.keys()){
            for(let i = 0; i < 1 + mono[1].length; i++){
                if((!mono[0].includes(k+i)) && ( i === 0 || Math.pow(p,k) <= mono[1][i-1])){
                    let q_mono = new Set(mono[0]);
                    q_mono = [...q_mono];
                    let ind;
                    if(q_mono.length > 0){
                        ind = q_mono.filter((x) => k+i<= x && x <= Math.max(...q_mono)).length;
                    } else {
                        ind = 0;
                    }
                    let coeff = -(2*(ind%2) - 1) * old_answer.get(mono); // (-1)^ind * old_answer[mono]
                    let lst = mono[0].slice();
                    if(ind === 0){
                        lst.push(k+i)
                    } else {
                        lst.splice(lst.length -ind, 0, k+i);
                    }
                    q_mono = lst;
                    let p_mono = mono[1].slice();
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
    if(s.length === 0){
        result = answer;
    } else {
        result = new MilnorBasis(p);
        for(const [e,r] of answer.keys()){
            let coeff = answer.get([e,r]);
            let prod = milnor_multiplication_even(r, s, p);
            for(let k of prod.keys()){
                result.addTerm([e,k], coeff*prod.get(k));
            }
        }
    }
    return result;
}

function milnor_multiplication(r, s, p){
    if(p === 2){
        return milnor_multiplication_even(r, s, p);
    } else {
        return milnor_multiplication_full(r, s, p);
    }
}

class MilnorBasis extends Vector {
    constructor(p){
        super(p);
    }

    static unit(p){
        let result = new MilnorBasis(p);
        if(p===2){
            result.set([],1);
        } else {
            result.set([[],[]],1)
        }
        return result;
    }

    static Q(n,p){
        let result = new MilnorBasis(p);
        if(Array.isArray(n)){
            for(let k of n){
                error.checkNonnegativeInteger("Q",n,k);
            }
            if(n.length > 1){
                return n.map( (i) => MilnorBasis.Q(i,p)).reduce((a,b)=>a.mult(b), MilnorBasis.unit(p));
            } else {
                n = n[0];
            }
        }
        error.checkNonnegativeInteger("Q", [n], n);
        if(p === 2){
            let m = Array(n+1).fill(0);
            m[n] = 1;
            result.set(m, 1);
        } else {
            result.set([[n],[]],1);
        }
        return result;
    }

    static Sq(n,p){
        let result = new MilnorBasis(p);
        if(!Array.isArray(n)){
            n = [n];
        }
        for(let k of n){
            error.checkNonnegativeInteger("Sq",n,k);
        }
        if(p === 2){
            result.set(n,1);
        } else {
            throw "Error";
            // result.set([[n%2],[Math.floor(n/2)]],1);
        }
        return result;
    }

    static P(n,p){
        let result = new MilnorBasis(p);
        if(!Array.isArray(n)){
            n = [n];
        }
        for(let k of n){
            error.checkNonnegativeInteger("P",n,k);
        }
        if(p === 2){
            result.set(n.map((k) => 2*k),1);
        } else {
            result.set([[],n],1);
        }
        return result;
    }

    static bP(n,p){
        let result = new MilnorBasis(p);
        if(!Array.isArray(n)){
            n = [n];
        }
        for(let k of n){
            error.checkNonnegativeInteger("bP",n,k);
        }
        if(p === 2){
            return milnor_multiplication_even([1], n.map((k) => 2*k), 2);
        } else {
            result.set([[1],n],1);
        }
        return result;
    }

    static b(p){
        let result = new MilnorBasis(p);
        if(p===2){
            result.set([1]);
        } else {
            result.set([[0],[]],1);
        }
        return result;
    }

    static pst(st, p) {
        let s = st[0];
        let t = st[1];
        error.checkNonnegativeInteger("pst",st,s);
        error.checkNonnegativeInteger("pst",st,t);
        let n = Array(t).fill(0);
        n[t-1] = Math.pow(p, s);
        if( p === 2 ){
            return MilnorBasis.Sq(n, p);
        } else {
            return MilnorBasis.P(n, p);
        }
    }

    mult(other_vector){
        if(other_vector.constructor === Number){
            return this.scale(other_vector);
        }
        return super.mult((mono1, mono2) => milnor_multiplication(mono1, mono2, this.p ), other_vector);
    }

    homogenousQ(){
        if(this.homogenous !== undefined){
            return this.homgenous;
        }
        let degs = new Set();
        let deg = 0;
        for(let [m,c] of this){
            deg = 0;
            if( this.p === 2) {
                for (let i = 0; i < m.length; i++) {
                    deg += (Math.pow(2, i + 1) - 1) * m[i];
                }
            } else {
                let Q = m[0];
                let P = m[1];
                for(let i = 0; i < Q.length; i++){
                    if(Q[i] > 0){
                        deg += 2*Math.pow(this.p,i) - 1
                    }
                }
                for(let i = 0; i < P.length; i++){
                    deg += (2*Math.pow(p, i + 1) - 2) * P[i];
                }
            }
            degs.add(deg);

        }
        this.homgenous = (degs.size <= 1);
        if(degs.size <= 1){
            this.deg = deg;
        }
        return this.homgenous;
    }

    degree(){
        if(this.homogenousQ()){
            return this.deg;
        } else {
            throw new Error("Inhomogenous element");
        }
    }

    excess(){
        let exc = Number.MAX_SAFE_INTEGER;
        for(let [m,c] of this){
            let monomial_excess;
            if(this.p === 2){
                monomial_excess = m.reduce((a,b)=>a+b);
            } else {
                monomial_excess = m[0].reduce((a,b)=>a+b) + 2*m[1].reduce((a,b)=>a+b);
            }
            if(monomial_excess < exc){
                exc = monomial_excess;
            }
        }
        return exc;
    }

    toLatex(){
        let result = [];
        if(this.p === 2){
            for(let [k,v] of this){
                let str = k.map( (n) => `Sq^{${n}}`).join("");
                result.push(str);
            }
        } else {
            for(let [k,v] of this){
                let str = k.map( (n,idx) => {
                    if(idx % 2 === 0){
                        if(n===1){
                            return "b";
                        } else {
                            return "";
                        }
                    }
                    return `P^{${n}}`;
                }).join("");
                if(v!== 1){
                    str = v + str;
                }
                result.push(str);
            }
        }
        return result.join(" + ");
    }

    toString(){
        let result = [];
        if(this.p === 2){
            for(let [k,v] of this){
                let str = `Sq(${k.join(", ")})`;
                result.push(str);
            }
        } else {
            for(let [k,v] of this){
                let qstr = k[0].map((n) => `Q${n}`).join("*");
                let pstr = k[1].length > 0 ? `P(${k[1].join(", ")})` : "";
                let str = [qstr, pstr].filter((s) => s!=="").join("*");
                if(v !== 1){
                    str = v + "*" + str;
                }
                result.push(str);
            }
        }
        result = result.join(" + ");
        if(result === ""){
            result = 0;
        }
        return result;
    }

    inspect(depth,opts){
        return this.toString();
    }

    toJSON(){
        return "1s@" + this.toString() + "@s1";
    }
}

module.exports = MilnorBasis;




let Infinity = Number.MAX_SAFE_INTEGER;

function milnor_basis_even(n, p, profile, trunc){
    if(n === 0) {
        return [[]];
    }
    let result = [];
    for (let mono of combinatorics.WeightedIntegerVectors(n, combinatorics.xi_degrees(n, p, true))) {
        mono.reverse();
        let exponents = mono.slice();
        while (exponents.length > 0 && exponents[exponents.length -1] === 0) {
            exponents.pop()
        }
        let okay = true;
        if (profile !== undefined && profile.length > 0) {
            for (let i = 0; i < exponents.length; i++) {
                if (profile.length > i && exponents[i] >= Math.pow(p, profile[i])) {
                    okay = false;
                    break;
                }
                if (profile.length <= i && trunc < Infinity && exponents[i] >= Math.pow(p, trunc)) {
                    okay = false;
                    break
                }
            }
        } else {
            okay = (! trunc || trunc === Infinity);
        }
        if (okay) {
            result.push(exponents);
        }
    }
    return result;
}



function milnor_basis(n, p=2, kwds = {}){
    /*
    Milnor basis in dimension `n` with profile function ``profile``.

    INPUT:

    - ``n`` - non-negative integer

    - ``p`` - positive prime number (optional, default 2)

    - ``profile`` - profile function (optional, default None).
      Together with ``truncation_type``, specify the profile function
      to be used; None means the profile function for the entire
      Steenrod algebra.  See
      :mod:`sage.algebras.steenrod.steenrod_algebra` and
      :func:`SteenrodAlgebra <sage.algebras.steenrod.steenrod_algebra.SteenrodAlgebra>`
      for information on profile functions.

    - ``truncation_type`` - truncation type, either 0 or Infinity
      (optional, default Infinity if no profile function is specified,
      0 otherwise)

    OUTPUT: tuple of mod p Milnor basis elements in dimension n

    At the prime 2, the Milnor basis consists of symbols of the form
    `\text{Sq}(m_1, m_2, ..., m_t)`, where each
    `m_i` is a non-negative integer and if `t>1`, then
    `m_t \neq 0`. At odd primes, it consists of symbols of the
    form `Q_{e_1} Q_{e_2} ... P(m_1, m_2, ..., m_t)`,
    where `0 \leq e_1 < e_2 < ...`, each `m_i` is a
    non-negative integer, and if `t>1`, then
    `m_t \neq 0`.

    EXAMPLES::

        sage: from sage.algebras.steenrod.steenrod_algebra_bases import milnor_basis
        sage: milnor_basis(7)
        ((0, 0, 1), (1, 2), (4, 1), (7,))
        sage: milnor_basis(7, 2)
        ((0, 0, 1), (1, 2), (4, 1), (7,))
        sage: milnor_basis(4, 2)
        ((1, 1), (4,))
        sage: milnor_basis(4, 2, profile=[2,1])
        ((1, 1),)
        sage: milnor_basis(4, 2, profile=(), truncation_type=0)
        ()
        sage: milnor_basis(4, 2, profile=(), truncation_type=Infinity)
        ((1, 1), (4,))
        sage: milnor_basis(9, 3)
        (((1,), (1,)), ((0,), (2,)))
        sage: milnor_basis(17, 3)
        (((2,), ()), ((1,), (3,)), ((0,), (0, 1)), ((0,), (4,)))
        sage: milnor_basis(48, p=5)
        (((), (0, 1)), ((), (6,)))
        sage: len(milnor_basis(100,3))
        13
        sage: len(milnor_basis(200,7))
        0
        sage: len(milnor_basis(240,7))
        3
        sage: len(milnor_basis(240,7, profile=((),()), truncation_type=Infinity))
        3
        sage: len(milnor_basis(240,7, profile=((),()), truncation_type=0))
        0
    */
    let generic = kwds['generic'] || p !== 2;

    if(n === 0){
        if(!generic){
            return [[]];
        } else {
            return [[[],[]]];
        }
    }
    let profile = kwds["profile"];
    let trunc = kwds["truncation_type"];
    if(trunc === undefined){
        if(profile !== undefined){
            trunc = 0;
        } else {
            trunc = Infinity;
        }
    }
    let result = [];
    if(!generic) {
        return milnor_basis_even(n, 2, profile, trunc);
    } else {  // p odd
        // first find the P part of each basis element.
        // in this part of the code (the P part), all dimensions are
        // divided by 2(p-1).
        for (let dim = 0; dim < n / (2 * (p - 1)) + 1; dim++) {
            let P_result = milnor_basis_even(dim, p, profile && profile[0], trunc);
            P_result = P_result.map( (x) => x.length === 0 ? [0] : x);
            // now find the Q part of the basis element.
            // dimensions here are back to normal.
            for (let p_mono of P_result) {
                let deg = n - 2 * dim * (p - 1);
                let q_degrees = combinatorics.xi_degrees(Math.floor((deg - 1)/(2*(p-1))), p).map( d =>  1+2*(p-1)*d);
                q_degrees.push(1);
                let q_degrees_decrease = q_degrees;
                q_degrees.reverse();
                if (deg % (2 * (p - 1)) <= q_degrees.length) {
                    // if this inequality fails, no way to have a partition
                    // with distinct parts.
                    for (let sigma of combinatorics.restricted_partitions(deg, q_degrees_decrease, true)) {
                        let q_mono = [];
                        for(let index = 0; index < q_degrees.length; index++){
                            if(sigma.includes(q_degrees[index])){
                                q_mono.push(index)
                            }
                        }
                        // check profile:
                        let okay = true;
                        if (profile !== undefined && profile[1].length > 0){
                            // check profile function for q_mono
                            for(let i of q_mono) {
                                if ((profile[1].length > i && profile[1][i] === 1) || (profile[1].length <= i && trunc === 0)) {
                                    okay = false;
                                    break
                                }
                            }
                        } else {
                            // profile is empty
                            okay = (trunc === Infinity);
                        }
                        if (okay) {
                            if (p_mono.length === 1 && p_mono[0] === 0) {
                                p_mono = [];
                            }
                            result.push([q_mono, p_mono]);
                        }
                    }
                }
            }
        }
    }
    return result;
}


/*
console.log(milnor_basis(7));
console.log([[0, 0, 1], [1, 2], [4, 1], [7]]);
console.log(milnor_basis(7, 2))
console.log([[0, 0, 1], [1, 2], [4, 1], [7]])
console.log( milnor_basis(4, 2))
console.log([[1, 1], [4]])
console.log( milnor_basis(4, 2, {profile: [2,1]}));
console.log([[1, 1]])
console.log( milnor_basis(4, 2, { profile : [], truncation_type : 0}))
console.log([])
console.log( milnor_basis(4, 2, { profile : [], truncation_type : Infinity}))
console.log([[1, 1], [4]]);*/
// console.log(milnor_basis_even(2,3));
//
// console.log( milnor_basis(9, 3));
// console.log([[[1], [1]], [[0], [2]]]);
// console.log( milnor_basis(17, 3));
// console.log([[[2], []], [[1], [3]], [[0], [0, 1]], [[0], [4]]]);
// console.log( milnor_basis(48, p=5));
// console.log([[[], [0, 1]], [[], [6]]]);
// console.log( milnor_basis(100,3).length);
// console.log(13);
// console.log( milnor_basis(200,7).length);
// console.log(0);
// console.log(milnor_basis(240,7).length);
// console.log(3);
// console.log(milnor_basis(240,7, { profile : [[],[]], truncation_type : Infinity}).length);
// console.log(3);
// console.log(milnor_basis(240,7, { profile : [[],[]], truncation_type : 0}).length);
// console.log(0);

// console.log(milnor_multiplication([],[1,2],2));
//
//
// console.log(milnor_multiplication_even([2], [2], 2));
// console.log(milnor_multiplication([[3],[]], [[1],[]], 5).m);

/*
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
//console.log(multinomial_odd([1,2,4], 107)) // 105*!/*/
