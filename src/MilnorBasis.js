let Vector = require("./vector.js");
let multinomial = require("./combinatorics.js").multinomial;


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
        if(Array.isArray(n) && n.length > 1){
            return n.map( (i) => Q(i,p)).reduce((a,b)=>a.mult(b), MilnorBasis.unit(p));
        }
        if(Array.isArray(n)){
            n = n[0];
        }
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
        if(p === 2){
            throw "Error";
            //result.set([2*n],1);
        } else {
            result.set([[],n],1);
        }
        return result;
    }

    static bP(n,p){
        let result = new MilnorBasis(p);
        if(p === 2){
            return milnor_multiplication_even([1],n,2);
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


//console.log(milnor_multiplication([],[1,2],2));


//console.log(milnor_multiplication_even([2], [2], 2));
//console.log(milnor_multiplication([[3],[]], [[1],[]], 5).m);

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
