let Vector = require("./vector.js");
let binomial = require("./combinatorics.js").binomial;


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
    let result = new SerreCartanBasis(p);
    if(mono.length === 1){
        result.set(mono,1);
        return result;
    }
    if(!generic && mono.length === 2){
        return adem(...mono);
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

class SerreCartanBasis extends Vector {
    constructor(p){
        super(p);
    }

    static Sq(n,p){
        let result = new SerreCartanBasis(p);
        if(p === 2){
            result.set(n,1);
        } else {
            result.set([n[0]%2,Math.floor(n[0]/2),0],1);
        }
        return result;
    }

    static P(n,p){
        let result = new SerreCartanBasis(p);
        if(p === 2){
            result.set([2*n[0]],1);
        } else {
            result.set([0,n[0],0],1);
        }
        return result;
    }

    static bP(n,p){
        let result = new SerreCartanBasis(p);
        if(p === 2){
            result.set([2*n[0]+1],1);
        } else {
            result.set([1,n[0],0],1);
        }
        return result;
    }

    static b(p){
        let result = new SerreCartanBasis(p);
        result.set([1],1);
        return result;
    }

    static Q(n,p){
        let result = new SerreCartanBasis(p);
        throw new Error("Not implemented");
        //return result;
    }


    mult(other_vector){
        if(other_vector.constructor === Number){
            return this.scale(other_vector);
        }
        if( this.p === 2 ){
            return super.mult((mono1,mono2) => {
                return make_mono_admissible(mono1.concat(mono2), this.p );
            }, other_vector);
        } else {
            return super.mult((mono1,mono2) => {
                if(mono1[mono1.length-1] === 1 && mono2[0] === 1){
                    return new SerreCartanBasis(p);
                } else {
                    let prod = mono1.slice(0,mono1.length-1).concat([mono1[mono1.length-1] + mono2[0]], mono2.slice(1));
                    return make_mono_admissible(prod, this.p);
                }
            }, other_vector);
        }
    }

    inspect(depth,opts){
        return this.toString();
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
                let str = k.map( (n) => `Sq${n}`).join("*");
                result.push(str);
            }
        } else {
            for(let [k,v] of this){
                let str = k.map( (n,idx) => {
                    if(idx % 2 === 0){
                        if(n === 1){
                            return "b";
                        } else {
                            return "";
                        }
                    }
                    return `P${n}`;
                }).filter((s) => s!=="").join("*");
                if(v!== 1){
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
}

module.exports = SerreCartanBasis;




//console.log(adem(2,2).m);
//console.log( make_mono_admissible([2,2]).m);

/*
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
console.log(adem(1,1,2, 5).m); // {(0, 3, 1): 1, (1, 3, 0): 2} */
//console.log(make_mono_admissible([0,1,0,1,0],3).m);
